compute_maf <- function(geno){
  f <- mean(geno,na.rm = TRUE)/2
  return(min(f, 1-f))
}

compute_missing <- function(geno){
  miss <- sum(is.na(geno))/length(geno)
  return(miss)
}

compute_non_missing_y <- function(y){
  nonmiss <- sum(!is.na(y))
  return(nonmiss)
}
  
compute_all_missing_y <- function(y){
  allmiss <- all(is.na(y))
  return(allmiss)
}

mean_impute <- function(geno){
  f <- apply(geno, 2, function(x) mean(x,na.rm = TRUE))
  for (i in 1:length(f)) geno[,i][which(is.na(geno[,i]))] <- f[i]
  return(geno)
}

is_zero_variance <- function(x) {
  if (length(unique(x))==1) return(T)
  else return(F)
}

#' @importFrom matrixStats colVars
filter_X <- function(X, missing_rate_thresh, maf_thresh, var_thresh=0) {
    rm_col <- which(apply(X, 2, compute_missing) > missing_rate_thresh)
    if (length(rm_col)) X <- X[, -rm_col]
    rm_col <- which(apply(X, 2, compute_maf) <= maf_thresh)
    if (length(rm_col)) X <- X[, -rm_col]
    rm_col <- which(apply(X, 2, is_zero_variance))
    if (length(rm_col)) X <- X[, -rm_col]
    X <- mean_impute(X)
    if (var_thresh>0) {
      rm_col <- which(matrixStats::colVars(X) < var_thresh)
      if (length(rm_col)) X <- X[, -rm_col]
    }
    return(X)
}

filter_Y <- function(Y, n_nonmiss){
  rm_col <- which(apply(Y, 2, compute_non_missing_y) < n_nonmiss)
  if (length(rm_col)) Y <- Y[, -rm_col]
  if(is.matrix(Y)){
    rm_rows <- which(apply(Y, 1, compute_all_missing_y))
    if (length(rm_rows)) Y <- Y[-rm_rows, ]  
  } else {
    Y <- Y[which(!is.na(Y))]
  }
  return(list(Y=Y, rm_rows = rm_rows))
}

load_genotype_data <- function(genotype, keep_indel = TRUE) {
  # Read genotype data using plink
  geno <- read_plink(genotype)
  # Process row names
  rownames(geno$bed) <- sapply(strsplit(rownames(geno$bed), ":"), `[`, 2)
  # Remove indels if specified
  if (!keep_indel) {
    is_indel <- with(geno$bim, grepl("[^ATCG]", V5) | grepl("[^ATCG]", V6) | nchar(V5) > 1 | nchar(V6) > 1)
    geno_bed <- geno$bed[, !is_indel]
  } else {
    geno_bed <- geno$bed
  }
  return(geno_bed)
}

filter_by_common_samples <- function(dat, common_samples) {
  dat[common_samples, , drop = FALSE] %>% .[order(rownames(.)), ]
}

prepare_data_list <- function(geno_bed, phenotype, covariate, region, imiss_cutoff, maf_cutoff, mac_cutoff, xvar_cutoff, keep_samples = NULL) {
  data_list <- tibble(
    covariate_path = covariate, 
    phenotype_path = phenotype
  ) %>%
    mutate(
      # Load covariates and transpose
      covar = map(covariate_path, ~ read_delim(.x, "\t", col_types = cols()) %>% select(-1) %>% na.omit() %>% t()),
      # Load phenotype data
      Y = map(phenotype_path, ~ {
        tabix_region(.x, region) %>% select(-4) %>% t() %>% as.matrix()
      }),
      # Determine common samples across Y, covar, and geno_bed
      common_samples = map2(covar, Y, ~ intersect(intersect(rownames(.x), rownames(.y)), rownames(geno_bed))),
      # Further intersect with keep_samples if provided
      common_samples = if (!is.null(keep_samples) && length(keep_samples) > 0) {
        map(common_samples, ~ intersect(.x, keep_samples))
      } else {
        common_samples
      },
      # Filter data based on common samples
      Y = map2(Y, common_samples, ~ filter_by_common_samples(.x, .y)),
      covar = map2(covar, common_samples, ~ filter_by_common_samples(.x, .y)),
      # Apply filter_X on the geno_bed data filtered by common samples
      X = map(common_samples, ~ {
        filtered_geno_bed <- filter_by_common_samples(geno_bed, .x)
        maf_val <- max(maf_cutoff, mac_cutoff / (2 * nrow(filtered_geno_bed)))
        filter_X(filtered_geno_bed, imiss_cutoff, maf_val, xvar_cutoff)
      }),
      # Track dropped samples
      dropped_samples_Y = map2(Y, common_samples, ~ setdiff(rownames(.x), .y)),
      dropped_samples_X = map2(X, common_samples, ~ setdiff(rownames(.x), .y)),
      dropped_samples_covar = map2(covar, common_samples, ~ setdiff(rownames(.x), .y))
    ) %>%
    select(covar, Y, X, dropped_samples_Y, dropped_samples_X, dropped_samples_covar)

  return(data_list)
}

prepare_X_matrix <- function(geno_bed, data_list, imiss_cutoff, maf_cutoff, mac_cutoff, xvar_cutoff) {
  # Calculate the union of all samples from data_list: any of X, covar and Y would do
  all_samples_union = map(data_list$covar, ~rownames(.x)) %>% unlist() %>% unique()
  # Find the intersection of these samples with the samples in geno_bed
  common_samples = intersect(all_samples_union, rownames(geno_bed))
  # Filter geno_bed using common_samples
  X_filtered = filter_by_common_samples(geno_bed, common_samples)
  # Calculate MAF cutoff considering the number of common samples
  maf_val = max(maf_cutoff, mac_cutoff / (2 * length(common_samples)))
  # Apply further filtering on X
  X_filtered = filter_X(X_filtered, imiss_cutoff, maf_val, xvar_cutoff)
  print(paste0("Dimension of input genotype data is row:", nrow(X_filtered), " column: ", ncol(X_filtered) ))
  return(X_filtered)
}


add_X_residuals <- function(data_list, scale_residuals = FALSE) {
  # Compute residuals for X and add them to data_list
  data_list <- data_list %>%
    mutate(
      lm_res_X = map2(X, covar, ~ .lm.fit(x = cbind(1, .y), y = .x)),
      X_resid_mean = map(lm_res_X, ~ apply(.x$residuals, 2, mean)),
      X_resid_sd = map(lm_res_X, ~ apply(.x$residuals, 2, sd)),
      X_resid = map(lm_res_X, ~ {
        residuals <- .x$residuals
        if (scale_residuals) {
          residuals <- scale(residuals)
        }
        residuals
      })
    )

  return(data_list)
}

add_Y_residuals <- function(data_list, conditions, y_as_matrix = FALSE, scale_residuals = FALSE) {
  # Compute residuals, their mean, and standard deviation, and add them to data_list
  data_list <- data_list %>%
    mutate(
      lm_res = map2(Y, covar, ~ .lm.fit(x = cbind(1, .y), y = .x)),
      Y_resid_mean = map(lm_res, ~ apply(.x$residuals, 2, mean)),
      Y_resid_sd = map(lm_res, ~ apply(.x$residuals, 2, sd)),
      Y_resid = map(lm_res, ~ {
        residuals <- .x$residuals
        if (scale_residuals) {
          residuals <- scale(residuals)
        }
        residuals
      })
    )

  if(y_as_matrix) {
    Y_resid_matrix = data_list %>%
                     select(Y_resid) %>%
                     unnest(Y_resid) %>%
                     t() %>%
                     as.matrix()
    colnames(Y_resid_matrix) <- conditions
    data_list$Y_resid <- Y_resid_matrix
  } else {
    Y_resid_list <- map(data_list$Y_resid, t) # Transpose back
    names(Y_resid_list) <- conditions
    data_list$Y_resid <- Y_resid_list
  }
  return(data_list)
}


#' @importFrom plink2R read_plink
#' @importFrom readr read_delim
#' @import purrr dplyr tibble
#' @importFrom utils read.table
#' @importFrom tidyr unnest
#' @importFrom stringr str_split
load_regional_association_data <- function(genotype, # PLINK file
                                           phenotype, # a vector of phenotype file names 
                                           covariate, # a vector of covariate file names corresponding to the phenotype file vector
                                           region, # a string of chr:start-end
                                           conditions, # a vector of strings
                                           maf_cutoff = 0,
                                           mac_cutoff = 0,
                                           xvar_cutoff = 0,
                                           imiss_cutoff = 0,
                                           y_as_matrix = FALSE,
                                           keep_indel = TRUE,
                                           keep_samples = c(),
                                           scale_residuals = FALSE) {
    ## Load genotype
    geno <- load_genotype_data(genotype, keep_indel, keep_samples)
    ## Load phenotype and covariates and perform some pre-processing
    ### including Y ( cov ) and specific X and covar match, filter X variants based on the overlapped samples.
    data_list <- prepare_data_list(geno, phenotype, covariate, region, imiss_cutoff,
                                    maf_cutoff, mac_cutoff, xvar_cutoff)
    maf_list <- lapply(data_list$X, function(x) apply(x, 2, compute_maf))
    ## Get residue Y for each of condition and its mean and sd
    data_list <- add_Y_residuals(data_list, conditions, y_as_matrix, scale_residuals)
    ## Get residue X for each of condition and its mean and sd
    data_list <- add_X_residuals(data_list, scale_residuals)
    # Get X matrix for union of samples
    X <- prepare_X_matrix(geno, data_list, imiss_cutoff, maf_cutoff, mac_cutoff, xvar_cutoff)
    region <- unlist(strsplit(region, ":", fixed = TRUE))
    ## residual_Y: if y_as_matrix is true, then return a matrix of R conditions, with column names being the names of the conditions (phenotypes) and row names being sample names. Even for one condition it has to be a matrix with just one column. if y_as_matrix is false, then return a list of y either vector or matrix (CpG for example), and they need to match with residual_X in terms of which samples are missing.
    ## residual_X: is a list of R conditions each is a matrix, with list names being the names of conditions, column names being SNP names and row names being sample names.
    ## X: is the somewhat original genotype matrix output from `filter_X`, with column names being SNP names and row names being sample names. Sample names of X should match example sample names of residual_Y matrix form (not list); but the matrices inside residual_X would be subsets of sample name of residual_Y matrix form (not list).
    return (list(
      residual_Y = data_list$Y_resid,
      residual_X = data_list$X_resid,
      residual_Y_scalar = if(scale_residuals) data_list$Y_resid_sd else 1,
      residual_X_scalar = if(scale_residuals) data_list$X_resid_sd else 1,
      dropped_sample = list(X=data_list$dropped_samples_X,Y=data_list$dropped_samples_Y,covar=data_list$dropped_samples_covar),
      covar = data_list$covar,
      Y = data_list$Y,
      X_data = data_list$X,
      X = X,
      maf = maf_list,
      chrom = region[1],
      grange = unlist(strsplit(region[2], "-", fixed = TRUE))
    ))
}

#' @return A list
#' @export
load_regional_univariate_data <- function(...) {
  dat <- load_regional_association_data(y_as_matrix = FALSE, ...)
  return (list(
          residual_Y = dat$residual_Y,
          residual_X = dat$residual_X,
          residual_Y_scalar = dat$residual_Y_scalar,
          residual_X_scalar = dat$residual_X_scalar,
          X = dat$X,
          dropped_sample = dat$dropped_sample,
          maf = dat$maf,
          chrom = dat$chrom,
          grange = dat$grange
          ))
}

#' @return A list
#' @export
load_regional_regression_data <- function(...) {
  dat <- load_regional_association_data(...)
  return (list(
          Y = dat$Y,
          X_data = dat$X_data,
          covar = dat$covar,
          dropped_sample = dat$dropped_sample,
          maf = dat$maf,
          chrom = dat$chrom,
          grange = dat$grange
          ))
}

#' @return A list
#' @export
load_regional_multivariate_data <- function(matrix_y_min_complete = NULL, # when Y is saved as matrix, remove those with non-missing counts less than this cutoff
                                            ...) {
  dat = load_regional_association_data(y_as_matrix = TRUE, ...)
  if (!is.null(matrix_y_min_complete)) {
    Y = filter_Y(dat$residual_Y, matrix_y_min_complete)
    if (length(Y$rm_rows)>0) {
      X =  dat$X[-Y$rm_rows, ]
      Y_scalar = dat$residual_Y_scalar[-Y$rm_rows]
      dropped_sample = rownames(dat$residual_Y)[Y$rm_rows]
    }else{
     X = dat$X
     Y_scalar = dat$residual_Y_scalar
     dropped_sample = dat$dropped_sample
    }   
  } else {
    Y = dat$residual_Y
    X = dat$X
    Y_scalar = dat$residual_Y_scalar
    dropped_sample = dat$dropped_sample
  }
  return (list(
        residual_Y = Y,
        residual_Y_scalar = Y_scalar,
        dropped_sample = dropped_sample,
        X = X,
        maf = dat$maf,
        chrom = dat$chrom,
        grange = dat$grange
        ))
}