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

#' @import dplyr
extract_tensorqtl_data <- function(path, region) {
        tabix_region(path, region) %>%
            mutate(variant = paste(`#CHROM`, POS, REF, ALT, sep = ":")) %>%
            select(-c(3, 6:9)) %>%
            distinct(variant, .keep_all = TRUE) %>%
            as.matrix
}

# This function extracts tensorQTL results for given region for multiple summary statistics files
#' @import dplyr
#' @export 
load_multitrait_tensorqtl_sumstat <- function(sumstats_paths, region, trait_names) {
        extract_component <- function(df, component_index) {
            df %>%
            select(6:ncol(df)) %>%
            mutate(across(everything(), ~as.numeric(strsplit(as.character(.), ":")[[1]][component_index]))) %>%
            as.matrix
        }
        Y <- lapply(sumstats_paths, extract_tensorqtl_data, region)

        # Combine matrices
        combined_matrix <- Reduce(function(x, y) merge(x, y, by = c("variant", "#CHROM", "POS", "REF", "ALT")), Y) %>%
            distinct(variant, .keep_all = TRUE)

        dat <- list(
            bhat = extract_component(combined_matrix, 1),
            sbhat = extract_component(combined_matrix, 2)
        )

        rownames(dat$bhat) <- rownames(dat$sbhat) <- combined_matrix$variant
        colnames(dat$bhat) <- colnames(dat$sbhat) <- trait_names
        return(dat)
}

load_genotype_data <- function(genotype, keep_indel = TRUE) {
  # Only geno_bed is used in the functions
  geno = read_plink(genotype)
  rownames(geno$bed) = read.table(text = rownames(geno$bed), sep= ":")$V2
  ## if indel is false, remove the indel in the genotype
  if (!keep_indel) {
    geno_bim = geno$bim %>%
      rename("chrom" = "V1","variant_id" = "V2","alt" = "V5","ref"="V6") %>%
      mutate(indel = ifelse(grepl("[^ATCG]",alt)=="TRUE"|grepl("[^ATCG]",ref)=="TRUE"|nchar(alt)>1|nchar(ref)>1,1, 0))
    geno_bed = geno$bed[,geno_bim$indel==0]
  } else {
    geno_bed = geno$bed
  }
  return(geno_bed)
}

prepare_data_list <- function(phenotype, covariate, geno_bed, imiss_cutoff, maf_cutoff, mac_cutoff, xvar_cutoff) {
  ## Load phenotype and covariates and perform some pre-processing
  ### including Y ( cov ) and specific X and covar match, filter X variants based on the overlapped samples.
  data_list = tibble(covariate_path = covariate, phenotype_path =phenotype) %>%
        mutate(
          covar = map(covariate_path, ~read_delim(.x,"\t")%>%select(-1)%>%na.omit%>%t()),
          Y = map2(phenotype_path,covar, ~{
            y_data <- tabix_region(.x, region)%>%select(-4)%>%select(rownames(.y))%>%t()%>%as.matrix
            return(y_data)
          }),
          Y = map(Y, ~.x%>%na.omit),    # remove na where Y raw data has na which block regression
          dropped_sample = map2(covar, Y , ~rownames(.x)[!rownames(.x) %in% rownames(.y)]),
          covar = map2(covar, Y , ~.x[intersect(.x%>%rownames,rownames(.y)),]), # remove the dropped samples from Y
          X_data = map(covar,~ filter_X( geno_bed[intersect(rownames(.x),rownames(geno_bed)),], imiss_cutoff, max(maf_cutoff, mac_cutoff/(2*length(intersect(rownames(.x),rownames(geno_bed))) ) ), xvar_cutoff)   ))
  return(data_list)
}

compute_residuals <- function(data_list) {
  ## Get residue Y for each of condition and its mean and sd
  data_list = data_list %>%
    mutate(
      Y_resid_mean = map2(Y, covar, ~.lm.fit(x = cbind(1,.y), y = .x)$residuals %>% apply(.,2,mean)),
      Y_resid_sd = map2(Y, covar, ~.lm.fit(x = cbind(1,.y), y = .x)$residuals %>% apply(.,2,sd)),
      Y_resid = map2(Y, covar, ~.lm.fit(x = cbind(1,.y), y = .x)$residuals %>% scale %>% t %>% as_tibble)) ## T so that it can be unnest
  return(data_list)
}

process_Y_residuals <- function(data_list, y_as_matrix, conditions) {
  if(y_as_matrix) {
      Y_resid = data_list %>% select(Y_resid) %>% unnest(Y_resid) %>% t %>% as.matrix
      colnames(Y_resid) = conditions
      print(paste("Dimension of Y matrix:", nrow(Y_resid), ncol(Y_resid)))
  } else {
      Y_resid = map(data_list$Y_resid,~.x %>% t) # Transpose back 
      names(Y_resid) = conditions
  }
  return(Y_resid)
}

prepare_X_matrix <- function(geno_bed, data_list, imiss_cutoff, maf_cutoff, xvar_cutoff) {
  # Get X matrix for union of samples
  all_samples = map(data_list$covar, ~rownames(.x)) %>% unlist %>% unique()
  maf_cutoff = max(maf_cutoff,mac_cutoff/(2*length(all_samples)))
  X = filter_X(geno_bed[all_samples,], imiss_cutoff, maf_cutoff, xvar_cutoff) ## Filter X for mvSuSiE
  return(X)
}

prepare_X_residuals <- function(data_list, X) {
  ## Get residue X for each of condition and its mean and sd
  print(paste0("Dimension of input genotype data is row:", nrow(X), " column: ", ncol(X) ))
  X_list = data_list %>%
    mutate(
      X_resid_mean= map2(X_data,covar,~.lm.fit(x = cbind(1,.y), y = .x)$residuals %>% data.frame() %>% apply(.,2,mean)),
      X_resid_sd= map2(X_data,covar,~.lm.fit(x = cbind(1,.y), y = .x)$residuals %>% data.frame() %>% apply(.,2,sd)),
      X_resid = map2(X_data,covar,~.lm.fit(x = cbind(1,.y), y = .x)$residuals %>% scale)) %>%
    select(X_resid_mean,X_resid_sd,X_resid)
  return(X_list)
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
                                           keep_indel = TRUE) {
    ## Load genotype
    geno <- load_genotype_data(genotype, keep_indel)
    ## Load phenotype and covariates and perform some pre-processing
    ### including Y ( cov ) and specific X and covar match, filter X variants based on the overlapped samples.
    data_list <- prepare_data_list(covariate, phenotype, region, geno$bed, imiss_cutoff,
                                    maf_cutoff, mac_cutoff, xvar_cutoff)
    ## Get residue Y for each of condition and its mean and sd
    data_list <- compute_residuals(data_list)
    Y_resid <- process_Y_residuals(data_list, y_as_matrix, conditions)
    # Get X matrix for union of samples
    X <- prepare_X_matrix(geno$bed, data_list$covar, imiss_cutoff, maf_cutoff, xvar_cutoff)
    X_list <- prepare_X_residuals(data_list, X)
    maf_list = lapply(data_list$X_data, function(x) apply(x, 2, compute_maf))
    region <- unlist(strsplit(region, ":", fixed = TRUE))
    ## Get residue X for each of condition and its mean and sd
    print(paste0("Dimension of input genotype data is row:", nrow(X), " column: ", ncol(X) ))
    ## residual_Y_scaled: if y_as_matrix is true, then return a matrix of R conditions, with column names being the names of the conditions (phenotypes) and row names being sample names. Even for one condition it has to be a matrix with just one column. if y_as_matrix is false, then return a list of y either vector or matrix (CpG for example), and they need to match with residual_X_scaled in terms of which samples are missing.
    ## residual_X_scaled: is a list of R conditions each is a matrix, with list names being the names of conditions, column names being SNP names and row names being sample names.
    ## X: is the somewhat original genotype matrix output from `filter_X`, with column names being SNP names and row names being sample names. Sample names of X should match example sample names of residual_Y_scaled matrix form (not list); but the matrices inside residual_X_scaled would be subsets of sample name of residual_Y_scaled matrix form (not list).
    return (list(
      residual_Y_scaled = Y_resid,
      residual_X_scaled = X_list$X_resid,
      residual_Y_sd = data_list$Y_resid_sd,
      residual_X_sd = X_list$X_resid_sd,
      dropped_sample = data_list$dropped_sample,
      covar = data_list$covar,
      Y = data_list$Y,
      X = X,
      X_data = data_list$X_data,
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
          residual_Y_scaled = dat$residual_Y_scaled,
          residual_X_scaled = dat$residual_X_scaled,
          residual_Y_sd = dat$residual_Y_sd,
          residual_X_sd = dat$residual_X_sd,
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
    Y = filter_Y(dat$residual_Y_scaled, matrix_y_min_complete)
    if (length(Y$rm_rows)>0) {
      X =  dat$X[-Y$rm_rows, ]
      Y_sd = dat$residual_Y_sd[-Y$rm_rows]
      dropped_sample = rownames(dat$residual_Y_scaled)[Y$rm_rows]
    }else{
     X = dat$X
     Y_sd = dat$residual_Y_sd
     dropped_sample = dat$dropped_sample
    }   
  } else {
    Y = dat$residual_Y_scaled
    X = dat$X
    Y_sd = dat$residual_Y_sd
    dropped_sample = dat$dropped_sample
  }
  return (list(
        residual_Y_scaled = Y,
        residual_Y_sd = Y_sd,
        dropped_sample = dropped_sample,
        X = X,
        maf = dat$maf,
        chrom = dat$chrom,
        grange = dat$grange
        ))
}