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
  rm_rows <- NULL
  if(is.matrix(Y)){
    rm_rows <- which(apply(Y, 1, compute_all_missing_y))
    if (length(rm_rows)) Y <- Y[-rm_rows, ]  
  } else {
    Y <- Y[which(!is.na(Y))]
  }
  return(list(Y=Y, rm_rows = rm_rows))
}

format_variant_id <- function(names_vector) {
    gsub("_", ":", names_vector)
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

parse_region <- function(region) {
  region_split <- unlist(strsplit(region, ":", fixed = TRUE))
  chrom <- gsub("^chr", "", region_split[1])  # Remove 'chr' prefix if present
  grange <- as.numeric(unlist(strsplit(region_split[2], "-", fixed = TRUE)))
  return(list(chrom = chrom, grange = grange))
}

NoSNPsError <- function(message) {
  structure(list(message = message), class = c("NoSNPsError", "error", "condition"))
}

#' Load genotype data for a specific region using data.table for efficiency
#' 
#' @param genotype Path to the genotype data file (without extension).
#' @param region The target region in the format "chr:start-end".
#' @param keep_indel Whether to keep indel SNPs.
#' @return A vector of SNP IDs in the specified region.
#' @importFrom snpStats read.plink
#' @importFrom data.table fread
#' @export
load_genotype_region <- function(genotype, region = NULL, keep_indel = TRUE) {
  if (!is.null(region)) {
    # Get SNP IDs from bim file
    parsed_region <- parse_region(region)
    chrom <- parsed_region$chrom
    start <- parsed_region$grange[1]
    end <- parsed_region$grange[2]
    # 6 columns for bim file
    col_types <- c("character", "character", "NULL", "integer", "NULL", "NULL")
    # Read a few lines of the bim file to check for 'chr' prefix
    bim_sample <- fread(paste0(genotype, ".bim"), nrows = 5, header = FALSE, colClasses = col_types)
    chr_prefix_present <- any(grepl("^chr", bim_sample$V1))
    # Read the bim file and remove 'chr' prefix if present
    bim_data <- fread(paste0(genotype, ".bim"), header = FALSE, colClasses = col_types)
    if (chr_prefix_present) {
      bim_data[, V1 := gsub("^chr", "", V1)]
    }
    snp_ids <- bim_data[V1 == chrom & start <= V4 & V4 <= end, V2]
    if (length(snp_ids) == 0) {
      stop(NoSNPsError(paste("No SNPs found in the specified region", region)))
    }
  } else {
    snp_ids <- NULL
  }
  # Read genotype data using snpStats read.plink
  geno <- read.plink(genotype, select.snps = snp_ids)

  # Remove indels if specified
  if (!keep_indel) {
    is_indel <- with(geno$map, grepl("[^ATCG]", allele.1) | grepl("[^ATCG]", allele.2) | nchar(allele.1) > 1 | nchar(allele.2) > 1)
    geno_bed <- geno$genotypes[, !is_indel]
  } else {
    geno_bed <- geno$genotypes
  }
  # By default, plink usage dosage of the *major* allele, since allele A1 is
  # usually the minor allele and the code "1" refers to the second allele A2,
  # so that "11" is A2/A2 or major/major.
 
  # We always use minor allele dosage, to be consistent with the output from
  # plink --recodeA which used minor allele dosage by default.
  return(2 - as(geno_bed, "numeric"))
}

load_covariate_data <- function(covariate_path) {
  return(map(covariate_path, ~ read_delim(.x, "\t", col_types = cols()) %>% select(-1) %>% t()))
}

load_phenotype_data <- function(phenotype_path, region) {
  return(map(phenotype_path, ~ tabix_region(.x, region) %>% select(-4) %>% t() %>% as.matrix()))
}

filter_by_common_samples <- function(dat, common_samples) {
  dat[common_samples, , drop = FALSE] %>% .[order(rownames(.)), ]
}

#' @importFrom readr read_delim cols
prepare_data_list <- function(geno_bed, phenotype, covariate, imiss_cutoff, maf_cutoff, mac_cutoff, xvar_cutoff, keep_samples = NULL) {
  data_list <- tibble(
    covar = covariate,
    Y = phenotype) %>%
    mutate(
      # Determine common complete samples across Y, covar, and geno_bed, considering missing values
      common_complete_samples = map2(covar, Y, ~ {
        covar_non_na <- rownames(.x)[!apply(.x, 1, function(row) all(is.na(row)))]
        y_non_na <- rownames(.y)[!apply(.y, 1, function(row) all(is.na(row)))]
        if (length(intersect(intersect(covar_non_na, y_non_na), rownames(geno_bed))) == 0) {
          stop("No common complete samples between genotype and phenotype/covariate data")
        } 
        intersect(intersect(covar_non_na, y_non_na), rownames(geno_bed))
      }),
      # Further intersect with keep_samples if provided
      common_complete_samples = if (!is.null(keep_samples) && length(keep_samples) > 0) {
        map(common_complete_samples, ~ intersect(.x, keep_samples))
      } else {
        common_complete_samples
      },
      # Determine dropped samples before filtering
      dropped_samples_covar = map2(covar, common_complete_samples, ~ setdiff(rownames(.x), .y)),
      dropped_samples_Y = map2(Y, common_complete_samples, ~ setdiff(rownames(.x), .y)),
      dropped_samples_X = map(common_complete_samples, ~ setdiff(rownames(geno_bed), .x)),
      # Filter data based on common complete samples
      Y = map2(Y, common_complete_samples, ~ filter_by_common_samples(.x, .y)),
      covar = map2(covar, common_complete_samples, ~ filter_by_common_samples(.x, .y)),
      # Apply filter_X on the geno_bed data filtered by common complete samples
      X = map(common_complete_samples, ~ {
        filtered_geno_bed <- filter_by_common_samples(geno_bed, .x)
        mac_val <- if (nrow(filtered_geno_bed) == 0) 0 else (mac_cutoff / (2 * nrow(filtered_geno_bed)))
        maf_val <- max(maf_cutoff, mac_val)
        filter_X(filtered_geno_bed, imiss_cutoff, maf_val, var_thresh = xvar_cutoff)
      })
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
  message(paste0("Dimension of input genotype data is row: ", nrow(X_filtered), " column: ", ncol(X_filtered) ))
  return(X_filtered)
}

add_X_residuals <- function(data_list, scale_residuals = FALSE) {
  # Compute residuals for X and add them to data_list
  data_list <- data_list %>%
    mutate(
      lm_res_X = map2(X, covar, ~ .lm.fit(x = cbind(1, .y), y = .x)$residuals %>% as.matrix()),
      X_resid_mean = map(lm_res_X, ~ apply(.x, 2, mean)),
      X_resid_sd = map(lm_res_X, ~ apply(.x, 2, sd)),
      X_resid = map(lm_res_X, ~ {
        if (scale_residuals) {
          scale(.x)
        } else {
          .x
        }
      })
    )

  return(data_list)
}

add_Y_residuals <- function(data_list, conditions, y_as_matrix = FALSE, scale_residuals = FALSE) {
  # Compute residuals, their mean, and standard deviation, and add them to data_list
  data_list <- data_list %>%
    mutate(
      lm_res = map2(Y, covar, ~ .lm.fit(x = cbind(1, .y), y = .x)$residuals %>% as.matrix()),
      Y_resid_mean = map(lm_res, ~ apply(.x, 2, mean)),
      Y_resid_sd = map(lm_res, ~ apply(.x, 2, sd)),
      Y_resid = map(lm_res, ~ {
        if (scale_residuals) {
          scale(.x)
        } else {
          .x
        }
      })
    )

  if(y_as_matrix) {
    # FIXME: double check the logic here see if NA is padded into it when there are missing data input
    Y_resid_matrix = data_list %>%
                     select(Y_resid) %>%
                     unnest(Y_resid) %>%
                     as.matrix()
    colnames(Y_resid_matrix) <- conditions
    data_list$Y_resid <- Y_resid_matrix
  } else {
    names(data_list$Y_resid) <- conditions
  }
  return(data_list)
}

#' @importFrom plink2R read_plink
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
                                           keep_samples = NULL,
                                           scale_residuals = FALSE) {
    ## Load genotype
    geno <- load_genotype_region(genotype, region, keep_indel)
    ## Load phenotype and covariates and perform some pre-processing
    covar <- load_covariate_data(covariate)
    pheno <- load_phenotype_data(phenotype, region)
    ### including Y ( cov ) and specific X and covar match, filter X variants based on the overlapped samples.
    data_list <- prepare_data_list(geno, pheno, covar, imiss_cutoff,
                                    maf_cutoff, mac_cutoff, xvar_cutoff, keep_samples)
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
      residual_Y_scalar = if(scale_residuals) data_list$Y_resid_sd else rep(1, length(data_list$Y_resid)),
      residual_X_scalar = if(scale_residuals) data_list$X_resid_sd else rep(1, length(data_list$X_resid)),
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

#' Load, Validate, and Consolidate TWAS Weights from Multiple RDS Files
#'
#' This function loads TWAS weight data from multiple RDS files, checks for the presence
#' of specified region and condition. If variable_name_obj is provided, it aligns and
#' consolidates weight matrices based on the object's variant names, filling missing data
#' with zeros. If variable_name_obj is NULL, it checks that all files have the same row
#' numbers for the condition and consolidates weights accordingly.
#'
#' @param weight_db_files Vector of file paths for RDS files containing TWAS weights. 
#' Each element organized as region/condition/weights
#' @param condition The specific condition to be checked and consolidated across all files.
#' @param region The specific region to be checked and consolidated across all files.
#' @param variable_name_obj The name of the variable/object to fetch from each file, if not NULL.
#' @return A consolidated matrix of weights for the specified condition.
#' @examples
#' # Example usage (replace with actual file paths, condition, region, and variable_name_obj):
#' weight_db_files <- c("path/to/file1.rds", "path/to/file2.rds")
#' condition <- "example_condition"
#' region <- "example_region"
#' variable_name_obj <- "example_variable" # or NULL for standard processing
#' consolidated_weights <- load_twas_weights(weight_db_files, condition, region, variable_name_obj)
#' print(consolidated_weights)
#' @import dplyr
#' @export
load_twas_weights <- function(weight_db_files, condition, region, 
                              variable_name_obj = NULL,
                              twas_weights_table = "twas_weights") {
  ## Internal function to load and validate data from RDS files
  load_and_validate_data <- function(weight_db_files, condition, region) {
    all_data <- lapply(weight_db_files, readRDS)

    ## Check if the specified region and condition are available in all files
    for (data in all_data) {
      if (!(region %in% names(data)) || !(condition %in% names(data[[region]]))) {
        stop("The specified region and/or condition is not available in all RDS files.")
      }
      if (!is.null(variable_name_obj) && !(variable_name_obj %in% names(data[[region]][[condition]]))) {
        stop("The specified variable_name_obj is not available in all RDS files.")
      }
    }

    return(all_data)
  }

  ## Internal function to align and merge weight matrices
  align_and_merge <- function(weights_list, variable_objs) {
    # Get the complete list of variant names across all files
    all_variants <- unique(unlist(lapply(variable_objs, names)))
    
    # Initialize the consolidated matrix with zeros
    consolidated_matrix <- matrix(0, nrow = length(all_variants), ncol = 0)
    existing_colnames <- character(0)

    # Fill the matrix with weights, aligning by variant names
    for (i in seq_along(weights_list)) {
      temp_matrix <- matrix(0, nrow = length(all_variants), ncol = ncol(weights_list[[i]]))
      rownames(temp_matrix) <- all_variants
      idx <- match(names(variable_objs[[i]]), all_variants)
      temp_matrix[idx, ] <- weights_list[[i]]

      # Ensure no duplicate column names
      new_colnames <- colnames(weights_list[[i]])
      if (any(duplicated(c(existing_colnames, new_colnames)))) {
        stop("Duplicate column names detected during merging process.")
      }
      existing_colnames <- c(existing_colnames, new_colnames)

      consolidated_matrix <- cbind(consolidated_matrix, temp_matrix)
    }

    colnames(consolidated_matrix) <- existing_colnames
    return(consolidated_matrix)
  }

  ## Internal function to consolidate weights for a given condition and region
  consolidate_weights <- function(all_data, condition, region, variable_name_obj) {
    weights <- lapply(all_data, function(data) sapply(data[[region]][[condition]][[twas_weights_table]], cbind))
    if (is.null(variable_name_obj)) {
      # Standard processing: Check for identical row numbers and consolidate
      row_numbers <- sapply(weights, function(data) nrow(data))
      if (length(unique(row_numbers)) > 1) {
        stop("Not all files have the same number of rows for the specified condition.")
      }
      weights <- sapply(weights, cbind)
    } else {
      # Processing with variable_name_obj: Align and merge data, fill missing with zeros
      variable_objs <- lapply(all_data, function(data) data[[region]][[condition]][[variable_name_obj]])
      weights <- align_and_merge(weights, variable_objs)
    }
    return(weights)
  }

  ## Load, validate, and consolidate data
  try({
    all_data <- load_and_validate_data(weight_db_files, condition, region)
    consolidated_weights <- consolidate_weights(all_data, condition, region, variable_name_obj)
    return(consolidated_weights)
  }, silent = TRUE)
}