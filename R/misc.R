matxMax <- function(mtx) {
  return(arrayInd(which.max(mtx), dim(mtx)))
}

compute_maf <- function(geno) {
  f <- mean(geno, na.rm = TRUE) / 2
  return(min(f, 1 - f))
}

compute_missing <- function(geno) {
  miss <- sum(is.na(geno)) / length(geno)
  return(miss)
}

compute_non_missing_y <- function(y) {
  nonmiss <- sum(!is.na(y))
  return(nonmiss)
}

compute_all_missing_y <- function(y) {
  allmiss <- all(is.na(y))
  return(allmiss)
}

mean_impute <- function(geno) {
  f <- apply(geno, 2, function(x) mean(x, na.rm = TRUE))
  for (i in 1:length(f)) geno[, i][which(is.na(geno[, i]))] <- f[i]
  return(geno)
}

is_zero_variance <- function(x) {
  if (length(unique(x)) == 1) {
    return(T)
  } else {
    return(F)
  }
}

#' @importFrom matrixStats colVars
filter_X <- function(X, missing_rate_thresh, maf_thresh, var_thresh = 0) {
  rm_col <- which(apply(X, 2, compute_missing) > missing_rate_thresh)
  if (length(rm_col)) X <- X[, -rm_col]
  rm_col <- which(apply(X, 2, compute_maf) <= maf_thresh)
  if (length(rm_col)) X <- X[, -rm_col]
  rm_col <- which(apply(X, 2, is_zero_variance))
  if (length(rm_col)) X <- X[, -rm_col]
  X <- mean_impute(X)
  if (var_thresh > 0) {
    rm_col <- which(matrixStats::colVars(X) < var_thresh)
    if (length(rm_col)) X <- X[, -rm_col]
  }
  return(X)
}

filter_Y <- function(Y, n_nonmiss) {
  rm_col <- which(apply(Y, 2, compute_non_missing_y) < n_nonmiss)
  if (length(rm_col)) Y <- Y[, -rm_col]
  rm_rows <- NULL
  if (is.matrix(Y)) {
    rm_rows <- which(apply(Y, 1, compute_all_missing_y))
    if (length(rm_rows)) Y <- Y[-rm_rows, ]
  } else {
    Y <- Y[which(!is.na(Y))]
  }
  return(list(Y = Y, rm_rows = rm_rows))
}

format_variant_id <- function(names_vector) {
  gsub("_", ":", names_vector)
}

#' Converted  Variant ID into a properly structured data frame
#' @param variant_id A data frame or character vector representing variant IDs.
#'   Expected formats are a data frame with columns "chrom", "pos", "A1", "A2",
#'   or a character vector in "chr:pos:A2:A1" or "chr:pos_A2_A1" format.
#' @return A data frame with columns "chrom", "pos", "A1", "A2", where 'chrom'
#'   and 'pos' are integers, and 'A1' and 'A2' are allele identifiers.
#' @noRd
variant_id_to_df <- function(variant_id) {
  # Check if target_variants is already a data.frame with the required columns
  if (is.data.frame(variant_id)) {
    if (!all(c("chrom", "pos", "A1", "A2") %in% names(variant_id))) {
      names(variant_id) <- c("chrom", "pos", "A2", "A1")
    }
    # Ensure that 'chrom' values are integers
    variant_id$chrom <- ifelse(grepl("^chr", variant_id$chrom),
      as.integer(sub("^chr", "", variant_id$chrom)), # Remove 'chr' and convert to integer
      as.integer(variant_id$chrom)
    ) # Convert to integer if not already
    variant_id$pos <- as.integer(variant_id$pos)
    return(variant_id)
  }
  # Function to split a string and create a data.frame
  create_dataframe <- function(string) {
    string <- gsub("_", ":", string)
    parts <- strsplit(string, ":", fixed = TRUE)
    data <- data.frame(do.call(rbind, parts), stringsAsFactors = FALSE)
    colnames(data) <- c("chrom", "pos", "A2", "A1")
    # Ensure that 'chrom' values are integers
    data$chrom <- ifelse(grepl("^chr", data$chrom),
      as.integer(sub("^chr", "", data$chrom)), # Remove 'chr' and convert to integer
      as.integer(data$chrom)
    ) # Convert to integer if not already
    data$pos <- as.integer(data$pos)
    return(data)
  }
  return(create_dataframe(variant_id))
}

load_genotype_data <- function(genotype, keep_indel = TRUE) {
  # Read genotype data using plink
  geno <- plink2R::read_plink(genotype)
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

#' @importFrom stringr str_split
#' @export
parse_region <- function(region) {
  if (!is.character(region) || length(region) != 1) {
    return(region)
  }

  if (!grepl("^chr[0-9XY]+:[0-9]+-[0-9]+$", region)) {
    stop("Input string format must be 'chr:start-end'.")
  }
  parts <- str_split(region, "[:-]")[[1]]
  df <- data.frame(
    chrom = gsub("^chr", "", parts[1]),
    start = as.integer(parts[2]),
    end = as.integer(parts[3])
  )

  return(df)
}

# Retrieve a nested element from a list structure
#' @export
get_nested_element <- function(nested_list, name_vector) {
  if (is.null(name_vector)) {
    return(NULL)
  }
  current_element <- nested_list
  for (name in name_vector) {
    if (is.null(current_element[[name]])) {
      stop("Element not found in the list")
    }
    current_element <- current_element[[name]]
  }
  return(current_element)
}

NoSNPsError <- function(message) {
  structure(list(message = message), class = c("NoSNPsError", "error", "condition"))
}

#' Load genotype data for a specific region using data.table for efficiency
#'
#' By default, plink usage dosage of the *major* allele, since "effect allele" A1 is
#' usually the minor allele and the code "1" refers to the "other allele" A2,
#' so that "11" is A2/A2 or major/major. We always use effect allele dosage, to
#' be more consistent with the minor allele based convention ie, plink --recodeA which used minor allele
#' dosage by default.
#'
#' @param genotype Path to the genotype data file (without extension).
#' @param region The target region in the format "chr:start-end".
#' @param keep_indel Whether to keep indel SNPs.
#' @return A vector of SNP IDs in the specified region.
#' 
#' @importFrom data.table fread
#' @importFrom magrittr %>%
#' @importFrom snpStats read.plink
#' @export
load_genotype_region <- function(genotype, region = NULL, keep_indel = TRUE) {
  if (!is.null(region)) {
    # Get SNP IDs from bim file
    parsed_region <- parse_region(region)
    chrom <- parsed_region$chrom
    start <- parsed_region$start
    end <- parsed_region$end
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
  return(2 - as(geno_bed, "numeric"))
}

#' @importFrom purrr map
#' @importFrom readr read_delim cols
#' @importFrom dplyr select mutate across everything
#' @importFrom magrittr %>%
#' @noRd
load_covariate_data <- function(covariate_path) {
  return(map(covariate_path, ~ read_delim(.x, "\t", col_types = cols()) %>%
    select(-1) %>%
    mutate(across(everything(), as.numeric)) %>%
    t()))
}

NoPhenotypeError <- function(message) {
  structure(list(message = message), class = c("NoPhenotypeError", "error", "condition"))
}

#' @importFrom purrr map2 compact
#' @importFrom readr read_delim cols
#' @importFrom dplyr filter select mutate across everything
#' @importFrom magrittr %>%
#' @noRd
load_phenotype_data <- function(phenotype_path, region, extract_region_name = NULL, region_name_col = NULL, tabix_header = TRUE) {
  if (is.null(extract_region_name)) {
    extract_region_name <- rep(list(NULL), length(phenotype_path))
  } else if (is.list(extract_region_name) && length(extract_region_name) != length(phenotype_path)) {
    stop("extract_region_name must be NULL or a list with the same length as phenotype_path.")
  } else if (!is.null(extract_region_name) && !is.list(extract_region_name)) {
    stop("extract_region_name must be NULL or a list.")
  }

  # Use `map2` to iterate over `phenotype_path` and `extract_region_name` simultaneously
  # `compact` should remove all NULL element
  phenotype_data <- compact(map2(phenotype_path, extract_region_name, ~ {
    tabix_data <- if (!is.null(region)) tabix_region(.x, region, tabix_header = tabix_header) else read_delim(.x, "\t", col_types = cols())
    if (nrow(tabix_data) == 0) {
      message(paste("Phenotype file ", .x, " is empty for the specified region", if (!is.null(region)) "" else region))
      return(NULL)
    }
    if (!is.null(.y) && is.vector(.y) && !is.null(region_name_col) && (region_name_col %% 1 == 0)) {
      if (region_name_col <= ncol(tabix_data)) {
        region_col_name <- colnames(tabix_data)[region_name_col]
        tabix_data <- tabix_data %>%
          filter(.data[[region_col_name]] %in% .y) %>%
          t()
        colnames(tabix_data) <- tabix_data[region_name_col, ]
        return(tabix_data)
      } else {
        stop("region_name_col is out of bounds for the number of columns in tabix_data.")
      }
    } else {
      return(tabix_data %>% t())
    }
  }))

  # Check if all phenotype files are empty
  if (length(phenotype_data) == 0) {
    stop(NoPhenotypeError(paste("All phenotype files are empty for the specified region", if (!is.null(region)) "" else region)))
  }
  return(phenotype_data)
}

#' @importFrom purrr map
#' @importFrom tibble as_tibble
#' @importFrom dplyr mutate
#' @importFrom magrittr %>%
#' @noRd
extract_phenotype_coordinates <- function(phenotype_list) {
  return(map(phenotype_list, ~ t(.x[1:3, ]) %>%
    as_tibble() %>%
    mutate(start = as.numeric(start), end = as.numeric(end))))
}

#' @importFrom magrittr %>%
#' @noRd
filter_by_common_samples <- function(dat, common_samples) {
  dat[common_samples, , drop = FALSE] %>% .[order(rownames(.)), ]
}

#' @importFrom tibble tibble
#' @importFrom dplyr mutate select
#' @importFrom purrr map map2
#' @importFrom magrittr %>%
#' @noRd
prepare_data_list <- function(geno_bed, phenotype, covariate, imiss_cutoff, maf_cutoff, mac_cutoff, xvar_cutoff, phenotype_header = 4, keep_samples = NULL) {
  data_list <- tibble(
    covar = covariate,
    Y = lapply(phenotype, function(x) apply(x[-c(1:phenotype_header), , drop = F], c(1, 2), as.numeric))
  ) %>%
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
      # Apply filter_X on the geno_bed data filtered by common complete samples and then format column names
      X = map(common_complete_samples, ~ {
        filtered_geno_bed <- filter_by_common_samples(geno_bed, .x)
        mac_val <- if (nrow(filtered_geno_bed) == 0) 0 else (mac_cutoff / (2 * nrow(filtered_geno_bed)))
        maf_val <- max(maf_cutoff, mac_val)
        filtered_data <- filter_X(filtered_geno_bed, imiss_cutoff, maf_val, var_thresh = xvar_cutoff)
        colnames(filtered_data) <- format_variant_id(colnames(filtered_data)) # Format column names right after filtering
        filtered_data
      })
    ) %>%
    select(covar, Y, X, dropped_samples_Y, dropped_samples_X, dropped_samples_covar)
  return(data_list)
}

#' @importFrom purrr map
#' @importFrom dplyr intersect
#' @importFrom magrittr %>%
#' @noRd
prepare_X_matrix <- function(geno_bed, data_list, imiss_cutoff, maf_cutoff, mac_cutoff, xvar_cutoff) {
  # Calculate the union of all samples from data_list: any of X, covar and Y would do
  all_samples_union <- map(data_list$covar, ~ rownames(.x)) %>%
    unlist() %>%
    unique()
  # Find the intersection of these samples with the samples in geno_bed
  common_samples <- intersect(all_samples_union, rownames(geno_bed))
  # Filter geno_bed using common_samples
  X_filtered <- filter_by_common_samples(geno_bed, common_samples)
  # Calculate MAF cutoff considering the number of common samples
  maf_val <- max(maf_cutoff, mac_cutoff / (2 * length(common_samples)))
  # Apply further filtering on X
  X_filtered <- filter_X(X_filtered, imiss_cutoff, maf_val, xvar_cutoff)
  colnames(X_filtered) <- format_variant_id(colnames(X_filtered))

  # To keep a log message
  variants <- as.data.frame(do.call(rbind, lapply(format_variant_id(colnames(X_filtered)), function(x) strsplit(x, ":")[[1]][1:2])), stringsAsFactors = FALSE)
  message(paste0("Dimension of input genotype data is ", nrow(X_filtered), " rows and ", ncol(X_filtered), " columns for genomic region of ", variants[1, 1], ":", min(as.integer(variants[, 2])), "-", max(as.integer(variants[, 2]))))
  return(X_filtered)
}

#' @importFrom purrr map map2
#' @importFrom dplyr mutate
#' @importFrom stats lm.fit sd
#' @importFrom magrittr %>%
#' @noRd
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

#' @importFrom purrr map map2
#' @importFrom dplyr mutate
#' @importFrom stats lm.fit sd
#' @importFrom magrittr %>%
#' @noRd
add_Y_residuals <- function(data_list, conditions, scale_residuals = FALSE) {
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

  names(data_list$Y_resid) <- conditions

  return(data_list)
}


#' Load regional association data
#'
#' This function loads genotype, phenotype, and covariate data for a specific region and performs data preprocessing.
#'
#' @param genotype PLINK bed file containing genotype data.
#' @param phenotype A vector of phenotype file names.
#' @param covariate A vector of covariate file names corresponding to the phenotype file vector.
#' @param region A string of chr:start-end for the phenotype region.
#' @param conditions A vector of strings representing different conditions or groups.
#' @param maf_cutoff Minimum minor allele frequency (MAF) cutoff. Default is 0.
#' @param mac_cutoff Minimum minor allele count (MAC) cutoff. Default is 0.
#' @param xvar_cutoff Maximum variant missingness cutoff. Default is 0.
#' @param imiss_cutoff Maximum individual missingness cutoff. Default is 0.
#' @param association_window A string of chr:start-end for the association analysis window (cis or trans). If not provided, all genotype data will be loaded.
#' @param extract_region_name A string (e.g., gene ID ENSG00000269699) to subset the information when there are multiple regions available. Default is NULL.
#' @param region_name_col Column name containing the region name. Default is NULL.
#' @param keep_indel Logical indicating whether to keep insertions/deletions (INDELs). Default is TRUE.
#' @param keep_samples A vector of sample names to keep. Default is NULL.
#' @param phenotype_header Number of rows to skip at the beginning of the transposed phenotype file (default is 4 for chr, start, end, and ID).
#' @param scale_residuals Logical indicating whether to scale residuals. Default is FALSE.
#' @param tabix_header Logical indicating whether the tabix file has a header. Default is TRUE.
#'
#' @return A list containing the following components:
#' \itemize{
#'   \item residual_Y: A list of residualized phenotype values (either a vector or a matrix).
#'   \item residual_X: A list of residualized genotype matrices for each condition.
#'   \item residual_Y_scalar: Scaling factor for residualized phenotype values.
#'   \item residual_X_scalar: Scaling factor for residualized genotype values.
#'   \item dropped_sample: A list of dropped samples for X, Y, and covariates.
#'   \item covar: Covariate data.
#'   \item Y: Original phenotype data.
#'   \item X_data: Original genotype data.
#'   \item X: Filtered genotype matrix.
#'   \item maf: Minor allele frequency (MAF) for each variant.
#'   \item chrom: Chromosome of the region.
#'   \item grange: Genomic range of the region (start and end positions).
#'   \item Y_coordinates: Phenotype coordinates if a region is specified.
#' }
#'
#' @export
load_regional_association_data <- function(genotype, # PLINK file
                                           phenotype, # a vector of phenotype file names
                                           covariate, # a vector of covariate file names corresponding to the phenotype file vector
                                           region, # a string of chr:start-end for phenotype region
                                           conditions, # a vector of strings
                                           maf_cutoff = 0,
                                           mac_cutoff = 0,
                                           xvar_cutoff = 0,
                                           imiss_cutoff = 0,
                                           association_window = NULL, #  a string of chr:start-end for association analysis window (cis or trans). If not provided all genotype data will be loaded
                                           extract_region_name = NULL, # a string of eg gene ID ENSG00000269699, this is helpful if we only want to keep a subset of the information when there are multiple regions available
                                           region_name_col = NULL,
                                           keep_indel = TRUE,
                                           keep_samples = NULL,
                                           phenotype_header = 4, # skip first 4 rows of transposed phenotype for chr, start, end and ID
                                           scale_residuals = FALSE,
                                           tabix_header = TRUE) {
  ## Load genotype
  geno <- load_genotype_region(genotype, association_window, keep_indel)
  ## Load phenotype and covariates and perform some pre-processing
  covar <- load_covariate_data(covariate)
  pheno <- load_phenotype_data(phenotype, region, extract_region_name = extract_region_name, region_name_col = region_name_col, tabix_header = tabix_header)
  ### including Y ( cov ) and specific X and covar match, filter X variants based on the overlapped samples.
  data_list <- prepare_data_list(geno, pheno, covar, imiss_cutoff,
    maf_cutoff, mac_cutoff, xvar_cutoff,
    phenotype_header = phenotype_header, keep_samples = keep_samples
  )
  maf_list <- lapply(data_list$X, function(x) apply(x, 2, compute_maf))
  ## Get residue Y for each of condition and its mean and sd
  data_list <- add_Y_residuals(data_list, conditions, scale_residuals)
  ## Get residue X for each of condition and its mean and sd
  data_list <- add_X_residuals(data_list, scale_residuals)
  # Get X matrix for union of samples
  X <- prepare_X_matrix(geno, data_list, imiss_cutoff, maf_cutoff, mac_cutoff, xvar_cutoff)
  region <- if (!is.null(region)) unlist(strsplit(region, ":", fixed = TRUE))
  ## residual_Y: a list of y either vector or matrix (CpG for example), and they need to match with residual_X in terms of which samples are missing.
  ## residual_X: is a list of R conditions each is a matrix, with list names being the names of conditions, column names being SNP names and row names being sample names.
  ## X: is the somewhat original genotype matrix output from `filter_X`, with column names being SNP names and row names being sample names. Sample names of X should match example sample names of residual_Y matrix form (not list); but the matrices inside residual_X would be subsets of sample name of residual_Y matrix form (not list).
  return(list(
    residual_Y = data_list$Y_resid,
    residual_X = data_list$X_resid,
    residual_Y_scalar = if (scale_residuals) data_list$Y_resid_sd else rep(1, length(data_list$Y_resid)),
    residual_X_scalar = if (scale_residuals) data_list$X_resid_sd else rep(1, length(data_list$X_resid)),
    dropped_sample = list(X = data_list$dropped_samples_X, Y = data_list$dropped_samples_Y, covar = data_list$dropped_samples_covar),
    covar = data_list$covar,
    Y = data_list$Y,
    X_data = data_list$X,
    X = X,
    maf = maf_list,
    chrom = region[1],
    grange = if (!is.null(region)) unlist(strsplit(region[2], "-", fixed = TRUE)) else NULL,
    Y_coordinates = if (!is.null(region)) extract_phenotype_coordinates(pheno) else NULL
  ))
}

#' @return A list
#' @export
load_regional_univariate_data <- function(...) {
  dat <- load_regional_association_data(...)
  return(list(
    residual_Y = dat$residual_Y,
    residual_X = dat$residual_X,
    residual_Y_scalar = dat$residual_Y_scalar,
    residual_X_scalar = dat$residual_X_scalar,
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
  return(list(
    Y = dat$Y,
    X_data = dat$X_data,
    covar = dat$covar,
    dropped_sample = dat$dropped_sample,
    maf = dat$maf,
    chrom = dat$chrom,
    grange = dat$grange
  ))
}

# return matrix of R conditions, with column names being the names of the conditions (phenotypes) and row names being sample names. Even for one condition it has to be a matrix with just one column.
#' @noRd
pheno_list_to_mat <- function(data_list) {
  all_row_names <- unique(unlist(lapply(data_list$residual_Y, rownames)))
  # Step 2: Align matrices and fill with NA where necessary
  aligned_mats <- lapply(data_list$residual_Y, function(mat) {
    expanded_mat <- matrix(NA, nrow = length(all_row_names), ncol = 1, dimnames = list(all_row_names, NULL))
    common_rows <- intersect(rownames(mat), all_row_names)
    expanded_mat[common_rows, ] <- mat[common_rows, ]
    return(expanded_mat)
  })
  Y_resid_matrix <- do.call(cbind, aligned_mats)
  colnames(Y_resid_matrix) <- names(data_list$residual_Y)
  data_list$residual_Y <- Y_resid_matrix
  return(data_list)
}

#' @return A list
#' @export
load_regional_multivariate_data <- function(matrix_y_min_complete = NULL, # when Y is saved as matrix, remove those with non-missing counts less than this cutoff
                                            ...) {
  dat <- pheno_list_to_mat(load_regional_association_data(...))
  if (!is.null(matrix_y_min_complete)) {
    Y <- filter_Y(dat$residual_Y, matrix_y_min_complete)
    if (length(Y$rm_rows) > 0) {
      X <- dat$X[-Y$rm_rows, ]
      Y_scalar <- dat$residual_Y_scalar[-Y$rm_rows]
      dropped_sample <- rownames(dat$residual_Y)[Y$rm_rows]
    } else {
      X <- dat$X
      Y_scalar <- dat$residual_Y_scalar
      dropped_sample <- dat$dropped_sample
    }
  } else {
    Y <- dat$residual_Y
    X <- dat$X
    Y_scalar <- dat$residual_Y_scalar
    dropped_sample <- dat$dropped_sample
  }
  return(list(
    residual_Y = Y,
    residual_Y_scalar = Y_scalar,
    dropped_sample = dropped_sample,
    X = X,
    maf = dat$maf,
    chrom = dat$chrom,
    grange = dat$grange
  ))
}

#' @return A list
#' @export
load_regional_functional_data <- function(...) {
  dat <- load_regional_association_data(...)
  return(dat)
}
