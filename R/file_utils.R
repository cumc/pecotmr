# This needs pgenlibr package
# devtools::install_github("chrchang/plink-ng", subdir="/2.0/pgenlibr")

# read PLINK files
#' @importFrom dplyr rename
#' @importFrom data.table fread
#' @importFrom tools file_path_sans_ext
read_pvar <- function(pgen) {
  pvarf <- paste0(file_path_sans_ext(pgen), ".pvar")
  pvardt <- fread(pvarf, skip = "#CHROM")
  pvardt <- rename(pvardt,
    "chrom" = "#CHROM", "pos" = "POS",
    "alt" = "ALT", "ref" = "REF", "id" = "ID"
  )
  pvardt <- pvardt[, c("chrom", "id", "pos", "alt", "ref")]
  return(pvardt)
}

#' @importFrom data.table fread
#' @importFrom tools file_path_sans_ext
read_bim <- function(bed) {
  bimf <- paste0(file_path_sans_ext(bed), ".bim")
  bim <- fread(bimf)
  colnames(bim) <- c("chrom", "id", "gpos", "pos", "a1", "a0")
  return(bim)
}

#' @importFrom data.table fread
#' @importFrom tools file_path_sans_ext
read_psam <- function(pgen) {
  psamf <- paste0(file_path_sans_ext(pgen), ".psam")
  psam <- fread(psamf, header = T)
  colnames(psam)[1:2] <- c("FID", "IID")
  return(psam)
}

#' @importFrom data.table fread
#' @importFrom tools file_path_sans_ext
read_fam <- function(bed) {
  famf <- paste0(file_path_sans_ext(bed), ".fam")
  return(fread(famf, header = F))
}

# open pgen/pvar PLINK 2 data format
#' @importFrom pgenlibr NewPgen
open_pgen <- function(pgenf) {
  return(NewPgen(pgenf))
}

# open bed/bim/fam: A PLINK 1 .bed is a valid .pgen
#' @importFrom pgenlibr NewPgen
open_bed <- function(bed) {
  raw_s_ct <- nrow(read_fam(bed))
  return(NewPgen(bed, raw_sample_ct = raw_s_ct))
}

#' @importFrom pgenlibr GetVariantCt ReadList
read_pgen <- function(pgen, variantidx = NULL, meanimpute = F) {
  if (is.null(variantidx)) {
    variantidx <- 1:GetVariantCt(pgen)
  }

  ReadList(pgen,
    variant_subset = variantidx,
    meanimpute = meanimpute
  )
}

#' @importFrom data.table fread
#' @importFrom dplyr as_tibble mutate filter
#' @importFrom tibble tibble
#' @importFrom magrittr %>%
#' @importFrom stringr str_detect

tabix_region <- function(file, region, tabix_header = "auto", target = "", target_column_index = "") {
  cmd_output <- tryCatch(
    {
      fread(cmd = paste0("tabix -h ", file, " ", region), sep = "auto", header = tabix_header)
    },
    error = function(e) NULL
  )

  if (!is.null(cmd_output) && target != "" && target_column_index != "") {
    cmd_output <- cmd_output %>%
      filter(str_detect(.[[target_column_index]], target))
  } else if (!is.null(cmd_output) && target != "") {
    cmd_output <- cmd_output %>%
      mutate(text = apply(., 1, function(row) paste(row, collapse = "_"))) %>%
      filter(str_detect(text, target)) %>%
      select(-text)
  }

  if (is.null(cmd_output) || nrow(cmd_output) == 0) {
    return(tibble())
  }

  cmd_output %>%
    as_tibble() %>%
    mutate(
      !!names(.)[1] := as.character(.[[1]]),
      !!names(.)[2] := as.numeric(.[[2]])
    )
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
#' @param xvar_cutoff Minimum variance cutoff. Default is 0.
#' @param imiss_cutoff Maximum individual missingness cutoff. Default is 0.
#' @param association_window A string of chr:start-end for the association analysis window (cis or trans). If not provided, all genotype data will be loaded.
#' @param extract_region_name A list of vectors of strings (e.g., gene ID ENSG00000269699) to subset the information when there are multiple regions available. Default is NULL.
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
                                           association_window = NULL,
                                           extract_region_name = NULL,
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
    ### change the ncol of each matrix
    expanded_mat <- matrix(NA, nrow = length(all_row_names), ncol = ncol(mat), dimnames = list(all_row_names, colnames(mat)))
    common_rows <- intersect(rownames(mat), all_row_names)
    expanded_mat[common_rows, ] <- mat[common_rows, ]
    return(expanded_mat)
  })
  Y_resid_matrix <- do.call(cbind, aligned_mats)
  if (!is.null(names(data_list$residual_Y))) {
    colnames(Y_resid_matrix) <- names(data_list$residual_Y)
  }
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
    maf = apply(X, 2, compute_maf),
    chrom = dat$chrom,
    grange = dat$grange,
    X_variance= matrixStats::colVars(X)
  ))
}

#' @return A list
#' @export
load_regional_functional_data <- function(...) {
  dat <- load_regional_association_data(...)
  return(dat)
}



# Function to remove gene name at the end of context name
#' @export
clean_context_names <- function(context, gene) {
  context_parts <- str_split(context, "_")[[1]]
  # Remove gene name if it matches the last part of the context
  if (tail(context_parts, n = 1) == gene) context_parts <- head(context_parts, -1)
  cleaned_context <- paste(context_parts, collapse = "_")
  return(cleaned_context)
}

#' Load, Validate, and Consolidate TWAS Weights from Multiple RDS Files
#'
#' This function loads TWAS weight data from multiple RDS files, checks for the presence
#' of specified region and condition. If variable_name_obj is provided, it aligns and
#' consolidates weight matrices based on the object's variant names, filling missing data
#' with zeros. If variable_name_obj is NULL, it checks that all files have the same row
#' numbers for the condition and consolidates weights accordingly.
#'
#' @param weight_db_file weight_db_files Vector of file paths for RDS files containing TWAS weights..
#' Each element organized as region/condition/weights
#' @param condition The specific condition to be checked and consolidated across all files.
#' @param variable_name_obj The name of the variable/object to fetch from each file, if not NULL.
#' @return A consolidated list of weights for the specified condition and a list of SuSiE results.
#' @examples
#' # Example usage (replace with actual file paths, condition, region, and variable_name_obj):
#' weight_db_files <- c("path/to/file1.rds", "path/to/file2.rds")
#' condition <- "example_condition"
#' region <- "example_region"
#' variable_name_obj <- "example_variable" # or NULL for standard processing
#' consolidated_weights <- load_twas_weights(weight_db_files, condition, region, variable_name_obj)
#' print(consolidated_weights)
#' @export
load_twas_weights <- function(weight_db_files, conditions = NULL,
                              variable_name_obj = c("preset_variants_result", "variant_names"),
                              susie_obj = c("preset_variants_result", "susie_result_trimmed"),
                              twas_weights_table = "twas_weights") {
  ## Internal function to load and validate data from RDS files
  load_and_validate_data <- function(weight_db_files, conditions, variable_name_obj) {
    all_data <- do.call(c, lapply(unname(weight_db_files), function(rds_file) {
      db <- readRDS(rds_file)
      if ("mnm_result" %in% names(db)) {
        db <- list(mnm_rs = db)
      } else {
        gene <- names(db)
        # Check if region from all RDS files are the same
        if (length(gene) != 1) {
          stop("More than one region provided in the RDS file. ")
        } else {
          names(db[[gene]]) <- sapply(names(db[[gene]]), function(context) clean_context_names(context, gene))
        }
      }
      return(db)
    }))
    # Check if region from all RDS files are the same
    unique_regions <- unique(names(all_data)[!names(all_data) %in% "mnm_rs"])
    if (length(unique_regions) != 1) {
      stop("The RDS files do not refer to the same region.")
    } else {
      # Combine the lists with the same region name
      combined_all_data <- lapply(split(all_data, names(all_data)), function(lst) {
        result <- list()
        for (item in lst) {
          for (name in names(item)) {
            if (name %in% names(result)) {
              if (result[[name]] != item[[name]]) {
                stop("Different twas weight data provided for identical context name. ")
              }
            } else {
              result[[name]] <- item[[name]]
            }
          }
        }
        return(result)
      })
    }
    # merge univariate and multivariate results for same gene-context pair
    if ("mnm_rs" %in% names(combined_all_data)) {
      gene <- names(combined_all_data)[!names(combined_all_data) %in% "mnm_rs"]
      overl_contexts <- names(combined_all_data[["mnm_rs"]])[names(combined_all_data[["mnm_rs"]]) %in% names(combined_all_data[[gene]])]
      multi_variants <- get_nested_element(combined_all_data$mnm_rs$mnm_result, variable_name_obj)
      for (context in overl_contexts) {
        uni_variants <- get_nested_element(combined_all_data[[gene]][[context]], variable_name_obj)
        multi_weights <- setNames(rep(0, length(uni_variants)), uni_variants)
        multi_weights <- lapply(combined_all_data[["mnm_rs"]][[context]]$twas_weights, function(weight_list) {
          aligned_weights <- setNames(rep(0, length(uni_variants)), uni_variants)
          aligned_weights[multi_variants[multi_variants %in% uni_variants]] <- unlist(weight_list)[multi_variants %in% uni_variants]
          aligned_weights <- as.matrix(aligned_weights)
        })
        combined_all_data[[gene]][[context]]$twas_weights <- c(combined_all_data[[gene]][[context]]$twas_weights, multi_weights)
        combined_all_data[[gene]][[context]]$twas_cv_result$performance <- c(
          combined_all_data[[gene]][[context]]$twas_cv_result$performance,
          combined_all_data[["mnm_rs"]][[context]]$twas_cv_result$performance
        )
      }
      combined_all_data[["mnm_rs"]] <- NULL
    }
    combined_all_data <- combined_all_data[[1]]
    # Set default for 'conditions' if they are not specified
    if (is.null(conditions)) {
      conditions <- names(combined_all_data)
    }
    ## Check if the specified condition and variable_name_obj are available in all files
    if (!all(conditions %in% names(combined_all_data))) {
      stop("The specified condition is not available in all RDS files.")
    }
    return(combined_all_data)
  }
  # Only extract the variant_names and SuSiE results
  extract_variants_and_susie_results <- function(combined_all_data, conditions) {
    # skip conditions do not have susie results available.
    conditions <- conditions[sapply(conditions, function(cond) {
      tryCatch(
        {
          !is.null(get_nested_element(combined_all_data, c(cond, susie_obj)))
        },
        error = function(e) {
          FALSE
        }
      )
    })]
    combined_susie_result <- lapply(conditions, function(condition) {
      result <- list(
        variant_names = get_nested_element(combined_all_data, c(condition, variable_name_obj)),
        susie_result = get_nested_element(combined_all_data, c(condition, susie_obj))
      )
      if ("top_loci" %in% names(get_nested_element(combined_all_data, c(condition, "preset_variants_result")))) {
        result$top_loci <- combined_all_data[[condition]]$preset_variants_result$top_loci
      }
      if ("target" %in% names(combined_all_data[[condition]])) {
        result$target <- get_nested_element(combined_all_data, c(condition, "target"))
      }
      if ("region_info" %in% names(combined_all_data[[condition]])) {
        result$region_info <- get_nested_element(combined_all_data, c(condition, "region_info"))
      }
      return(result)
    })
    names(combined_susie_result) <- conditions
    return(combined_susie_result)
  }
  # Internal function to align and merge weight matrices
  align_and_merge <- function(weights_list, variable_objs) {
    consolidated_list <- list()
    # Fill the matrix with weights, aligning by variant names
    for (i in seq_along(weights_list)) {
      # get conditon specific variant names
      all_variants <- unique(unlist(variable_objs[[i]]))
      # Initialize the temp matrix with zeros
      existing_colnames <- character(0)
      temp_matrix <- matrix(0, nrow = length(all_variants), ncol = ncol(weights_list[[i]]))
      rownames(temp_matrix) <- all_variants
      if (!length(all_variants) == nrow(weights_list[[i]])) {
        stop(paste0(
          "Variant number mismatch in twas weights: ", nrow(weights_list[[i]]), " and variant number in susie result: ", length(all_variants),
          " for context ", names(weights_list)[i], ". "
        ))
      }
      idx <- match(variable_objs[[i]], all_variants)
      temp_matrix[idx, ] <- weights_list[[i]]
      # Ensure no duplicate column names
      new_colnames <- colnames(weights_list[[i]])
      dups <- duplicated(c(existing_colnames, new_colnames))
      if (any(dups)) {
        duplicated_names <- paste(c(existing_colnames, new_colnames)[dups], collapse = ", ")
        stop("Duplicate column names detected during merging process: ", duplicated_names, ".")
      }
      existing_colnames <- c(existing_colnames, new_colnames)

      # consolidated_list[[i]] <- matrix(as.numeric(temp_matrix), nrow = nrow(temp_matrix), byrow = TRUE)
      consolidated_list[[i]] <- temp_matrix
      colnames(consolidated_list[[i]]) <- existing_colnames
    }
    return(consolidated_list)
  }

  # Internal function to consolidate weights for given condition
  consolidate_weights_list <- function(combined_all_data, conditions, variable_name_obj, twas_weights_table) {
    combined_weights_by_condition <- lapply(conditions, function(condition) {
      temp_list <- get_nested_element(combined_all_data, c(condition, twas_weights_table))
      temp_list <- temp_list[!names(temp_list) %in% "variant_names"]
      sapply(temp_list, cbind)
    })
    names(combined_weights_by_condition) <- conditions
    if (is.null(variable_name_obj)) {
      # Standard processing: Check for identical row numbers and consolidate
      row_numbers <- sapply(combined_weights_by_condition, function(data) nrow(data))
      if (length(unique(row_numbers)) > 1) {
        stop("Not all files have the same number of rows for the specified condition.")
      }
      weights <- combined_weights_by_condition
    } else {
      # Processing with variable_name_obj: Align and merge data, fill missing with zeros
      variable_objs <- lapply(conditions, function(condition) {
        get_nested_element(combined_all_data, c(condition, variable_name_obj))
      })
      weights <- align_and_merge(combined_weights_by_condition, variable_objs)
    }
    names(weights) <- conditions
    return(weights)
  }

  ## Load, validate, and consolidate data
  try(
    {
      combined_all_data <- load_and_validate_data(weight_db_files, conditions, variable_name_obj)
      # update condition in case of merging rds files
      conditions <- names(combined_all_data)
      combined_susie_result <- extract_variants_and_susie_results(combined_all_data, conditions)
      conditions <- names(combined_susie_result)
      weights <- consolidate_weights_list(combined_all_data, conditions, variable_name_obj, twas_weights_table)
      performance_tables <- lapply(conditions, function(condition) {
        get_nested_element(combined_all_data, c(condition, "twas_cv_result", "performance"))
      })
      names(performance_tables) <- conditions
      return(list(susie_results = combined_susie_result, weights = weights, twas_cv_performance = performance_tables))
    },
    silent = TRUE
  )
}

#' Load summary statistic data
#'
#' This function formats the input summary statistics dataframe with uniform column names
#' to fit into the SuSiE pipeline. The mapping is performed through the specified column file.
#' Additionally, it extracts sample size, case number, control number, and variance of Y.
#' Missing values in n_sample, n_case, and n_control are backfilled with median values.
#'
#' @param sumstat_path File path to the summary statistics.
#' @param column_file_path File path to the column file for mapping.
#' @param n_sample User-specified sample size. If unknown, set as 0 to retrieve from the sumstat file.
#' @param n_case User-specified number of cases.
#' @param n_control User-specified number of controls.
#'
#' @param region The region where tabix use to subset the input dataset.
#' @param target User-specified gene/phenotype name used to further subset the phenotype data.
#' @param target_column_index Filter this specific column for the target.
#' @return A list of rss_input, including the column-name-formatted summary statistics,
#' sample size (n), and var_y.
#'
#' @importFrom dplyr mutate group_by summarise
#' @importFrom magrittr %>%
#' @importFrom data.table fread
#' @export
load_rss_data <- function(sumstat_path, column_file_path, subset = TRUE, n_sample = 0, n_case = 0, n_control = 0, target = "", region = "", target_column_index = "") {
  # Read and preprocess column mapping
  column_data <- read.table(column_file_path, header = FALSE, sep = ":", stringsAsFactors = FALSE) %>%
    rename(standard = V1, original = V2)

  # Initialize sumstats variable
  sumstats <- NULL
  var_y <- NULL

  sumstats <- load_tsv_region(sumstat_path = sumstat_path, region = region, target = target, target_column_index = target_column_index)

  # Standardize column names based on mapping
  for (name in colnames(sumstats)) {
    if (name %in% column_data$original) {
      index <- which(column_data$original == name)
      colnames(sumstats)[colnames(sumstats) == name] <- column_data$standard[index]
    }
  }

  if (!"z" %in% colnames(sumstats) && all(c("beta", "se") %in%
    colnames(sumstats))) {
    sumstats$z <- sumstats$beta / sumstats$se
  }
  if (!"beta" %in% colnames(sumstats) && "z" %in% colnames(sumstats)) {
    sumstats$beta <- sumstats$z
    sumstats$se <- 1
  }
  for (col in c("n_sample", "n_case", "n_control")) {
    if (col %in% colnames(sumstats)) {
      sumstats[[col]][is.na(sumstats[[col]])] <- median(sumstats[[col]],
        na.rm = TRUE
      )
    }
  }
  if (n_sample != 0 && (n_case + n_control) != 0) {
    stop("Please provide sample size, or case number with control number, but not both")
  } else if (n_sample != 0) {
    n <- n_sample
  } else if ((n_case + n_control) != 0) {
    n <- n_case + n_control
    phi <- n_case / n
    var_y <- 1 / (phi * (1 - phi))
  } else {
    if ("n_sample" %in% colnames(sumstats)) {
      n <- median(sumstats$n_sample)
    } else if (all(c("n_case", "n_control") %in% colnames(sumstats))) {
      n <- median(sumstats$n_case + sumstats$n_control)
      phi <- median(sumstats$n_case / n)
      var_y <- 1 / (phi * (1 - phi))
    } else {
      warning("Sample size and variance of Y could not be determined from the summary statistics.")
      n <- NULL
    }
  }
  return(list(sumstats = sumstats, n = n, var_y = var_y))
}

#' Load customized tsv data
#'
#' This function load the input data. If the input sumstat data is .gz and tabixed, then can use the region parameter to subset the data
#' and filter by target column
#' Otherwise, it will only filter by target column since tabix command won't function (this apply to .tsv, .txt files)
#'
#'
#' @param sumstat_path File path to the summary statistics.
#' @param region The region where tabix use to subset the input dataset. Format: chr:start-end (eg: 9:10000-50000)
#' @param target User-specified gene/phenotype name used to further subset the phenotype data.
#' @param target_column_index Filter this specific column for the target.
#'
#' @return A dataframe of the subsetted summary statistics,
#'
#' @importFrom data.table fread
#' @export

load_tsv_region <- function(sumstat_path, region = "", target = "", target_column_index = "") {
  sumstats <- NULL
  cmd <- NULL
  if (grepl(".gz$", sumstat_path)) {
    if (is.null(sumstats) || nrow(sumstats) == 0) {
      if (target != "" && region != "" && target_column_index != "") {
        # region specified, target specified
        cmd <- paste0("zcat ", sumstat_path, " | head -1 && tabix ", sumstat_path, " ", region, " | awk '$", target_column_index, " ~ /", target, "/'")
      } else if (target != "" && region == "" && target_column_index != "") {
        # region not specified, target specified
        cmd <- paste0("zcat ", sumstat_path, " | awk '$", target_column_index, " ~ /", target, "/'")
      } else if (region != "" && (target_column_index == "" || target == "")) {
        # region specified, target not specified
        cmd <- paste0("zcat ", sumstat_path, " | head -1 && tabix ", sumstat_path, " ", region)
      } else {
        # both not specified. Instead of reading a command, read the file path instead
        cmd <- sumstat_path
      }
      sumstats <- tryCatch(
        {
          fread(cmd)
        },
        error = function(e) {
          # the gz file cannot be processed by tabix / target column missing
          stop("Data read error. Please make sure this gz file is tabixed and the target column exists.")
        }
      )
    }
  } else {
    # Non-gz files, cannot be tabixed. Load the whole dataset and only apply target filter
    warning("Not a tabixed gz file, loading the whole data.")
    if (target != "" && target_column_index != "") {
      # target specified
      sumstats <- fread(sumstat_path)
      keep_index <- which(str_detect(sumstats[[target_column_index]], target))
      sumstats <- sumstats[keep_index, ]
    } else {
      # target not specified, the whole dataset
      sumstats <- fread(sumstat_path)
    }
  }

  return(sumstats)
}
