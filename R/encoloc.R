#' xQTL GWAS Enrichment Analysis
#'
#' This function processes GWAS and xQTL finemapped data files and then computes QTL enrichment.
#' For details on the parameters `pi_gwas`, `pi_qtl`, `lambda`, `ImpN`, and `num_threads`,
#' refer to the documentation of the `compute_qtl_enrichment` function.
#'
#' @param xqtl_files Vector of xQTL RDS file paths.
#' @param gwas_files Vector of GWAS RDS file paths.
#' @param xqtl_finemapping_obj Optional table name in xQTL RDS files (default 'susie_fit').
#' @param gwas_finemapping_obj Optional table name in GWAS RDS files (default 'susie_fit').
#' @param pi_gwas Optional parameter for GWAS enrichment estimation (see `compute_qtl_enrichment`).
#' @param pi_qtl Optional parameter for xQTL enrichment estimation (see `compute_qtl_enrichment`).
#' @param lambda Shrinkage parameter for enrichment computation (see `compute_qtl_enrichment`).
#' @param ImpN Importance parameter for enrichment computation (see `compute_qtl_enrichment`).
#' @param num_threads Number of threads for parallel processing (see `compute_qtl_enrichment`).
#' @return The output from the compute_qtl_enrichment function.
#' @examples
#' gwas_files <- c("gwas_file1.rds", "gwas_file2.rds")
#' xqtl_files <- c("xqtl_file1.rds", "xqtl_file2.rds")
#' result <- xqtl_enrichment_wrapper(gwas_files, xqtl_files)
#' @export
xqtl_enrichment_wrapper <- function(xqtl_files, gwas_files, 
                                    gwas_finemapping_obj = NULL, xqtl_finemapping_obj = NULL,
                                    gwas_varname_obj = NULL, xqtl_varname_obj = NULL,
                                    pi_gwas = NULL, pi_qtl = NULL, 
                                    lambda = 1.0, ImpN = 25,
                                    num_threads = 1) {

  process_finemapped_data <- function(xqtl_files, gwas_files,
                                    gwas_finemapping_obj = NULL, xqtl_finemapping_obj = NULL,
                                    gwas_varname_obj = NULL, xqtl_varname_obj = NULL) {
    # Process GWAS data
    gwas_pip <- list()
    for (file in gwas_files) {
        raw_data <- readRDS(file)
        gwas_data <- if (!is.null(gwas_finemapping_obj)) get_nested_element(raw_data, gwas_finemapping_obj) else raw_data 
        pip <- gwas_data$pip
        if (!is.null(gwas_varname_obj)) names(pip) <- get_nested_element(raw_data, gwas_varname_obj)
        gwas_pip <- c(gwas_pip, list(pip))
    }

    # Check for unique variant names in GWAS pip vectors
    all_variant_names <- unique(unlist(lapply(gwas_pip, names)))
    if(length(unique(all_variant_names)) != length(all_variant_names)) {
        stop("Non-unique variant names found in GWAS data with different pip values.")
    }
    gwas_pip <- unlist(gwas_pip)

    # Process xQTL data
    xqtl_data <- lapply(xqtl_files, function(file) {
        raw_data <- readRDS(file)
        xqtl_data <- if (!is.null(xqtl_finemapping_obj)) get_nested_element(raw_data, xqtl_finemapping_obj) else raw_data
        list(alpha = xqtl_data$alpha, pip = setNames(xqtl_data$pip, get_nested_element(raw_data, xqtl_varname_obj)), 
            prior_variance = xqtl_data$V)
    })
    
    # Return results as a list
    return(list(gwas_pip = gwas_pip, xqtl_data = xqtl_data))
  } 
  
  # Load data
  dat <- process_finemapped_data(xqtl_files, gwas_files, gwas_finemapping_obj, xqtl_finemapping_obj, gwas_varname_obj, xqtl_varname_obj)
  # Compute QTL enrichment
  return(compute_qtl_enrichment(gwas_pip = dat$gwas_pip, susie_qtl_regions = dat$xqtl_data,
                                pi_gwas = pi_gwas, pi_qtl = pi_qtl,
                                lambda = lambda, ImpN = ImpN,
                                num_threads = num_threads))
} 



#' Functions for  Colocalization Analysis Wrapper

# Function to filter and order colocalization results
filter_and_order_coloc_results <- function(coloc_results_fil) {
  # Ensure the input has more than one column
  if (ncol(coloc_results_fil) <= 1) {
    stop("Insufficient number of columns in colocalization results")
  }

  cs_num <- ncol(coloc_results_fil) - 1
  ordered_results <- list()

  for (n in 1:cs_num) {
    # Selecting relevant columns and ordering
    tmp_coloc_results_fil <- coloc_results_fil[, c(1, n+1)] %>% .[order(.[, 2], decreasing = TRUE), ]
    ordered_results[[n]] <- tmp_coloc_results_fil
  }

  return(ordered_results)
}

# Function to calculate cumulative sum
calculate_cumsum <- function(coloc_results) {
  cumsum(coloc_results[, 2])
}

# FIXME: better way to get analysis region?
# Function to extract variants and analysis region
extract_variants_and_region <- function(analysis_script_obj, variants) {
  variants <- variants %>% gsub("chr", "", .)
  analysis_region <- str_extract(analysis_script_obj, 'cis_window = "chr\\d+:[0-9]+-[0-9]+') %>% gsub('cis_window = "', '', .)
  list(variants = variants, analysis_region = analysis_region)
}

# Function to load and extract LD matrix
load_and_extract_ld_matrix <- function(ld_meta_file_path, analysis_region, variants) {
  # This is a placeholder for loading LD matrix, adjust as per your actual function
  ld_ref <- load_LD_matrix(LD_meta_file_path = ld_meta_file_path, region = analysis_region)
  ext_ld <- ld_ref[[2]][variants, variants]
  ext_ld
}

#' @importFrom susieR get_purity
# Function to calculate purity
calculate_purity <- function(variants, ext_ld, squared) {
  # This is a placeholder for calculating purity, adjust as per your actual function
  purity <- matrix(get_purity(variants, Xcorr = ext_ld, squared), 1, 3)
  purity
}

# Main processing function
process_coloc_results <- function(coloc_results, LD_meta_file_path,analysis_script_obj, PPH4_thres = 0.8, coloc_pip_thres = 0.95, squared = FALSE, min_abs_corr = 0.5, null_index = 0, coloc_index = "PP.H4.abf") {
  # Extract PIP values from coloc_result summary
  coloc_summary <- as.data.frame(coloc_result$summary)
  coloc_pip <- coloc_summary[, grepl("PP", colnames(coloc_summary))]

  # Filter and extract relevant columns from coloc_result results
  # PP.H4 is highest and > 0.8
  coloc_results_df <- as.data.frame(coloc_result$results)
  coloc_filter <- apply(coloc_pip, 1, function(row) {
    max_index <- which.max(row)
    max_value <- row[max_index]
    return(max_value > PPH4_thres && colnames(coloc_pip)[max_index] == coloc_index)
  })
  coloc_results_fil <- coloc_results_df[, c(1, which(coloc_filter) + 1), drop = FALSE]
  coloc_summary_fil <- coloc_summary[which(coloc_filter),, drop = FALSE]

  #prepare to calculate purity
  ordered_results <- filter_and_order_coloc_results(coloc_results_fil)
  cs <- list()
  purity <- NULL

  for (n in 1:length(ordered_results)) {
    tmp_coloc_results_fil <- ordered_results[[n]]
    tmp_coloc_results_fil_csm <- calculate_cumsum(tmp_coloc_results_fil)
    cs[[n]] <- tmp_coloc_results_fil[, 1][1:(which(tmp_coloc_results_fil_csm > coloc_pip_thres) %>% min)]

    # Extract variants and analysis region
    extraction_result <- extract_variants_and_region(analysis_script_obj, cs[[n]] )
    variants <- extraction_result$variants
    analysis_region <- extraction_result$analysis_region

    # Load and extract LD matrix
    ext_ld <- load_and_extract_ld_matrix(LD_meta_file_path, analysis_region, variants)

    # Calculate purity
    if (null_index > 0 && null_index %in% variants) {
      purity <- rbind(purity, c(-9, -9, -9))
    } else {
      current_purity <- calculate_purity(variants, ext_ld, squared)
      purity <- rbind(purity, current_purity)
    }
  }

  # Process purity data
  purity <- as.data.frame(purity)
  if (squared) {
    colnames(purity) <- c("min.sq.corr", "mean.sq.corr", "median.sq.corr")
  } else {
    colnames(purity) <- c("min.abs.corr", "mean.abs.corr", "median.abs.corr")
  }

  threshold <- ifelse(squared, min_abs_corr^2, min_abs_corr)
  is_pure <- which(purity[, 1] >= threshold)

  # Finalize the result
  coloc_res <- list()
  if (length(is_pure) > 0) {
    cs <- cs[is_pure]
    purity <- purity[is_pure, ]
    true_summary <- coloc_summary_fil[is_pure, ]
    coloc_res$sets <- list(cs = cs, purity = purity, true_summary = true_summary)
  }

  return(coloc_res)
}




#' Colocalization Analysis Wrapper
#'
#' This function processes xQTL and multiple GWAS finemapped data files for colocalization analysis.
#'
#' @param xqtl_file Path to the xQTL RDS file.
#' @param gwas_files Vector of paths to GWAS RDS files.
#' @param gwas_finemapping_obj Optional table name in GWAS RDS files (default 'susie_fit').
#' @param xqtl_finemapping_obj Optional table name in xQTL RDS files (default 'susie_fit').
#' @param gwas_varname_obj Optional table name in GWAS RDS files (default 'susie_fit').
#' @param xqtl_varname_obj Optional table name in xQTL RDS files (default 'susie_fit').
#' @param xqtl_script_obj Optional table name in xQTL RDS files (default 'susie_fit').
#' @param LD_meta_file_path Path to the metadata of LD reference.
#' @param p1, p2, and p12 are results from xqtl_enrichment_wrapper (default 'p1=1e-4, p2=1e-4, p12=5e-6', same as coloc.bf_bf).
#' @param prior_tol When the prior variance is estimated, compare the estimated value to \code{prior_tol} at the end of the computation, 
#'   and exclude a single effect from PIP computation if the estimated prior variance is smaller than this tolerance value.
#' @return A list containing the processed xQTL and GWAS logBF matrices for colocalization analysis, coloc results, output from the compute_qtl_enrichment function
#' @examples
#' xqtl_file <- "xqtl_file.rds"
#' gwas_files <- c("gwas_file1.rds", "gwas_file2.rds")
#' result <- coloc_wrapper(xqtl_file, gwas_files, LD_meta_file_path)
#' @importFrom dplyr bind_rows
#' @importFrom tidyr replace_na
#' @importFrom coloc coloc.bf_bf
#' @export
coloc_wrapper <- function(xqtl_file, gwas_files, 
                          gwas_finemapping_obj = NULL, xqtl_finemapping_obj = NULL,
                          gwas_varname_obj = NULL, xqtl_varname_obj = NULL, 
                          xqtl_script_obj = NULL, 
                          LD_meta_file_path, prior_tol = 1e-9,
                          p1=1e-4, p2=1e-4, p12=5e-6, ...) {
    # Load and process GWAS data
    gwas_lbf_matrices <- lapply(gwas_files, function(file) {
        raw_data <- readRDS(file)
        gwas_data <- if (!is.null(gwas_finemapping_obj)) get_nested_element(raw_data, gwas_finemapping_obj) else raw_data 
        gwas_lbf_matrix <- as.data.frame(gwas_data$lbf_variable)
        gwas_lbf_matrix <- gwas_lbf_matrix[gwas_data$V > prior_tol,]
        if (!is.null(gwas_varname_obj)) names(pip) <- get_nested_element(raw_data, gwas_varname_obj)
        return(gwas_lbf_matrix)
    })

    # Validate uniqueness of column names across GWAS matrices
    all_gwas_colnames <- unique(unlist(lapply(gwas_lbf_matrices, colnames)))
    if(length(all_gwas_colnames) != sum(sapply(gwas_lbf_matrices, ncol))) {
        stop("Duplicate variant names found across GWAS regions analyzed. This is not expected.")
    }

    # Combine GWAS matrices and replace NAs with zeros
    combined_gwas_lbf_matrix <- bind_rows(gwas_lbf_matrices) %>%
                                mutate(across(everything(), ~replace_na(., 0)))


    # Process xQTL data
    raw_data <- readRDS(xqtl_file)
    xqtl_data <- if (!is.null(xqtl_finemapping_obj)) get_nested_element(raw_data, xqtl_finemapping_obj) else raw_data
    xqtl_lbf_matrix <- as.data.frame(xqtl_data$lbf_variable)
    xqtl_lbf_matrix <- xqtl_lbf_matrix[xqtl_data$V > prior_tol,]
    if (!is.null(xqtl_varname_obj)) colnames(xqtl_lbf_matrix) <- get_nested_element(raw_data, xqtl_varname_obj)

    #add 'chr' in colnames 
    add_chr_prefix <- function(df) {
      colnames(df) <- ifelse(grepl("chr", colnames(df)), colnames(df), paste0("chr", colnames(df)))
      return(df)
    }

    combined_gwas_lbf_matrix <- add_chr_prefix(combined_gwas_lbf_matrix)
    xqtl_lbf_matrix <- add_chr_prefix(xqtl_lbf_matrix)
    

    # Match column names and reorder matrices
    common_colnames <- intersect(colnames(xqtl_lbf_matrix), colnames(combined_gwas_lbf_matrix))
    xqtl_lbf_matrix <- xqtl_lbf_matrix[, common_colnames, drop = FALSE] %>% as.matrix
    combined_gwas_lbf_matrix <- combined_gwas_lbf_matrix[, common_colnames, drop = FALSE] %>% as.matrix

    # Report the number of dropped columns from xQTL matrix
    num_dropped_cols <- length(setdiff(colnames(xqtl_lbf_matrix), common_colnames))
    message("Number of columns dropped from xQTL matrix: ", num_dropped_cols)

    # COLOC function 
    coloc_res <- coloc.bf_bf(xqtl_lbf_matrix, combined_gwas_lbf_matrix, p1 = p1, p2 = p2, p12 = p12, ...)
    coloc_res <- c(coloc_res,process_coloc_results(coloc_result, LD_meta_file_path, get_nested_element(raw_data, xqtl_script_obj)))
    # post processing for coloc results
    return(coloc_res)
}