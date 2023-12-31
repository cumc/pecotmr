#' xQTL GWAS Enrichment Analysis
#'
#' This function processes GWAS and xQTL finemapped data files and then computes QTL enrichment.
#' For details on the parameters `pi_gwas`, `pi_qtl`, `lambda`, `ImpN`, and `num_threads`,
#' refer to the documentation of the `compute_qtl_enrichment` function.
#'
#' @param gwas_finemapped_data Vector of GWAS RDS file paths.
#' @param xqtl_finemapped_data Vector of xQTL RDS file paths.
#' @param gwas_finemapping_obj Optional table name in GWAS RDS files (default 'susie_fit').
#' @param xqtl_finemapping_obj Optional table name in xQTL RDS files (default 'susie_fit').
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
xqtl_enrichment_wrapper <- function(gwas_finemapped_data, xqtl_finemapped_data,
                                    gwas_finemapping_obj = "susie_fit",
                                    xqtl_finemapping_obj = "susie_fit",
                                    pi_gwas = NULL, pi_qtl = NULL, 
                                    lambda = 1.0, ImpN = 25,
                                    num_threads = 1) {
  process_finemapped_data <- function(gwas_finemapped_data, xqtl_finemapped_data,
                                    gwas_finemapping_obj = "susie_fit",
                                    xqtl_finemapping_obj = "susie_fit") {
    # Process GWAS data
    gwas_pip <- list()
    for (file in gwas_finemapped_data) {
        gwas_data <- readRDS(file)[[gwas_finemapping_obj]]
        pip <- gwas_data$pip
        names(pip) <- gwas_data$variant_names
        gwas_pip <- c(gwas_pip, list(pip))
    }

    # Check for unique variant names in GWAS pip vectors
    all_variant_names <- unique(unlist(lapply(gwas_pip, names)))
    if(length(unique(all_variant_names)) != length(all_variant_names)) {
        stop("Non-unique variant names found in GWAS data with different pip values.")
    }
    gwas_pip <- unlist(gwas_pip)

    # Process xQTL data
    xqtl_data <- lapply(xqtl_finemapped_data, function(file) {
        xqtl_data <- readRDS(file)[[xqtl_finemapping_obj]]
        list(alpha = xqtl_data$alpha,
             pip = setNames(xqtl_data$pip, xqtl_data$variant_names),
             prior_variance = xqtl_data$prior_variance)
    })

    # Return results as a list
    return(list(gwas_pip = gwas_pip, xqtl_data = xqtl_data))
  } 
  
  # Load data
  dat <- process_finemapped_data(gwas_finemapped_data, xqtl_finemapped_data,
                                                 gwas_finemapping_obj, xqtl_finemapping_obj)

  # Compute QTL enrichment
  return(compute_qtl_enrichment(gwas_pip = dat$gwas_pip, susie_qtl_regions = dat$xqtl_data,
                                pi_gwas = pi_gwas, pi_qtl = pi_qtl,
                                lambda = lambda, ImpN = ImpN,
                                num_threads = num_threads))
} 


#' Colocalization Analysis Wrapper
#'
#' This function processes xQTL and multiple GWAS finemapped data files for colocalization analysis.
#'
#' @param xqtl_file Path to the xQTL RDS file.
#' @param gwas_files Vector of paths to GWAS RDS files.
#' @param gwas_finemapping_obj Optional table name in GWAS RDS files (default 'susie_fit').
#' @param xqtl_finemapping_obj Optional table name in xQTL RDS files (default 'susie_fit').
#' @return A list containing the processed xQTL and GWAS logBF matrices for colocalization analysis.
#' @examples
#' xqtl_file <- "xqtl_file.rds"
#' gwas_files <- c("gwas_file1.rds", "gwas_file2.rds")
#' result <- coloc_wrapper(xqtl_file, gwas_files)
#' @importFrom dplyr bind_rows
#' @importFrom tidyr replace_na
#' @export
coloc_wrapper <- function(xqtl_file, gwas_files, 
                          gwas_finemapping_obj = "susie_fit", 
                          xqtl_finemapping_obj = "susie_fit") {
    # Load and process GWAS data
    gwas_lbf_matrices <- lapply(gwas_files, function(file_path) {
        gwas_data <- readRDS(file_path)[[gwas_finemapping_obj]]
        gwas_lbf_matrix <- as.data.frame(gwas_data$lbf_variable)
        colnames(gwas_lbf_matrix) <- gwas_data$variant_names
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
    xqtl_data <- readRDS(xqtl_file)[[xqtl_finemapping_obj]]
    xqtl_lbf_matrix <- as.data.frame(xqtl_data$lbf_variable)
    colnames(xqtl_lbf_matrix) <- xqtl_data$variant_names

    # Match column names and reorder matrices
    common_colnames <- intersect(colnames(xqtl_lbf_matrix), colnames(combined_gwas_lbf_matrix))
    xqtl_lbf_matrix <- xqtl_lbf_matrix[, common_colnames, drop = FALSE]
    combined_gwas_lbf_matrix <- combined_gwas_lbf_matrix[, common_colnames, drop = FALSE]

    # Report the number of dropped columns from xQTL matrix
    num_dropped_cols <- length(setdiff(colnames(xqtl_lbf_matrix), common_colnames))
    message("Number of columns dropped from xQTL matrix: ", num_dropped_cols)

    # Return the processed data for now
    # FIXME: To complete this, we need to apply colocalization analysis on these single effect logBF and the enrichment analysis results computed by `xqtl_enrichment_wrapper`
    return(list(xqtl_lbf_matrix = xqtl_lbf_matrix, combined_gwas_lbf_matrix = combined_gwas_lbf_matrix))
}