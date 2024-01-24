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

  get_nested_element <- function(nested_list, name_vector) {
    if (is.null(name_vector)) return (NULL)
    current_element <- nested_list
    for (name in name_vector) {
      if (is.null(current_element[[name]])) {
        stop("Element not found in the list")
      }
      current_element <- current_element[[name]]
    }
    return(current_element)
  }

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


#' Colocalization Analysis Wrapper
#'
#' This function processes xQTL and multiple GWAS finemapped data files for colocalization analysis.
#'
#' @param xqtl_file Path to the xQTL RDS file.
#' @param gwas_files Vector of paths to GWAS RDS files.
#' @param gwas_finemapping_obj Optional table name in GWAS RDS files (default 'susie_fit').
#' @param xqtl_finemapping_obj Optional table name in xQTL RDS files (default 'susie_fit').
#' @param p1, p2, and p12 are results from xqtl_enrichment_wrapper (default 'p1=1e-4, p2=1e-4, p12=5e-6', same as coloc.bf_bf)

#' @return A list containing the processed xQTL and GWAS logBF matrices for colocalization analysis, coloc results, output from the compute_qtl_enrichment function
#' @examples
#' xqtl_file <- "xqtl_file.rds"
#' gwas_files <- c("gwas_file1.rds", "gwas_file2.rds")
#' result <- coloc_wrapper(xqtl_file, gwas_files)
#' @importFrom dplyr bind_rows
#' @importFrom tidyr replace_na
#' @importFrom coloc coloc.bf_bf
#' @export
coloc_wrapper <- function(xqtl_file, gwas_files, 
                          gwas_finemapping_obj = NULL, xqtl_finemapping_obj = NULL,
                          gwas_varname_obj = NULL, xqtl_varname_obj = NULL,
                          p1=1e-4, p2=1e-4, p12=5e-6, ...) {
    # Load and process GWAS data
    gwas_lbf_matrices <- lapply(gwas_files, function(file) {
        raw_data <- readRDS(file)
        gwas_data <- if (!is.null(gwas_finemapping_obj)) get_nested_element(raw_data, gwas_finemapping_obj) else raw_data 
        gwas_lbf_matrix <- as.data.frame(gwas_data$lbf_variable)
        gwas_lbf_matrix <- gwas_lbf_matrix[gwas_data$V > 0,]
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
    raw_data <- readRDS(file)
    xqtl_data <- if (!is.null(xqtl_finemapping_obj)) get_nested_element(raw_data, xqtl_finemapping_obj) else raw_data
    xqtl_lbf_matrix <- as.data.frame(xqtl_data$lbf_variable)
    xqtl_lbf_matrix <- xqtl_lbf_matrix[xqtl_data$V > 0,]
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
    return(coloc_res)
}