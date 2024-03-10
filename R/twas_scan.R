#' TWAS Analysis
#'
#' Performs TWAS analysis using the provided weights matrix, GWAS summary statistics database,
#' and LD matrix. It extracts the necessary GWAS summary statistics and LD matrix based on the
#' specified variants and computes the z-score and p-value for each gene.
#'
#' @param weights_matrix A matrix containing weights for all methods.
#' @param gwas_sumstats_db A data frame containing the GWAS summary statistics.
#' @param LD_matrix A matrix representing linkage disequilibrium between variants.
#' @param extract_variants_objs A vector of variant identifiers to extract from the GWAS and LD matrix.
#'
#' @return A list with TWAS z-scores and p-values across four methods for each gene.
#' @export
twas_analysis <- function(weights_matrix, gwas_sumstats_db, LD_matrix, extract_variants_objs) {
  #
  # Extract gwas_sumstats
  gwas_sumstats_subset <- gwas_sumstats_db[match(extract_variants_objs, gwas_sumstats_db$variant_id), ]
  # Validate that the GWAS subset is not empty
  if (nrow(gwas_sumstats_subset) == 0 | all(is.na(gwas_sumstats_subset))) {
    stop("No GWAS summary statistics found for the specified variants.")
  }
  # Check if extract_variants_objs are in the rownames of LD_matrix
  valid_indices <- extract_variants_objs %in% rownames(LD_matrix)
  if (!any(valid_indices)) {
    stop("None of the specified variants are present in the LD matrix.")
  }
  # Extract only the valid indices from extract_variants_objs
  valid_variants_objs <- extract_variants_objs[valid_indices]
  # Extract LD_matrix subset using valid indices
  LD_matrix_subset <- LD_matrix[valid_variants_objs, valid_variants_objs]
  # Caculate the z score and pvalue of each gene
  twas_z_pval <- apply(
    as.matrix(weights_matrix), 2,
    function(x) twas_z(x, gwas_sumstats_subset$z, R = LD_matrix_subset)
  )
  return(twas_z_pval)
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
#' @return A consolidated list of weights for the specified condition and a list of susie_trimmed_results.
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
load_twas_weights <- function(weight_db_files, conditions = NULL,
                              variable_name_obj = "variant_names",
                              twas_weights_table = "twas_weights") {
  ## Internal function to load and validate data from RDS files
  load_and_validate_data <- function(weight_db_files, conditions, variable_name_obj) {
    all_data <- lapply(weight_db_files, readRDS)
    unique_regions <- unique(unlist(lapply(all_data, function(data) names(data))))
    # Check if region from all RDS files are the same
    if (length(unique_regions) != 1) {
      stop("The RDS files do not refer to the same region.")
    } else {
      # Assuming all data refer to the same region, now combine data by conditions
      combined_all_data <- do.call("c", lapply(all_data, function(data) data[[1]]))
    }
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
  # Only extract the variant_names and susie_result_trimmed
  extract_variants_and_susie_results <- function(combined_all_data, conditions) {
    combined_susie_result_trimmed <- lapply(conditions, function(condition) {
      list(
        variant_names = get_nested_element(combined_all_data, c(condition, "preset_variants_result", "variant_names")),
        susie_result_trimmed = get_nested_element(combined_all_data, c(condition, "preset_variants_result", "susie_result_trimmed")),
        top_loci = get_nested_element(combined_all_data, c(condition, "preset_variants_result", "top_loci")),
        region_info = get_nested_element(combined_all_data, c(condition, "region_info"))
      )
    })
    names(combined_susie_result_trimmed) <- conditions
    return(combined_susie_result_trimmed)
  }
  # Internal function to align and merge weight matrices
  align_and_merge <- function(weights_list, variable_objs) {
    # Get the complete list of variant names across all files
    all_variants <- unique(unlist(variable_objs))
    consolidated_list <- list()
    # Fill the matrix with weights, aligning by variant names
    for (i in seq_along(weights_list)) {
      # Initialize the temp matrix with zeros
      existing_colnames <- character(0)
      temp_matrix <- matrix(0, nrow = length(all_variants), ncol = ncol(weights_list[[i]]))
      rownames(temp_matrix) <- all_variants
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
    # Set default for 'conditions' if they are not specified
    if (is.null(conditions)) {
      conditions <- names(combined_all_data)
    }
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
        get_nested_element(combined_all_data, c(condition, twas_weights_table, variable_name_obj))
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
      combined_susie_result_trimmed <- extract_variants_and_susie_results(combined_all_data, conditions)
      weights <- consolidate_weights_list(combined_all_data, conditions, variable_name_obj, twas_weights_table)
      return(list(susie_results = combined_susie_result_trimmed, weights = weights))
    },
    silent = TRUE
  )
}
