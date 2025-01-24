#' Univariate Analysis Pipeline
#'
#' This function performs univariate analysis for fine-mapping and Transcriptome-Wide Association Study (TWAS)
#' with optional cross-validation.
#'
#' @param X A matrix of genotype data where rows represent samples and columns represent genetic variants.
#' @param Y A vector of phenotype measurements.
#' @param X_scalar A scalar or vector to rescale X to its original scale.
#' @param Y_scalar A scalar to rescale Y to its original scale.
#' @param maf A vector of minor allele frequencies for each variant in X.
#' @param X_variance Optional variance of X. Default is NULL.
#' @param other_quantities A list of other quantities to be passed to susie_post_processor. Default is an empty list.
#' @param imiss_cutoff Individual missingness cutoff. Default is 1.0.
#' @param maf_cutoff Minor allele frequency cutoff. Default is NULL.
#' @param xvar_cutoff Variance cutoff for X. Default is 0.05.
#' @param ld_reference_meta_file An optional path to a file containing linkage disequilibrium reference data. Default is NULL.
#' @param pip_cutoff_to_skip Cutoff value for skipping analysis based on PIP values. Default is 0.
#' @param init_L Initial number of components for SuSiE model optimization. Default is 5.
#' @param max_L The maximum number of components in SuSiE. Default is 20.
#' @param l_step Step size for increasing the number of components during SuSiE optimization. Default is 5.
#' @param signal_cutoff Cutoff value for signal identification in PIP values. Default is 0.025.
#' @param coverage A vector of coverage probabilities for credible sets. Default is c(0.95, 0.7, 0.5).
#' @param twas_weights Whether to compute TWAS weights. Default is TRUE.
#' @param sample_partition Sample partition for cross-validation. Default is NULL.
#' @param max_cv_variants The maximum number of variants to be included in cross-validation. Default is -1 (no limit).
#' @param cv_folds The number of folds to use for cross-validation. Default is 5.
#' @param cv_threads The number of threads to use for parallel computation in cross-validation. Default is 1.
#' @param verbose Verbosity level. Default is 0.
#'
#' @return A list containing the univariate analysis results.
#' @importFrom susieR susie
#' @export
univariate_analysis_pipeline <- function(
    # input data
    X,
    Y,
    maf,
    X_scalar = 1,
    Y_scalar = 1,
    X_variance = NULL,
    other_quantities = list(),
    # filters
    imiss_cutoff = 1.0,
    maf_cutoff = NULL,
    xvar_cutoff = 0,
    ld_reference_meta_file = NULL,
    pip_cutoff_to_skip = 0,
    # methods parameter configuration
    init_L = 5,
    max_L = 20,
    l_step = 5,
    # fine-mapping results summary
    signal_cutoff = 0.025,
    coverage = c(0.95, 0.7, 0.5),
    # TWAS weights and CV for TWAS weights
    twas_weights = TRUE,
    sample_partition = NULL,
    max_cv_variants = -1,
    cv_folds = 5,
    cv_threads = 1,
    verbose = 0) {
  # Input validation
  if (!is.matrix(X) || !is.numeric(X)) stop("X must be a numeric matrix")
  if (!is.vector(Y) && !(is.matrix(Y) && ncol(Y) == 1) || !is.numeric(Y)) stop("Y must be a numeric vector or a single column matrix")
  if (nrow(X) != length(Y)) stop("X and Y must have the same number of rows/length")
  if (!is.numeric(maf) || length(maf) != ncol(X)) stop("maf must be a numeric vector with length equal to the number of columns in X")
  if (any(maf < 0 | maf > 1)) stop("maf values must be between 0 and 1")
  if (!is.numeric(X_scalar) || (length(X_scalar) != 1 && length(X_scalar) != ncol(X))) stop("X_scalar must be a numeric scalar or vector with length equal to the number of columns in X")
  if (!is.numeric(Y_scalar) || length(Y_scalar) != 1) stop("Y_scalar must be a numeric scalar")
  if (!is.numeric(init_L) || init_L <= 0) stop("init_L must be a positive integer")
  if (!is.numeric(max_L) || max_L <= 0) stop("max_L must be a positive integer")
  if (!is.numeric(l_step) || l_step <= 0) stop("l_step must be a positive integer")

  # Initial PIP check
  if (pip_cutoff_to_skip != 0) {
    if (pip_cutoff_to_skip < 0) {
      # automatically determine the cutoff to use
      pip_cutoff_to_skip <- 3 * 1 / ncol(X)
    }
    top_model_pip <- susie(X, Y, L = 1)$pip
    if (!any(top_model_pip > pip_cutoff_to_skip)) {
      message(paste("Skipping follow-up analysis: No signals above PIP threshold", pip_cutoff_to_skip, "in initial model screening."))
      return(list())
    } else {
      message(paste("Follow-up on region because signals above PIP threshold", pip_cutoff_to_skip, "were detected in initial model screening."))
    }
  }

  # Filter variants if LD reference is provided
  if (!is.null(ld_reference_meta_file)) {
    variants_kept <- filter_variants_by_ld_reference(colnames(X), ld_reference_meta_file)
    X <- X[, variants_kept$data, drop = FALSE]
    maf <- maf[variants_kept$idx]
    if (length(X_scalar) > 1) X_scalar <- X_scalar[variants_kept$idx]
  }

  # Filter X based on missingness, MAF, and variance
  if (!is.null(imiss_cutoff) || !is.null(maf_cutoff)) {
    X_filtered <- filter_X(X, imiss_cutoff, maf_cutoff, var_thresh = xvar_cutoff, maf = maf, X_variance = X_variance)
    kept_indices <- match(colnames(X_filtered), colnames(X))
    maf <- maf[kept_indices]
    if (length(X_scalar) > 1) X_scalar <- X_scalar[kept_indices]
    X <- X_filtered
  }

  # Main analysis
  st <- proc.time()
  res <- list()

  # SuSiE analysis with optimization
  message("Fitting SuSiE model on input data with L optimization...")
  res$susie_fitted <- susie_wrapper(X, Y,
    init_L = init_L, max_L = max_L, l_step = l_step,
    refine = TRUE, coverage = coverage[1]
  )

  # Process SuSiE results
  susie_result_trimmed <- susie_post_processor(
    res$susie_fitted, X, Y, X_scalar, Y_scalar, maf,
    secondary_coverage = if (length(coverage) > 1) coverage[-1] else NULL,
    signal_cutoff = signal_cutoff,
    other_quantities = other_quantities
  )
  res <- c(res, susie_result_trimmed)
  res$total_time_elapsed <- proc.time() - st

  # TWAS weights and cross-validation
  if (twas_weights) {
    res$twas_weights_result <- twas_weights_pipeline(
      X, Y, res$susie_fitted,
      cv_folds = cv_folds,
      max_cv_variants = max_cv_variants,
      cv_threads = cv_threads,
      sample_partition = sample_partition
    )
    if ("top_loci" %in% names(res)) res$twas_weights_result$susie_weights_intermediate$top_loci <- res$top_loci
  }

  return(res)
}

#' RSS Analysis Pipeline
#'
#' This function performs an end-to-end RSS analysis pipeline, including data loading,
#' preprocessing, quality control, imputation, and SuSiE RSS analysis. It provides flexibility
#' in specifying various analysis options and parameters.
#'
#' @param sumstat_path File path to the summary statistics.
#' @param column_file_path File path to the column file for mapping.
#' @param LD_data A list containing combined LD variants data that is generated by load_LD_matrix.
#' @param n_sample User-specified sample size. If unknown, set as 0 to retrieve from the sumstat file.
#' @param n_case User-specified number of cases.
#' @param n_control User-specified number of controls.
#' @param skip_region A character vector specifying regions to be skipped in the analysis (optional).
#'                    Each region should be in the format "chrom:start-end" (e.g., "1:1000000-2000000").
#' @param L Initial number of causal configurations to consider in the analysis (default: 8).
#' @param max_L Maximum number of causal configurations to consider when dynamically adjusting L (default: 20).
#' @param l_step Step size for increasing L when the limit is reached during dynamic adjustment (default: 5).
#' @param qc_method Quality control method to use. Options are "rss_qc", "dentist", or "slalom" (default: "rss_qc").
#' @param analysis_method Analysis method to use. Options are "susie_rss", "single_effect", or "bayesian_conditional_regression" (default: "susie_rss").
#' @param impute Logical; if TRUE, performs imputation for outliers identified in the analysis (default: TRUE).
#' @param impute_opts A list of imputation options including rcond, R2_threshold, and minimum_ld (default: list(rcond = 0.01, R2_threshold = 0.6, minimum_ld = 5)).
#' @param coverage Coverage levels for SuSiE RSS analysis (default: c(0.95, 0.7, 0.5)).
#' @param pip_cutoff_to_skip PIP cutoff to skip imputation (default: 0).
#' @param signal_cutoff Signal cutoff for susie_post_processor (default: 0.025).
#'
#' @return A list containing the final_result and input_rss_data.
#'   - final_result: A list containing the results of various SuSiE RSS analyses.
#'   - input_rss_data: A processed data frame containing summary statistics after preprocessing.
#'
#' @importFrom magrittr %>%
#' @export
rss_analysis_pipeline <- function(
    sumstat_path, column_file_path, LD_data, n_sample = 0, n_case = 0, n_control = 0, target = "", region = "", target_column_index = "", skip_region = NULL,
    qc_method = c("rss_qc", "dentist", "slalom"),
    finemapping_method = c("susie_rss", "single_effect", "bayesian_conditional_regression"),
    finemapping_opts = list(
      init_L = 5, max_L = 20, l_step = 5,
      coverage = c(0.95, 0.7, 0.5), signal_cutoff = 0.025
    ),
    impute = TRUE, impute_opts = list(rcond = 0.01, R2_threshold = 0.6, minimum_ld = 5, lamb = 0.01),
    pip_cutoff_to_skip = 0, remove_indels = FALSE) {
  res <- list()
  rss_input <- load_rss_data(
    sumstat_path = sumstat_path, column_file_path = column_file_path,
    n_sample = n_sample, n_case = n_case, n_control = n_control, target = target, region = region, target_column_index = target_column_index
  )

  sumstats <- rss_input$sumstats
  n <- rss_input$n
  var_y <- rss_input$var_y

  # Preprocess the input data
  preprocess_results <- rss_basic_qc(sumstats, LD_data, skip_region = skip_region, remove_indels = remove_indels)
  sumstats <- preprocess_results$sumstats
  LD_mat <- preprocess_results$LD_mat

  if (pip_cutoff_to_skip != 0) {
    if (pip_cutoff_to_skip < 0) {
      # automatically determine the cutoff to use
      pip_cutoff_to_skip <- 3 * 1 / ncol(X)
    }
    top_model_pip <- susie_rss_wrapper(z = sumstats$z, R = LD_mat, L = 1, n = n, var_y = var_y)$pip
    if (!any(top_model_pip > pip_cutoff_to_skip)) {
      message(paste("Skipping follow-up analysis: No signals above PIP threshold", pip_cutoff_to_skip, "in initial model screening."))
      return(list(rss_data_analyzed = sumstats))
    } else {
      message(paste("Follow-up on region because signals above PIP threshold", pip_cutoff_to_skip, "were detected in initial model screening."))
    }
  }

  # Perform quality control
  if (!is.null(qc_method)) {
    qc_results <- summary_stats_qc(sumstats, LD_data, n = n, var_y = var_y, method = qc_method)
    sumstats <- qc_results$sumstats
    LD_mat <- qc_results$LD_mat
  }

  # Perform imputation
  if (impute) {
    impute_results <- raiss(LD_data$ref_panel, sumstats, LD_data$combined_LD_matrix, rcond = impute_opts$rcond, R2_threshold = impute_opts$R2_threshold, minimum_ld = impute_opts$minimum_ld, lamb = impute_opts$lamb)
    sumstats <- impute_results$result_filter
    LD_mat <- impute_results$LD_mat
  }

  # Perform fine-mapping
  if (!is.null(finemapping_method)) {
    pri_coverage <- finemapping_opts$coverage[1]
    sec_coverage <- if (length(finemapping_opts$coverage) > 1) finemapping_opts$coverage[-1] else NULL
    res <- susie_rss_pipeline(sumstats, LD_mat,
      n = n, var_y = var_y,
      L = finemapping_opts$init_L, max_L = finemapping_opts$max_L, l_step = finemapping_opts$l_step,
      analysis_method = finemapping_method,
      coverage = pri_coverage,
      secondary_coverage = sec_coverage,
      signal_cutoff = finemapping_opts$signal_cutoff
    )
    if (!is.null(qc_method)) {
      res$outlier_number <- qc_results$outlier_number
    }
  }
  if (impute & !is.null(qc_method)) {
    method_name <- paste0(toupper(qc_method), "_RAISS_imputed")
  } else if (!impute & !is.null(qc_method)) {
    method_name <- toupper(qc_method)
  } else {
    method_name <- "NO_QC"
  }
  result_list <- list()
  result_list[[method_name]] <- res
  result_list[["rss_data_analyzed"]] <- sumstats
  return(result_list)
}
