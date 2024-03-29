#' Univariate analysis pipeline
#' @importFrom susieR susie
#' @export
univariate_analysis_pipeline <- function(X, Y, X_scalar, Y_scalar, maf, other_quantities = list(),
                                         pip_cutoff_to_skip = 0, region_info = NULL,
                                         finemapping = TRUE, twas_weights = TRUE,
                                         finemapping_opts = list(
                                           init_L = 5, max_L = 20, l_step = 5,
                                           coverage = c(0.95, 0.7, 0.5), signal_cutoff = 0.025
                                         ),
                                         twas_weights_opts = list(
                                           cv_folds = 5, min_cv_maf = 0, max_cv_variants = -1, seed = 999, cv_threads = 1,
                                           ld_reference_meta_file = NULL, region_info = NULL
                                         )) {
  if (pip_cutoff_to_skip > 0) {
    top_model_pip <- susie(X, Y, L = 1)$pip
    if (!any(top_model_pip > pip_cutoff_to_skip)) {
      message(paste("Skipping follow-up analysis: No signals above PIP threshold", pip_cutoff_to_skip, "in initial model screening."))
      return(list())
    } else {
      message(paste("Follow-up on region because signals above PIP threshold", pip_cutoff_to_skip, "were detected in initial model screening."))
    }
  }
  pri_coverage <- finemapping_opts$coverage[1]
  sec_coverage <- if (length(finemapping_opts$coverage) > 1) finemapping_opts$coverage[-1] else NULL
  st <- proc.time()
  if (finemapping) {
    res <- susie_wrapper(X, Y, init_L = finemapping_opts$init_L, max_L = finemapping_opts$max_L, l_step = finemapping_opts$l_step, refine = TRUE, coverage = pri_coverage)
    res <- susie_post_processor(res, X, Y, X_scalar, Y_scalar, maf,
      secondary_coverage = sec_coverage, signal_cutoff = finemapping_opts$signal_cutoff,
      other_quantities = other_quantities
    )
  } else {
    res <- list()
  }
  if (twas_weights) {
    twas_weights_output <- twas_weights_pipeline(X, Y, maf,
      susie_fit = res$susie_result_trimmed,
      ld_reference_meta_file = twas_weights_opts$ld_reference_meta_file,
      X_scalar = X_scalar, y_scalar = Y_scalar,
      cv_folds = twas_weights_opts$cv_folds, coverage = pri_coverage, secondary_coverage = sec_coverage, signal_cutoff = finemapping_opts$signal_cutoff,
      min_cv_maf = twas_weights_opts$min_cv_maf, max_cv_variants = twas_weights_opts$max_cv_variants, cv_seed = twas_weights_opts$seed, cv_threads = twas_weights_opts$cv_threads
    )
    # clean up the output database
    res <- c(res, twas_weights_output)
    res$twas_weights <- lapply(res$twas_weights, function(x) {
      rownames(x) <- NULL
      return(x)
    })
  }
  if (!is.null(region_info)) res$region_info <- region_info
  res$total_time_elapsed <- proc.time() - st
  return(res)
}

#' TWAS Weights Pipeline
#'
#' This function performs weights computation for Transcriptome-Wide Association Study (TWAS)
#' incorporating various steps such as filtering variants by linkage disequilibrium reference panel variants,
#' fitting models using SuSiE and other methods, and calculating TWAS weights and predictions.
#' Optionally, it can perform cross-validation for TWAS weights.
#'
#' @param X A matrix of genotype data where rows represent samples and columns represent genetic variants.
#' @param y A vector of phenotype measurements for each sample.
#' @param susie_fit An object returned by the SuSiE function, containing the SuSiE model fit.
#' @param maf A vector of minor allele frequencies for each variant in X.
#' @param ld_reference_meta_file An optional path to a file containing linkage disequilibrium reference data. If provided, variants in X are filtered based on this reference.
#' @param cv_folds The number of folds to use for cross-validation. Set to 0 to skip cross-validation.
#' @param X_scalar A scalar or vector to scale the genotype data. Defaults to 1 (no scaling).
#' @param y_scalar A scalar to scale the phenotype data. Defaults to 1 (no scaling).
#' @param coverage The coverage probability used in SuSiE for credible set construction. Defaults to 0.95.
#' @param secondary_coverage A vector of secondary coverage probabilities for credible set refinement. Defaults to c(0.7, 0.5).
#' @param signal_cutoff A threshold for determining significant signals in the SuSiE output. Defaults to 0.05.
#' @param mr_ash_max_iter The maximum number of iterations for the MR-ASH method. Defaults to 100.
#' @param min_cv_maf The minimum minor allele frequency for variants to be included in cross-validation. Defaults to 0.05.
#' @param max_cv_variants The maximum number of variants to be included in cross-validation. Defaults to -1 which means no limit.
#' @param cv_seed The seed for random number generation in cross-validation. Defaults to 999.
#' @param cv_threads The number of threads to use for parallel computation in cross-validation. Defaults to 1.
#' @return A list containing results from the TWAS pipeline, including TWAS weights, predictions, and optionally cross-validation results.
#' @export
#' @examples
#' # Example usage (assuming appropriate objects for X, y, susie_fit, and maf are available):
#' twas_results <- twas_weights_pipeline(X, y, maf, susie_fit)
twas_weights_pipeline <- function(X, y, maf, susie_fit, ld_reference_meta_file = NULL, cv_folds = 5, X_scalar = 1, y_scalar = 1,
                                  coverage = 0.95, secondary_coverage = c(0.7, 0.5), signal_cutoff = 0.05,
                                  mr_ash_max_iter = 100, min_cv_maf = 0.05, max_cv_variants = -1, cv_seed = 999, cv_threads = 1) {
  res <- list()
  if (!is.null(susie_fit)) {
    L <- length(which(susie_fit$V > 1E-9))
    init_L <- max(1, L - 2)
    max_L <- L + 3
  } else {
    # susie_fit did not detect anything significant
    init_L <- 2
    max_L <- 2
  }
  if (!is.null(ld_reference_meta_file)) {
    variants_kept <- filter_variants_by_ld_reference(colnames(X), ld_reference_meta_file)
    X <- X[, variants_kept$data, drop = FALSE]
    maf <- maf[variants_kept$idx]
    res$preset_variants_result <- susie_wrapper(X, y, init_L = init_L, max_L = max_L, refine = TRUE, coverage = coverage)
    res$preset_variants_result <- susie_post_processor(res$preset_variants_result, X, y, if (X_scalar == 1) 1 else X_scalar[variants_kept$idx], y_scalar, maf, secondary_coverage = secondary_coverage, signal_cutoff = signal_cutoff)
    res$preset_variants_result$analysis_script <- NULL
    res$preset_variants_result$sumstats <- NULL
    susie_fit <- res$preset_variants_result$susie_result_trimmed
  }
  weight_methods <- list(enet_weights = list(), lasso_weights = list(), mrash_weights = list(init_prior_sd = TRUE, max.iter = mr_ash_max_iter))
  if (!is.null(susie_fit)) weight_methods$susie_weights <- list(susie_fit = susie_fit)
  # get TWAS weights
  res$twas_weights <- twas_weights(X, y, weight_methods = weight_methods)
  # get TWAS predictions for possible next steps such as computing correlations between predicted expression values
  res$twas_predictions <- twas_predict(X, res$twas_weights)
  if (cv_folds > 1) {
    # A few cutting corners to run CV faster at the disadvantage of SuSiE and mr.ash:
    # 1. reset SuSiE to not using refine or adaptive L (more or less the default SuSiE)
    # 2. at most 100 iterations for mr.ash allowed
    # 3. only use a subset of top signals and common variants
    if (!is.null(susie_fit)) weight_methods$susie_weights <- list(refine = FALSE, init_L = max_L, max_L = max_L)
    variants_for_cv <- c()
    if (!is.null(res$preset_variants_result$top_loci) && nrow(res$preset_variants_result$top_loci) > 0) {
      variants_for_cv <- res$preset_variants_result$top_loci[, 1]
    }
    common_var <- colnames(X)[which(maf > min_cv_maf)]
    if (max_cv_variants < 0) max_cv_variants <- Inf
    if (length(common_var) + length(variants_for_cv) > max_cv_variants) {
      common_var <- sample(common_var, max_cv_variants - length(variants_for_cv), replace = FALSE)
    }
    variants_for_cv <- unique(c(variants_for_cv, common_var))
    res$twas_cv_result <- twas_weights_cv(X, y, fold = cv_folds, weight_methods = weight_methods, seed = cv_seed, max_num_variants = max_cv_variants, num_threads = cv_threads, variants_to_keep = if (length(variants_for_cv) > 0) variants_for_cv else NULL)
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
    sumstat_path, column_file_path, LD_data, n_sample = 0, n_case = 0, n_control = 0, skip_region = NULL,
    qc_method = c("rss_qc", "dentist", "slalom"),
    finemapping_method = c("susie_rss", "single_effect", "bayesian_conditional_regression"),
    finemapping_opts = list(
      init_L = 5, max_L = 20, l_step = 5,
      coverage = c(0.95, 0.7, 0.5), signal_cutoff = 0.025
    ),
    impute = TRUE, impute_opts = list(rcond = 0.01, R2_threshold = 0.6, minimum_ld = 5, lamb = 0.01),
    pip_cutoff_to_skip = 0) {
  rss_input <- load_rss_data(
    sumstat_path = sumstat_path, column_file_path = column_file_path,
    n_sample = n_sample, n_case = n_case, n_control = n_control
  )

  sumstats <- rss_input$sumstats
  n <- rss_input$n
  var_y <- rss_input$var_y
    
  # Preprocess the input data
  preprocess_results <- rss_basic_qc(sumstats, LD_data, skip_region = skip_region)
  sumstats <- preprocess_results$sumstats
  LD_mat <- preprocess_results$LD_mat

  if (pip_cutoff_to_skip > 0) {
    top_model <- susie_rss_wrapper(z = sumstats$z, R = LD_mat, L = 1, n = n, var_y = var_y)
    top_model_pip <- top_model$pip
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
  res <- list()
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
  }
  return(list(result = res, rss_data_analyzed = sumstats))
}