#' Univariate analysis pipeline
#' @importFrom susieR susie
#' @export
univariate_analysis_pipeline <- function(X, Y, X_scalar, Y_scalar, maf, other_quantities = list(),
                                         pip_cutoff_to_skip = 0, region_info = NULL,
                                         finemapping = TRUE, twas_weights = TRUE,
                                         finemapping_opts = list(
                                           init_L = 5, max_L = 20,
                                           coverage = c(0.95, 0.7, 0.5), signal_cutoff = 0.025
                                         ),
                                         twas_weights_opts = list(
                                           cv_folds = 5, min_cv_maf = 0, max_cv_variants = -1, seed = 999, cv_threads = 1,
                                           ld_reference_meta_file = NULL, region_info = NULL
                                         )) {
  if (pip_cutoff_to_skip > 0) {
    # return a NULL set if the top loci model does not show any potentially significant variants
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
    res <- susie_wrapper(X, Y, init_L = finemapping_opts$init_L, max_L = finemapping_opts$max_L, refine = TRUE, coverage = pri_coverage)
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
  if (cv_folds > 0) {
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
