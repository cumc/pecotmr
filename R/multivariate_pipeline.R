#' Multivariate Analysis Pipeline
#'
#' This function performs weights computation for Transcriptome-Wide Association Study (TWAS) with fitting
#' models using mvSuSiE and mr.mash with the option of using a limited number of variants selected from
#' mvSuSiE fine-mapping for computing TWAS weights with cross-validation.
#'
#' @param X A matrix of genotype data where rows represent samples and columns represent genetic variants.
#' @param Y A matrix of phenotype measurements, representing samples and columns represent conditions.
#' @param maf A list of vectors for minor allele frequencies for each variant in X.
#' @param max_L The maximum number of components in mvSuSiE. Default is 30.
#' @param ld_reference_meta_file An optional path to a file containing linkage disequilibrium reference data. If provided, variants in X are filtered based on this reference.
#' @param pip_cutoff_to_skip Cutoff value for skipping conditions based on PIP values. Default is 0.
#' @param signal_cutoff Cutoff value for signal identification in PIP values for susie_post_processor. Default is 0.025.
#' @param coverage A vector of coverage probabilities, with the first element being the primary coverage and the rest being secondary coverage probabilities for credible set refinement. Defaults to c(0.95, 0.7, 0.5).
#' @param data_driven_prior_matrices A list of data-driven covariance matrices for mr.mash weights.
#' @param data_driven_prior_matrices_cv A list of data-driven covariance matrices for mr.mash weights in cross-validation.
#' @param canonical_prior_matrices If set to TRUE, will compute canonical covariance matrices and add them into the prior covariance matrix list in mrmash_wrapper. Default is TRUE.
#' @param sample_partition Sample partition for cross-validation.
#' @param mrmash_max_iter The maximum number of iterations for mr.mash. Default is 5000.
#' @param mvsusie_max_iter The maximum number of iterations for mvSuSiE. Default is 200.
#' @param min_cv_maf The minimum minor allele frequency for variants to be included in cross-validation. Default is 0.05.
#' @param max_cv_variants The maximum number of variants to be included in cross-validation. Defaults to -1 which means no limit.
#' @param cv_folds The number of folds to use for cross-validation. Set to 0 to skip cross-validation. Default is 5.
#' @param cv_threads The number of threads to use for parallel computation in cross-validation. Defaults to 1.
#' @param cv_seed The seed for random number generation in cross-validation. Defaults to 999.
#' @param prior_weights_min The minimum weight for prior covariance matrices. Default is 1e-4.
#' @param twas_weights If set to TRUE, computes TWAS weights. Default is FALSE.
#' @param verbose Verbosity level. Default is 0.
#'
#' @return A list containing the multivariate analysis results.
#' @examples
#' library(pecotmr)
#'
#' data(multitrait_data)
#' attach(multitrait_data)
#'
#' data_driven_prior_matrices <- list(
#'   U = prior_matrices,
#'   w = rep(1 / length(prior_matrices), length(prior_matrices))
#' )
#'
#' data_driven_prior_matrices_cv <- lapply(prior_matrices_cv, function(x) {
#'   list(U = x, w = rep(1 / length(x), length(x)))
#' })
#'
#' result <- multivariate_analysis_pipeline(
#'   X = multitrait_data$X[, 1:100],
#'   Y = multitrait_data$Y,
#'   maf = colMeans(multitrait_data$X[, 1:100]),
#'   max_L = 10,
#'   ld_reference_meta_file = NULL,
#'   max_cv_variants = -1,
#'   pip_cutoff_to_skip = 0.025,
#'   signal_cutoff = 0.025,
#'   data_driven_prior_matrices = data_driven_prior_matrices,
#'   data_driven_prior_matrices_cv = data_driven_prior_matrices_cv,
#'   canonical_prior_matrices = TRUE,
#'   sample_partition = NULL,
#'   cv_folds = 5,
#'   cv_seed = 999,
#'   cv_threads = 2,
#'   prior_weights_min = 1e-4,
#'   twas_weights = TRUE
#' )
#' @importFrom mvsusieR mvsusie create_mixture_prior
#' @export
multivariate_analysis_pipeline <- function(
    X, Y, maf, max_L = 30, ld_reference_meta_file = NULL, pip_cutoff_to_skip = 0, signal_cutoff = 0.025, coverage = c(0.95, 0.7, 0.5), data_driven_prior_matrices = NULL,
    data_driven_prior_matrices_cv = NULL, canonical_prior_matrices = TRUE, sample_partition = NULL,
    mrmash_max_iter = 5000, mvsusie_max_iter = 200, min_cv_maf = 0.05, max_cv_variants = -1, cv_folds = 5,
    cv_threads = 1, cv_seed = 999, prior_weights_min = 1e-4, twas_weights = FALSE, verbose = 0) {
  skip_conditions <- function(X, Y, pip_cutoff_to_skip) {
    if (length(pip_cutoff_to_skip) == 1 && is.numeric(pip_cutoff_to_skip)) {
      pip_cutoff_to_skip <- rep(pip_cutoff_to_skip, ncol(Y))
    } else if (length(pip_cutoff_to_skip) != ncol(Y)) {
      stop("pip_cutoff_to_skip must be a single number or a vector of the same length as ncol(Y).")
    }
    cols_to_keep <- logical(ncol(Y))
    for (r in 1:ncol(Y)) {
      if (pip_cutoff_to_skip[r] > 0) {
        non_missing_indices <- which(!is.na(Y[, r]))
        X_non_missing <- X[non_missing_indices, ]
        Y_non_missing <- Y[non_missing_indices, r]

        top_model_pip <- susie(X_non_missing, Y_non_missing, L = 1)$pip

        if (any(top_model_pip > pip_cutoff_to_skip[r])) {
          cols_to_keep[r] <- TRUE
        } else {
          message(paste0(
            "Skipping condition ", colnames(Y)[r], ", because all top_model_pip < pip_cutoff_to_skip = ",
            pip_cutoff_to_skip[r], ". Top loci model does not show any potentially significant variants."
          ))
        }
      }
    }

    Y_filtered <- Y[, cols_to_keep, drop = FALSE]

    if (ncol(Y_filtered) <= 1) {
      warning("After filtering, Y has ", ncol(Y_filtered), " column(s) left. Returning NULL.")
      return(NULL)
    } else {
      return(Y_filtered)
    }
  }

  initialize_multivariate_prior <- function(condition_names, data_driven_prior_matrices,
                                            data_driven_prior_matrices_cv, cv_folds, prior_weights_min) {
    if (!is.null(data_driven_prior_matrices)) {
      data_driven_prior_matrices <- list(matrices = data_driven_prior_matrices$U, weights = data_driven_prior_matrices$w)
      data_driven_prior_matrices <- create_mixture_prior(mixture_prior = data_driven_prior_matrices, weights_tol = prior_weights_min, include_indices = condition_names)
    }

    if (!is.null(data_driven_prior_matrices_cv)) {
      data_driven_prior_matrices_cv <- lapply(
        data_driven_prior_matrices_cv,
        function(x) {
          x <- list(matrices = x$U, weights = x$w)
          create_mixture_prior(mixture_prior = x, weights_tol = prior_weights_min, include_indices = condition_names)
        }
      )
    } else {
      if (!is.null(data_driven_prior_matrices)) {
        data_driven_prior_matrices_cv <- lapply(1:cv_folds, function(x) {
          return(data_driven_prior_matrices)
        })
      }
    }

    return(list(
      data_driven_prior_matrices = data_driven_prior_matrices, data_driven_prior_matrices_cv = data_driven_prior_matrices_cv
    ))
  }

  # Skip conditions based on PIP values
  Y <- skip_conditions(X, Y, pip_cutoff_to_skip)

  # Return empty list if all conditions are skipped
  if (is.null(Y)) {
    return(list())
  }

  # Filter data based on remaining conditions
  filtered_data <- initialize_multivariate_prior(colnames(Y), data_driven_prior_matrices,
    data_driven_prior_matrices_cv, cv_folds,
    prior_weights_min = prior_weights_min
  )

  data_driven_prior_matrices <- filtered_data$data_driven_prior_matrices
  data_driven_prior_matrices_cv <- filtered_data$data_driven_prior_matrices_cv

  if (twas_weights) {
    message("Fitting mr.mash model on input data ...")
    mrmash_fitted <- mrmash_wrapper(
      X = X, Y = Y, data_driven_prior_matrices = data_driven_prior_matrices,
      canonical_prior_matrices = canonical_prior_matrices, max_iter = mrmash_max_iter
    )
    resid_Y <- mrmash_fitted$V
  } else {
    resid_Y <- mr.mash.alpha:::compute_cov_flash(Y)
    mrmash_fitted <- NULL
  }

  pri_coverage <- coverage[1]
  sec_coverage <- if (length(coverage) > 1) coverage[-1] else NULL

  res <- list()
  message("Fitting mvSuSiE model on input data ...")
  mvsusie_fitted <- mvsusie(X,
    Y = Y, L = max_L, prior_variance = data_driven_prior_matrices,
    residual_variance = resid_Y, precompute_covariances = TRUE, compute_objective = TRUE,
    estimate_residual_variance = FALSE, estimate_prior_variance = TRUE, estimate_prior_method = "EM",
    max_iter = mvsusie_max_iter, n_thread = 1, approximate = FALSE, verbosity = verbose, coverage = pri_coverage
  )

  # Process mvSuSiE results
  res$mnm_result <- susie_post_processor(
    mvsusie_fitted, X, NULL, 1, 1,
    maf = maf, secondary_coverage = sec_coverage, signal_cutoff = signal_cutoff, mode = "mvsusie"
  )
  res$mnm_result$mrmash_result <- mrmash_fitted

  # Run TWAS pipeline
  if (twas_weights) {
    res <- twas_multivariate_weights_pipeline(X, Y, maf, res,
      resid_Y = resid_Y, cv_folds = cv_folds, sample_partition = sample_partition,
      ld_reference_meta_file = ld_reference_meta_file, max_cv_variants = max_cv_variants,
      mvsusie_max_iter = mvsusie_max_iter, mrmash_max_iter = mrmash_max_iter, signal_cutoff = signal_cutoff,
      coverage = pri_coverage, secondary_coverage = sec_coverage,
      canonical_prior_matrices = canonical_prior_matrices, data_driven_prior_matrices = data_driven_prior_matrices,
      data_driven_prior_matrices_cv = data_driven_prior_matrices_cv, cv_seed = cv_seed,
      min_cv_maf = min_cv_maf, cv_threads = cv_threads
    )
  }
  return(res)
}

#' TWAS Weights Multivariate Pipeline
#'
#' @importFrom mvsusieR mvsusie
#' @export
twas_multivariate_weights_pipeline <- function(
    X, Y, maf, res, ld_reference_meta_file = NULL, cv_folds = 5, sample_partition = NULL,
    data_driven_prior_matrices = NULL, data_driven_prior_matrices_cv = NULL, canonical_prior_matrices = FALSE, resid_Y,
    mvsusie_max_iter = 200, mrmash_max_iter = 5000,
    signal_cutoff = 0.05, coverage = 0.95, secondary_coverage = c(0.7, 0.5),
    min_cv_maf = 0.05, max_cv_variants = -1, cv_seed = 999, cv_threads = 1) {
  determine_max_L <- function(mvsusie_prefit) {
    if (!is.null(mvsusie_prefit)) {
      L <- length(which(mvsusie_prefit$V > 1E-9)) + 2
    } else {
      L <- 2
    }
    return(L)
  }

  copy_twas_results <- function(res, twas_weight, twas_predictions) {
    for (i in names(res)) {
      if (i == "mnm_result") next
      res[[i]]$twas_weights <- lapply(twas_weight, function(wgts) {
        wgts[, i]
      })
      res[[i]]$twas_predictions <- lapply(twas_predictions, function(pred) {
        pred[, i]
      })
    }
    return(res)
  }

  copy_twas_cv_results <- function(res, twas_cv_result) {
    for (i in names(res)) {
      if (i == "mnm_result") next
      res[[i]]$twas_cv_result$sample_partition <- twas_cv_result$sample_partition
      res[[i]]$twas_cv_result$prediction <- lapply(
        twas_cv_result$prediction,
        function(predicted) {
          predicted[, i]
        }
      )
      res[[i]]$twas_cv_result$performance <- lapply(
        twas_cv_result$performance,
        function(perform) {
          perform[i, ]
        }
      )
      res[[i]]$twas_cv_result$time_elapsed <- twas_cv_result$time_elapsed
    }
    return(res)
  }

  # Main analysis starts
  mvsusie_fitted <- res$mnm_result$susie_result_trimmed
  max_L <- determine_max_L(mvsusie_fitted)
  if (!is.null(ld_reference_meta_file)) {
    # Filter variants and update mvsusie_fitted
    variants_kept <- filter_variants_by_ld_reference(colnames(X), ld_reference_meta_file)
    X <- X[, variants_kept$data, drop = FALSE]
    maf <- maf[variants_kept$idx]
    message("Fitting mvSuSiE model on preset variants ...")
    res$mnm_result$preset_variants_result <- mvsusie(
      X = X, Y = Y, L = max_L, prior_variance = data_driven_prior_matrices,
      residual_variance = resid_Y, precompute_covariances = T, compute_objective = T,
      estimate_residual_variance = F, estimate_prior_variance = T, estimate_prior_method = "EM",
      max_iter = mvsusie_max_iter, n_thread = 1, approximate = F, verbosity = verbose, coverage = coverage
    )
    res$mnm_result$preset_variants_result <- susie_post_processor(
      res$mnm_result$preset_variants_result, X, NULL, 1, 1,
      maf = maf, secondary_coverage = secondary_coverage, signal_cutoff = signal_cutoff, mode = "mvsusie"
    )
    res$mnm_result$preset_variants_result$analysis_script <- NULL
    res$mnm_result$preset_variants_result$sumstats <- NULL
    mvsusie_fitted <- res$mnm_result$preset_variants_result$susie_result_trimmed
  }

  # Compute TWAS weights and predictions
  weight_methods <- list(
    mrmash_weights = list(
      mrmash_fit = res$mnm_result$mrmash_result,
      data_driven_prior_matrices = data_driven_prior_matrices,
      canonical_prior_matrices = canonical_prior_matrices, max_iter = mrmash_max_iter,
      verbose = verbose
    ),
    mvsusie_weights = list(
      mvsusie_fit = mvsusie_fitted,
      prior_variance = data_driven_prior_matrices,
      residual_variance = resid_Y, L = max_L, max_iter = mvsusie_max_iter,
      verbosity = verbose
    )
  )
  message("Computing TWAS weights for multivariate analysis methods ...")
  # get TWAS weights
  twas_weights_res <- twas_weights(X = X, Y = Y, weight_methods = weight_methods)
  # get TWAS predictions for possible next steps such as computing correlations between predicted expression values
  twas_predictions <- twas_predict(X, twas_weights_res)

  # copy TWAS results by condition
  res <- copy_twas_results(res, twas_weights_res, twas_predictions)

  # Perform cross-validation if specified
  if (cv_folds > 1) {
    weight_methods <- list(
      mrmash_weights = list(
        data_driven_prior_matrices = data_driven_prior_matrices,
        data_driven_prior_matrices_cv = data_driven_prior_matrices_cv,
        canonical_prior_matrices = canonical_prior_matrices,
        max_iter = mrmash_max_iter,
        verbose = verbose
      ),
      mvsusie_weights = list(
        prior_variance = data_driven_prior_matrices,
        data_driven_prior_matrices_cv = data_driven_prior_matrices_cv,
        residual_variance = resid_Y, L = max_L,
        max_iter = mvsusie_max_iter,
        verbosity = verbose
      )
    )
    variants_for_cv <- c()
    if (!is.null(res$mnm_result$preset_variants_result$top_loci) && nrow(res$mnm_result$preset_variants_result$top_loci) > 0) {
      variants_for_cv <- res$mnm_result$preset_variants_result$top_loci[, 1]
    }
    common_var <- colnames(X)[which(maf > min_cv_maf)]
    if (max_cv_variants < 0) max_cv_variants <- Inf
    if (length(common_var) + length(variants_for_cv) > max_cv_variants) {
      common_var <- sample(common_var, max_cv_variants - length(variants_for_cv), replace = FALSE)
    }
    variants_for_cv <- unique(c(variants_for_cv, common_var))
    message("Performing cross-validation to assess TWAS weights ...")
    twas_cv_result <- twas_weights_cv(
      X = X, Y = Y, fold = cv_folds,
      weight_methods = weight_methods,
      sample_partition = sample_partition,
      num_threads = cv_threads, seed = cv_seed,
      max_num_variants = max_cv_variants,
      variants_to_keep = if (length(variants_for_cv) > 0) variants_for_cv else NULL,
      data_driven_prior_matrices_cv = data_driven_prior_matrices_cv
    )
    res <- copy_twas_cv_results(res, twas_cv_result)
  }
  return(res)
}
