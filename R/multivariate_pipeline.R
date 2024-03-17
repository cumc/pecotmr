#' Multivariate Pipeline
#'
#' This function performs weights computation for Transcriptome-Wide Association Study (TWAS) with fitting
#' models using mvSuSiE and mr.mash with the option of using a limited number of variants selected from
#' mvSuSiE fine-mapping for computing TWAS weights with cross-validation.
#'
#' @param X A matrix of genotype data where rows represent samples and columns represent genetic variants.
#' @param y A matrix of phenotype measurements, representing samples and columns represent conditions.
#' @param maf A list of vectors for minor allele frequencies for each variant in X.
#' @param dropped_sample A list of lists with X, Y, and covar, each list contains a list of dropped samples as vectors for each condition. Can be obtained from the output of load_regional_multivariate_data.
#' @param max_L The maximum number of components in mvSuSiE. Default is 30.
#' @param ld_reference_meta_file An optional path to a file containing linkage disequilibrium reference data. If provided, variants in X are filtered based on this reference.
#' @param X_scalar Scalar for the genotype data, used in residual scaling at the step of mvsusie_post_processor. Defaults to 1 (no scaling).
#' @param y_scalar Scalar for the phenotype data, used in residual scaling at the step of mvsusie_post_processor. Defaults to 1 (no scaling).
#' @param pip_cutoff_to_skip Cutoff value for skipping conditions based on PIP values. Default is 0.
#' @param signal_cutoff Cutoff value for signal identification in PIP values for mvsusie_post_processor. Default is 0.025.
#' @param secondary_coverage A vector of secondary coverage probabilities for credible set refinement. Defaults to c(0.7, 0.5).
#' @param mrmash_weights_prior_matrices A list of data-driven covariance matrices for mr.mash weights.
#' @param mrmash_weights_prior_matrices_cv A list of data-driven covariance matrices for mr.mash weights in cross-validation.
#' @param prior_canonical_matrices If set to TRUE, will compute canonical covariance matrices and add them into the prior covariance matrix list in mrmash_wrapper. Default is TRUE.
#' @param sample_partition Sample partition for cross-validation.
#' @param mrmash_max_iter The maximum number of iterations for mr.mash. Default is 5000.
#' @param mvsusie_max_iter The maximum number of iterations for mvSuSiE. Default is 200.
#' @param max_cv_variants The maximum number of variants to be included in cross-validation. Defaults to 5000.
#' @param cv_folds The number of folds to use for cross-validation. Set to 0 to skip cross-validation. Default is 5.
#' @param cv_threads The number of threads to use for parallel computation in cross-validation. Defaults to 1.
#' @param cv_seed The seed for random number generation in cross-validation. Defaults to 999.
#'
#' @importFrom mvsusieR mvsusie create_mixture_prior
#' @export
run_multivariate_pipeline <- function(
    X, y, maf, max_L = 30, ld_reference_meta_file = NULL, pip_cutoff_to_skip = rep(
      0,
      ncol(y)
    ), signal_cutoff = 0.025, secondary_coverage = c(0.5, 0.7), mrmash_weights_prior_matrices = NULL,
    mrmash_weights_prior_matrices_cv = NULL, prior_canonical_matrices = TRUE, sample_partition = NULL,
    mrmash_max_iter = 5000, mvsusie_max_iter = 200, max_cv_variants = 5000, cv_folds = 5,
    cv_threads = 1, cv_seed = 999, weights_tol = 1e-3, twas_weights = FALSE) {

  skip_conditions <- function(y, pip_cutoff_to_skip) {
    for (r in 1:ncol(y)) {
      if (pip_cutoff_to_skip[r] > 0) {
        top_model_pip <- susie(X, y[, r], L = 1)$pip
        if (!any(top_model_pip > pip_cutoff_to_skip[r])) {
          message(paste0(
            "Skipping condition ", colnames(y)[r], ", because all top_model_pip < pip_cutoff_to_skip = ",
            pip_cutoff_to_skip[r], ", top loci model does not show any potentially significant variants."
          ))
          y[, r] <- NA
        }
      }
    }
    return(y)
  }

   filter_prior_matrices <- function(
     condition_names , mrmash_weights_prior_matrices,
      mrmash_weights_prior_matrices_cv, weights_tol = 10^-3) {
    
      if (!is.null(mrmash_weights_prior_matrices)) {
         names(mrmash_weights_prior_matrices)[1:2] <- c("matrices","weights")  
         mrmash_weights_prior_matrices <- create_mixture_prior(mixture_prior = mrmash_weights_prior_matrices, weights_tol = weights_tol, include_indices = condition_names)
      }
    
      if (!is.null(mrmash_weights_prior_matrices_cv)) {
        mrmash_weights_prior_matrices_cv <- lapply(
          mrmash_weights_prior_matrices_cv,
          function(x) {
            names(x)[1:2] <- c("matrices","weights")
            create_mixture_prior(mixture_prior = x, weights_tol = weights_tol, include_indices = condition_names)
          }
        )
      }

    return(list(
      mrmash_weights_prior_matrices = mrmash_weights_prior_matrices, mrmash_weights_prior_matrices_cv = mrmash_weights_prior_matrices_cv
    ))
  }

  # Skip conditions based on PIP values
  y <- skip_conditions(y, pip_cutoff_to_skip)

  # Return empty list if all conditions are skipped
  if (all(is.na(y))) {
    return(list())
  }

  # Filter data based on remaining conditions
  filtered_data <- filter_prior_matrices(colnames(y), mrmash_weights_prior_matrices,
    mrmash_weights_prior_matrices_cv
  )

  filtered_mrmash_weights_prior_matrices <- filtered_data$mrmash_weights_prior_matrices
  filtered_mrmash_weights_prior_matrices_cv <- filtered_data$mrmash_weights_prior_matrices_cv


  # Fit mr.mash model
  mrmash_output <- mrmash_wrapper(
    X = X, Y = y, prior_data_driven = filtered_mrmash_weights_prior_matrices$prior_variance,
    prior_grid = NULL, prior_canonical_matrices = prior_canonical_matrices, max_iter = mrmash_max_iter
  )
  resid_Y <- mrmash_output$V

  # Fit mvSuSiE model
  mvsusie_fitted <- mvsusie(X,
    Y = y, L = max_L, prior_variance = filtered_mrmash_weights_prior_matrices, prior_weights = filtered_mrmash_weights_prior_matrices$prior_variance$weights,
    residual_variance = resid_Y, precompute_covariances = F, compute_objective = T, estimate_residual_variance = F,
    estimate_prior_variance = T, estimate_prior_method = "EM", max_iter = mvsusie_max_iter,
    n_thread = 1, approximate = F
  )
  return(mvsusie_fitted)
  # Process mvSuSiE results
  # mvsusie_prefit <- mvsusie_post_processor(
  #   mvsusie_fitted, X, y, X_scalar, y_scalar,
  #   maf, secondary_coverage, signal_cutoff, dropped_samples
  # )
    
  # Run TWAS pipeline
  if (!twas_weights) {
     result <- twas_multivariate_weights_pipeline(X, y, maf, mvsusie_prefit, mvsusie_fitted,
               prior = prior, resid_Y = resid_Y, dropped_samples = dropped_samples, cv_folds = cv_folds,
               ld_reference_meta_file = ld_reference_meta_file, max_cv_variants = max_cv_variants,
               mvsusie_max_iter = mvsusie_max_iter, mrmash_max_iter = mrmash_max_iter, signal_cutoff = signal_cutoff,
               pip_cutoff_to_skip = pip_cutoff_to_skip, secondary_coverage = secondary_coverage,
               prior_canonical_matrices = prior_canonical_matrices, mrmash_weights_prior_matrices = mrmash_weights_prior_matrices,
               mrmash_weights_prior_matrices_cv = mrmash_weights_prior_matrices_cv, cv_seed = cv_seed
              )
     return(result)
   }
}

#' TWAS Weights Multivariate Pipeline
#'
#' @importFrom mvsusieR mvsusie
#' @export
twas_multivariate_weights_pipeline <- function(
    X, y, maf, mvsusie_prefit, mvsusie_fitted,
    dropped_samples, prior, resid_Y, X_scalar = rep(1, ncol(y)), y_scalar = rep(
      1,
      ncol(y)
    ), max_cv_variants = 5000, mvsusie_max_iter = 200, mrmash_max_iter = 5000,
    pip_cutoff_to_skip = 0, signal_cutoff = 0.025, secondary_coverage = c(0.5, 0.7),
    ld_reference_meta_file = NULL, cv_folds = 5, sample_partition = NULL, mrmash_weights_prior_matrices = NULL,
    mrmash_weights_prior_matrices_cv = NULL, prior_canonical_matrices = FALSE, cv_seed = 999,
    cv_threads = 1) {
  run_twas_cv <- function(res, X, y, top_sig_idx, cv_folds, sample_partition, mrmash_weights_prior_matrices_cv,
                          mrmash_max_iter, prior_canonical_matrices, mvsusie_fitted, resid_Y, max_L,
                          mvsusie_max_iter, cv_threads, cv_seed) {
    weight_methods <- list(mrmash_weights = list(), mvsusie_weights = list())
    message(paste0(
      "Performing cross-validation with ", length(top_sig_idx),
      " variants."
    ))

    twas_cv_result <- twas_weights_cv(
      X = X[, top_sig_idx], Y = y, fold = cv_folds,
      weight_methods = weight_methods, sample_partition = sample_partition,
      mrmash_weights_prior_matrices = mrmash_weights_prior_matrices_cv, mrmash_max_iter = mrmash_max_iter,
      prior_canonical_matrices = prior_canonical_matrices, mvsusie_fit = mvsusie_fitted,
      residual_variance = resid_Y, L = max_L, mvsusie_max_iter = mvsusie_max_iter,
      num_threads = cv_threads, seed = cv_seed
    )

    for (i in names(res)) {
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

  select_top_variants <- function(mvsusie_fitted, max_cv_variants) {
    if (max_cv_variants > length(mvsusie_fitted$pip)) {
      max_cv_variants <- length(mvsusie_fitted$pip)
    }
    top_sig_idx <- which(mvsusie_fitted$pip %in% mvsusie_fitted$pip[order(-mvsusie_fitted$pip)[1:max_cv_variants]])
    mvsusie_fitted$coef <- mvsusie_fitted$coef[c(1, top_sig_idx + 1), ]
    return(top_sig_idx)
  }

  split_twas_results <- function(res, twas_weight, twas_predictions) {
    for (i in names(res)) {
      res[[i]]$twas_weights <- lapply(twas_weight, function(wgts) {
        wgts[, i]
      })
      res[[i]]$twas_predictions <- lapply(twas_predictions, function(pred) {
        pred[, i]
      })
    }
    return(res)
  }

  compute_twas_weights <- function(X, y, mrmash_weights_prior_matrices, prior_canonical_matrices,
                                   mrmash_max_iter, mvsusie_fitted, prior, resid_Y, max_L, mvsusie_max_iter) {
    weight_methods <- list(
      mrmash_weights = list(
        prior_data_driven_matrices = mrmash_weights_prior_matrices,
        prior_canonical_matrices = prior_canonical_matrices, max_iter = mrmash_max_iter
      ),
      mvsusie_weights = list(
        mvsusie_fit = mvsusie_fitted, prior_variance = prior,
        residual_variance = resid_Y, L = max_L, max_iter = mvsusie_max_iter
      )
    )
    twas_weight <- twas_weights(X = X, Y = y, weight_methods = weight_methods)
    return(twas_weight)
  }

  mvsusie_preset_variants <- function(X, y, maf, ld_reference_meta_file, max_L,
                                      prior, resid_Y, mvsusie_max_iter) {
    variants_kept <- filter_variants_by_ld_reference(colnames(X), ld_reference_meta_file)
    X <- X[, variants_kept$data, drop = FALSE]
    maf <- lapply(maf, function(x, idx) {
      x[idx]
    }, idx = variants_kept$idx)

    mvsusie_fitted <- mvsusie(
      X = X, Y = y, L = max_L, prior_variance = prior,
      residual_variance = resid_Y, precompute_covariances = F, compute_objective = T,
      estimate_residual_variance = F, estimate_prior_variance = T, estimate_prior_method = "EM",
      max_iter = mvsusie_max_iter, n_thread = 1, approximate = F
    )

    return(list(X = X, maf = maf, mvsusie_fitted = mvsusie_fitted))
  }

  determine_max_L <- function(mvsusie_prefit) {
    max_L <- c()
    for (r in 1:length(mvsusie_prefit)) {
      if (!is.null(mvsusie_prefit[[r]])) {
        L <- length(which(mvsusie_prefit[[r]]$V > 1e-09))
        max_L <- c(max_L, L + 3)
      } else {
        max_L <- c(max_L, 2)
      }
    }
    max_L <- unique(max_L)
    if (length(max_L) >= 2) {
      max_L <- max(max_L)
    }
    return(max_L)
  }

  max_L <- determine_max_L(mvsusie_prefit)

  res <- list()

  # Filter variants
  if (!is.null(ld_reference_meta_file)) {
    filtered_variants <- mvsusie_preset_variants(
      X, y, maf, ld_reference_meta_file,
      max_L, prior, resid_Y, mvsusie_max_iter
    )
    X <- filtered_variants$X
    maf <- filtered_variants$maf
    mvsusie_fitted <- filtered_variants$mvsusie_fitted
  }

  # Process mvSuSiE results with filtered variants
  res <- mvsusie_post_processor(
    mvsusie_fitted, X, y, X_scalar, y_scalar, maf,
    secondary_coverage, signal_cutoff, dropped_samples
  )

  # Compute TWAS weights and predictions
  twas_weight <- compute_twas_weights(
    X, y, mrmash_weights_prior_matrices, prior_canonical_matrices,
    mrmash_max_iter, mvsusie_fitted, prior, resid_Y, max_L, mvsusie_max_iter
  )
  twas_predictions <- twas_predict(X, twas_weight)

  # Split TWAS results by condition
  res <- split_twas_results(res, twas_weight, twas_predictions)

  # Perform cross-validation if specified
  if (cv_folds > 0) {
    top_sig_idx <- select_top_variants(mvsusie_fitted, max_cv_variants)
    res <- run_twas_cv(
      res, X, y, top_sig_idx, cv_folds, sample_partition, mrmash_weights_prior_matrices_cv,
      mrmash_max_iter, prior_canonical_matrices, mvsusie_fitted, resid_Y, max_L,
      mvsusie_max_iter, cv_threads, cv_seed
    )
  }

  return(res)
}
