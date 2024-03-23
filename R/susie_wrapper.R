#' Convert Log Bayes Factors to Single Effects PIP
#'
#' This function converts log Bayes factors (LBF) to alpha values, optionally
#' using prior weights. It handles numerical stability by adjusting with the
#' maximum LBF value.
#'
#' @param lbf Numeric vector of log Bayes factors.
#' @param prior_weights Optional numeric vector of prior weights for each element in lbf.
#' @return A named numeric vector of alpha values corresponding to the input LBF.
#' @examples
#' lbf <- c(-0.5, 1.2, 0.3)
#' alpha <- lbf_to_alpha_vector(lbf)
#' print(alpha)
lbf_to_alpha_vector <- function(lbf, prior_weights = NULL) {
  if (is.null(prior_weights)) prior_weights <- rep(1 / length(lbf), length(lbf))
  maxlbf <- max(lbf)

  # If maxlbf is 0, return a vector of zeros
  if (maxlbf == 0) {
    return(setNames(rep(0, length(lbf)), names(lbf)))
  }

  # w is proportional to BF, subtract max for numerical stability
  w <- exp(lbf - maxlbf)

  # Posterior prob for each SNP
  w_weighted <- w * prior_weights
  weighted_sum_w <- sum(w_weighted)
  alpha <- w_weighted / weighted_sum_w

  return(alpha)
}

#' Applies the 'lbf_to_alpha_vector' function row-wise to a matrix of log Bayes factors
#' to convert them to Single Effect PIP values.
#'
#' @param lbf Matrix of log Bayes factors.
#' @return A matrix of alpha values with the same dimensions as the input LBF matrix.
#' @examples
#' lbf_matrix <- matrix(c(-0.5, 1.2, 0.3, 0.7, -1.1, 0.4), nrow = 2)
#' alpha_matrix <- lbf_to_alpha(lbf_matrix)
#' print(alpha_matrix)
#' @export
lbf_to_alpha <- function(lbf) t(apply(lbf, 1, lbf_to_alpha_vector))

# FIXME: this function does not work properly. Need to fix it in multiple aspects
#' Adjust SuSiE Weights
#'
#' This function adjusts the SuSiE weights based on a set of intersected variants.
#' It subsets various components like lbf_matrix, mu, and scale factors based on these variants.
#'
#' @param weight_db_file A RDS file containing TWAS weights.
#' @param condition specific condition.
#' @param keep_variants Vector of variant names to keep.
#' @param allele_qc Optional
#' @return A list of adjusted xQTL coefficients and remained variants ids
#' @export
adjust_susie_weights <- function(twas_weights_results, condition, keep_variants, allele_qc = TRUE) {
  # Intersect the rownames of weights with keep_variants
  twas_weights_variants <- get_nested_element(twas_weights_results, c("susie_results", condition, "variant_names"))
  # allele flip twas weights matrix variants name
  if (allele_qc) {
    weights_matrix <- get_nested_element(twas_weights_results, c("weights", condition))
    weights_matrix_qced <- allele_qc(twas_weights_variants, gwas_LD_list$combined_LD_variants, weights_matrix, 1:ncol(weights_matrix))
    intersected_indices <- which(weights_matrix_qced$qc_summary$keep == TRUE)
  } else {
    keep_variants_transformed <- ifelse(!startsWith(keep_variants, "chr"), paste0("chr", keep_variants), keep_variants)
    intersected_variants <- intersect(twas_weights_variants, keep_variants_transformed)
    intersected_indices <- match(intersected_variants, twas_weights_variants)
  }
  if (length(intersected_indices) == 0) {
    stop("Error: No intersected variants found. Please check 'twas_weights' and 'keep_variants' inputs to make sure there are variants left to use.")
  }
  # Subset lbf_matrix, mu, and x_column_scale_factors
  lbf_matrix <- get_nested_element(twas_weights_results, c("susie_results", condition, "susie_result_trimmed", "lbf_variable"))
  mu <- get_nested_element(twas_weights_results, c("susie_results", condition, "susie_result_trimmed", "mu"))
  x_column_scal_factors <- get_nested_element(twas_weights_results, c("susie_results", condition, "susie_result_trimmed", "X_column_scale_factors"))

  lbf_matrix_subset <- lbf_matrix[, intersected_indices]
  mu_subset <- mu[, intersected_indices]
  x_column_scal_factors_subset <- x_column_scal_factors[intersected_indices]

  # Convert lbf_matrix to alpha and calculate adjusted xQTL coefficients
  adjusted_xqtl_alpha <- lbf_to_alpha(lbf_matrix_subset)
  adjusted_xqtl_coef <- colSums(adjusted_xqtl_alpha * mu_subset) / x_column_scal_factors_subset

  return(list(adjusted_susie_weights = adjusted_xqtl_coef, remained_variants_ids = names(adjusted_xqtl_coef)))
}

#' @importFrom susieR susie
#' @export
susie_wrapper <- function(X, y, init_L = 10, max_L = 30, l_step = 5, ...) {
  if (init_L == max_L) {
    return(susie(X, y, L = init_L, median_abs_corr = 0.8, ...))
  }
  L <- init_L
  # Perform SuSiE by dynamically increasing L
  gst <- proc.time()
  while (TRUE) {
    st <- proc.time()
    res <- susie(X, y,
      L = L,
      median_abs_corr = 0.8, ...
    )
    res$time_elapsed <- proc.time() - st
    if (!is.null(res$sets$cs)) {
      if (length(res$sets$cs) >= L && L <= max_L) {
        L <- L + l_step
      } else {
        break
      }
    } else {
      break
    }
  }
  message(paste("Total time elapsed for susie_wrapper:", (proc.time() - gst)[3]))
  return(res)
}

#' Wrapper Function for SuSiE RSS with Dynamic L Adjustment
#'
#' This function performs SuSiE RSS analysis, dynamically adjusting the number of causal configurations (L)
#' and applying quality control and imputation as necessary. It includes the total phenotypic variance `var_y`
#' as one of its parameters to align with the `susie_rss` function's interface.
#'
#' @param z Z score vector.
#' @param R LD matrix.
#' @param bhat Vector of effect size estimates.
#' @param shat Vector of standard errors for effect size estimates.
#' @param var_y Total phenotypic variance.
#' @param n Sample size; if NULL, certain functionalities that require sample size will be skipped.
#' @param L Initial number of causal configurations to consider.
#' @param max_L Maximum number of causal configurations to consider.
#' @param l_step Step size for increasing L when the limit is reached.
#' @param zR_discrepancy_correction Logical indicating if z-score and R matrix discrepancy correction should be performed.
#' @param ... Extra parameters to pass to the susie_rss function.
#' @return SuSiE RSS fit object after dynamic L adjustment
#' @importFrom susieR susie_rss
#' @export
susie_rss_wrapper <- function(z, R, bhat, shat, n = NULL, var_y = NULL, L = 10, max_L = 30, l_step = 5,
                              zR_discrepancy_correction = FALSE, coverage = 0.95, ...) {
  if (L == 1) {
    return(susie_rss(
      z = z, R = R, bhat = bhat, shat = shat, var_y = var_y, n = n,
      L = 1, max_iter = 1, median_abs_corr = 0.8, correct_zR_discrepancy = FALSE, coverage = coverage, ...
    ))
  }
  if (L == max_L) {
    return(susie_rss(
      z = z, R = R, bhat = bhat, shat = shat, var_y = var_y, n = n, L = L, median_abs_corr = 0.8,
      correct_zR_discrepancy = zR_discrepancy_correction, coverage = coverage, ...
    ))
  }
  while (TRUE) {
    st <- proc.time()
    susie_rss_result <- susie_rss(
      z = z, R = R, bhat = bhat, shat = shat, var_y = var_y, n = n, L = L, median_abs_corr = 0.8,
      correct_zR_discrepancy = zR_discrepancy_correction, coverage = coverage, ...
    )
    susie_rss_result$time_elapsed <- proc.time() - st
    # Check for convergence and adjust L if necessary
    if (!is.null(susie_rss_result$sets$cs)) {
      if (length(susie_rss_result$sets$cs) >= L && L <= max_L) {
        L <- L + l_step # Increase L for the next iteration
      } else {
        break
      }
    } else {
      break # Break the loop if no credible sets are found
    }
  }

  return(susie_rss_result)
}

#' Run the SuSiE RSS pipeline
#'
#' This function runs the SuSiE RSS pipeline, including single-effect regression, no QC analysis,
#' and optional QC and Bayesian conditional analysis. It processes the input summary statistics and LD data
#' to perform various SuSiE analyses, providing results in a structured output.
#'
#' @param sumstat A list or data frame containing summary statistics with necessary columns.
#' @param R The LD matrix.
#' @param ref_panel Reference panel for QC and imputation.
#' @param n Sample size.
#' @param L Initial number of causal configurations to consider in the analysis.
#' @param var_y Variance of Y.
#' @param QC Perform quality control (default: TRUE).
#' @param impute Perform imputation (default: TRUE).
#' @param bayesian_conditional_analysis Perform Bayesian conditional analysis (default: TRUE).
#' @param lamb Regularization parameter for the RAiSS imputation method.
#' @param rcond Condition number for the RAiSS imputation method.
#' @param R2_threshold R-squared threshold for the RAiSS imputation method.
#' @param max_L Maximum number of components for QC.
#' @param l_step Step size for increasing L when the limit is reached during dynamic adjustment.
#' @param minimum_ld Minimum LD for QC.
#' @param coverage Coverage level for susie_rss analysis (default: 0.95).
#' @param secondary_coverage Secondary coverage levels for susie_rss analysis (default: c(0.7, 0.5)).
#' @param pip_cutoff_to_skip PIP cutoff to skip imputation (default: 0.025).
#' @param signal_cutoff Signal cutoff for susie_post_processor (default: 0.1).
#'
#' @return A list containing the results of various SuSiE RSS analyses.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr arrange select
#' @export
susie_rss_pipeline <- function(sumstat, R, ref_panel, n, L, var_y, QC = TRUE, impute = TRUE, bayesian_conditional_analysis = TRUE, lamb = 0.01, rcond = 0.01, R2_threshold = 0.6,
                               max_L = 20, l_step = 5, minimum_ld = 5, coverage = 0.95,
                               secondary_coverage = c(0.7, 0.5), pip_cutoff_to_skip = 0.025, signal_cutoff = 0.1) {
  if (!is.null(sumstat$z)) {
    z <- sumstat$z
    bhat <- NULL
    shat <- NULL
  } else if ((!is.null(sumstat$beta)) && (!is.null(sumstat$se))) {
    z <- sumstat$beta / sumstat$se
    bhat <- NULL
    shat <- NULL
  } else {
    stop("Sumstat should have z or (bhat and shat)")
  }

  final_result <- list()

  LD_extract <- R[sumstat$variant_id, sumstat$variant_id, drop = FALSE]
  single_effect_res <- susie_rss_wrapper(z = z, R = LD_extract, bhat = bhat, shat = shat, L = 1, n = n, var_y = var_y, coverage = coverage)
  if (max(single_effect_res$pip) < pip_cutoff_to_skip) {
    cat(paste0("No PIP larger than ", pip_cutoff_to_skip, " in this region."))
    if (impute) {
      z <- sumstat$z
      known_zscores <- sumstat %>% arrange(pos)
      final_result$sumstats_qc_impute <- raiss(ref_panel, known_zscores, R, lamb = lamb, rcond = rcond, R2_threshold = R2_threshold, minimum_ld = minimum_ld)$result_nofilter
      final_result$sumstats_qc_impute$chrom <- as.numeric(final_result$sumstats_qc_impute$chrom)
    }
  } else {
    single_effect_post <- susie_post_processor(single_effect_res, data_x = LD_extract, data_y = list(z = z), signal_cutoff = signal_cutoff, secondary_coverage = secondary_coverage, mode = "susie_rss")
    final_result$single_effect_regression <- single_effect_post
    result_noqc <- susie_rss_wrapper(z = z, R = LD_extract, bhat = bhat, shat = shat, n = n, L = L, var_y = var_y, coverage = coverage)
    result_noqc_post <- susie_post_processor(result_noqc, data_x = LD_extract, data_y = list(z = z), signal_cutoff = signal_cutoff, secondary_coverage = secondary_coverage, mode = "susie_rss")
    final_result$noqc <- result_noqc_post

    if (QC) {
      result_qced <- susie_rss_qc(sumstat, ref_panel = ref_panel, R = R, n = n, L = L, impute = impute, lamb = lamb, rcond = rcond, R2_threshold = R2_threshold, max_L = max_L, minimum_ld = minimum_ld, l_step = l_step, var_y = var_y, coverage = coverage)
      var_impute_kept <- names(result_qced$qc_impute_result$pip)
      result_qced_impute_post <- susie_post_processor(result_qced$qc_impute_result, data_x = R[var_impute_kept, var_impute_kept, drop = FALSE], data_y = list(z = result_qced$qc_impute_result$z), signal_cutoff = signal_cutoff, secondary_coverage = secondary_coverage, mode = "susie_rss")
      result_qced_only_post <- susie_post_processor(result_qced$qc_only_result, data_x = LD_extract, data_y = list(z = z), signal_cutoff = signal_cutoff, secondary_coverage = secondary_coverage, mode = "susie_rss")

      final_result$qc_impute <- result_qced_impute_post
      final_result$qc_only <- result_qced_only_post
      final_result$qc_only$outlier <- result_qced$qc_only_result$zR_outliers

      result_qced$sumstats_qc_impute_filtered$chrom <- as.numeric(result_qced$sumstats_qc_impute_filtered$chrom)
      result_qced$sumstats_qc_impute$chrom <- as.numeric(result_qced$sumstats_qc_impute$chrom)
      final_result$sumstats_qc_impute_filtered <- result_qced$sumstats_qc_impute_filtered %>% select(-beta, -se)
      final_result$sumstats_qc_impute <- result_qced$sumstats_qc_impute %>% select(-beta, -se)
    }

    if (bayesian_conditional_analysis) {
      conditional_noqc <- susie_rss_wrapper(z = z, R = LD_extract, bhat = bhat, shat = shat, n = n, L = L, max_iter = 1, var_y = var_y, coverage = coverage)
      conditional_noqc_post <- susie_post_processor(conditional_noqc, data_x = LD_extract, data_y = list(z = z), signal_cutoff = signal_cutoff, secondary_coverage = secondary_coverage, mode = "susie_rss")
      final_result$conditional_regression_noqc <- conditional_noqc_post

      if (QC) {
        outlier <- result_qced$qc_only_result$zR_outliers
        if (!is.null(outlier) & length(outlier) != 0) {
          z_rmoutlier <- z[-outlier]
          LD_rmoutlier <- LD_extract[-outlier, -outlier, drop = FALSE]
        } else {
          z_rmoutlier <- z
          LD_rmoutlier <- LD_extract
        }
        conditional_qc_only <- susie_rss_wrapper(z = z_rmoutlier, R = LD_rmoutlier, bhat = bhat, shat = shat, n = n, L = L, max_iter = 1, var_y = var_y, coverage = coverage)
        conditional_qc_only_post <- susie_post_processor(conditional_qc_only, data_x = LD_rmoutlier, data_y = list(z = z_rmoutlier), signal_cutoff = signal_cutoff, secondary_coverage = secondary_coverage, mode = "susie_rss")

        if (impute) {
          conditional_qced_impute <- susie_rss_wrapper(z = result_qced$qc_impute_result$z, R = R[var_impute_kept, var_impute_kept, drop = FALSE], bhat = bhat, shat = shat, n = n, L = L, max_iter = 1, var_y = var_y, coverage = coverage)
          conditional_qced_impute_post <- susie_post_processor(conditional_qced_impute, data_x = R[var_impute_kept, var_impute_kept, drop = FALSE], data_y = list(z = result_qced$qc_impute_result$z), signal_cutoff = signal_cutoff, secondary_coverage = secondary_coverage, mode = "susie_rss")
          final_result$conditional_regression_qc_impute <- conditional_qced_impute_post
        }

        final_result$conditional_regression_qc_only <- conditional_qc_only_post
      }
    }
  }

  return(final_result)
}

#' Post-process SuSiE Analysis Results
#'
#' This function processes the results from SuSiE (Sum of Single Effects) genetic analysis.
#' It extracts and processes various statistics and indices based on the provided SuSiE object and other parameters.
#' The function can operate in 3 modes: 'susie', 'susie_rss', 'mvsusie', based on the method used for the SuSiE analysis.
#'
#' @param susie_output Output from running susieR::susie() or susieR::susie_rss() or mvsusieR::mvsusie()
#' @param data_x Genotype data matrix for 'susie' or Xcorr matrix for 'susie_rss'.
#' @param data_y Phenotype data vector for 'susie' or summary stats object for 'susie_rss' (a list contain attribute betahat and sebetahat AND/OR z). i.e. data_y = list(betahat = ..., sebetahat = ...), or NULL for mvsusie
#' @param X_scalar Scalar for the genotype data, used in residual scaling.
#' @param y_scalar Scalar for the phenotype data, used in residual scaling.
#' @param maf Minor Allele Frequencies vector.
#' @param secondary_coverage Vector of coverage thresholds for secondary conditional analysis.
#' @param signal_cutoff Cutoff value for signal identification in PIP values.
#' @param other_quantities A list of other quantities to be added to the final object.
#' @param prior_eff_tol Prior effective tolerance.
#' @param mode Specify the analysis mode: 'susie' or 'susie_rss'.
#' @return A list containing modified SuSiE object along with additional post-processing information.
#' @examples
#' # Example usage for SuSiE
#' # result <- susie_post_processor(susie_output, X_data, y_data, maf, mode = "susie")
#' # Example usage for SuSiE RSS
#' # result <- susie_post_processor(susie_output, Xcorr, z, maf, mode = "susie_rss")
#' @importFrom dplyr full_join
#' @importFrom purrr map_int pmap
#' @importFrom susieR get_cs_correlation susie_get_cs
#' @importFrom stringr str_replace
#' @export
susie_post_processor <- function(susie_output, data_x, data_y, X_scalar, y_scalar, maf = NULL,
                                 secondary_coverage = c(0.5, 0.7), signal_cutoff = 0.1,
                                 other_quantities = NULL, prior_eff_tol = 1e-9, min_abs_corr = 0.5,
                                 median_abs_corr = 0.8,
                                 mode = c("susie", "susie_rss", "mvsusie")) {
  mode <- match.arg(mode)
  get_cs_index <- function(snps_idx, susie_cs) {
    # Use pmap to iterate over each vector in susie_cs
    idx_lengths <- tryCatch(
      {
        pmap(list(x = susie_cs), function(x) {
          # Check if snps_idx is in the CS and return the length of the CS if it is
          if (snps_idx %in% x) {
            return(length(x))
          } else {
            return(NA_integer_)
          }
        }) %>% unlist()
      },
      error = function(e) NA_integer_
    )
    idx <- which(!is.na(idx_lengths))
    # idx length should be either 1 or 0
    # But in some rare cases there will be a convergence issue resulting in a variant belong to multiple CS
    # In which case we will keep one of them, with a warning
    if (length(idx) > 0) {
      if (length(idx) > 1) {
        smallest_cs_idx <- which.min(idx_lengths[idx])
        selected_cs <- idx[smallest_cs_idx]
        selected_length <- idx_lengths[selected_cs]

        warning(sprintf(
          "Variable %d found in multiple CS: %s. Keeping smallest: CS %d (length %d).",
          snps_idx, paste(idx, collapse = ", "), selected_cs, selected_length
        ))
        idx <- selected_cs # Keep index with smallest length
      }
      return(idx)
    } else {
      return(NA_integer_)
    }
  }
  get_top_variants_idx <- function(susie_output, signal_cutoff) {
    c(which(susie_output$pip >= signal_cutoff), unlist(susie_output$sets$cs)) %>%
      unique() %>%
      sort()
  }
  get_cs_info <- function(susie_output_sets_cs, top_variants_idx) {
    cs_info_pri <- map_int(top_variants_idx, ~ get_cs_index(.x, susie_output_sets_cs))
    ifelse(is.na(cs_info_pri), 0, as.numeric(str_replace(names(susie_output_sets_cs)[cs_info_pri], "L", "")))
  }
  get_cs_and_corr <- function(susie_output, coverage, data_x, mode = c("susie", "susie_rss", "mvsusie")) {
    if (mode %in% c("susie", "mvsusie")) {
      susie_output_secondary <- list(sets = susie_get_cs(susie_output, X = data_x, coverage = coverage, min_abs_corr = min_abs_corr, median_abs_corr = median_abs_corr), pip = susie_output$pip)
      susie_output_secondary$cs_corr <- get_cs_correlation(susie_output_secondary, X = data_x)
      susie_output_secondary
    } else {
      susie_output_secondary <- list(sets = susie_get_cs(susie_output, Xcorr = data_x, coverage = coverage, min_abs_corr = min_abs_corr, median_abs_corr = median_abs_corr), pip = susie_output$pip)
      susie_output_secondary$cs_corr <- get_cs_correlation(susie_output_secondary, Xcorr = data_x)
      susie_output_secondary
    }
  }

  # Initialize result list
  res <- list(
    variant_names = format_variant_id(names(susie_output$pip))
  )
  analysis_script <- load_script()
  if (analysis_script != "") res$analysis_script <- analysis_script
  if (!is.null(other_quantities)) res$other_quantities <- other_quantities
  if (mode == "mvsusie") {
    res$condition_names <- susie_output$condition_names
  }
  if (!is.null(data_y)) {
    # Mode-specific processing
    if (mode == "susie") {
      # Processing specific to susie_post_processor
      res$sumstats <- univariate_regression(data_x, data_y)
      y_scalar <- if (is.null(y_scalar) || all(y_scalar == 1)) 1 else y_scalar
      X_scalar <- if (is.null(X_scalar) || all(X_scalar == 1)) 1 else X_scalar
      res$sumstats$betahat <- res$sumstats$betahat * y_scalar / X_scalar
      res$sumstats$sebetahat <- res$sumstats$sebetahat * y_scalar / X_scalar
      res$sample_names <- rownames(data_y)
    } else if (mode == "susie_rss") {
      # Processing specific to susie_rss_post_processor
      res$sumstats <- data_y
    }
  }
  if (!is.null(susie_output$V)) {
    # for fSuSiE there is no V for now
    eff_idx <- which(susie_output$V > prior_eff_tol)
  } else {
    eff_idx <- 1:nrow(susie_output$alpha)
  }
  if (length(eff_idx) > 0) {
    # Prepare for top loci table
    top_variants_idx_pri <- get_top_variants_idx(susie_output, signal_cutoff)
    cs_pri <- get_cs_info(susie_output$sets$cs, top_variants_idx_pri)
    susie_output$cs_corr <- if (mode %in% c("susie", "mvsusie")) get_cs_correlation(susie_output, X = data_x) else get_cs_correlation(susie_output, Xcorr = data_x)
    top_loci_list <- list("coverage_0.95" = data.frame(variant_idx = top_variants_idx_pri, cs_idx = cs_pri, stringsAsFactors = FALSE))

    ## Loop over each secondary coverage value
    sets_secondary <- list()
    for (sec_cov in secondary_coverage) {
      sets_secondary[[paste0("coverage_", sec_cov)]] <- get_cs_and_corr(susie_output, sec_cov, data_x, mode)
      top_variants_idx_sec <- get_top_variants_idx(sets_secondary[[paste0("coverage_", sec_cov)]], signal_cutoff)
      cs_sec <- get_cs_info(sets_secondary[[paste0("coverage_", sec_cov)]]$sets$cs, top_variants_idx_sec)
      top_loci_list[[paste0("coverage_", sec_cov)]] <- data.frame(variant_idx = top_variants_idx_sec, cs_idx = cs_sec, stringsAsFactors = FALSE)
    }

    # Iterate over the remaining tables, rename and merge them
    names(top_loci_list[[1]])[2] <- paste0("cs_", names(top_loci_list)[1])
    top_loci <- top_loci_list[[1]]
    for (i in 2:length(top_loci_list)) {
      names(top_loci_list[[i]])[2] <- paste0("cs_", names(top_loci_list)[i])
      top_loci <- full_join(top_loci, top_loci_list[[i]], by = "variant_idx")
    }
    if (nrow(top_loci) > 0) {
      top_loci[is.na(top_loci)] <- 0
      variants <- res$variant_names[top_loci$variant_idx]
      pip <- susie_output$pip[top_loci$variant_idx]
      top_loci_cols <- c("variant_id", if (!is.null(res$sumstats$betahat)) "betahat", if (!is.null(res$sumstats$sebetahat)) "sebetahat", if (!is.null(res$sumstats$z)) "z", if (!is.null(maf)) "maf", "pip", colnames(top_loci)[-1])
      res$top_loci <- data.frame(variants, stringsAsFactors = FALSE)
      res$top_loci$betahat <- if (!is.null(res$sumstats$betahat)) res$sumstats$betahat[top_loci$variant_idx] else NULL
      res$top_loci$sebetahat <- if (!is.null(res$sumstats$sebetahat)) res$sumstats$sebetahat[top_loci$variant_idx] else NULL
      res$top_loci$z <- if (!is.null(res$sumstats$z)) res$sumstats$z[top_loci$variant_idx] else NULL
      res$top_loci$maf <- if (!is.null(maf)) maf[top_loci$variant_idx] else NULL
      res$top_loci$pip <- pip
      res$top_loci <- cbind(res$top_loci, top_loci[, -1])
      colnames(res$top_loci) <- top_loci_cols
      rownames(res$top_loci) <- NULL
    }
    names(susie_output$pip) <- NULL
    res$susie_result_trimmed <- list(
      pip = susie_output$pip,
      sets = susie_output$sets,
      cs_corr = susie_output$cs_corr,
      sets_secondary = lapply(sets_secondary, function(x) x[names(x) != "pip"]),
      alpha = susie_output$alpha[eff_idx, , drop = FALSE],
      lbf_variable = susie_output$lbf_variable[eff_idx, , drop = FALSE],
      V = if (!is.null(susie_output$V)) susie_output$V[eff_idx] else NULL,
      niter = susie_output$niter
    )
    if (mode == "susie") {
      res$susie_result_trimmed$X_column_scale_factors <- susie_output$X_column_scale_factors
      res$susie_result_trimmed$mu <- susie_output$mu[eff_idx, , drop = FALSE]
      res$susie_result_trimmed$mu2 <- susie_output$mu2[eff_idx, , drop = FALSE]
    }
    if (mode == "mvsusie") {
      # res$susie_result_trimmed$b1 = susie_output$b1[eff_idx, , , drop = FALSE]
      # res$susie_result_trimmed$b2 = susie_output$b2[eff_idx, , , drop = FALSE]
      res$susie_result_trimmed$b1_rescaled <- susie_output$b1_rescaled[eff_idx, , , drop = FALSE]
      res$susie_result_trimmed$coef <- susie_output$coef
      res$susie_result_trimmed$clfsr <- susie_output$conditional_lfsr[eff_idx, , , drop = FALSE]
      # other lfsr can be computed:
      # se_lfsr <- mvsusie_single_effect_lfsr(clfsr, alpha)
      # lfsr <- mvsusie_get_lfsr(clfsr, alpha)
    }
    class(res$susie_result_trimmed) <- "susie"
  }
  return(res)
}
