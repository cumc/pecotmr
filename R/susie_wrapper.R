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
lbf_to_alpha_vector = function(lbf, prior_weights = NULL) {
      if (is.null(prior_weights)) prior_weights = rep(1/length(lbf), length(lbf))
      maxlbf = max(lbf)
      # If maxlbf is 0, return a vector of zeros
      if (maxlbf == 0) {
        return(setNames(rep(0, length(lbf)), names(lbf)))
      }
      # w is proportional to BF, subtract max for numerical stability.
      w = exp(lbf - maxlbf)
      # Posterior prob for each SNP.
      w_weighted = w * prior_weights
      weighted_sum_w = sum(w_weighted)
      alpha = w_weighted / weighted_sum_w
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
lbf_to_alpha = function(lbf) t(apply(lbf, 1, lbf_to_alpha_vector))

#' @importFrom susieR susie
#' @export 
susie_wrapper = function(X, y, init_L = 10, max_L = 30, coverage = 0.95, max_iter=500, l_step = 5) {
        L = init_L
        # Perform SuSiE by dynamically increase L
        while (TRUE) {
            st = proc.time()
            res <- susie(X,y, L=L,
                             max_iter=500,
                             estimate_residual_variance=TRUE,
                             estimate_prior_variance=TRUE,
                             refine=TRUE,
                             compute_univariate_zscore=FALSE,
                             min_abs_corr=0.5,
                             median_abs_corr=0.8,
                             coverage=coverage)
            res$time_elapsed <- proc.time() - st
            if (!is.null(res$sets$cs)) {
                if (length(res$sets$cs)>=L && L<=max_L) {
                  L = L + l_step
                } else {
                  break
                }
            } else {
              break
           }
        }
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
                              zR_discrepancy_correction = FALSE, ...) {
  if (L == 1) {
    return(susie_rss(z = z, R = R, bhat = bhat, shat = shat, var_y = var_y, n = n, 
                    L = 1, max_iter = 1, correct_zR_discrepancy = FALSE, ...))
  }
  while (TRUE) {
      st = proc.time()
      susie_rss_result <- susie_rss(z = z, R = R, bhat = bhat, shat = shat, var_y = var_y, n = n, L = L,
                                    correct_zR_discrepancy = zR_discrepancy_correction, ...)
      susie_rss_result$time_elapsed <- proc.time() - st
    # Check for convergence and adjust L if necessary
    if (!is.null(susie_rss_result$sets$cs)) {
      if (length(susie_rss_result$sets$cs) >= L && L <= max_L) {
        L <- L + l_step  # Increase L for the next iteration
      } else {
        break
      }
    } else {
      break  # Break the loop if no credible sets are found
    }
  }

  # Error handling if the model fails to converge


  return(susie_rss_result)
}

#' SuSiE RSS Analysis with Quality Control and Imputation
#'
#' Performs SuSiE RSS analysis with optional quality control steps that include
#' z-score and LD matrix discrepancy correction and imputation for outliers. It leverages
#' the `susie_rss` function for the core analysis and provides additional functionality
#' for handling data discrepancies and missing values.
#'
#' @param z Numeric vector of z-scores corresponding to the effect size estimates, with names matching the reference panel's variant IDs.
#' @param R Numeric matrix representing the LD (linkage disequilibrium) matrix.
#' @param ref_panel Data frame with 'chrom', 'pos', 'variant', 'A1', 'A2' column that matches the names of z.
#' @param bhat Optional numeric vector of effect size estimates.
#' @param shat Optional numeric vector of standard errors associated with the effect size estimates.
#' @param var_y Optional numeric value representing the total phenotypic variance.
#' @param n Optional numeric value representing the sample size used in the analysis. 
#' @param L Initial number of causal configurations to consider in the analysis.
#' @param max_L Maximum number of causal configurations to consider when dynamically adjusting L.
#' @param l_step Step size for increasing L when the limit is reached during dynamic adjustment.
#' @param lamb Regularization parameter for the RAiSS imputation method.
#' @param rcond Condition number for the RAiSS imputation method.
#' @param R2_threshold R-squared threshold for the RAiSS imputation method.
#' @param minimum_ld Minimum number of LD values for the RAiSS imputation method.
#' @param impute Logical; if TRUE, performs imputation for outliers identified in the analysis.
#' @param output_qc Logical; if TRUE, includes QC-only results in the output.
#' @return A list containing the results of the SuSiE RSS analysis after applying quality control measures and optional imputation.
#' @importFrom susieR susie_rss
#' @import dplyr
#' @export
susie_rss_qc <- function(z, R, ref_panel, bhat=NULL, shat=NULL, var_y=NULL, n = NULL, L = 10, max_L = 20, l_step = 5, 
                        lamb = 0.01, rcond = 0.01, R2_threshold = 0.6, minimum_ld = 5, impute = TRUE, output_qc = TRUE, ...) {
  
    
    
  ## Input validation for z-scores and reference panel
  if (length(z) != nrow(ref_panel)) {
    stop("The length of z-scores does not match the number of rows in the reference panel.")
  }
  
  if (is.null(names(z))) {
    warning("Z-score names are NULL. Assuming names match reference panel variant_id.")
  } else if (!all(names(z) %in% ref_panel$variant_id)) {
    stop("Names of z-scores do not match the reference panel variant_id.")
  }

  ## Perform initial SuSiE RSS analysis with discrepancy correction
  result <- susie_rss(z=z, R=R, bhat=bhat, shat=shat, var_y=var_y, n=n, L=max_L,  
                     correct_zR_discrepancy=TRUE, track_fit = TRUE, ...)

  ## Initialize result_final
  result_final <- NULL

  ## Imputation for outliers if enabled and required
  if (impute && !is.null(result$zR_outliers) && length(result$zR_outliers) > 0) {
    ## Extracting known z-scores excluding outliers
    outlier = result$zR_outliers
    known_zscores = ref_panel[-outlier]
    known_zscorse$Z = z[-outlier]
    known_zscores = known_zscore %>% arrange(pos)
    
    ## Imputation logic using RAiSS or other methods
    imputation_result <- raiss(ref_panel, known_zscores, R, lamb = lamb, rcond = rcond, 
                               R2_threshold = R2_threshold, minimum_ld = minimum_ld)
    
    ## Filter out variants not included in the imputation result
    filtered_out_variant <- setdiff(ref_panel$variant_id, imputation_result$variant_id)
    
    ## Update the LD matrix excluding filtered variants
    LD_extract_filtered <- if (length(filtered_out_variant) > 0) {
      filtered_out_id <- match(filtered_out_variant, ref_panel$variant_id)
      as.matrix(R)[-filtered_out_id, -filtered_out_id]
    } else {
      as.matrix(R)
    }

    ## Re-run SuSiE RSS with imputed z-scores and updated LD matrix
    result_final$qc_impute_result <- susie_rss_wrapper(z=imputation_result$Z, R=LD_extract_filtered, bhat=bhat, shat=shat, var_y=var_y, 
                                      n=n, L=L, max_L=max_L, l_step=l_step, zR_discrepancy_correction=FALSE, ...)
    result_final$qc_impute_result$z = imputation_result$Z
  }
    if(output_qc){
    result_final$qc_only_result <- result
    }
    return(result_final)
}



#' Post-process SuSiE or SuSiE_rss Analysis Results
#'
#' This function processes the results from SuSiE or SuSiE_rss (Sum of Single Effects) genetic analysis.
#' It extracts and processes various statistics and indices based on the provided SuSiE object and other parameters.
#' The function can operate in two modes: 'susie' and 'susie_rss', based on the method used for the SuSiE analysis.
#'
#' @param susie_output Output from running susieR::susie() or susieR::susie_rss()
#' @param data_x Genotype data matrix for 'susie' or Xcorr matrix for 'susie_rss'.
#' @param data_y Phenotype data vector for 'susie' or summary stats object for 'susie_rss' (a list contain attribute betahat and sebetahat AND/OR z). i.e. data_y = list(betahat = ..., sebetahat = ...)
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
#' # result <- combined_susie_post_processor(susie_output, X_data, y_data, maf, mode = "susie")
#' # Example usage for SuSiE_rss
#' # result <- combined_susie_post_processor(susie_output, Xcorr, z, maf, mode = "susie_rss")
#' @importFrom dplyr full_join
#' @importFrom purrr map_int pmap
#' @importFrom susieR get_cs_correlation susie_get_cs
#' @importFrom stringr str_replace
#' @export
susie_post_processor <- function(susie_output, data_x, data_y, X_scalar, y_scalar, maf, 
                                secondary_coverage = c(0.5, 0.7), signal_cutoff = 0.1, 
                                other_quantities = NULL, prior_eff_tol = 1e-9, 
                                mode = c("susie", "susie_rss")) {
    mode <- match.arg(mode)
    get_cs_index <- function(snps_idx, susie_cs) {
        idx <- tryCatch(
            which(
                pmap(list(a = susie_cs), function(a) snps_idx %in% a) %>% unlist()
            ),
            error = function(e) NA_integer_
        )
        if(length(idx) == 0) return(NA_integer_)
        return(idx)
    }
    get_top_variants_idx <- function(susie_output, signal_cutoff) {
        c(which(susie_output$pip >= signal_cutoff), unlist(susie_output$sets$cs)) %>% unique %>% sort
    }
    get_cs_info <- function(susie_output_sets_cs, top_variants_idx) {
        cs_info_pri <- map_int(top_variants_idx, ~get_cs_index(.x, susie_output_sets_cs))
        ifelse(is.na(cs_info_pri), 0, as.numeric(str_replace(names(susie_output_sets_cs)[cs_info_pri], "L", "")))
    }
    get_cs_and_corr <- function(susie_output, coverage, data_x, mode = c("susie", "susie_rss")) {
        if (mode == "susie") {
            susie_output_secondary <- list(sets = susie_get_cs(susie_output, X = data_x, coverage = coverage), pip = susie_output$pip)
            susie_output_secondary$cs_corr <- get_cs_correlation(susie_output_secondary, X = data_x)
            susie_output_secondary
        } else {
            susie_output_secondary <- list(sets = susie_get_cs(susie_output, Xcorr = data_x, coverage = coverage), pip = susie_output$pip)
            susie_output_secondary$cs_corr <- get_cs_correlation(susie_output_secondary, Xcorr = data_x)
            susie_output_secondary
        }
    }

    # Initialize result list
    res <- list(other_quantities = other_quantities,
                susie_result_trimmed = list(),
                analysis_script = load_script(),
                variant_names = format_variant_id(names(susie_output$pip)))

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

    eff_idx <- which(susie_output$V > prior_eff_tol)
    if (length(eff_idx) > 0) {
        # Prepare for top loci table
        top_variants_idx_pri <- get_top_variants_idx(susie_output, signal_cutoff)
        cs_pri <- get_cs_info(susie_output$sets$cs, top_variants_idx_pri)
        susie_output$cs_corr <- if (mode=="susie") get_cs_correlation(susie_output, X = data_x) else get_cs_correlation(susie_output, Xcorr = data_x) 
        top_loci_list <- list("coverage_0.95" = data.frame(variant_idx = top_variants_idx_pri, cs_idx = cs_pri, stringsAsFactors=FALSE))
        
        ## Loop over each secondary coverage value
        sets_secondary <- list()
        for (sec_cov in secondary_coverage) {
            sets_secondary[[paste0("coverage_", sec_cov)]] <- get_cs_and_corr(susie_output, sec_cov, data_x, mode)
            top_variants_idx_sec <- get_top_variants_idx(sets_secondary[[paste0("coverage_", sec_cov)]], signal_cutoff)
            cs_sec <- get_cs_info(sets_secondary[[paste0("coverage_", sec_cov)]]$sets$cs, top_variants_idx_sec)
            top_loci_list[[paste0("coverage_", sec_cov)]] <- data.frame(variant_idx = top_variants_idx_sec, cs_idx = cs_sec, stringsAsFactors=FALSE)
        }

        # Iterate over the remaining tables, rename and merge them
        names(top_loci_list[[1]])[2] <- paste0("cs_", names(top_loci_list)[1])
        top_loci <- top_loci_list[[1]]
        for (i in 2:length(top_loci_list)) {
            names(top_loci_list[[i]])[2] <- paste0("cs_", names(top_loci_list)[i])
            top_loci <- full_join(top_loci, top_loci_list[[i]], by = "variant_idx")
        }
        top_loci[is.na(top_loci)] <- 0
        variants <- res$variant_names[top_loci$variant_idx]
        pip <- susie_output$pip[top_loci$variant_idx]
        top_loci_cols <- c("variant_id" , if (!is.null(res$sumstats$betahat)) "betahat", if (!is.null(res$sumstats$sebetahat)) "sebetahat", if (!is.null(res$sumstats$z)) "z", if (!is.null(maf)) "maf", "pip" , colnames(top_loci)[-1])
        res$top_loci <- data.frame(variants, stringsAsFactors = FALSE)
        res$top_loci$betahat = if (!is.null(res$sumstats$betahat)) res$sumstats$betahat[top_loci$variant_idx] else NULL
        res$top_loci$sebetahat = if (!is.null(res$sumstats$sebetahat)) res$sumstats$sebetahat[top_loci$variant_idx] else NULL
        res$top_loci$z = if (!is.null(res$sumstats$z)) res$sumstats$z[top_loci$variant_idx] else NULL
        res$top_loci$maf = if (!is.null(maf)) maf[top_loci$variant_idx] else NULL
        res$top_loci$pip = pip
        res$top_loci = cbind(res$top_loci, top_loci[,-1])
        colnames(res$top_loci) <- top_loci_cols
        rownames(res$top_loci) <- NULL
        names(susie_output$pip) <- NULL
        res$susie_result_trimmed <- list(
            pip = susie_output$pip,
            sets = susie_output$sets,
            cs_corr = susie_output$cs_corr,
            sets_secondary = sets_secondary,
            alpha = susie_output$alpha[eff_idx, , drop = FALSE],
            lbf_variable = susie_output$lbf_variable[eff_idx, , drop = FALSE],
            mu = susie_output$mu[eff_idx, , drop = FALSE],
            mu2 = susie_output$mu2[eff_idx, , drop = FALSE],
            V = susie_output$V[eff_idx],
            X_column_scale_factors = if (mode == "susie") susie_output$X_column_scale_factors else NULL
        )
        class(res$susie_result_trimmed) <- "susie"
    } 
    return(res)
}