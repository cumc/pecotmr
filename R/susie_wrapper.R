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
  
  # w is proportional to BF, subtract max for numerical stability
  w = exp(lbf - maxlbf)
  
  # Posterior prob for each SNP
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
susie_wrapper = function(X, y, init_L = 10, max_L = 30, l_step = 5, ...) {
  if (init_L == max_L) {
    return(susie(X, y, L = init_L, median_abs_corr = 0.8, ...))
  }
  L = init_L
  # Perform SuSiE by dynamically increasing L
  gst = proc.time()
  while (TRUE) {
    st = proc.time()
    res <- susie(X, y, L = L,
                 median_abs_corr = 0.8, ...)
    res$time_elapsed <- proc.time() - st
    if (!is.null(res$sets$cs)) {
      if (length(res$sets$cs) >= L && L <= max_L) {
        L = L + l_step
      } else {
        break
      }
    } else {
      break
    }
  }
  message(paste("Total time elapsed for susie_wrapper:", (proc.time()-gst)[3]))
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
    return(susie_rss(z = z, R = R, bhat = bhat, shat = shat, var_y = var_y, n = n, 
                    L = 1, max_iter = 1, median_abs_corr = 0.8, correct_zR_discrepancy = FALSE, coverage = coverage, ...))
  }
  if (L == max_L) {
    return(susie_rss(z = z, R = R, bhat = bhat, shat = shat, var_y = var_y, n = n, L = L, median_abs_corr = 0.8,
                                    correct_zR_discrepancy = zR_discrepancy_correction, coverage = coverage, ...))
  }
  while (TRUE) {
      st = proc.time()
      susie_rss_result <- susie_rss(z = z, R = R, bhat = bhat, shat = shat, var_y = var_y, n = n, L = L, median_abs_corr = 0.8,
                                    correct_zR_discrepancy = zR_discrepancy_correction,coverage = coverage, ...)
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

  return(susie_rss_result)
}

#' Get input parameters for the SuSiE_rss function
#'
#' This function formats the input summary statistics dataframe with uniform column names
#' to fit into the SuSiE pipeline. The mapping is performed through the specified column file.
#' Additionally, it extracts sample size, case number, control number, and variance of Y.
#'
#' @param sumstat_path File path to the summary statistics.
#' @param column_file_path File path to the column file for mapping.
#' @param n_sample User-specified sample size. If unknown, set as 0 to retrieve from the sumstat file.
#' @param n_case User-specified number of cases.
#' @param n_control User-specified number of controls.
#'
#' @return A list of rss_input, including the column-name-formatted summary statistics,
#' sample size (n), and var_y.
#'
#' @importFrom data.table fread
#' @export
get_rss_input = function(sumstat_path, column_file_path, n_sample, n_case, n_control) {
    var_y = NULL
    sumstats = fread(sumstat_path)
    column_data <- read.table(column_file_path, header = FALSE, sep = ":", stringsAsFactors = FALSE)
    colnames(column_data) = c("standard", "original")
    count = 1
    for (name in colnames(sumstats)) {
        if (name %in% column_data$original) {
            index = which(column_data$original == name)
            colnames(sumstats)[count] = column_data$standard[index]
        }
        count = count + 1
    }

    if (length(sumstats$z) == 0) {
        sumstats$z = sumstats$beta / sumstats$se
    }

    if (length(sumstats$beta) == 0) {
        sumstats$beta = sumstats$z
        sumstats$se = 1
    }

    if (n_sample != 0 & (n_case + n_control) != 0) {
        stop("Please provide sample size, or case number with control number, but not both")
    } else if (n_sample != 0) {
        n = n_sample
    } else if ((n_case + n_control) != 0) {
        n = n_case + n_control
        phi = n_case/n
        var_y = residual_variance = 1 / (phi * (1 - phi))
    } else {
        if (length(sumstats$n_sample) != 0) {
            n = median(sumstats$n_sample)
        } else if (length(sumstats$n_case) != 0 & length(sumstats$n_control) != 0) {
            n = median(sumstats$n_case + sumstats$n_control)
            phi = median(sumstats$n_case / n)
            var_y = residual_variance = 1 / (phi * (1 - phi))
        } else {
            n = NULL
        }
    }
    return(rss_input = list(sumstats = sumstats, n = n, var_y = var_y))
}

#' Preprocess input data for RSS analysis
#'
#' This function preprocesses summary statistics and LD data for RSS analysis.
#' It performs allele quality control, flipping alleles as necessary, and removes
#' specified regions from the analysis.
#'
#' @param sumstats A data frame containing summary statistics with columns "chrom", "pos", "A1", and "A2".
#' @param LD_data A list containing combined LD variants data that is generated by load_LD_matrix.
#' @param skip_region A character vector specifying regions to be skipped in the analysis (optional).
#'
#' @return A processed data frame containing summary statistics after preprocessing.
#' @import dplyr
#' @import tibble
#'
#' @export
rss_input_preprocess = function(sumstats, LD_data, skip_region = NULL) {
    target_variants = sumstats[, c("chrom", "pos", "A1", "A2")]
    ref_variants = LD_data$combined_LD_variants
    allele_flip = allele_qc(target_variants, ref_variants, sumstats, col_to_flip = c("beta", "z"), match.min.prop = 0.2, remove_dups = FALSE, flip = TRUE, remove = TRUE)

    if (length(skip_region) != 0) {
        skip_table = tibble(region = skip_region) %>% separate(region, into = c("chrom", "start", "end"), sep = "[:-]")
        skip_variant = c()
        for (i in 1:nrow(skip_table)) {
            variant = allele_flip$target_data_qced %>% filter(chrom == skip_table$chrom[i] & pos > skip_table$start[i] & pos < skip_table$end[i]) %>% pull(variant_id)
            skip_variant = c(skip_variant, variant)
        }
        allele_flip$target_data_qced = allele_flip$target_data_qced %>% filter(!(variant_id) %in% skip_variant)
    }

    sumstat_processed = allele_flip$target_data_qced %>% arrange(pos)
    return(sumstat_processed)
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
#' @import dplyr 
#' @import susieR 
#' @import tibble
#' @export
susie_rss_pipeline = function(sumstat, R, ref_panel, n, L, var_y, QC = TRUE, impute = TRUE, bayesian_conditional_analysis = TRUE, lamb = 0.01, rcond = 0.01, R2_threshold = 0.6,
                              max_L = 20, l_step = 5, minimum_ld = 5, coverage = 0.95,
                              secondary_coverage = c(0.7, 0.5), pip_cutoff_to_skip = 0.025, signal_cutoff = 0.1) {
  
  if (!is.null(sumstat$z)) {
    z = sumstat$z 
    bhat = NULL
    shat = NULL
  } else if ((!is.null(sumstat$beta)) && (!is.null(sumstat$se))) {
    z = sumstat$beta / sumstat$se
    bhat = NULL
    shat = NULL
  } else {
    stop("Sumstat should have z or (bhat and shat)")
  }
  
  final_result = list()
  
  LD_extract = R[sumstat$variant_id, sumstat$variant_id, drop = FALSE]
  single_effect_res = susie_rss_wrapper(z = z, R = LD_extract, bhat = bhat, shat = shat, L = 1, n = n, var_y = var_y, coverage = coverage)
  if (max(single_effect_res$pip) < pip_cutoff_to_skip) {
    cat(paste0("No PIP larger than ", pip_cutoff_to_skip, " in this region."))
    if (impute) {
      z = sumstat$z
      known_zscores = sumstat %>% arrange(pos)
      final_result$sumstats_qc_impute = raiss(ref_panel, known_zscores, R, lamb = lamb, rcond = rcond, R2_threshold = R2_threshold, minimum_ld = minimum_ld)$result_nofilter
      final_result$sumstats_qc_impute$chrom = as.numeric(final_result$sumstats_qc_impute$chrom)
    }
  } else {
    single_effect_post = susie_post_processor(single_effect_res, data_x = LD_extract, data_y = list(z = z), signal_cutoff = signal_cutoff, secondary_coverage = secondary_coverage, mode = "susie_rss")
    final_result$single_effect_regression = single_effect_post
    result_noqc = susie_rss_wrapper(z = z, R = LD_extract, bhat = bhat, shat = shat, n = n, L = L, var_y = var_y, coverage = coverage)
    result_noqc_post = susie_post_processor(result_noqc, data_x = LD_extract, data_y = list(z = z), signal_cutoff = signal_cutoff, secondary_coverage = secondary_coverage, mode = "susie_rss")
    final_result$noqc = result_noqc_post
    
    if (QC) {
      result_qced = susie_rss_qc(sumstat, ref_panel = ref_panel, R = R, n = n, L = L, impute = impute, lamb = lamb, rcond = rcond, R2_threshold = R2_threshold, max_L = max_L, minimum_ld = minimum_ld, l_step = l_step, var_y = var_y, coverage = coverage)  
      var_impute_kept = names(result_qced$qc_impute_result$pip)
      result_qced_impute_post = susie_post_processor(result_qced$qc_impute_result, data_x = R[var_impute_kept, var_impute_kept, drop = FALSE], data_y = list(z = result_qced$qc_impute_result$z), signal_cutoff = signal_cutoff, secondary_coverage = secondary_coverage, mode = "susie_rss")
      result_qced_only_post = susie_post_processor(result_qced$qc_only_result, data_x = LD_extract, data_y = list(z = z), signal_cutoff = signal_cutoff, secondary_coverage = secondary_coverage, mode = "susie_rss")
      
      final_result$qc_impute = result_qced_impute_post
      final_result$qc_only = result_qced_only_post
      final_result$qc_only$outlier = result_qced$qc_only_result$zR_outliers
      
      result_qced$sumstats_qc_impute_filtered$chrom = as.numeric(result_qced$sumstats_qc_impute_filtered$chrom)
      result_qced$sumstats_qc_impute$chrom = as.numeric(result_qced$sumstats_qc_impute$chrom)
      final_result$sumstats_qc_impute_filtered = result_qced$sumstats_qc_impute_filtered %>% select(-beta, -se)
      final_result$sumstats_qc_impute = result_qced$sumstats_qc_impute %>% select(-beta, -se)
    }
    
    if (bayesian_conditional_analysis) {
      conditional_noqc = susie_rss_wrapper(z = z, R = LD_extract, bhat = bhat, shat = shat, n = n, L = L, max_iter = 1, var_y = var_y, coverage = coverage)
      conditional_noqc_post = susie_post_processor(conditional_noqc, data_x = LD_extract, data_y = list(z = z), signal_cutoff = signal_cutoff, secondary_coverage = secondary_coverage, mode = "susie_rss")
      final_result$conditional_regression_noqc = conditional_noqc_post
      
      if (QC) {
        outlier = result_qced$qc_only_result$zR_outliers
        if(!is.null(outlier) & length(outlier) != 0){
            z_rmoutlier = z[-outlier]
            LD_rmoutlier = LD_extract[-outlier, -outlier, drop = FALSE]
        }else{
            z_rmoutlier = z
            LD_rmoutlier = LD_extract
        }
        conditional_qc_only = susie_rss_wrapper(z = z_rmoutlier , R = LD_rmoutlier, bhat = bhat, shat = shat, n = n, L = L, max_iter = 1, var_y = var_y, coverage = coverage)
        conditional_qc_only_post = susie_post_processor(conditional_qc_only, data_x = LD_rmoutlier, data_y = list(z = z_rmoutlier), signal_cutoff = signal_cutoff, secondary_coverage = secondary_coverage, mode = "susie_rss")
        
        if (impute) {
          conditional_qced_impute = susie_rss_wrapper(z = result_qced$qc_impute_result$z, R = R[var_impute_kept, var_impute_kept, drop = FALSE], bhat = bhat, shat = shat, n = n, L = L, max_iter = 1, var_y = var_y, coverage = coverage)
          conditional_qced_impute_post = susie_post_processor(conditional_qced_impute, data_x = R[var_impute_kept, var_impute_kept, drop = FALSE], data_y = list(z = result_qced$qc_impute_result$z), signal_cutoff = signal_cutoff, secondary_coverage = secondary_coverage, mode = "susie_rss")
          final_result$conditional_regression_qc_impute = conditional_qced_impute_post
        }
        
        final_result$conditional_regression_qc_only =  conditional_qc_only_post
      }
    }
  }
  
  return(final_result)
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
#' @param ref_panel Data frame with 'chrom', 'pos', 'variant_id', 'A1', 'A2' column that matches the names of z.
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
susie_rss_qc <- function(sumstat, R, ref_panel, bhat=NULL, shat=NULL, var_y=NULL, n = NULL, L = 10, max_L = 20, l_step = 5, 
                        lamb = 0.01, rcond = 0.01, R2_threshold = 0.6, minimum_ld = 5, impute = TRUE, output_qc = TRUE, coverage = 0.95, ...) {
  z = sumstat$z
  ## Input validation for z-scores and reference panel
  
  if (is.null(names(z))) {
    warning("Z-score names are NULL. Assuming names match reference panel variant_id.")
  } else if (!all(names(z) %in% ref_panel$variant_id)) {
    stop("Names of z-scores do not match the reference panel variant_id.")
  }

  ## Perform initial SuSiE RSS analysis with discrepancy correction
  ## FIXME: should parse ...  directly and use do.call  for function calls
  LD_extract = R[sumstat$variant_id, sumstat$variant_id, drop = FALSE]
  result <- susie_rss(z=z, R=LD_extract, bhat=bhat, shat=shat, var_y=var_y, n=n, L=max_L,  
                     correct_zR_discrepancy=TRUE, track_fit = TRUE, max_iter = 100, coverage = coverage)

  ## Initialize result_final
  result_final <- NULL

  ## Imputation for outliers if enabled and required
  if (impute) {
    ## Extracting known z-scores excluding outliers
    if(!is.null(result$zR_outliers) & length(result$zR_outliers) != 0){
        outlier = result$zR_outliers
      
        known_zscores = sumstat[-outlier,]
        known_zscores = known_zscores %>% arrange(pos)
    }else{
        known_zscores = sumstat %>% arrange(pos)
    }
    
    ## Imputation logic using RAiSS or other methods
    imputation_res <- raiss(ref_panel, known_zscores, R, lamb = lamb, rcond = rcond, 
                               R2_threshold = R2_threshold, minimum_ld = minimum_ld)
      
    imputation_result_nofilter = imputation_res$result_nofilter
    imputation_result_filter = imputation_res$result_filter 
    ## Filter out variants not included in the imputation result
    filtered_out_variant <- setdiff(ref_panel$variant_id, imputation_result_filter$variant_id)
    
    ## Update the LD matrix excluding filtered variants
    LD_extract_filtered <- if (length(filtered_out_variant) > 0) {
      filtered_out_id <- match(filtered_out_variant, ref_panel$variant_id)
      as.matrix(R)[-filtered_out_id, -filtered_out_id]
    } else {
      as.matrix(R)
    }

    ## Re-run SuSiE RSS with imputed z-scores and updated LD matrix
    result_final$qc_impute_result <- susie_rss_wrapper(z=imputation_result_filter$z, R=LD_extract_filtered, bhat=bhat, shat=shat, var_y=var_y, 
                                      n=n, L=L, max_L=max_L, l_step=l_step, zR_discrepancy_correction=FALSE, coverage = coverage, ...)
    result_final$qc_impute_result$z = imputation_result_filter$z
    result_final$sumstats_qc_impute_filtered = imputation_result_filter
    result_final$sumstats_qc_impute = imputation_result_nofilter
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
susie_post_processor <- function(susie_output, data_x, data_y, X_scalar, y_scalar, maf = NULL, 
                                secondary_coverage = c(0.5, 0.7), signal_cutoff = 0.1, 
                                other_quantities = NULL, prior_eff_tol = 1e-9, min_abs_corr= 0.5, 
                                median_abs_corr = 0.8,
                                mode = c("susie", "susie_rss")) {
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
                }, error = function(e) NA_integer_ 
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

                warning(sprintf("Variable %d found in multiple CS: %s. Keeping smallest: CS %d (length %d).",
                                snps_idx, paste(idx, collapse = ', '), selected_cs, selected_length))
                idx <- selected_cs  # Keep index with smallest length
            }
            return(idx)
        } else {
            return(NA_integer_)
        }
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
            susie_output_secondary <- list(sets = susie_get_cs(susie_output, X = data_x, coverage = coverage, min_abs_corr=min_abs_corr, median_abs_corr=median_abs_corr), pip = susie_output$pip)
            susie_output_secondary$cs_corr <- get_cs_correlation(susie_output_secondary, X = data_x)
            susie_output_secondary
        } else {
            susie_output_secondary <- list(sets = susie_get_cs(susie_output, Xcorr = data_x, coverage = coverage, min_abs_corr=min_abs_corr, median_abs_corr=median_abs_corr), pip = susie_output$pip)
            susie_output_secondary$cs_corr <- get_cs_correlation(susie_output_secondary, Xcorr = data_x)
            susie_output_secondary
        }
    }

    # Initialize result list
    res <- list(other_quantities = other_quantities,
                analysis_script = load_script(),
                variant_names = format_variant_id(names(susie_output$pip)))
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
        if (nrow(top_loci) > 0) {
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
        }
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
            V = if (!is.null(susie_output$V)) susie_output$V[eff_idx] else NULL,
            niter = susie_output$niter,
            X_column_scale_factors = if (mode == "susie") susie_output$X_column_scale_factors else NULL
        )
        class(res$susie_result_trimmed) <- "susie"
    } 
    return(res)
}