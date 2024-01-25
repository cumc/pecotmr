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

#' Post-process SuSiE Analysis Results
#'
#' This function processes the results from SuSiE (Sum of Single Effects) genetic analysis. 
#' It extracts and processes various statistics and indices based on the provided SuSiE object and other parameters.
#' It has been updated to handle multiple secondary coverage values, performing secondary analysis for each and combining the top loci.
#' Also, it tracks cs_secondary correlation for each coverage value and includes this information in the top_loci table.
#'
#' @param susie_output Output from running susieR::susie()
#' @param X_data The genotype data matrix.
#' @param y_data The phenotype data vector.
#' @param X_scalar Scalar for the genotype data, used in residual scaling.
#' @param y_scalar Scalar for the phenotype data, used in residual scaling.
#' @param maf Minor Allele Frequencies vector.
#' @param secondary_coverage Vector of coverage thresholds for secondary conditional analysis.
#' @param signal_cutoff Cutoff value for signal identification in PIP values. Default is 0.1.
#' @param other_quantities A list of other quantities to be added to the final object.
#' @return A list containing modified SuSiE object along with additional post-processing information.
#' @examples
#' # Example usage
#' # susie_result <- susie_post_processor(susie_output, X_data, y_data, X_scalar, y_scalar, maf, c(0.5, 0.7))
#' @importFrom susieR get_cs_correlation susie_get_cs
#' @importFrom stringr str_replace
#' @export
susie_post_processor <- function(susie_output, X_data, y_data, X_scalar, y_scalar, maf, 
                                 secondary_coverage = c(0.5, 0.7), signal_cutoff = 0.1, 
                                 other_quantities = list(), prior_eff_tol = 0) {
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
        # it is okay to have PIP = NULL here
        c(which(susie_output$pip >= signal_cutoff), unlist(susie_output$sets$cs)) %>% unique %>% sort
    }
    get_cs_info <- function(susie_output_sets_cs, top_variants_idx) {
        cs_info_pri <- map_int(top_variants_idx, ~get_cs_index(.x, susie_output_sets_cs))
        ifelse(is.na(cs_info_pri), 0, as.numeric(str_replace(names(susie_output_sets_cs)[cs_info_pri], "L", "")))
    }
    get_cs_and_corr <- function(susie_output, coverage, X_data) {
        susie_output_secondary <- list(sets = susie_get_cs(susie_output, X_data, coverage = coverage), pip = susie_output$pip)
        susie_output_secondary$cs_corr <- get_cs_correlation(susie_output_secondary, X = X_data)
        susie_output_secondary
    }
    # Compute univariate regression results 
    res <- list(sumstats = univariate_regression(X_data, y_data), 
                analysis_script = load_script(), 
                other_quantities = other_quantities,
                sample_names = rownames(y_data),
                variant_names = format_variant_id(names(susie_output$pip)),
                susie_result_trimmed = list()
            )
    eff_idx <- which(susie_output$V > prior_eff_tol)
    if (length(eff_idx) > 0) {
        # Prepare for top loci table
        top_variants_idx_pri <- get_top_variants_idx(susie_output, signal_cutoff)
        cs_pri <- get_cs_info(susie_output$sets$cs, top_variants_idx_pri)
        susie_output$cs_corr <- get_cs_correlation(susie_output, X = X_data)
        top_loci_list <- list("coverage_0.95" = data.frame(variant_idx = top_variants_idx_pri, cs_idx = cs_pri, stringsAsFactors=F))

        ## Loop over each secondary coverage value
        sets_secondary <- list()
        for (sec_cov in secondary_coverage) {
            sets_secondary[[paste0("coverage_", sec_cov)]] <- get_cs_and_corr(susie_output, sec_cov, X_data)
            top_variants_idx_sec <- get_top_variants_idx(sets_secondary[[paste0("coverage_", sec_cov)]], signal_cutoff)
            cs_sec <- get_cs_info(sets_secondary[[paste0("coverage_", sec_cov)]]$sets$cs, top_variants_idx_sec)
            top_loci_list[[paste0("coverage_", sec_cov)]] <-data.frame(variant_idx = top_variants_idx_sec, cs_idx = cs_sec, stringsAsFactors=F) 
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
        y_scalar <- if (is.null(y_scalar) || all(y_scalar == 1)) 1 else y_scalar[top_loci$variant_idx]
        X_scalar <- if (is.null(X_scalar) || all(X_scalar == 1)) 1 else X_scalar[top_loci$variant_idx]
        top_loci_cols <- c("variant_id", if (!is.null(maf)) "maf", "bhat", "sbhat", "pip", colnames(top_loci)[-1])
        res$top_loci <- data.frame(variants,
                                    res$sumstats$betahat[top_loci$variant_idx] * y_scalar / X_scalar, 
                                    res$sumstats$sebetahat[top_loci$variant_idx] * y_scalar / X_scalar, 
                                    pip, stringsAsFactors = FALSE)
        res$top_loci$maf = if (!is.null(maf)) maf[top_loci$variant_idx] else NULL
        res$top_loci = cbind(res$top_loci , top_loci[,-1])
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
            X_column_scale_factors = susie_output$X_column_scale_factors
        )
        class(res$susie_result_trimmed) <- "susie"
    } 
    return(res)
}
                           
                           
#' Post-process SuSiE_rss Analysis Results
#'
#' This function processes the results from SuSiE_rss. 
#' It extracts and processes various statistics and indices based on the provided SuSiE object and other parameters.
#' It has been updated to handle multiple secondary coverage values, performing secondary analysis for each and combining the top loci.
#' Also, it tracks cs_secondary correlation for each coverage value and includes this information in the top_loci table.
#'
#' @param susie_output Output from running susieR:susie_rss()
#' @param z Z score vector. (can be generated by beta/se)
#' @param maf Minor Allele Frequencies vector.
#' @param secondary_coverage Vector of coverage thresholds for secondary conditional analysis.
#' @param signal_cutoff Cutoff value for signal identification in PIP values. Default is 0.1.
#' @param other_quantities A list of other quantities to be added to the final object.
#' @return A list containing modified SuSiE object along with additional post-processing information.
#' @examples
#' # Example usage
#' # susie_result <- susie_post_processor(susie_output, z, Xcorr, maf, c(0.5, 0.7))
#' @importFrom susieR get_cs_correlation susie_get_cs
#' @importFrom stringr str_replace
#' @export                           
susie_rss_post_processor <- function(susie_output, z ,Xcorr, secondary_coverage = c(0.5, 0.7), maf,
                                    signal_cutoff = 0.1, other_quantities = list(), prior_eff_tol = 0){
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
        # it is okay to have PIP = NULL here
        c(which(susie_output$pip >= signal_cutoff), unlist(susie_output$sets$cs)) %>% unique %>% sort
    }
            
    get_cs_info <- function(susie_output_sets_cs, top_variants_idx) {
        cs_info_pri <- map_int(top_variants_idx, ~get_cs_index(.x, susie_output_sets_cs))
        ifelse(is.na(cs_info_pri), 0, as.numeric(str_replace(names(susie_output_sets_cs)[cs_info_pri], "L", "")))
    }
            
    get_cs_and_corr <- function(susie_output, coverage, Xcorr) {
        susie_output_secondary <- list(sets = susie_get_cs(susie_output, Xcorr = Xcorr, coverage = coverage), pip = susie_output$pip)
        susie_output_secondary$cs_corr <- get_cs_correlation(susie_output_secondary, Xcorr = Xcorr)
        susie_output_secondary
    }
            
    # Compute univariate regression results 
    res <- list(sumstats = z, 
                analysis_script = load_script(),
                other_quantities = other_quantities,
                variant_names = format_variant_id(names(susie_output$pip)),
                susie_result_trimmed = list()
            )
            
     eff_idx <- which(susie_output$V > prior_eff_tol)
            
    if (length(eff_idx) > 0) {
        # Prepare for top loci table
        top_variants_idx_pri <- get_top_variants_idx(susie_output, signal_cutoff)
        cs_pri <- get_cs_info(susie_output$sets$cs, top_variants_idx_pri)
        susie_output$cs_corr <- get_cs_correlation(susie_output, Xcorr = Xcorr)
        top_loci_list <- list("coverage_0.95" = data.frame(variant_idx = top_variants_idx_pri, cs_idx = cs_pri, stringsAsFactors=F))

         ## Loop over each secondary coverage value
        sets_secondary <- list()
        for (sec_cov in secondary_coverage) {
            sets_secondary[[paste0("coverage_", sec_cov)]] <- get_cs_and_corr(susie_output, sec_cov, Xcorr)
            top_variants_idx_sec <- get_top_variants_idx(sets_secondary[[paste0("coverage_", sec_cov)]], signal_cutoff)
            cs_sec <- get_cs_info(sets_secondary[[paste0("coverage_", sec_cov)]]$sets$cs, top_variants_idx_sec)
            top_loci_list[[paste0("coverage_", sec_cov)]] <-data.frame(variant_idx = top_variants_idx_sec, cs_idx = cs_sec, stringsAsFactors=F) 
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
        top_loci_cols <- c("variant_id", if (!is.null(maf)) "maf", "z", "pip", colnames(top_loci)[-1])
        res$top_loci <- data.frame(variants, 
                        z = z[top_loci$variant_idx],
                        pip, stringsAsFactors = FALSE)
        res$top_loci$maf = if (!is.null(maf)) maf[top_loci$variant_idx] else NULL
        res$top_loci  = cbind(res$top_loci,top_loci[,-1])
                           
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
            X_column_scale_factors = susie_output$X_column_scale_factors
        )
        class(res$susie_result_trimmed) <- "susie"
    } 
    return(res)
    
}