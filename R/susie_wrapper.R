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
            res$analysis_time <- proc.time() - st
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
#' @param fobj A list representing the SuSiE analysis object.
#' @param X_data The genotype data matrix.
#' @param y_data The phenotype data vector.
#' @param X_scalar Scalar for the genotype data, used in residual scaling.
#' @param y_scalar Scalar for the phenotype data, used in residual scaling.
#' @param maf Minor Allele Frequencies vector.
#' @param secondary_coverage Vector of coverage thresholds for secondary conditional analysis.
#' @param signal_cutoff Cutoff value for signal identification in PIP values. Default is 0.1.
#' @param other_quantities A list of other quantities to be added to the final object.
#' @return A modified SuSiE object with additional post-processing information.
#' @examples
#' # Example usage
#' # susie_result <- susie_post_processor(fobj, X_data, y_data, X_scalar, y_scalar, maf, c(0.5, 0.7))
#' @export
susie_post_processor <- function(fobj, X_data, y_data, X_scalar, y_scalar, maf, 
                                 secondary_coverage = c(0.5, 0.7), signal_cutoff = 0.1, 
                                 other_quantities = list()) {
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

    # Compute univariate regression results 
    res <- list(sumstats = univariate_regression(X_data, y_data), other_quantities = other_quantities)
    eff_idx <- which(fobj$V > 0)
    if (length(eff_idx) > 0) {
        fobj$analysis_script <- load_script()
        fobj$cs_corr <- get_cs_correlation(fobj, X = X_data)
        top_variants_idx <- c(which(fobj$pip >= signal_cutoff), unlist(fobj$sets$cs)) %>% unique %>% sort

        fobj$cs_secondary_corr <- vector("list", length(secondary_coverage))
        names(fobj$cs_secondary_corr) <- as.character(secondary_coverage)
        fobj$sets_secondary <- vector("list", length(secondary_coverage))
        names(fobj$sets_secondary) <- as.character(secondary_coverage)

        ## Loop over each secondary coverage value
        for (sec_cov in secondary_coverage) {
            fobj_secondary <- list(sets = susie_get_cs(fobj, X_data, coverage = sec_cov))
            fobj$cs_secondary_corr[[as.character(sec_cov)]] <- get_cs_correlation(fobj_secondary, X = X_data)
            ## Store secondary sets for each coverage value and merge top loci
            fobj$sets_secondary[[as.character(sec_cov)]] <- fobj_secondary$sets
            variants_index_secondary <- c(which(fobj$pip >= signal_cutoff), unlist(fobj$sets_secondary[[as.character(sec_cov)]]$cs)) %>% unique %>% sort
            top_variants_idx <- union(top_variants_idx, variants_index_secondary)
        }

        fobj$phenotype_name <- colnames(y_data)
        fobj$sample_names <- rownames(y_data)
        fobj$variant_names <- format_variant_id(names(fobj$pip))
        ## Prepare for the top loci table
        variants <- format_variant_id(names(fobj$pip)[top_variants_idx])
        pip <- fobj$pip[top_variants_idx]

        cs_info_pri <- map_int(top_variants_idx, ~get_cs_index(.x, fobj$sets$cs))
        cs_pri <- if (is.na(cs_info_pri)) 0 else as.numeric(str_replace(names(fobj$sets$cs)[cs_info_pri], "L", ""))

        ## Compute secondary CS information
        cs_secondary_info <- matrix(NA_integer_, nrow = length(top_variants_idx), ncol = length(secondary_coverage))
        colnames(cs_secondary_info) <- paste0("cs_secondary_", as.character(secondary_coverage))

        for (i in seq_along(secondary_coverage)) {
            sec_cov_name <- as.character(secondary_coverage[i])
            for (variant_idx in top_variants_idx) {
                cs_secondary_info[variant_idx, i] <- get_cs_index(variant_idx, fobj$sets_secondary[[sec_cov_name]]$cs)
            }
        }

        Y_resid_scalar <- if (!is.null(y_scalar)) y_scalar else 1
        X_resid_scalar <- if (!is.null(X_scalar)) X_scalar[top_variants_idx] else 1

        top_loci_cols <- c("variant_id", "bhat", "sbhat", "pip", "cs_index_primary", colnames(cs_secondary_info))
        fobj$top_loci <- data.frame(variants, maf = if (!is.null(maf)) maf[top_variants_idx] else NULL, fobj$sumstats$betahat[top_variants_idx] * Y_resid_scalar / X_resid_scalar, fobj$sumstats$sebetahat[top_variants_idx] * Y_resid_scalar / X_resid_scalar, pip, cs_info_pri, cs_secondary_info, stringsAsFactors = FALSE)
        colnames(fobj$top_loci) <- top_loci_cols

        res$susie_result_trimmed <- list(
            phenotype_name = fobj$phenotype_name,
            sample_names = fobj$sample_names,
            variant_names = fobj$variant_names,
            pip = fobj$pip
            sets = fobj$sets,
            cs_corr = fobj$cs_corr,
            cs_secondary_corr = fobj$cs_secondary_corr,
            sets_secondary = fobj$sets_secondary,
            alpha = fobj$alpha[eff_idx, , drop = FALSE],
            mu = fobj$mu[eff_idx, , drop = FALSE],
            mu2 = fobj$mu2[eff_idx, , drop = FALSE],
            V = fobj$V[eff_idx],
            X_column_scale_factors = fobj$X_column_scale_factors,
        )
        class(res$susie_result_trimmed) <- "susie"
    } else {
        res <- list(analysis_script = load_script(), pip = fobj$pip, variant_names = format_variant_id(names(fobj$pip)), sumstats = fobj$sumstats)
        names(res$pip) <- NULL
    }
    return(fobj)
}