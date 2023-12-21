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
#'
#' @param fobj A list representing the SuSiE analysis object.
#' @param X_data The genotype data matrix.
#' @param y_data The phenotype data vector.
#' @param X_scalar Scalar for the genotype data, used in residual scaling.
#' @param y_scalar Scalar for the phenotype data, used in residual scaling.
#' @param maf Minor Allele Frequencies vector.
#' @param secondary_coverage Coverage threshold for secondary conditional analysis. Default is 0.5.
#' @param signal_cutoff Cutoff value for signal identification in PIP values. Default is 0.1.
#' @param other_quantities A list of other quantities to be added to the final object.
#' @return A modified SuSiE object with additional post-processing information.
#' @examples
#' # Example usage
#' # susie_result <- susie_post_processor(fobj, X_data, y_data, X_scalar, y_scalar, maf)
#' @export
susie_post_processor <- function(fobj, X_data, y_data, X_scalar, y_scalar, maf, 
                                 secondary_coverage = 0.5, signal_cutoff = 0.1, other_quantities = list()) {
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
    eff_idx <- which(fobj$V > 0)
    if (length(eff_idx) > 0) {
        fobj$analysis_script <- load_script()
        fobj$cs_corr <- get_cs_correlation(fobj, X = X_data)
        fobj$cs_snps <- format_variant_id(names(fobj$pip[unlist(fobj$sets$cs)]))
        fobj_secondary <- list(sets = susie_get_cs(fobj, X_data, coverage = secondary_coverage))
        fobj$cs_secondary_corr <- get_cs_correlation(fobj_secondary, X = X_data)
        fobj$sets_secondary <- fobj_secondary$sets
        fobj_secondary <- NULL
        fobj$phenotype_name <- colnames(y_data)
        fobj$sample_names <- rownames(y_data)
        fobj$variant_names <- format_variant_id(names(fobj$pip))
        variants_index <- c(which(fobj$pip >= signal_cutoff), unlist(fobj$sets$cs)) %>% unique %>% sort
        variants_index <- if (length(variants_index) == 0) which.max(fobj$pip) else variants_index
        variants_index_secondary <- c(which(fobj$pip >= signal_cutoff), unlist(fobj$sets_secondary$cs)) %>% unique %>% sort
        variants_index_secondary <- if (length(variants_index_secondary) == 0) which.max(fobj$pip) else variants_index_secondary
        variants_merge <- unique(c(variants_index, variants_index_secondary)) %>% sort
        variants <- format_variant_id(names(fobj$pip)[variants_merge])
        pip <- fobj$pip[variants_merge]
        cs_info_pri <- map_int(variants_index, ~get_cs_index(.x, fobj$sets$cs))
        cs_info_sec <- map_int(variants_index_secondary, ~get_cs_index(.x, fobj$sets_secondary$cs))
        cs_pri <- if (is.na(cs_info_pri)) 0 else as.numeric(str_replace(names(fobj$sets$cs)[cs_info_pri], "L", ""))
        cs_sec <- if (is.na(cs_info_sec)) 0 else as.numeric(str_replace(names(fobj$sets_secondary$cs)[cs_info_sec], "L", ""))
        cs_index_primary <- cs_index_secondary <- rep(NA, length(variants_merge))
        cs_index_primary[match(variants_index, variants_merge)] <- cs_pri
        cs_index_secondary[match(variants_index_secondary, variants_merge)] <- cs_sec
        Y_resid_scalar <- if (!is.null(y_scalar)) y_scalar else 1
        X_resid_scalar <- if (!is.null(X_scalar)) X_scalar[variants_merge] else 1
        fobj$sumstats <- univariate_regression(X_data, y_data)
        top_loci_cols <- c("variant_id", "bhat", "sbhat", "pip", "cs_index_primary", "cs_index_secondary")
        fobj$top_loci <- data.frame(variants, maf = if (!is.null(maf)) maf[variants_merge] else NULL, fobj$sumstats$betahat[variants_merge] * Y_resid_scalar / X_resid_scalar, fobj$sumstats$sebetahat[variants_merge] * Y_resid_scalar / X_resid_scalar, pip, cs_index_primary, cs_index_secondary, stringsAsFactors=F)
        colnames(fobj$top_loci) <- c("variant_id", if (!is.null(maf)) "maf" else NULL, "bhat", "sbhat", "pip", "cs_index_primary", "cs_index_secondary")
        rownames(fobj$top_loci) <- NULL
        fobj$alpha <- fobj$alpha[eff_idx, , drop = FALSE]
        fobj$mu <- fobj$mu[eff_idx, , drop = FALSE]
        fobj$mu2 <- fobj$mu2[eff_idx, , drop = FALSE]
        fobj$V <- fobj$V[eff_idx]
        fobj$Xr <- NULL
        fobj$fitted <- NULL
        fobj$lbf_variable <- NULL
        fobj$alpha <- NULL
        fobj$mu <- NULL
        fobj$mu2 <- NULL
        fobj$X_column_scale_factors <- NULL
        fobj$pip <- NULL
        for (item in names(other_quantities)) {
            fobj[[item]] <- other_quantities[[item]]
        }
    } else {
        fobj <- list(analysis_script = load_script(), pip = fobj$pip, variant_names = format_variant_id(names(fobj$pip)))
        names(fobj$pip) <- NULL
    }
    return(fobj)
}

