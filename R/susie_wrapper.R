#' @importFrom susieR get_cs_correlation univariate_regression susie_get_cs
#' @importFrom stringr str_replace
#' @export 
susie_post_processor <- function(fobj, X_data, y_data, X_sd, y_sd, maf, secondary_coverage = 0.5, signal_cutoff = 0.1, other_quantities=list()) {
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
    eff_idx = which(fobj$V>0)
    if (length(eff_idx)>0) {
        fobj$sets_secondary = susie_get_cs(fobj, X_data, coverage=secondary_coverage)
        fobj$analysis_script = load_script()
        fobj$cs_corr = get_cs_correlation(fobj, X=X_data)
        fobj$cs_snps = gsub("_",":",names(fobj$pip[unlist(fobj$sets$cs)]))
        fobj$phenotype_name = colnames(y_data)
        fobj$sample_names = rownames(y_data)
        fobj$variant_names = gsub("_",":",names(fobj$pip))
        variants_index = c(which(fobj$pip >= signal_cutoff), unlist(fobj$sets$cs)) %>% unique %>% sort
        if (length(variants_index)==0) {
            variants_index = which.max(fobj$pip)
        }
        variants_index_secondary = c(which(fobj$pip >= signal_cutoff), unlist(fobj$sets_secondary$cs)) %>% unique %>% sort
        if (length(variants_index_secondary)==0) {
            variants_index_secondary = which.max(fobj$pip)
        }
        variants_merge = unique(c(variants_index, variants_index_secondary))%>%sort
        variants = gsub("_",":",names(fobj$pip)[variants_merge])
        pip = fobj$pip[variants_merge]
        cs_info_pri = map_int(variants_index, ~get_cs_index(.x, fobj$sets$cs))
        cs_info_sec = map_int(variants_index_secondary, ~get_cs_index(.x, fobj$sets_secondary$cs))
        cs_pri= ifelse(is.na(cs_info_pri), 0, str_replace(names(fobj$sets$cs)[cs_info_pri], "L", "") %>% as.numeric)
        cs_sec= ifelse(is.na(cs_info_sec), 0, str_replace(names(fobj$sets_secondary$cs)[cs_info_sec], "L", "") %>% as.numeric)
        cs_index_primary = cs_index_secondary = rep(NA, length(variants_merge))
        cs_index_primary[match(variants_index,variants_merge)]=cs_pri
        cs_index_secondary[match(variants_index_secondary,variants_merge)]=cs_sec
        if (!is.null(y_sd)) {
            Y_resid_sd = fdat$residual_Y_sd[[r]]
        } else {
            Y_resid_sd = 1
        }
        if (!is.null(X_sd)) {
            X_resid_sd = fdat$residual_X_sd[[r]][variants_merge]
        } else {
            X_resid_sd = 1
        }
        univariate_res = univariate_regression(X_data[, variants_merge, drop=F], y_data)
        if (!is.null(maf)) {
            fobj$top_loci = data.frame(variants, maf[variants_merge], univariate_res$betahat*Y_resid_sd/X_resid_sd, univariate_res$sebetahat*Y_resid_sd/X_resid_sd, pip, cs_index_primary,cs_index_secondary)
            colnames(fobj$top_loci) = c("variant_id", "maf", "bhat", "sbhat", "pip", "cs_index_primary","cs_index_secondary")
        } else {
            fobj$top_loci = data.frame(variants, univariate_res$betahat*Y_resid_sd/X_resid_sd, univariate_res$sebetahat*Y_resid_sd/X_resid_sd, pip, cs_index_primary,cs_index_secondary)
            colnames(fobj$top_loci) = c("variant_id", "bhat", "sbhat", "pip", "cs_index_primary","cs_index_secondary")
        }
        rownames(fobj$top_loci) = NULL
        fobj$alpha = fobj$alpha[eff_idx,,drop=F]
        fobj$mu = fobj$mu[eff_idx,,drop=F]
        fobj$mu2 = fobj$mu2[eff_idx,,drop=F]
        fobj$V = fobj$V[eff_idx]
        fobj$Xr = NULL
        fobj$fitted = NULL
        colnames(fobj$lbf_variable) = NULL
        colnames(fobj$alpha) = NULL
        colnames(fobj$mu) = NULL
        colnames(fobj$mu2) = NULL
        names(fobj$X_column_scale_factors) = NULL
        names(fobj$pip) = NULL
        # copy other useful information
        for (item in names(other_quantities)) fobj[[item]] = other_quantities[[item]]
        class(fobj) = "list"
    } else {
        fobj = list(analysis_script = load_script(), pip = fobj$pip, variant_names = gsub("_",":",names(fobj$pip)))
        names(fobj$pip) = NULL
    }
    return(fobj)
}

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