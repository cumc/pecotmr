#' @importFrom susieR get_cs_correlation univariate_regression susie_get_cs
post_process_susie <- function(fobj, fdat, r, secondary_coverage = 0.7, signal_cutoff = 0.1) {
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
        fobj$sets_secondary = susie_get_cs(fobj, fdat$residual_X_scaled[[r]], coverage=secondary_coverage)
        fobj$analysis_script = load_script()
        fobj$cs_corr = get_cs_correlation(fobj, X=fdat$residual_X_scaled[[r]])
        fobj$cs_snps = gsub("_",":",names(fobj$pip[unlist(fobj$sets$cs)]))
        fobj$phenotype_name = colnames(fdat$residual_Y_scaled[[r]])
        fobj$dropped_samples = fdat$dropped_sample[[r]]
        fobj$sample_names = rownames(fdat$residual_Y_scaled[[r]])
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
        maf = fdat$maf[[r]][variants_merge]
        variants = gsub("_",":",names(fobj$pip)[variants_merge])
        pip = fobj$pip[variants_merge]
        cs_info_pri = map_int(variants_index, ~get_cs_index(.x, fobj$sets$cs))
        cs_info_sec = map_int(variants_index_secondary, ~get_cs_index(.x, fobj$sets_secondary$cs))
        cs_pri= ifelse(is.na(cs_info_pri), 0, str_replace(names(fobj$sets$cs)[cs_info_pri], "L", "") %>% as.numeric)
        cs_sec= ifelse(is.na(cs_info_sec), 0, str_replace(names(fobj$sets_secondary$cs)[cs_info_sec], "L", "") %>% as.numeric)
        cs_index_primary = cs_index_secondary = rep(NA, length(variants_merge))
        cs_index_primary[match(variants_index,variants_merge)]=cs_pri
        cs_index_secondary[match(variants_index_secondary,variants_merge)]=cs_sec
        Y_resid_sd = fdat$residual_Y_sd[[r]]
        X_resid_sd = fdat$residual_X_sd[[r]][variants_merge]
        univariate_res = univariate_regression(fdat$residual_X_scaled[[r]][, variants_merge, drop=F], fdat$residual_Y_scaled[[r]])
        fobj$top_loci = data.frame(variants, maf, univariate_res$betahat*Y_resid_sd/X_resid_sd, univariate_res$sebetahat*Y_resid_sd/X_resid_sd, pip, cs_index_primary,cs_index_secondary)
        colnames(fobj$top_loci) = c("variant_id", "maf", "bhat", "sbhat", "pip", "cs_index_primary","cs_index_secondary")
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
        class(fobj) = "list"
        # generate weights for TWAS using some alternative approaches --- this is not exactly fine-mapping but makes sense to do it here for data production
      fobj$susie_weights = susie_weights(fobj)
      fobj$susie_r2 = cor(fdat$residual_X_scaled[[r]] %*% fobj$susie_weights, fdat$residual_Y_scaled[[r]])^2
      fobj$enet_weights = glmnet_weights(fdat$residual_X_scaled[[r]], fdat$residual_Y_scaled[[r]])
      fobj$enet_r2 = cor(fdat$residual_X_scaled[[r]] %*% fobj$enet_weights, fdat$residual_Y_scaled[[r]])^2
      fobj$lasso_weights = glmnet_weights(fdat$residual_X_scaled[[r]], fdat$residual_Y_scaled[[r]], alpha = 1)
      fobj$lasso_r2 = cor(fdat$residual_X_scaled[[r]] %*% fobj$lasso_weights, fdat$residual_Y_scaled[[r]])^2
      fobj$mr_ash_weights = mr_ash_weights(fdat$residual_X_scaled[[r]], fdat$residual_Y_scaled[[r]], beta.init=fobj$lasso_weights)
      fobj$mr_ash_r2 = cor(fdat$residual_X_scaled[[r]] %*% fobj$mr_ash_weights, fdat$residual_Y_scaled[[r]])^2

    } else {
        fobj = list(analysis_script = load_script(), pip = fobj$pip, variant_names = gsub("_",":",names(fobj$pip)))
        names(fobj$pip) = NULL
    }
    return(fobj)
}