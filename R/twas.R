fine_twas <- function(weight_ids, gwas_ids, uber_ids, genotypeMatrix, ld_uber_ids, weights, zscores, return_stat = FALSE) {
        modifiers <- handle_weights(gwas_ids, weight_ids)
        weights <- modifiers * weights

        common_variants <- intersect(
            ld_uber_ids,
            uber_ids)
      
        if(length(common_variants) > 1){
      
        genotypeMatrix <- genotypeMatrix[, which(ld_uber_ids %in% common_variants)]
        genotypeMatrix <- data.matrix(genotypeMatrix)
        # pick the variants with MAC > 0, maybe change to MAC > 10 in the future
        genotypeMatrix <- genotypeMatrix[,which(colSums(genotypeMatrix,na.rm = T)>0)] 

        # mean impute X
        genetype_data_imputed <- apply(genotypeMatrix, 2, function(x){
            pos <- which(is.na(x))
            if (length(pos) != 0){
                x[pos] <- mean(x,na.rm = TRUE)
                }
                return(x)
            })
        # need to add `large = T` in Rfast::cora, or could get wrong results 
        ld_matrix <- Rfast::cora(genetype_data_imputed,large = T)
        colnames(ld_matrix) <- rownames(ld_matrix) <- colnames(genetype_data_imputed)
        } else {
        ld_matrix <- matrix(1)
        }
  
        stat <- t(weights) %*% zscores
        denom <- t(weights) %*% ld_matrix %*% weights
        zscore <- stat/sqrt(denom)

        if (length(zscores) == 1) {
            zscore <- zscores[[1]]
            if (weights[[1]] < 0) {
                zscore <- zscore * -1
            }
        }
        
        pval <- pchisq( zscore * zscore, 1, lower.tail = FALSE)
    
        #debug_print(gwas_ids, modifiers, weights, zscores, stat, denom, paste0("${_output:d}", "/ptwas-scan.debug"))

        return(if (return_stat) zscore else pval[1])
    }

pval_acat <- function(pvals) {
    if (length(pvals) == 1) {
        return(pvals[0])
    }
    stat <- 0.00
    pval_min <- 1.00

    stat <- sum(qcauchy(pvals))
    pval_min <- min(pval_min, min(qcauchy(pvals)))

    return(pcauchy(stat/length(pvals), lower.tail = FALSE))
}

library(harmonicmeanp)
pval_hmp <- function(pvals) {
	# https://search.r-project.org/CRAN/refmans/harmonicmeanp/html/pLandau.html
    pvalues <- unique(pvals)
	L <- length(pvalues)
	HMP <- L/sum(pvalues^-1)

	LOC_L1 <- 0.874367040387922
	SCALE <- 1.5707963267949

	return(pLandau(1/HMP, mu = log(L) + LOC_L1, sigma = SCALE, lower.tail = FALSE))
}

pval_global <- function(pvals, comb_method = "HMP", naive=FALSE) {
    # assuming sstats has tissues as columns and rows as pvals
    min_pval <- min(pvals)
    n_total_tests <- pvals %>% unique() %>% length() # There should be one unique pval per tissue
    global_pval <- if (comb_method == "HMP") pval_HMP(pvals) else pval_ACAT(pvals) # pval vector
    naive_pval <- min(n_total_tests*min_pval, 1.0)
    return(if (naive) naive_pval else global_pval) # global_pval and naive_pval
}