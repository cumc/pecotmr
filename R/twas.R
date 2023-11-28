# FIXME: check that weights have the same length as z-scores
#' @importFrom Rfast cora
#' @export
twas_z <- function(weights, z, R=NULL, X=NULL) {
    if (is.null(R)) {
        # mean impute X
        genetype_data_imputed <- apply(X, 2, function(x){
            pos <- which(is.na(x))
            if (length(pos) != 0){
                x[pos] <- mean(x,na.rm = TRUE)
                }
                return(x)
            })
        # need to add `large = T` in Rfast::cora, or could get wrong results 
        # R <- cora(genetype_data_imputed, large = T)
        # FIXME: test and enable using cora if Rfast is installed
        # See how susieR did it
        R <- cor(genetype_data_imputed) 
        colnames(R) <- rownames(R) <- colnames(genetype_data_imputed)
    }
    stat <- t(weights) %*% z
    denom <- t(weights) %*% R %*% weights
    zscore <- stat/sqrt(denom)
    pval <- pchisq( zscore * zscore, 1, lower.tail = FALSE)
    return(list(z=zscore, pval=pval))
}

#' @importFrom susieR coef.susie
#' @export
susie_weights <- function(susie_fit) {
    coef.susie(susie_fit)[-1]
}

# Get a reasonable setting for the standard deviations of the mixture
# components in the mixture-of-normals prior based on the data (X, y).
# Input se is an estimate of the residual *variance*, and n is the
# number of standard deviations to return. This code is adapted from
# the autoselect.mixsd function in the ashr package.
#' @importFrom susieR univariate_regression
init_prior_sd <- function (X, y, n = 30) {
  res <- univariate_regression(X, y)
  smax <- 3*max(res$betahat)
  seq(0, smax, length.out = n)
}

#' @importFrom glmnet cv.glmnet
#' @importFrom stats coef
#' @export
glmnet_weights <- function(X, y, alpha=0.5) {
	eff.wgt = matrix(0, ncol=1, nrow=ncol(X))
	sds = apply(X, 2, sd)
	keep = sds != 0 & !is.na(sds)
	enet = cv.glmnet(x=X[,keep], y=y, alpha=alpha, nfold=5, intercept=T, standardize=F)
	eff.wgt[keep] = coef( enet , s = "lambda.min")[2:(sum(keep)+1)]
	return(eff.wgt)
}

#' @examples 
#' wgt.lasso = glmnet_weights(X, y, alpha=1)
#' wgt.mr.ash = mr_ash_weights(eqtl$X, eqtl$y_res, beta.init=wgt.lasso)
#' @importFrom mr.ash.alpha mr.ash
#' @importFrom stats predict
#' @export
mr_ash_weights <- function(X, y, init_prior_sd=TRUE, ...) {
    sa2 = NULL
    if (init_prior_sd) sa2 = init_prior_sd(X,y)^2
    fit.mr.ash = mr.ash(X, y, sa2=sa2, ...)
    predict(fit.mr.ash, type = "coefficients")[-1]
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

#' @importFrom harmonicmeanp pLandau
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
                           
           
#utmost_twas_z
## the utmost paper assumes X are not standardized, in the formula below - the input of X is assumed to be standardized 
library(GBJ)
           
twas_joint_z <- function(ld, Bhat, gwas_z){
    idx <- which(rownames(ld) %in% rownames(Bhat))
    D <- ld[idx,idx]
    
    cov_y <- crossprod(Bhat, D) %*% Bhat
    y_sd <- sqrt(diag(cov_y))
    x_sd <- rep(1, nrow(bhat))   #we assume X is standardized
    
    #get gamma matrix MxM (snp x snp) 
    g <- lapply(colnames(Bhat), function(x){
        gm <- diag(x_sd/y_sd[x], length(x_sd), length(x_sd))
        return(gm)
        })
        names(g) <- colnames(Bhat)

    ######### Get TWAS - Z statistics & P-value, GBJ test ########  
    z <- do.call(rbind, lapply(colnames(Bhat), function(x){
            Zi <- crossprod(Bhat[,x], g[[x]]) %*% as.numeric(gwas_z[,"Z"])
            pval <- 2*pnorm(abs(Zi), lower.tail=FALSE)
            Zp <- c(Zi, pval)
            names(Zp) <- c("Z", "pval")
            return(Zp)}))
    rownames(z) = colnames(Bhat) 

    # GBJ test 
    lam <- matrix(rep(NA,ncol(Bhat)*nrow(Bhat)), nrow = ncol(Bhat))
        rownames(lam) <- colnames(Bhat)
        for (p in colnames(Bhat)) {
            la <- as.matrix(Bhat[,p] %*% g[[p]])
            lam[p, ] <-  la
            }
    
    sig <- tcrossprod((lam %*% D), lam)
    gbj <- GBJ(test_stats=z[,1], cor_mat=sig)
    rs <- list("Z" =z, "GBJ"=gbj)
    return(rs)
}