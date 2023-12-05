# requires mr.mash.alpha package.
# remotes::install_github("stephenslab/mr.mash.alpha")
#
#' @title mr.mash wrapper
#' 
#' @description Take precomputed prior grid and mixture prior to compute weights with mr.mash
#' 
#' @param X n x p matrix of genotype, where n is the total number of individuals and p denotes
#' the number of SNPs.
#' 
#' @param Y n x r matrix of residual for expression, where n is the total number of individuals 
#' and r is the total number of conditions (tissue/cell-types). 
#' 
#' @param prior_data_driven_matrices A list of data-driven covariance matrices.
#' 
#' @param prior_grid a vector of scaling factors to be used in fitting mr.mash model. 
#'
#' 
#' @return A mr.mash fit, stored as a list with some or all of the
#' following elements:
#' 
#' \item{mu1}{p x r matrix of posterior means for the regression
#'   coeffcients. }
#' 
#' \item{S1}{r x r x p array of posterior covariances for the
#'   regression coeffcients. }
#' 
#' \item{w1}{p x K matrix of posterior assignment probabilities to the
#'   mixture components. }
#'   
#' \item{V}{r x r residual covariance matrix. }
#' 
#' \item{w0}{K-vector with (updated, if \code{update_w0=TRUE}) prior mixture weights, each 
#'  associated with the respective covariance matrix in \code{S0}}.
#'   
#' \item{S0}{r x r x K array of prior covariance matrices
#'   on the regression coefficients. }
#' 
#' \item{intercept}{r-vector containing posterior mean estimate of the
#'   intercept.}
#' 
#' \item{fitted}{n x r matrix of fitted values.}
#' 
#' \item{G}{r x r covariance matrix of fitted values.}
#' 
#' \item{pve}{r-vector of proportion of variance explained by the covariates.}
#' 
#' \item{ELBO}{Evidence Lower Bound (ELBO) at last iteration.}
#' 
#' \item{progress}{A data frame including information regarding
#'   convergence criteria at each iteration.}
#'   
#' \item{converged}{\code{TRUE} or \code{FALSE}, indicating whether
#'   the optimization algorithm converged to a solution within the chosen tolerance
#'   level.}
#' 
#' \item{elapsed_time}{computation runtime for fitting mr.mash. }
#'   
#' \item{Y}{n x r matrix of responses at last iteration (only relevant when missing values
#'   are present in the input Y).}
#' 
#' 
#' @examples 
#' ###Set seed
#' set.seed(123)
#' prior_grid <- runif(17,0.00005,0.05)
#'
#' #simulate X and Y, sample id, snp id
#' sample_id <- paste0("P000", str_pad(1:400, 3, pad = "0"))
#' X<- matrix(sample(0:2, size=n*p, replace=TRUE, prob=c(.65,.30,.05)), nrow=n)
#' rownames(X) <- sample_id
#' colnames(X) <- paste0("rs", sample(10000:100000, p)) #snp names
#' tissues <- c('Adipose Tissue','Muscle Tissue', 'Brain Tissue','Liver Tissue','Kidney Tissue',
#' 'Heart Tissue','Lung Tissue')
#' Y <- matrix(runif(n*r,-2,2), nrow=n)
#' Y <- scale(Y)
#' colnames(Y) <- tissues
#' rownames(Y) <- sample_id 
#' 
#' 
#' set.seed(Sys.time())
#' components <- c('XtX','tFLASH_default','FLASH_default','tFLASH_nonneg','FLASH_nonneg','PCA')
#' 
#' prior_data_driven_matrices <- list()
#' for (i in components){
#'    A <- matrix(runif(r^2)*2-1, ncol=r) 
#'    cov <- t(A) %*% A
#'    colnames(cov) <- tissues
#'    rownames(cov) <- tissues
#'    prior_data_driven_matrices[[i]]<- cov
#'    }
#' 
#' #fit mr.mash without cross validation
#' res <- mrmash_wrapper(X=X, Y=Y, prior_data_driven_matrices=prior_data_driven_matrices, prior_grid=prior_grid)
#' 
#' 
#' @importFrom mr.mash.alpha compute_canonical_covs mr.mash expand_covs
#' 
#' @export
#'

mrmash_wrapper <- function(X, 
                           Y, 
                           prior_data_driven_matrices, 
                           prior_grid, 
                           nthreads=2,
                           prior_canonical_matrices = FALSE,
                           standardize=TRUE, 
                           update_w0 = TRUE,
                           update_V=TRUE, 
                           update_V_method = "full", 
                           B_init_method="enet",
                           max_iter=5000, 
                           tol = 0.01, 
                           verbose = FALSE){

  if (is.null(prior_grid)) {
    stop("Please provide prior grid.")
  }

  if (is.null(prior_data_driven_matrices) && !isTRUE(prior_canonical_matrices)) {
    stop("Please provide prior_data_driven_matrices or set prior_canonical_matrices=TRUE.")
  }  
    
    
  ### Compute canonical matrices, if requested
  if (isTRUE(prior_canonical_matrices)) {
    prior_canonical_matrices <- compute_canonical_covs(ncol(Y), singletons=TRUE, hetgrid=c(0, 0.25, 0.5, 0.75, 1))
    if (!is.null(prior_data_driven_matrices)) {
      S0_raw <- c(prior_canonical_matrices, prior_data_driven_matrices)
    } else {
      S0_raw <- prior_canonical_matrices 
    }
  } else {
    S0_raw <- prior_data_driven_matrices
  }
    
    
  ### Compute prior covariance
  S0 <- expand_covs(S0_raw, prior_grid, zeromat=TRUE)
  time1 <- proc.time()


  if (B_init_method == "enet") {
      out <- compute_coefficients_univ_glmnet(X, Y, alpha=0.5, standardize=standardize, nthreads=nthreads, Xnew=NULL)
  } else if (B_init_method == "glasso") {
      out <- compute_coefficients_glasso(X, Y, standardize=standardize, nthreads=nthreads, Xnew=NULL)
  }

    
  B_init <- out$Bhat
  w0 <- compute_w0(B_init, length(S0))

    
  ### Filter prior components based on weights
  comps_filtered <- filter_S0_w0(S0=S0, w0=w0)
  S0 <- comps_filtered$S0
  w0 <- comps_filtered$w0
  rm(comps_filtered)

  ### Fit mr.mash
  fit_mrmash <- mr.mash(X=X, Y=Y, S0=S0, w0=w0, update_w0=update_w0, tol=tol,
                        max_iter=max_iter, convergence_criterion="ELBO", compute_ELBO=TRUE,
                        standardize=standardize, verbose=verbose, update_V=update_V, update_V_method=update_V_method,
                        nthreads=nthreads, mu1_init=B_init)

  time2 <- proc.time()
  elapsed_time <- time2["elapsed"] - time1["elapsed"]
  fit_mrmash$analysis_time <- elapsed_time

  return(fit_mrmash)
    
}





########## utility functions for mr.mash#########

#' @importFrom doMC registerDoMC
#' @importFrom mr.mash.alpha extract_missing_Y_pattern compute_V_init impute_missing_Y
#' @importFrom glmnet cv.glmnet
#'
### Function to compute initial estimates of the coefficients from group-lasso
compute_coefficients_glasso <- function(X, Y, standardize, nthreads, Xnew=NULL, version=c("Rcpp", "R")) {
  version <- match.arg(version)
  n <- nrow(X)
  p <- ncol(X)
  r <- ncol(Y)
  Y_has_missing <- any(is.na(Y))
  tissue_names <- colnames(Y)

  if (Y_has_missing) {
    ### Extract per-individual Y missingness patterns
    Y_miss_patterns <- extract_missing_Y_pattern(Y)

    ### Compute V and its inverse
    V <- compute_V_init(X, Y, matrix(0, p, r), method="flash")
    Vinv <- chol2inv(chol(V))

    ### Initialize missing Ys
    muy <- colMeans(Y, na.rm=TRUE)
    for (l in 1:r) {
      Y[is.na(Y[, l]), l] <- muy[l]
    }

    ### Compute expected Y (assuming B=0)
    mu <- matrix(rep(muy, each=n), n, r)

    ### Impute missing Ys
    Y <- impute_missing_Y(Y=Y, mu=mu, Vinv=Vinv, miss=Y_miss_patterns$miss, non_miss=Y_miss_patterns$non_miss,
                          version=version)$Y
  }

  ## Fit group-lasso
  if (nthreads > 1) {
    registerDoMC(nthreads)
    paral <- TRUE
  } else {
    paral <- FALSE
  }

  cvfit_glmnet <- cv.glmnet(x=X, y=Y, family="mgaussian", alpha=1, standardize=standardize, parallel=paral)
  coeff_glmnet <- coef(cvfit_glmnet, s="lambda.min")

  ## Build matrix of initial estimates for mr.mash
  B <- matrix(as.numeric(NA), nrow=p, ncol=r)

  for (i in 1:length(coeff_glmnet)) {
    B[, i] <- as.vector(coeff_glmnet[[i]])[-1]
  }

  ## Make predictions if requested. (In compute_coefficients_glasso() cuntion)
  if (!is.null(Xnew)) {
    Yhat_glmnet <- drop(predict(cvfit_glmnet, newx=Xnew, s="lambda.min"))
    colnames(Yhat_glmnet) <- tissue_names
    res <- list(Bhat=B, Ytrain=Y, Yhat_new=Yhat_glmnet)
  } else {
    res <- list(Bhat=B, Ytrain=Y)
  }
  return(res)
}


#' @importFrom doMC registerDoMC
#' @importFrom glmnet cv.glmnet
#'
### Function to compute coefficients for univariate glmnet

compute_coefficients_univ_glmnet <- function(X, Y, alpha, standardize, nthreads, Xnew=NULL) {
  r <- ncol(Y)

  linreg <- function(i, X, Y, alpha, standardize, nthreads, Xnew) {
    if (nthreads > 1) {
      registerDoMC(nthreads)
      paral <- TRUE
    } else {
      paral <- FALSE
    }

    samples_kept <- which(!is.na(Y[, i]))
    Ynomiss <- Y[samples_kept, i]
    Xnomiss <- X[samples_kept, ]

    cvfit <- cv.glmnet(x=Xnomiss, y=Ynomiss, family="gaussian", alpha=alpha, standardize=standardize, parallel=paral)
    coeffic <- as.vector(coef(cvfit, s="lambda.min"))
    lambda_seq <- cvfit$lambda

    ## Make predictions if requested
    if (!is.null(Xnew)) {
      yhat_glmnet <- drop(predict(cvfit, newx=Xnew, s="lambda.min"))
      res <- list(bhat=coeffic, lambda_seq=lambda_seq, yhat_new=yhat_glmnet)
    } else {
      res <- list(bhat=coeffic, lambda_seq=lambda_seq)
    }

    return(res)
  }

  out <- lapply(1:r, linreg, X, Y, alpha, standardize, nthreads, Xnew)

  Bhat <- sapply(out, "[[", "bhat")

  if (!is.null(Xnew)) {
    Yhat_new <- sapply(out, "[[", "yhat_new")
    colnames(Yhat_new) <- colnames(Y)
    results <- list(Bhat=Bhat[-1, ], intercept=Bhat[1, ], Yhat_new=Yhat_new)
  } else {
    results <- list(Bhat=Bhat[-1, ], intercept=Bhat[1, ])
  }
  return(results)
}



### Filter S0 and w0Drop mixture components with weight equal to 0
filter_S0_w0 <- function(S0, w0, thresh=.Machine$double.eps) {
  comps_to_keep <- which(w0 > thresh)
  S0 <- S0[comps_to_keep]
  w0 <- w0[comps_to_keep]

  return(list(S0=S0, w0=w0))
}



### Compute prior weights from coefficients estimates
compute_w0 <- function(Bhat, ncomps) {
  prop_nonzero <- sum(rowSums(abs(Bhat)) > 0) / nrow(Bhat)
  w0 <- c((1 - prop_nonzero), rep(prop_nonzero / (ncomps - 1), (ncomps - 1)))

  if (sum(w0 != 0) < 2) {
    w0 <- rep(1 / ncomps, ncomps)
  }

  return(w0)
}
