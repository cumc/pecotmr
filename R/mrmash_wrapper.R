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
#' @importFrom mr.mash.alpha compute_canonical_covs mr.mash expand_covs compute_univariate_sumstats
#' @importFrom foreach foreach %do%
#' @export
#'

mrmash_wrapper <- function(X, 
                           Y, 
                           prior_data_driven_matrices=NULL, 
                           prior_grid=NULL, 
                           nthreads=2,
                           prior_canonical_matrices = FALSE,
                           standardize=FALSE, 
                           update_w0 = TRUE,
                           w0_threshold = 1e-8,
                           update_V=TRUE, 
                           update_V_method = "full", 
                           B_init_method="enet",
                           max_iter=5000, 
                           tol = 0.01, 
                           verbose = FALSE){

  if (is.null(prior_data_driven_matrices) && !isTRUE(prior_canonical_matrices)) {
    stop("Please provide prior_data_driven_matrices or set prior_canonical_matrices=TRUE.")
  }  
  
    
  Y_has_missing <- any(is.na(Y))
    
  if(Y_has_missing && B_init_method=="glasso"){
      stop("B_init_method=glasso can only be used without missing values in Y")
  }
    
    
  # Compute summary statistics and prior_grids
  # in the case of cross validation, X input becomes Xtrain, Y becomes Ytrain
  sumstats <- compute_univariate_sumstats(X, Y, standardize=standardize, standardize.response=FALSE, mc.cores=nthreads)
  prior_grid <- compute_grid(bhat=sumstats$Bhat, sbhat=sumstats$Shat)

    
  ### Compute canonical matrices, if requested
  if (isTRUE(prior_canonical_matrices)) {
    prior_canonical_matrices <- compute_canonical_covs(ncol(Y), singletons=TRUE, hetgrid=c(0, 0.25, 0.5, 0.75, 1))
      
    if (!is.null(prior_data_driven_matrices)) {
      prior_data_driven_matrices <- filter_datadriven_mats(Y, prior_data_driven_matrices)
      S0_raw <- c(prior_canonical_matrices, prior_data_driven_matrices)  
    } else {
      S0_raw <- prior_canonical_matrices 
    }
      
  } else {
    S0_raw <- filter_datadriven_mats(Y, prior_data_driven_matrices)
  }
 
    
  ### Compute prior covariance 
  S0 <- expand_covs(S0_raw, prior_grid, zeromat=TRUE)
  time1 <- proc.time()
    

  if (B_init_method == "enet") {
      out <- compute_coefficients_univ_glmnet(X, Y, alpha=0.5, standardize=standardize, nthreads=nthreads, Xnew=NULL)
  } else if (B_init_method == "glasso") {
      out <- compute_coefficients_glasso(X, Y, standardize=standardize, nthreads=nthreads, Xnew=NULL)
  }
    
  B_init <- as.matrix(out$Bhat)
  w0 <- compute_w0(B_init, length(S0))


  ### Fit mr.mash
  fit_mrmash <- mr.mash(X=X, Y=Y, S0=S0, w0=w0, update_w0=update_w0, tol=tol,
                        max_iter=max_iter, convergence_criterion="ELBO", compute_ELBO=TRUE,
                        standardize=standardize, verbose=verbose, update_V=update_V, update_V_method=update_V_method,
                        w0_threshold = w0_threshold, nthreads=nthreads, mu1_init=B_init)

    
  time2 <- proc.time()
  elapsed_time <- time2["elapsed"] - time1["elapsed"]
  fit_mrmash$analysis_time <- elapsed_time

  return(fit_mrmash)
      
    
}





########## utility functions for mr.mash#########

#' @importFrom doMC registerDoMC
#' @importFrom glmnet cv.glmnet
#'
### Function to compute initial estimates of the coefficients from group-lasso
compute_coefficients_glasso <- function(X, Y, standardize, nthreads, Xnew=NULL, version=c("Rcpp", "R")) {
  n <- nrow(X)
  p <- ncol(X)
  r <- ncol(Y)
  tissue_names <- colnames(Y)

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

  ## Make predictions if requested. 
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




### Compute prior weights from coefficients estimates
compute_w0 <- function(Bhat, ncomps) {
  prop_nonzero <- sum(rowSums(abs(Bhat)) > 0) / nrow(Bhat)
  w0 <- c((1 - prop_nonzero), rep(prop_nonzero / (ncomps - 1), (ncomps - 1)))

  if (sum(w0 != 0) < 2) {
    w0 <- rep(1 / ncomps, ncomps)
  }

  return(w0)
}


###Filter data-driven matrices 
filter_datadriven_mats<- function(Y, datadriven_mats){
  tissues_to_keep <- colnames(Y)
  datadriven_mats_filt <- lapply(datadriven_mats, function(x, to_keep){x[to_keep, to_keep]}, tissues_to_keep)
  return(datadriven_mats_filt)
}



######### Function to compute grids##########
#compute grid
compute_grid <- function(bhat, sbhat){
    grid_mins = c()
    grid_maxs = c()

    #endpoints = compute_grid_endpoints(sumstat)
    include = !(sbhat==0 | !is.finite(sbhat) | is.na(bhat))
    gmax = grid_max(bhat[include], sbhat[include])
    gmin = grid_min(bhat[include], sbhat[include])

    grid_mins = c(grid_mins, gmin)
    grid_maxs = c(grid_maxs, gmax)

    gmin_tot = min(grid_mins)
    gmax_tot = max(grid_maxs)
    grid = autoselect_mixsd(gmin_tot, gmax_tot, mult=sqrt(2))^2  

    return(grid)        
}


###Compute the minimum value for the grid
grid_min = function(bhat,sbhat){
  min(sbhat)
}

###Compute the maximum value for the grid
grid_max = function(bhat,sbhat){
  if (all(bhat^2 <= sbhat^2)) {
    8 * grid_min(bhat,sbhat) # the unusual case where we don't need much grid
  } else {
    2 * sqrt(max(bhat^2 - sbhat^2))
  }
}
   

###Function to compute the grid
autoselect_mixsd <- function(gmin, gmax, mult=2){
  if (mult == 0) {
    return(c(0, gmax/2))
  }
  else {
    npoint = ceiling(log2(gmax/gmin)/log2(mult))
    return(mult^((-npoint):0) * gmax)
  }
}

