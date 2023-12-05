# requires mr.mash.alpha package.
# remotes::install_github("stephenslab/mr.mash.alpha")
#
#' @title mr.mash pipeline 
#' 
#' @description Take precomputed prior grid and mixture prior to compute weights with mr.mash
#' 
#' @param X n x p matrix of genotype, where n is the total number of individuals and p denotes the number of SNPs.
#' 
#' @param Y n x r standardized residual for expression, where n is the total number of individuals and 
#' r is the total number of conditions (tissue/cell-types). 
#' 
#' @param prior_matrices mixture prior computed with ultimate deconvolution model. It is a list object
#' with the following elements:
#' 
#' \item{U}{A list containing the estimated prior covariance matrices. }
#' 
#' \item{w}{A vector containing the estimated mixture weights. }
#' 
#' @param prior_grid a vector of scaling factors to be used in fitting mr.mash model. 
#' 
#' @param sample_partition a data frame that contains sample partition information for sample IDs and fold number as integer value
#' for a specific fold in cross validation, with column name "id" and "fold".  
#' 
#' @return A mr.mash fit, stored as a list with some or all of the
#' following elements:
#' 
#' \item{w1_colsums}{posterior weight sums for each condition. }
#'
#' \item{Bhat}{p x r weight matrix. }
#'
#' \item{Ytest}{subset of Y matrix as the testing data set for k fold cross validation. }
#'
#' \item{Ybar_train}{mean expression level across all conditions. }
#'
#' \item{elapsed_time}{computation runtime for fitting mr.mash. }
#'
#' \item{Yhat_test_glmnet}{predicted expression from a cross-validated glmnet model. }
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
#' \item{w0}{K-vector with (updated, if \code{update_w0=TRUE}) prior mixture weights, each associated with
#'   the respective covariance matrix in \code{S0}}.
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
#' \item{Y}{n x r matrix of responses at last iteration (only relevant when missing values
#'   are present in the input Y).}
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
#' Y <- matrix(runif(n*r,-2,2), nrow=n)
#' Y <- scale(Y)
#' colnames(Y) <- tissues
#' rownames(Y) <- sample_id 
#'
#' #sample partition for cross validation
#' id <- sample(sample_id, 200)
#' fold <- 1
#' sample_partition <- data.frame(id, fold)
#' 
#' 
#' set.seed(Sys.time())
#' components <- c('XtX','tFLASH_default','FLASH_default_1','FLASH_default_2','FLASH_default_3','FLASH_default_4','
#' FLASH_default_5','tFLASH_nonneg','FLASH_nonneg_1','FLASH_nonneg_2','FLASH_nonneg_3','PCA_1','PCA_2','PCA_3','tPCA')
#' tissues <- c('Adipose Tissue','Muscle Tissue', 'Brain Tissue','Liver Tissue','Kidney Tissue','Heart Tissue','Lung Tissue')
#' 
#' prior_matrices <- list()
#' for (i in components){
#'    A <- matrix(runif(r^2)*2-1, ncol=r) 
#'    cov <- t(A) %*% A
#'    colnames(cov) <- tissues
#'    rownames(cov) <- tissues
#'    prior_matrices[["U"]][[i]]<- list(mat=cov)
#'    }
#' prior_matrices$w <- runif(length(components),0.01,0.2)
#' 
#' #fit mr.mash with cross validation
#' res <- mr.mash.pipeline(X=X, Y=Y, sample_partition=sample_partition, prior_matrices=prior_matrices, prior_grid=grid))
#' 
#' #fit mr.mash without cross validation
#' res <- mr.mash.pipeline(X=X, Y=Y, prior_matrices=prior_matrices, prior_grid=prior_grid))
#' 
#' 
#' @importFrom mr.mash.alpha compute_V_init extract_missing_Y_pattern impute_missing_Y compute_canonical_covs mr.mash expand_covs
#' @importFrom matrixStats colVars rowSds
#' @importFrom glmnet cv.glmnet
#' @importFrom doMC registerDoMC





########## utility functions for mr.mash#########

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

  ## Make predictions if requested
  if (!is.null(Xnew)) {
    Yhat_glmnet <- drop(predict(cvfit_glmnet, newx=Xnew, s="lambda.min"))
    colnames(Yhat_glmnet) <- tissue_names
    res <- list(Bhat=B, Ytrain=Y, Yhat_new=Yhat_glmnet)
  } else {
    res <- list(Bhat=B, Ytrain=Y)
  }
  return(res)
}



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


### Filter data-driven matrices and summary stats based on tissues used
filter_datadriven_mats_and_sumstats <- function(Y, datadriven_mats, sumstats) {
  tissues_to_keep <- colnames(Y)

  if (!is.list(datadriven_mats$U[[1]])) {
    datadriven_mats_filt <- lapply(datadriven_mats$U, function(x, to_keep) {x[to_keep, to_keep]}, tissues_to_keep)
  } else {
    datadriven_mats_filt <- lapply(datadriven_mats$U, function(x, to_keep) {x$mat[to_keep, to_keep]}, tissues_to_keep)
  }

  if (!is.null(sumstats)) {
    sumstats_filt <- lapply(sumstats[[1]], function(x, to_keep) {x[, to_keep]}, tissues_to_keep)
  } else {
    sumstats_filt <- sumstats
  }

  return(list(datadriven_mats_filt=datadriven_mats_filt, sumstats_filt=sumstats_filt))
}

### Split the data into training and test
split_data <- function(X, Y, sample_folds, fold) {
  test_ids <- sample_folds[which(sample_folds$fold == fold), "id"]
  Xtrain <- X[!(rownames(X) %in% test_ids), ]
  Ytrain <- Y[!(rownames(Y) %in% test_ids), ]
  Xtest <- X[rownames(X) %in% test_ids, ]
  Ytest <- Y[rownames(Y) %in% test_ids, ]

  return(list(Xtrain=Xtrain, Ytrain=Ytrain, Xtest=Xtest, Ytest=Ytest))
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

### Compute column sums of the posterior assignment probabilities
compute_posterior_weight_colsum <- function(w, S0_labels) {
  w1_colsums <- colSums(w)
  if (!is.null(S0_labels)) {
    tmp <- rep(0, length(S0_labels))
    names(tmp) <- S0_labels
    tmp[which(S0_labels %in% names(w1_colsums))] <- w1_colsums
    w1_colsums <- tmp
  }
  return(w1_colsums)
}




# fit mr.mash
#' @importFrom mr.mash.alpha compute_V_init extract_missing_Y_pattern impute_missing_Y compute_canonical_covs mr.mash expand_covs
#' @importFrom matrixStats colVars rowSds
#' @importFrom glmnet cv.glmnet
#' @importFrom doMC registerDoMC
#' @export

mr.mash.pipeline <- function(X, 
                             Y, 
                             prior_matrices, 
                             prior_grid,
                             sample_partition=NULL, 
                             nthreads=2,
                             prior_weights=NULL,
                             canonical_mats = FALSE,
                             standardize=TRUE, 
                             update_w0 = TRUE,
                             w0_threshold = 1e-8,
                             update_V=TRUE, 
                             update_V_method = "full", 
                             B_init_method="enet",
                             max_iter=5000, 
                             tol = 0.01, 
                             verbose = FALSE, 
                             save_model=FALSE, 
                             glmnet_pred = TRUE,
                             seed=999){
    
  set.seed(seed)

  datadriven_mats <- prior_matrices
    
  # check if cross validation
  if(is.null(sample_partition)) {                
    message("sample_partition not provided, fitting mr.mash without cross validation. ")
    } else {
      fold <- unique(sample_partition$fold)
      if (length(fold) != 1) {
        message("Please provide sample IDs for one fold.")
      }
    }

  # computing prior weights
  w0_init <- prior_weights
  if (is.null(w0_init)) {
    message("Prior weights not provided. Computing them from initial estimates of the coefficients.")
  }

  ### Drop tissues with < n_nonmiss_Y in data-driven matrices and sumstats
  datadriven_mats_sumstats_filt <- filter_datadriven_mats_and_sumstats(Y=Y, datadriven_mats=datadriven_mats, sumstats=NULL)
  S0_data <- datadriven_mats_sumstats_filt$datadriven_mats_filt
  sumstats <- datadriven_mats_sumstats_filt$sumstats_filt
  rm(datadriven_mats_sumstats_filt)

  if (is.null(prior_grid)) {
    stop("Please provide prior grid.")
  }

  if (is.null(prior_matrices)) {
    stop("Please provide prior_matrices.")
  }

  if (!is.null(sample_partition)) {
    ### Split the data into training and test sets
    dat_split <- split_data(X=X, Y=Y, sample_folds=sample_partition, fold=fold)
    Xtrain <- dat_split$Xtrain
    Ytrain <- dat_split$Ytrain
    Xtest <- dat_split$Xtest
    Ytest <- dat_split$Ytest
    rm(dat_split)
  }

  ### Compute canonical matrices, if requested
  if (canonical_mats) {
    S0_can <- compute_canonical_covs(ncol(Ytrain), singletons=TRUE, hetgrid=c(0, 0.25, 0.5, 0.75, 1))
    S0_raw <- c(S0_can, S0_data)
  } else {
    S0_raw <- S0_data
  } 
  ### Compute prior covariance
  S0 <- expand_covs(S0_raw, prior_grid, zeromat=TRUE)
  time1 <- proc.time()

  ### Compute initial estimates of regression coefficients and prior weights (if not provided)
  if (glmnet_pred) { 
    if (!is.null(sample_partition)) {
      Xnew <- Xtest
    } else {
      Xnew <- X
    } # end glmnet_pred = TRUE
  } else {
    Xnew <- NULL
  }

  if (B_init_method == "enet") {
    if (!is.null(sample_partition)) {
      out <- compute_coefficients_univ_glmnet(Xtrain, Ytrain, alpha=0.5, standardize=standardize, nthreads=nthreads, Xnew=Xnew)
    } else {
      out <- compute_coefficients_univ_glmnet(X, Y, alpha=0.5, standardize=standardize, nthreads=nthreads, Xnew=Xnew)
    }
  } else if (B_init_method == "glasso") {
    if (!is.null(sample_partition)) {
      out <- compute_coefficients_glasso(Xtrain, Ytrain, standardize=standardize, nthreads=nthreads, Xnew=Xnew)
    } else {
      out <- compute_coefficients_glasso(X, Y, standardize=standardize, nthreads=nthreads, Xnew=Xnew)
    }
  }

  B_init <- out$Bhat

  if (is.null(w0_init)) {
    external_w0 <- FALSE
    w0 <- compute_w0(B_init, length(S0))
  } else {
    external_w0 <- TRUE
    w0 <- w0_init
  }

  ### Filter prior components based on weights
  comps_filtered <- filter_S0_w0(S0=S0, w0=w0)
  S0 <- comps_filtered$S0
  w0 <- comps_filtered$w0
  rm(comps_filtered)

  ### Fit mr.mash
  if (is.null(sample_partition)) {
    Xtrain = X 
    Ytrain = Y
  } # For no CV

  fit_mrmash <- mr.mash(X=Xtrain, Y=Ytrain, S0=S0, w0=w0, update_w0=update_w0, tol=tol,
                        max_iter=max_iter, convergence_criterion="ELBO", compute_ELBO=TRUE,
                        standardize=standardize, verbose=verbose, update_V=update_V, update_V_method=update_V_method,
                        w0_threshold=w0_threshold, nthreads=nthreads, mu1_init=B_init)

  time2 <- proc.time()
  elapsed_time <- time2["elapsed"] - time1["elapsed"]
  
  ### Make predictions
  if (!is.null(sample_partition)) {
    Yhat_test <- predict(fit_mrmash, Xtest)
  }
  
  ### Save results
  if (external_w0) {
    if (!is.null(sample_partition)) {
      resu <- list(Ytest=Ytest, Yhat_test=Yhat_test, analysis_time=elapsed_time)
    }
  } else {
    if (w0_threshold > 0) {
      S0_labels = names(S0)
    } else {
      S0_labels = NULL
    } 
    w1_colsums <- compute_posterior_weight_colsum(fit_mrmash$w1, S0_labels)
    if (!is.null(sample_partition)) {
      resu <- list(w1_colsums=w1_colsums, Bhat=fit_mrmash$mu1, Ytest=Ytest, Yhat_test=Yhat_test, 
                   Ybar_train=colMeans(Ytrain, na.rm=TRUE), analysis_time=elapsed_time)
    } else {
      resu <- fit_mrmash
      resu$w1_colsums <- w1_colsums
      resu$analysis_time <- elapsed_time
    }
  }  

  if (save_model) {
    resu$model <-  fit_mrmash
  }

  if (!is.null(sample_partition)) {
    if (glmnet_pred) {
      resu$Yhat_test_glmnet <- out$Yhat_new
    }
  }
  return(resu)
}
