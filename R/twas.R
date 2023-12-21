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
    pval <- pchisq(zscore * zscore, 1, lower.tail = FALSE)
    return(list(z=zscore, pval=pval))
}

#' Cross-Validation for Transcriptome-Wide Association Studies (TWAS)
#'
#' Performs cross-validation for TWAS, supporting both univariate and multivariate methods. 
#' It can either create folds for cross-validation or use pre-defined sample partitions. 
#' For multivariate methods, it applies the method to the entire Y matrix for each fold.
#'
#' @param X A matrix of samples by features, where each row represents a sample and each column a feature.
#' @param Y A matrix (or vector, which will be converted to a matrix) of samples by outcomes, where each row corresponds to a sample.
#' @param fold An optional integer specifying the number of folds for cross-validation. 
#' If NULL, 'sample_partitions' must be provided.
#' @param sample_partitions An optional dataframe with predefined sample partitions, 
#' containing columns 'Sample' (sample names) and 'Fold' (fold number). If NULL, 'fold' must be provided.
#' @param weight_methods A list of methods and their specific arguments, formatted as list(method1 = method1_args, method2 = method2_args). 
#' methods in the list can be either univariate (applied to each column of Y) or multivariate (applied to the entire Y matrix).
#' @param seed An optional integer to set the random seed for reproducibility of sample splitting.
#' @param num_threads The number of threads to use for parallel processing.
#'        If set to -1, the function uses all available cores.
#'        If set to 0 or 1, no parallel processing is performed.
#'        If set to 2 or more, parallel processing is enabled with that many threads.
#' @return A list with the following components:
#' \itemize{
#'   \item `sample_partition`: A dataframe showing the sample partitioning used in the cross-validation.
#'   \item `results`: A list of results for each method. Each list element contains:
#'     \itemize{
#'       \item `predictions`: A matrix of predicted Y values for each fold.
#'       \item `rsq_pval`: A matrix with rows representing methods and two columns for R-squared and p-value, calculated for each column of Y.
#'     }
#' }
#' @importFrom parallel makeCluster
#' @importFrom parallel detectCores
#' @importFrom parallel stopCluster
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @importFrom doParallel registerDoParallel
#' @export
twas_weights_cv <- function(X, Y, fold = NULL, sample_partitions = NULL, weight_methods = NULL, seed = NULL, num_threads = 1) {
    # Validation checks
    if (is.null(fold) || !is.numeric(fold) || fold <= 0) {
        stop("Invalid value for 'fold'. It must be a positive integer.")
    }
    if (!is.matrix(X) || (!is.matrix(Y) && !is.vector(Y))) {
        stop("X must be a matrix and Y must be a matrix or a vector.")
    }
    if (is.vector(Y)) {
        Y <- matrix(Y, ncol = 1)
    }
    if (nrow(X) != nrow(Y)) {
        stop("The number of rows in X and Y must be the same.")
    }

    # Get sample names
    if (!is.null(rownames(X))) {
        sample_names <- rownames(X)
    } else if (!is.null(rownames(Y))) {
        sample_names <- rownames(Y)
    } else {
        sample_names <- 1:nrow(X)
    }

    # Create or use provided folds
    if (!is.null(fold)) {
        if (!is.null(seed)) set.seed(seed)
        sample_indices <- sample(nrow(X))
        folds <- cut(seq(1, nrow(X)), breaks = fold, labels = FALSE)
        sample_partition <- data.frame(Sample = sample_names[sample_indices], Fold = folds, stringsAsFactors = FALSE)
    } else if (!is.null(sample_partitions)) {
        if (!all(sample_partitions$Sample %in% sample_names)) {
            stop("Some samples in 'sample_partitions' do not match the samples in 'X' and 'Y'.")
        }
        folds <- sample_partitions$Fold
        sample_partition <- sample_partitions
    } else {
        stop("Either 'fold' or 'sample_partitions' must be provided.")
    }

    if (is.null(weight_methods)) {
        return(list(sample_partition = sample_partition))
    } else {
        # Hardcoded vector of multivariate weight_methods
        multivariate_weight_methods <- c('mrmash_weights')

        # Determine the number of cores to use
        num_cores <- ifelse(num_threads == -1, detectCores(), num_threads)
        # Perform CV with parallel processing
        process_method <- function(j) {
            fold_indices <- which(folds == j)
            training_indices <- setdiff(1:nrow(X), fold_indices)
            X_train <- X[training_indices, ]
            Y_train <- Y[training_indices, ,drop=F]
            X_test <- X[fold_indices, ]

            # Remove columns with zero standard error
            valid_columns <- apply(X_train, 2, function(col) sd(col) != 0)
            X_train <- X_train[, valid_columns, drop=F]

            setNames(lapply(names(weight_methods), function(method) {
                args <- weight_methods[[method]]
                if (method %in% multivariate_weight_methods) {
                    # Apply multivariate method to entire Y for this fold
                    weights_matrix <- do.call(method, c(list(X = X_train, Y = Y_train), args))
                    # Adjust the weights matrix to include zeros for invalid columns
                    full_weights_matrix <- matrix(0, nrow = ncol(X), ncol = ncol(Y))
                    rownames(full_weights_matrix) <- rownames(weights_matrix)
                    colnames(full_weights_matrix) <- colnames(weights_matrix)
                    full_weights_matrix[valid_columns, ] <- weights_matrix[valid_columns, ]
                    return(X_test %*% full_weights_matrix)
                } else {
                    sapply(1:ncol(Y_train), function(k) {
                        
                        weights <- do.call(method, c(list(X = X_train, y = Y_train[, k]), args))
                        full_weights <- rep(0, ncol(X))
                        full_weights[valid_columns] <- weights
                        # Handle NAs in weights
                        full_weights[is.na(full_weights)] <- 0
                        return(X_test %*% full_weights)
                    })
                }
            }), names(weight_methods))
        }

        if (num_cores >= 2) {
            cl <- makeCluster(num_cores)
            registerDoParallel(cl)
            fold_results <- foreach(j = 1:fold, .combine = 'c') %dopar% {
                process_method(j)
            }
            stopCluster(cl)
        } else {
            fold_results <- lapply(1:fold, process_method)
        }
        # Reorganize into Y_pred
        # After cross validation, each sample should have been in
        # test set at some point, and therefore has predicted value.
        # The prediction matrix is therefore exactly the same dimension as input Y
        Y_pred <- setNames(lapply(weight_methods, function(x) matrix(NA, nrow = nrow(Y), ncol = ncol(Y))), names(weight_methods))
        for (j in 1:length(fold_results)) {
            fold_indices <- which(folds == j)
            for (method in names(weight_methods)) {
                Y_pred[[method]][fold_indices, ] <- fold_results[[j]][[method]]
            }
        }
        names(Y_pred) <- gsub("_weights", "_predicted", names(Y_pred))
        # Compute rsq and p-value for each method
        rsq_pval_table <- matrix(NA, nrow = length(weight_methods), ncol = 2)
        rownames(rsq_pval_table) <- names(Y_pred)
        colnames(rsq_pval_table) <- c("rsq", "pval")
        for (r in 1:ncol(Y)) {
            for (m in 1:length(weight_methods)) {
                if ( !is.na(sd(Y_pred[[m]][, r])) && sd(Y_pred[[m]][, r]) != 0 ) {
                    lm_fit <- lm(Y[, r] ~ Y_pred[[m]][, r])
                    rsq_pval_table[m, 1] <- summary(lm_fit)$adj.r.squared
                    rsq_pval_table[m, 2] <- summary(lm_fit)$coefficients[2,4]
                } else {
                    rsq_pval_table[m, 1] <- NA
                    rsq_pval_table[m, 2] <- NA
                }
            }
        }
        return(list(sample_partition = sample_partition, prediction = Y_pred, rsq_pval = rsq_pval_table))
    }
} 

#' Run multiple TWAS weight methods
#'
#' Applies specified weight methods to the datasets X and Y, returning weight matrices for each method.
#' Handles both univariate and multivariate methods, and filters out columns in X with zero standard error.
#' This function utilizes parallel processing to handle multiple methods.
#'
#' @param X A matrix of samples by features, where each row represents a sample and each column a feature.
#' @param Y A matrix (or vector, which will be converted to a matrix) of samples by outcomes, where each row corresponds to a sample.
#' @param weight_methods A list of methods and their specific arguments, formatted as list(method1 = method1_args, method2 = method2_args). 
#' methods in the list are applied to the datasets X and Y.
#' @param num_threads The number of threads to use for parallel processing.
#'        If set to -1, the function uses all available cores.
#'        If set to 0 or 1, no parallel processing is performed.
#'        If set to 2 or more, parallel processing is enabled with that many threads.
#' @return A list where each element is named after a method and contains the weight matrix produced by that method.
#'
#' @importFrom parallel makeCluster
#' @importFrom parallel detectCores
#' @importFrom parallel stopCluster
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @importFrom doParallel registerDoParallel
#' @export
twas_weights <- function(X, Y, weight_methods, num_threads = 1) {
    if (!is.matrix(X) || (!is.matrix(Y) && !is.vector(Y))) {
        stop("X must be a matrix and Y must be a matrix or a vector.")
    }

    if (is.vector(Y)) {
        Y <- matrix(Y, ncol = 1)
    }

    if (nrow(X) != nrow(Y)) {
        stop("The number of rows in X and Y must be the same.")
    }

    # Determine number of cores to use
    num_cores <- ifelse(num_threads == -1, parallel::detectCores(), num_threads)

    process_method <- function(method_name) {
        # Hardcoded vector of multivariate methods
        multivariate_weight_methods <- c('mrmash_weights')
        args <- weight_methods[[method_name]]
        # Remove columns with zero standard error
        valid_columns <- apply(X, 2, function(col) sd(col) != 0)
        X_filtered <- X[, valid_columns]

        if (method_name %in% multivariate_weight_methods) {
            # Apply multivariate method
            weights_matrix <- do.call(method_name, c(list(X = X_filtered, Y = Y), args))
            # Adjust the weights matrix to include zeros for invalid columns
            full_weights_matrix <- matrix(0, nrow = ncol(X), ncol = ncol(Y))
            rownames(full_weights_matrix) <- rownames(weights_matrix)
            colnames(full_weights_matrix) <- colnames(weights_matrix)
            full_weights_matrix[valid_columns, ] <- weights_matrix[valid_columns, ]
            return(full_weights_matrix)
        } else {
            # Apply univariate method to each column of Y
            # Initialize it with zeros to avoid NA
            weights_matrix <- matrix(0, nrow = ncol(X_filtered), ncol = ncol(Y))
            for (k in 1:ncol(Y)) {
                weights_vector <- do.call(method_name, c(list(X = X_filtered, y = Y[, k]), args))
                weights_matrix[, k] <- weights_vector
            }
            return(weights_matrix)
        }
    }

    if (num_cores >= 2) {
        # Set up parallel backend to use multiple cores
        cl <- makeCluster(num_cores)
        registerDoParallel(cl)
        weights_list <- foreach(method_name = names(weight_methods), .combine = 'c') %dopar% {
            process_method(method_name)
        }
        stopCluster(cl)
    } else {
        weights_list <- lapply(names(weight_methods), process_method)
    }
    names(weights_list) <- names(weight_methods)

    return(weights_list)
}

#' @importFrom susieR coef.susie
#' @export
susie_weights <- function(X=NULL, y=NULL, susie_fit=NULL, ...) {
    if (is.null(susie_fit)) {
        # get susie_fit object
        susie_fit = susie_wrapper(X,y,...)
    }
    if ("intercept" %in% names(susie_fit) && "alpha" %in% names(susie_fit) && "mu" %in% names(susie_fit) && "X_column_scale_factors" %in% names(susie_fit)) {
        # This is designed to cope with output from pecotmr::susie_post_processor()
        return(coef.susie(susie_fit)[-1])
    } else {
        return(rep(0, length(susie_fit$pip)))
    }
}

#' @importFrom mr.mash.alpha coef.mr.mash
#' @export
mrmash_weights <- function(...) {
    res <- mrmash_wrapper(...)
    return(coef.mr.mash(res)[-1,])
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
glmnet_weights <- function(X, y, alpha) {
	eff.wgt = matrix(0, ncol=1, nrow=ncol(X))
	sds = apply(X, 2, sd)
	keep = sds != 0 & !is.na(sds)
	enet = cv.glmnet(x=X[,keep], y=y, alpha=alpha, nfold=5, intercept=T, standardize=F)
	eff.wgt[keep] = coef( enet , s = "lambda.min")[2:(sum(keep)+1)]
	return(eff.wgt)
}

#' @export 
enet_weights <- function(X, y) glmnet_weights(X,y,0.5) 

#' @export
lasso_weights <- function(X, y) glmnet_weights(X,y,1) 

#' @examples 
#' wgt.mr.ash = mrash_weights(eqtl$X, eqtl$y_res, beta.init=lasso_weights(X,y))
#' @importFrom mr.ash.alpha mr.ash
#' @importFrom stats predict
#' @export
mrash_weights <- function(X, y, init_prior_sd=TRUE, ...) {
    args_list <- list(...)
    if (!"beta.init" %in% names(args_list)) {
        args_list$beta.init <- lasso_weights(X, y)
    }
    fit.mr.ash <- do.call("mr.ash", c(list(X = X, y = y, sa2 = ifelse(init_prior_sd, init_prior_sd(X, y)^2, NULL)), args_list))
    predict(fit.mr.ash, type = "coefficients")[-1]
}

#' @examples 
#' wgt.lasso = glmnet_weights(X, y, alpha=1)
#' wgt.mr.ash = mrash_weights(eqtl$X, eqtl$y_res, beta.init=wgt.lasso)
#' @importFrom mr.ash.alpha mr.ash
#' @importFrom stats predict
#' @export
mrash_weights <- function(X, y, init_prior_sd=TRUE, ...) {
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
           

           
## TWAS with UTMOST framework
## the utmost paper assumes X are not standardized, in the formula below - the input of X is assumed to be standardized 
#' @importFrom GBJ GBJ
#' @export
           
twas_joint_z <- function(ld, Bhat, gwas_z){
    idx <- which(rownames(ld) %in% rownames(Bhat))
    D <- ld[idx,idx]
    
    cov_y <- crossprod(Bhat, D) %*% Bhat
    y_sd <- sqrt(diag(cov_y))
    x_sd <- rep(1, nrow(Bhat))   #we assume X is standardized
    
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
