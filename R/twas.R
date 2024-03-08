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
#' @param max_num_variants An optional integer to set the randomly selected maximum number of variants to use for CV purpose, to save computing time.
#' @param variants_to_keep An optional integer to ensure that the listed variants are kept in the CV when there is a limit on the max_num_variants to use.
#' @param num_threads The number of threads to use for parallel processing.
#'        If set to -1, the function uses all available cores.
#'        If set to 0 or 1, no parallel processing is performed.
#'        If set to 2 or more, parallel processing is enabled with that many threads.
#' @return A list with the following components:
#' \itemize{
#'   \item `sample_partition`: A dataframe showing the sample partitioning used in the cross-validation.
#'   \item `prediction`: A list of matrices with predicted Y values for each method and fold.
#'   \item `metrics`: A matrix with rows representing methods and columns for various metrics:
#'     \itemize{
#'       \item `corr`: Pearson's correlation between predicated and observed values.
#'       \item `adj_rsq`: Adjusted R-squared value (which indicates the proportion of variance explained by the model) that accounts for the number of predictors in the model.
#'       \item `pval`: P-value assessing the significance of the model's predictions.
#'       \item `RMSE`: Root Mean Squared Error, a measure of the model's prediction error.
#'       \item `MAE`: Mean Absolute Error, a measure of the average magnitude of errors in a set of predictions.
#'     }
#'   \item `time_elapsed`: The time taken to complete the cross-validation process.
#' }
#' @importFrom future plan multisession availableCores
#' @importFrom furrr future_map furrr_options
#' @importFrom purrr map
#' @export
twas_weights_cv <- function(X, Y, fold = NULL, sample_partitions = NULL, weight_methods = NULL, seed = NULL, max_num_variants = NULL, variants_to_keep = NULL, num_threads = 1, ...) {
    split_data <- function(X, Y, sample_partition, fold){
      if (is.null(rownames(X))) {
        warning("Row names in X are missing. Using row indices.")
        rownames(X) <- 1:nrow(X)
      }
      if (is.null(rownames(Y))) {
        warning("Row names in Y are missing. Using row indices.")
        rownames(Y) <- 1:nrow(Y)
      }
      test_ids <- sample_partition[which(sample_partition$Fold == fold), "Sample"]
      Xtrain <- X[!(rownames(X) %in% test_ids), ,drop=FALSE]
      Ytrain <- Y[!(rownames(Y) %in% test_ids), ,drop=FALSE]
      Xtest <- X[rownames(X) %in% test_ids, ,drop=FALSE]
      Ytest <- Y[rownames(Y) %in% test_ids, ,drop=FALSE]
      if (nrow(Xtrain) == 0 || nrow(Ytrain) == 0 || nrow(Xtest) == 0 || nrow(Ytest) == 0) {
        stop("Error: One of the datasets (train or test) has zero rows.")
      }
      return(list(Xtrain=Xtrain, Ytrain=Ytrain, Xtest=Xtest, Ytest=Ytest))
    }

    # Validation checks
    if (!is.null(fold) && (!is.numeric(fold) || fold <= 0)) {
        stop("Invalid value for 'fold'. It must be a positive integer.")
    }
    
    if (!is.matrix(X) || (!is.matrix(Y) && !is.vector(Y))) {
        stop("X must be a matrix and Y must be a matrix or a vector.")
    }
    
    if (is.vector(Y)) {
        Y <- matrix(Y, ncol = 1)
        message(paste("Y converted to matrix of", nrow(Y), "rows and", ncol(Y), "columns."))
    }
    
    if (nrow(X) != nrow(Y)) {
        stop("The number of rows in X and Y must be the same.")
    }

    # Check if X has row names and Y does not
    if (!is.null(rownames(X)) && is.null(rownames(Y))) {
        rownames(Y) <- rownames(X)
    }

    # Check if Y has row names and X does not
    if (!is.null(rownames(Y)) && is.null(rownames(X))) {
        rownames(X) <- rownames(Y)
    }

    # Get sample names
    if (!is.null(rownames(X))) {
        sample_names <- rownames(X)
    } else if (!is.null(rownames(Y))) {
        sample_names <- rownames(Y)
    } else {
        sample_names <- 1:nrow(X)
    }

    # Select variants if necessary
    if (!is.null(max_num_variants) && ncol(X)> max_num_variants) {
        if (!is.null(variants_to_keep) && length(variants_to_keep) > 0) {
            variants_to_keep <- intersect(variants_to_keep, colnames(X))
            remaining_columns <- setdiff(colnames(X), variants_to_keep)
            if (length(variants_to_keep) < max_num_variants) {
                additional_columns <- sample(remaining_columns, max_num_variants - length(variants_to_keep), replace = FALSE)
                selected_columns <- union(variants_to_keep, additional_columns)
                message(sprintf("Including %d specified variants and randomly selecting %d additional variants, for a total of %d variants out of %d for cross-validation purpose.",
                    length(variants_to_keep), length(additional_columns), length(selected_columns), ncol(X)))
            } else {
                selected_columns <- sample(variants_to_keep, max_num_variants, replace = FALSE)
                message(paste("Randomly selecting", length(selected_columns), "out of", length(variants_to_keep), "input variants for cross validation purpose."))
            }
        } else {
            selected_columns <- sort(sample(ncol(X), max_num_variants, replace = FALSE))
            message(paste("Randomly selecting", length(selected_columns), "out of", ncol(X), "variants for cross validation purpose."))
        }
        X <- X[, selected_columns]
    }
    
    arg <-list(...)
    
    # Create or use provided folds
    if (!is.null(fold)) {
        if (!is.null(seed)) set.seed(seed)
        
        if (!is.null(sample_partitions)) {
            if(fold != length(unique(sample_partition$Fold))){
                message(paste0("fold number provided does not match with sample partition, performing ", length(unique(sample_partition$Fold)),
                       " fold cross validation based on provided sample partition. "))
                }
            
            folds <- sample_partitions$Fold
            sample_partition <- sample_partitions
            
        } else {
            sample_indices <- sample(nrow(X))
            folds <- cut(seq(1, nrow(X)), breaks = fold, labels = FALSE)
            sample_partition <- data.frame(Sample = sample_names[sample_indices], Fold = folds, stringsAsFactors = FALSE)
        }
    } else if (!is.null(sample_partitions)) {
        if (!all(sample_partitions$Sample %in% sample_names)) {
            stop("Some samples in 'sample_partitions' do not match the samples in 'X' and 'Y'.")
        }
        sample_partition <- sample_partitions
        fold <- length(unique(sample_partition$Fold))
    } else {
        stop("Either 'fold' or 'sample_partitions' must be provided.")
    }

    st = proc.time()
    if (is.null(weight_methods)) {
        return(list(sample_partition = sample_partition))
    } else {
        # Hardcoded vector of multivariate weight_methods
        multivariate_weight_methods <- c('mrmash_weights', 'mvsusie_weights')

        # Determine the number of cores to use
        num_cores <- ifelse(num_threads == -1, availableCores(), num_threads)
        num_cores <- min(num_cores, availableCores())
        
        # Perform CV with parallel processing
        compute_method_predictions <- function(j) { 
            dat_split <- split_data(X, Y, sample_partition=sample_partition, fold=j)
            X_train <- dat_split$Xtrain
            Y_train <- dat_split$Ytrain
            X_test <- dat_split$Xtest
            Y_test <- dat_split$Ytest

            # Remove columns with zero standard error
            valid_columns <- apply(X_train, 2, function(col) sd(col) != 0)
            X_train <- X_train[, valid_columns, drop=F]

            setNames(lapply(names(weight_methods), function(method) {
                args <- weight_methods[[method]]
                if (method %in% multivariate_weight_methods) {
                    # Apply multivariate method to entire Y for this fold
                    prior_matrices <- arg$mrmash_weights_prior_matrices
                    prior_matrices <- prior_matrices[[paste0("fold_", j)]]
                    weights_matrix <- do.call(method, c(list(X = X_train, Y = Y_train, prior_data_driven_matrices=prior_matrices, 
                                                             cannonical_matrices=arg$prior_canonical_matrices,  
                                                             max_iter=arg$mrmash_max_iter), arg))
                    # Adjust the weights matrix to include zeros for invalid columns
                    full_weights_matrix <- matrix(0, nrow = ncol(X), ncol = ncol(Y))
                    rownames(full_weights_matrix) <- rownames(weights_matrix)
                    colnames(full_weights_matrix) <- colnames(weights_matrix)
                    full_weights_matrix[valid_columns, ] <- weights_matrix[valid_columns, ]
                    Y_pred <- X_test %*% full_weights_matrix
                    rownames(Y_pred) <- rownames(X_test)
                    return(Y_pred)
                } else {
                    Y_pred <- sapply(1:ncol(Y_train), function(k) {
                        if (!is.null(seed)) set.seed(seed)
                        weights <- do.call(method, c(list(X = X_train, y = Y_train[, k]), args))
                        full_weights <- rep(0, ncol(X))
                        full_weights[valid_columns] <- weights
                        # Handle NAs in weights
                        full_weights[is.na(full_weights)] <- 0
                        X_test %*% full_weights
                    })
                    rownames(Y_pred) <- rownames(X_test)
                    return(Y_pred)
                }
            }), names(weight_methods))
        }

        if (num_cores >= 2) {
            plan(multisession, workers = num_cores)
            fold_results <- future_map(1:fold, compute_method_predictions, .options = furrr_options(seed = seed))
        } else { 
            fold_results <- map(1:fold, compute_method_predictions)
        }
                                   
        # Reorganize into Y_pred
        # After cross validation, each sample should have been in
        # test set at some point, and therefore has predicted value.
        # The prediction matrix is therefore exactly the same dimension as input Y
        Y_pred <- setNames(lapply(weight_methods, function(x) `dimnames<-`(matrix(NA, nrow(Y), ncol(Y)), dimnames(Y))), names(weight_methods))
        for (j in 1:length(fold_results)) {
            for (method in names(weight_methods)) {
                Y_pred[[method]][rownames(fold_results[[j]][[method]]), ] <- fold_results[[j]][[method]]
            }
        }
                                  
        names(Y_pred) <- gsub("_weights", "_predicted", names(Y_pred))
     
        # Compute rsq, adj rsq, p-value, RMSE, and MAE for each method
        # metrics_table <- matrix(NA, nrow = length(weight_methods), ncol = 5)
        metrics_table <- list()
                                  
        for (m in names(weight_methods)){
            
            metrics_table[[m]] <-  matrix(NA, nrow = ncol(Y), ncol = 5)
            colnames(metrics_table[[m]]) <- c("corr", "adj_rsq", "adj_rsq_pval", "RMSE", "MAE")
            rownames(metrics_table[[m]]) <- colnames(Y)
            
            for (r in 1:ncol(Y)){
                method_predictions <- Y_pred[[gsub("_weights", "_predicted", m)]][, r]
                actual_values <- Y[, r]
                # Remove missing values in the first place
                na_indx <- which(is.na(actual_values))
                if (length(na_indx)!=0) {
                    method_predictions <- method_predictions[-na_indx] 
                    actual_values <- actual_values[-na_indx] 
                }
                if ( sd(method_predictions) != 0 ) {
                    lm_fit <- lm(actual_values ~ method_predictions)
                    
                    # Calculate raw correlation and and adjusted R-squared
                    metrics_table[[m]][r, "corr"] <- cor(actual_values, method_predictions)

                    
                    metrics_table[[m]][r, "adj_rsq"] <- summary(lm_fit)$adj.r.squared

                    # Calculate p-value
                    metrics_table[[m]][r, "adj_rsq_pval"] <- summary(lm_fit)$coefficients[2, 4]

                    # Calculate RMSE
                    residuals <- actual_values - method_predictions
                    metrics_table[[m]][r, "RMSE"] <- sqrt(mean(residuals^2))

                    # Calculate MAE
                    metrics_table[[m]][r, "MAE"] <- mean(abs(residuals))
                } else {
                    metrics_table[[m]][r, ] <- NA
                    message(paste0("Predicted values for condition ", r , " using ", m , 
                                   " have zero variance. Filling performance metric with NAs"))
                }
            }
        }
        names(metrics_table) <- gsub("_weights", "_performance", names(metrics_table))
        return(list(sample_partition = sample_partition, prediction = Y_pred, performance = metrics_table, time_elapsed = proc.time() - st))
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
#' @export
#' @importFrom future plan multisession availableCores
#' @importFrom furrr future_map furrr_options
#' @importFrom purrr map
twas_weights <- function(X, Y, weight_methods, num_threads = 1, seed = NULL) {
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
    num_cores <- ifelse(num_threads == -1, availableCores(), num_threads)
    num_cores <- min(num_cores, availableCores())

    compute_method_weights <- function(method_name) {
        # Hardcoded vector of multivariate methods
        multivariate_weight_methods <- c('mrmash_weights', 'mvsusie_weights')
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
            if(!is.null(seed)) set.seed(seed)
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
        plan(multisession, workers = num_cores)
        weights_list <- names(weight_methods) %>% future_map(compute_method_weights, .options = furrr_options(seed = seed))
    } else {
        weights_list <- names(weight_methods) %>% map(compute_method_weights)
    }
    names(weights_list) <- names(weight_methods)

    if (!is.null(colnames(X))) {
        weights_list <- lapply(weights_list, function(x) { rownames(x) <- colnames(X); return(x) })
    }
    return(weights_list)
}

#' @export
twas_predict <- function(X, weights_list) {
   setNames(lapply(weights_list, function(w) X %*% w), gsub("_weights", "_predicted", names(weights_list)))
}

#' @importFrom susieR coef.susie
#' @export
susie_weights <- function(X=NULL, y=NULL, susie_fit=NULL, ...) {
    if (is.null(susie_fit)) {
        # get susie_fit object
        susie_fit = susie_wrapper(X,y,...)
    }
    if ("alpha" %in% names(susie_fit) && "mu" %in% names(susie_fit) && "X_column_scale_factors" %in% names(susie_fit)) {
        # This is designed to cope with output from pecotmr::susie_post_processor()
        # We set intercept to 0 and later trim it off anyways
        susie_fit$intercept = 0 
        return(coef.susie(susie_fit)[-1])
    } else {
        return(rep(0, length(susie_fit$pip)))
    }
}
                               
#' Adjust SuSiE Weights
#'
#' This function adjusts the SuSiE weights based on a set of intersected variants.
#' It subsets various components like lbf_matrix, mu, and scale factors based on these variants.
#'
#' @param weight_db_file A RDS file containing TWAS weights.
#' @param condition specific condition.
#' @param keep_variants Vector of variant names to keep.
#' @param allele_qc Optional                   
#' @return A list of adjusted xQTL coefficients and remained variants ids
#' @export

adjust_susie_weights <- function(twas_weights_results, condition, keep_variants, allele_qc = TRUE) {
  # Intersect the rownames of weights with keep_variants
  twas_weights_variants <- get_nested_element(twas_weights_results, c("susie_results",condition, "variant_names"))
  # allele flip twas weights matrix variants name
  if(allele_qc == TRUE){
            weights_matrix <-  get_nested_element(twas_weights_results, c("weights",condition))
            weights_matrix_qced <- allele_qc(twas_weights_variants, gwas_LD_list$combined_LD_variants, weights_matrix, 1:ncol(weights_matrix))
            intersected_indices <- which(weights_matrix_qced$qc_summary$keep == TRUE)
  } else {
      keep_variants_transformed <- ifelse(!startsWith(keep_variants, "chr"), paste0("chr", keep_variants), keep_variants)
      intersected_variants <- intersect(twas_weights_variants, keep_variants_transformed)
      intersected_indices <- match(intersected_variants, twas_weights_variants)
  }
  if (length(intersected_indices) == 0) {
    stop("Error: No intersected variants found. Please check 'twas_weights' and 'keep_variants' inputs to make sure there are variants left to use.")
  }
  # Subset lbf_matrix, mu, and x_column_scale_factors
  lbf_matrix <- get_nested_element(twas_weights_results,c("susie_results",condition,"susie_result_trimmed","lbf_variable"))
  mu <- get_nested_element(twas_weights_results,c("susie_results",condition,"susie_result_trimmed","mu"))
  x_column_scal_factors <- get_nested_element(twas_weights_results,c("susie_results",condition,"susie_result_trimmed","X_column_scale_factors"))

  lbf_matrix_subset <- lbf_matrix[, intersected_indices]
  mu_subset <- mu[, intersected_indices]
  x_column_scal_factors_subset <- x_column_scal_factors[intersected_indices]

  # Convert lbf_matrix to alpha and calculate adjusted xQTL coefficients
  adjusted_xqtl_alpha <- lbf_to_alpha(lbf_matrix_subset)
  adjusted_xqtl_coef <- colSums(adjusted_xqtl_alpha * mu_subset) / x_column_scal_factors_subset

  return(list(adjusted_susie_weights = adjusted_xqtl_coef,remained_variants_ids = names(adjusted_xqtl_coef)))
}

#' @importFrom mr.mash.alpha coef.mr.mash
#' @export
mrmash_weights <- function(...) {
    res <- mrmash_wrapper(...)
    return(coef.mr.mash(res)[-1,])
}
                   
#' @importFrom mvsusieR coef.mvsusie mvsusie create_mixture_prior
#' @export               
mvsusie_weights <- function(mvsusie_fit=NULL, X=NULL, Y=NULL, prior_variance=NULL, residual_variance=NULL, L=30, mvsusie_max_iter=200, ...) {
    if (is.null(mvsusie_fit)) {
        message("Did not provide mvsusie_fit, fitting mvSuSiE now")
        if (is.null(X) || is.null(Y)) {
            stop("Both X and Y must be provided if mvsusie_fit is NULL.")
        }
        if (is.null(prior_variance)) prior_variance = create_mixture_prior(R=ncol(Y)) 
        if (is.null(residual_variance)) residual_variance = mr.mash.alpha:::compute_cov_flash(Y)
        
        mvsusie_fit = mvsusie(X=X, Y=Y, L=L, prior_variance=prior_variance, 
                            residual_variance=residual_variance, precompute_covariances=F, 
                            compute_objective=T, estimate_residual_variance=F, estimate_prior_variance=T, 
                            estimate_prior_method='EM', max_iter=mvsusie_max_iter, 
                            n_thread=1, approximate=F)
    }
    return(coef.mvsusie(mvsusie_fit)[-1,])
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
    fit.mr.ash <- do.call("mr.ash", c(list(X = X, y = y, sa2 = if (init_prior_sd) init_prior_sd(X, y)^2 else NULL), args_list))
    predict(fit.mr.ash, type = "coefficients")[-1]
}

pval_acat <- function(pvals) {
    if (length(pvals) == 1) {
        return(pvals[1])
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
    global_pval <- if (comb_method == "HMP") pval_hmp(pvals) else pval_acat(pvals) # pval vector
    naive_pval <- min(n_total_tests*min_pval, 1.0)
    return(if (naive) naive_pval else global_pval) # global_pval and naive_pval
}
           
## Multi-condition TWAS joint test
## the UTMOST paper assumes X are not standardized, in the formula below - the input of X is assumed to be standardized 
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

#' TWAS Weights Pipeline 
#'
#' This function performs weights computation for Transcriptome-Wide Association Study (TWAS)
#' incorporating various steps such as filtering variants by linkage disequilibrium reference panel variants,
#' fitting models using SuSiE and other methods, and calculating TWAS weights and predictions.
#' Optionally, it can perform cross-validation for TWAS weights.
#'
#' @param X A matrix of genotype data where rows represent samples and columns represent genetic variants.
#' @param y A vector of phenotype measurements for each sample.
#' @param susie_fit An object returned by the SuSiE function, containing the SuSiE model fit.
#' @param maf A vector of minor allele frequencies for each variant in X.
#' @param ld_reference_meta_file An optional path to a file containing linkage disequilibrium reference data. If provided, variants in X are filtered based on this reference.
#' @param cv_folds The number of folds to use for cross-validation. Set to 0 to skip cross-validation.
#' @param X_scalar A scalar or vector to scale the genotype data. Defaults to 1 (no scaling).
#' @param y_scalar A scalar to scale the phenotype data. Defaults to 1 (no scaling).
#' @param coverage The coverage probability used in SuSiE for credible set construction. Defaults to 0.95.
#' @param secondary_coverage A vector of secondary coverage probabilities for credible set refinement. Defaults to c(0.7, 0.5).
#' @param signal_cutoff A threshold for determining significant signals in the SuSiE output. Defaults to 0.05.
#' @param mr_ash_max_iter The maximum number of iterations for the MR-ASH method. Defaults to 100.
#' @param min_cv_maf The minimum minor allele frequency for variants to be included in cross-validation. Defaults to 0.05.
#' @param max_cv_variants The maximum number of variants to be included in cross-validation. Defaults to 5000.
#' @param cv_seed The seed for random number generation in cross-validation. Defaults to 999.
#' @param cv_threads The number of threads to use for parallel computation in cross-validation. Defaults to 1.
#' @return A list containing results from the TWAS pipeline, including TWAS weights, predictions, and optionally cross-validation results.
#' @export
#' @examples
#' # Example usage (assuming appropriate objects for X, y, susie_fit, and maf are available):
#' twas_results <- twas_weights_pipeline(X, y, maf, susie_fit)
twas_weights_pipeline <- function(X, y, maf, susie_fit, ld_reference_meta_file = NULL, cv_folds = 5, X_scalar = 1, y_scalar = 1, 
                          coverage = 0.95, secondary_coverage = c(0.7, 0.5), signal_cutoff = 0.05, 
                          mr_ash_max_iter = 100, min_cv_maf = 0.05, max_cv_variants = 5000, cv_seed = 999, cv_threads = 1) {
    res <- list()
    if (!is.null(susie_fit)) {
        L <- length(which(susie_fit$V > 1E-9))
        init_L <- max(1, L - 2)
        max_L <- L + 3
    } else {
        # susie_fit did not detect anything significant
        init_L <- 2
        max_L <- 2
    }
    if (!is.null(ld_reference_meta_file)) {
        variants_kept <- filter_variants_by_ld_reference(colnames(X), ld_reference_meta_file)
        X <- X[, variants_kept$data, drop = FALSE]
        maf <- maf[variants_kept$idx]
        res$preset_variants_result <- susie_wrapper(X, y, init_L = init_L, max_L = max_L, refine = TRUE, coverage = coverage)
        res$preset_variants_result <- susie_post_processor(res$preset_variants_result, X, y, if (X_scalar == 1) 1 else X_scalar[variants_kept$idx], y_scalar, maf, secondary_coverage = secondary_coverage, signal_cutoff = signal_cutoff)
        res$preset_variants_result$analysis_script <- NULL
        res$preset_variants_result$sumstats <- NULL
        susie_fit <- res$preset_variants_result$susie_result_trimmed
    }
    weight_methods <- list(enet_weights = list(), lasso_weights = list(), mrash_weights = list(init_prior_sd = TRUE, max.iter = mr_ash_max_iter))
    if (!is.null(susie_fit)) weight_methods$susie_weights <- list(susie_fit = susie_fit)
    # get TWAS weights
    res$twas_weights <- twas_weights(X, y, weight_methods = weight_methods)
    # get TWAS predictions for possible next steps such as computing correlations between predicted expression values
    res$twas_predictions <- twas_predict(X, res$twas_weights)
    if (cv_folds > 0) {
        # A few cutting corners to run CV faster at the disadvantage of SuSiE and mr.ash:
        # 1. reset SuSiE to not using refine or adaptive L (more or less the default SuSiE)
        # 2. at most 100 iterations for mr.ash allowed
        # 3. only use a subset of top signals and common variants
        if (!is.null(susie_fit)) weight_methods$susie_weights <- list(refine = FALSE, init_L = max_L, max_L = max_L)
        variants_for_cv <- c()
        if (!is.null(res$preset_variants_result$top_loci) && nrow(res$preset_variants_result$top_loci) > 0) {
            variants_for_cv <- res$preset_variants_result$top_loci[, 1]
        }
        common_var <- colnames(X)[which(maf > min_cv_maf)]
        if (length(common_var) + length(variants_for_cv) > max_cv_variants) {
            common_var <- sample(common_var, max_cv_variants - length(variants_for_cv), replace = FALSE)
        }
        variants_for_cv <- unique(c(variants_for_cv, common_var))
        res$twas_cv_result <- twas_weights_cv(X, y, fold = cv_folds, weight_methods = weight_methods, seed = cv_seed, max_num_variants = max_cv_variants, num_threads = cv_threads, variants_to_keep = if (length(variants_for_cv) > 0) variants_for_cv else NULL)
    }
    return(res)
}
