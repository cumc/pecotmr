#' Calculate TWAS z-score and p-value
#'
#' This function calculates the TWAS z-score and p-value given the weights, z-scores,
#' and optionally the correlation matrix (R) or the genotype matrix (X).
#'
#' @param weights A numeric vector of weights.
#' @param z A numeric vector of z-scores.
#' @param R An optional correlation matrix. If not provided, it will be calculated from the genotype matrix X.
#' @param X An optional genotype matrix. If R is not provided, X must be supplied to calculate the correlation matrix.
#'
#' @return A list containing the following elements:
#' \itemize{
#'   \item z: The TWAS z-score.
#'   \item pval: The corresponding p-value.
#' }
#'
#' @importFrom Rfast cora
#' @importFrom stats cor pchisq
#'
#' @export
twas_z <- function(weights, z, R = NULL, X = NULL) {
  # Check that weights and z-scores have the same length
  if (length(weights) != length(z)) {
    stop("Weights and z-scores must have the same length.")
  }

  if (is.null(R)) R <- compute_LD(X)

  stat <- t(weights) %*% z
  denom <- t(weights) %*% R %*% weights
  zscore <- stat / sqrt(denom)
  pval <- pchisq(zscore * zscore, 1, lower.tail = FALSE)

  return(list(z = zscore, pval = pval))
}

#' Multi-condition TWAS joint test
#'
#' This function performs a multi-condition TWAS joint test using the GBJ method.
#' It assumes that the input genotype matrix (X) is standardized.
#'
#' @param R An optional correlation matrix. If not provided, it will be calculated from the genotype matrix X.
#' @param X An optional genotype matrix. If R is not provided, X must be supplied to calculate the correlation matrix.
#' @param weights A matrix of weights, where each column corresponds to a different condition.
#' @param z A vector of GWAS z-scores.
#'
#' @return A list containing the following elements:
#' \itemize{
#'   \item Z: A matrix of TWAS z-scores and p-values for each condition.
#'   \item GBJ: The result of the GBJ test.
#' }
#'
#' @importFrom GBJ GBJ
#' @importFrom stats cor pnorm
#'
#' @export
twas_joint_z <- function(R = NULL, X = NULL, weights, z) {
  # Check that weights and z-scores have the same number of rows
  if (nrow(weights) != length(z)) {
    stop("Number of rows in weights must match the length of z-scores.")
  }


  if (is.null(R)) R <- compute_LD(X)

  idx <- which(rownames(R) %in% rownames(weights))
  D <- R[idx, idx]

  cov_y <- crossprod(weights, D) %*% weights
  y_sd <- sqrt(diag(cov_y))
  x_sd <- rep(1, nrow(weights)) # Assuming X is standardized

  # Get gamma matrix MxM (snp x snp)
  g <- lapply(colnames(weights), function(x) {
    gm <- diag(x_sd / y_sd[x], length(x_sd), length(x_sd))
    return(gm)
  })
  names(g) <- colnames(weights)

  ######### Get TWAS - Z statistics & P-value, GBJ test ########
  z_matrix <- do.call(rbind, lapply(colnames(weights), function(x) {
    Zi <- crossprod(weights[, x], g[[x]]) %*% as.numeric(z)
    pval <- 2 * pnorm(abs(Zi), lower.tail = FALSE)
    Zp <- c(Zi, pval)
    names(Zp) <- c("Z", "pval")
    return(Zp)
  }))
  rownames(z_matrix) <- colnames(weights)

  # GBJ test
  lam <- matrix(rep(NA, ncol(weights) * nrow(weights)), nrow = ncol(weights))
  rownames(lam) <- colnames(weights)

  for (p in colnames(weights)) {
    la <- as.matrix(weights[, p] %*% g[[p]])
    lam[p, ] <- la
  }

  sig <- tcrossprod((lam %*% D), lam)
  gbj <- GBJ(test_stats = z_matrix[, 1], cor_mat = sig)

  rs <- list("Z" = z_matrix, "GBJ" = gbj)
  return(rs)
}

#' TWAS Analysis
#'
#' Performs TWAS analysis using the provided weights matrix, GWAS summary statistics database,
#' and LD matrix. It extracts the necessary GWAS summary statistics and LD matrix based on the
#' specified variants and computes the z-score and p-value for each gene.
#'
#' @param weights_matrix A matrix containing weights for all methods.
#' @param gwas_sumstats_db A data frame containing the GWAS summary statistics.
#' @param LD_matrix A matrix representing linkage disequilibrium between variants.
#' @param extract_variants_objs A vector of variant identifiers to extract from the GWAS and LD matrix.
#'
#' @return A list with TWAS z-scores and p-values across four methods for each gene.
#' @export
twas_analysis <- function(weights_matrix, gwas_sumstats_db, LD_matrix, extract_variants_objs) {
  #
  # Extract gwas_sumstats
  gwas_sumstats_subset <- gwas_sumstats_db[match(extract_variants_objs, gwas_sumstats_db$variant_id), ]
  # Validate that the GWAS subset is not empty
  if (nrow(gwas_sumstats_subset) == 0 | all(is.na(gwas_sumstats_subset))) {
    stop("No GWAS summary statistics found for the specified variants.")
  }
  # Check if extract_variants_objs are in the rownames of LD_matrix
  valid_indices <- extract_variants_objs %in% rownames(LD_matrix)
  if (!any(valid_indices)) {
    stop("None of the specified variants are present in the LD matrix.")
  }
  # Extract only the valid indices from extract_variants_objs
  valid_variants_objs <- extract_variants_objs[valid_indices]
  # Extract LD_matrix subset using valid indices
  LD_matrix_subset <- LD_matrix[valid_variants_objs, valid_variants_objs]
  # Caculate the z score and pvalue of each gene
  twas_z_pval <- apply(
    as.matrix(weights_matrix), 2,
    function(x) twas_z(x, gwas_sumstats_subset$z, R = LD_matrix_subset)
  )
  return(twas_z_pval)
}

#' Cross-Validation for weights selection in Transcriptome-Wide Association Studies (TWAS)
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
  split_data <- function(X, Y, sample_partition, fold) {
    if (is.null(rownames(X))) {
      warning("Row names in X are missing. Using row indices.")
      rownames(X) <- 1:nrow(X)
    }
    if (is.null(rownames(Y))) {
      warning("Row names in Y are missing. Using row indices.")
      rownames(Y) <- 1:nrow(Y)
    }
    test_ids <- sample_partition[which(sample_partition$Fold == fold), "Sample"]
    Xtrain <- X[!(rownames(X) %in% test_ids), , drop = FALSE]
    Ytrain <- Y[!(rownames(Y) %in% test_ids), , drop = FALSE]
    Xtest <- X[rownames(X) %in% test_ids, , drop = FALSE]
    Ytest <- Y[rownames(Y) %in% test_ids, , drop = FALSE]
    if (nrow(Xtrain) == 0 || nrow(Ytrain) == 0 || nrow(Xtest) == 0 || nrow(Ytest) == 0) {
      stop("Error: One of the datasets (train or test) has zero rows.")
    }
    return(list(Xtrain = Xtrain, Ytrain = Ytrain, Xtest = Xtest, Ytest = Ytest))
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
  if (!is.null(max_num_variants) && ncol(X) > max_num_variants) {
    if (!is.null(variants_to_keep) && length(variants_to_keep) > 0) {
      variants_to_keep <- intersect(variants_to_keep, colnames(X))
      remaining_columns <- setdiff(colnames(X), variants_to_keep)
      if (length(variants_to_keep) < max_num_variants) {
        additional_columns <- sample(remaining_columns, max_num_variants - length(variants_to_keep), replace = FALSE)
        selected_columns <- union(variants_to_keep, additional_columns)
        message(sprintf(
          "Including %d specified variants and randomly selecting %d additional variants, for a total of %d variants out of %d for cross-validation purpose.",
          length(variants_to_keep), length(additional_columns), length(selected_columns), ncol(X)
        ))
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

  arg <- list(...)

  # Create or use provided folds
  if (!is.null(fold)) {
    if (!is.null(seed)) set.seed(seed)

    if (!is.null(sample_partitions)) {
      if (fold != length(unique(sample_partition$Fold))) {
        message(paste0(
          "fold number provided does not match with sample partition, performing ", length(unique(sample_partition$Fold)),
          " fold cross validation based on provided sample partition. "
        ))
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

  st <- proc.time()
  if (is.null(weight_methods)) {
    return(list(sample_partition = sample_partition))
  } else {
    # Hardcoded vector of multivariate weight_methods
    multivariate_weight_methods <- c("mrmash_weights", "mvsusie_weights")

    # Determine the number of cores to use
    num_cores <- ifelse(num_threads == -1, availableCores(), num_threads)
    num_cores <- min(num_cores, availableCores())

    # Perform CV with parallel processing
    compute_method_predictions <- function(j) {
      dat_split <- split_data(X, Y, sample_partition = sample_partition, fold = j)
      X_train <- dat_split$Xtrain
      Y_train <- dat_split$Ytrain
      X_test <- dat_split$Xtest
      Y_test <- dat_split$Ytest

      # Remove columns with zero standard error
      valid_columns <- apply(X_train, 2, function(col) sd(col) != 0)
      X_train <- X_train[, valid_columns, drop = F]

      setNames(lapply(names(weight_methods), function(method) {
        args <- weight_methods[[method]]
        if (method %in% multivariate_weight_methods) {
          # Apply multivariate method to entire Y for this fold
          prior_matrices <- arg$mrmash_weights_prior_matrices
          prior_matrices <- prior_matrices[[paste0("fold_", j)]]
          weights_matrix <- do.call(method, c(list(
            X = X_train, Y = Y_train, prior_data_driven_matrices = prior_matrices,
            cannonical_matrices = arg$prior_canonical_matrices,
            max_iter = arg$mrmash_max_iter
          ), arg))
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

    for (m in names(weight_methods)) {
      metrics_table[[m]] <- matrix(NA, nrow = ncol(Y), ncol = 5)
      colnames(metrics_table[[m]]) <- c("corr", "adj_rsq", "adj_rsq_pval", "RMSE", "MAE")
      rownames(metrics_table[[m]]) <- colnames(Y)

      for (r in 1:ncol(Y)) {
        method_predictions <- Y_pred[[gsub("_weights", "_predicted", m)]][, r]
        actual_values <- Y[, r]
        # Remove missing values in the first place
        na_indx <- which(is.na(actual_values))
        if (length(na_indx) != 0) {
          method_predictions <- method_predictions[-na_indx]
          actual_values <- actual_values[-na_indx]
        }
        if (sd(method_predictions) != 0) {
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
          message(paste0(
            "Predicted values for condition ", r, " using ", m,
            " have zero variance. Filling performance metric with NAs"
          ))
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
    multivariate_weight_methods <- c("mrmash_weights", "mvsusie_weights")
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
      if (!is.null(seed)) set.seed(seed)
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
    weights_list <- lapply(weights_list, function(x) {
      rownames(x) <- colnames(X)
      return(x)
    })
  }
  return(weights_list)
}

#' Predict outcomes using TWAS weights
#'
#' This function takes a matrix of predictors (\code{X}) and a list of TWAS (transcriptome-wide
#' association studies) weights (\code{weights_list}), and calculates the predicted outcomes by
#' multiplying \code{X} by each set of weights in \code{weights_list}. The names of the elements
#' in the output list are derived from the names in \code{weights_list}, with "_weights" replaced
#' by "_predicted".
#'
#' @param X A matrix or data frame of predictors where each row is an observation and each
#' column is a variable.
#' @param weights_list A list of numeric vectors representing the weights for each predictor.
#' The names of the list elements should follow the pattern \code{[outcome]_weights}, where
#' \code{[outcome]} is the name of the outcome variable that the weights are associated with.
#'
#' @return A named list of numeric vectors, where each vector is the predicted outcome for the
#' corresponding set of weights in \code{weights_list}. The names of the list elements are
#' derived from the names in \code{weights_list} by replacing "_weights" with "_predicted".
#'
#' @export
#' @examples
#' # Assuming `X` is your matrix of predictors and `weights_list` is your list of weights:
#' predicted_outcomes <- twas_predict(X, weights_list)
#' print(predicted_outcomes)
twas_predict <- function(X, weights_list) {
  setNames(lapply(weights_list, function(w) X %*% w), gsub("_weights", "_predicted", names(weights_list)))
}
