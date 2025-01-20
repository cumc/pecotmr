#' Calculate Quantile Rank Score
#' @param pheno numeric vector of phenotype values
#' @param covariates matrix/data.frame of covariates
#' @param tau quantile level
#' @return vector of rank scores
#' @import quantreg
#' @importFrom stats rnorm
#' @noRd
calculate_rank_score <- function(pheno, covariates, tau, seed = 123) {
  set.seed(seed) 
  # Add random variable to covariates
  covar_data <- data.frame(covariates)
  covar_data$d_rv <- rnorm(nrow(covar_data))
  
  # Fit quantile regression
  Qreg <- suppressWarnings(rq(pheno ~ ., data = covar_data, tau = tau, method = "fn"))
  
  # Get rank score
  coeff <- summary(Qreg, se = "ker")$coefficients
  SE_d_rv <- coeff[nrow(coeff), 2]
  a_i_tau <- (tau - ifelse(residuals(Qreg) < 0, 1, 0)) * SE_d_rv / sqrt(-tau^2 + tau)
  
  return(a_i_tau)
}

#' Fit Rank Scores for All Quantile Levels
#' @param pheno numeric vector of phenotype values
#' @param covariates matrix/data.frame of covariates
#' @param num_tau_levels integer number of quantile levels
#' @param num_cores integer number of cores for parallel processing
#' @return list of rank scores for each quantile level
#' @importFrom parallel mclapply
#' @noRd
fit_rank_scores <- function(pheno, covariates, num_tau_levels, num_cores = 1) {
  if (num_cores > 1 && requireNamespace("parallel", quietly = TRUE)) {
    mclapply(1:num_tau_levels, function(i) {
      tau <- i / (num_tau_levels + 1)
      calculate_rank_score(pheno, covariates, tau)
    }, mc.cores = num_cores)
  } else {
    lapply(1:num_tau_levels, function(i) {
      tau <- i / (num_tau_levels + 1)
      calculate_rank_score(pheno, covariates, tau)
    })
  }
}

#' Calculate Integrated Rank Score
#' @param rank_scores list of rank scores
#' @param method character "equal" or "ivw"
#' @param num_tau_levels integer number of quantile levels
#' @return vector of integrated rank scores
#' @importFrom quadprog solve.QP
#' @importFrom stats cov
#' @noRd
### calculate_integrated_score considering odds and even number tau
calculate_integrated_score <- function(rank_scores, method = "equal", num_tau_levels) {
  if (num_tau_levels %% 2 == 0) {  # even
    mid_point <- num_tau_levels / 2
    lower_half <- 1:mid_point
    upper_half <- (mid_point + 1):num_tau_levels
    n_pairs <- mid_point
  } else {  # odds
    mid_point <- ceiling(num_tau_levels / 2)
    lower_half <- 1:(mid_point-1)       
    middle <- mid_point                  
    upper_half <- (mid_point+1):num_tau_levels  
    n_pairs <- length(lower_half)  
  }
  
  if (method == "equal") {
    int_rank_score <- 0
    if (num_tau_levels %% 2 == 0) {
      # even number tau
      for (i in 1:num_tau_levels) {
        weight <- ifelse(i > mid_point, 1, -1)
        int_rank_score <- int_rank_score + weight * rank_scores[[i]]
      }
    } else {
      # odds number tau
      for (i in lower_half) int_rank_score <- int_rank_score - rank_scores[[i]]
      # int_rank_score <- int_rank_score + 0 * rank_scores[[middle]]  
      for (i in upper_half) int_rank_score <- int_rank_score + rank_scores[[i]]
    }
  } else if (method == "ivw") {
    a_i_tau_diff_matrix <- matrix(0, nrow = length(rank_scores[[1]]), ncol = n_pairs)
    
    for (i in 1:n_pairs) {
      if (num_tau_levels %% 2 == 0) {
        a_i_tau_diff_matrix[, i] <- rank_scores[[upper_half[i]]] - rank_scores[[lower_half[i]]]
      } else {
        a_i_tau_diff_matrix[, i] <- rank_scores[[upper_half[i]]] - rank_scores[[lower_half[i]]]
      }
    }
    
    # Optimization steps
    Y_QI_Var_Cov <- cov(a_i_tau_diff_matrix)
    Dmat <- Y_QI_Var_Cov
    dvec <- rep(0, n_pairs)
    Amat <- cbind(rep(1, n_pairs), diag(n_pairs))
    bvec <- c(1, rep(0, n_pairs))
    
    qp <- solve.QP(Dmat, dvec, t(Amat), bvec, meq = 1)
    weights_vec <- qp$solution
    
    int_rank_score <- 0
    for (i in 1:n_pairs) {
      int_rank_score <- int_rank_score + weights_vec[i] * 
        (rank_scores[[upper_half[i]]] - rank_scores[[lower_half[i]]])
    }
  }
  
  return(int_rank_score / n_pairs)  
}

#' Main QUAIL Rank Score Pipeline
#' @param phenotype numeric vector of phenotype values
#' @param covariates matrix/data.frame of covariates
#' @param num_tau_levels integer number of quantile levels
#' @param method character "equal" or "ivw"
#' @param num_cores integer number of cores for parallel processing
#' @return data.frame with integrated rank scores
#' @export
QUAIL_rank_score_pipeline <- function(phenotype, covariates, 
                                    num_tau_levels = 19,
                                    method = "equal",
                                    num_cores = 1) {
  # Input validation
  start_time <- Sys.time()   
  # Calculate rank scores for all quantile levels
  if (!is.null(phenotype)) {
    if (!is.numeric(phenotype)) {
      if (is.data.frame(phenotype) || is.matrix(phenotype)) {
        phenotype <- as.numeric(phenotype[[1]])
      } else {
        stop("phenotype must be a numeric vector.")
      }
    }
  }  

  rank_scores <- fit_rank_scores(phenotype, covariates, num_tau_levels, num_cores)
  
  # Calculate integrated rank score
  int_rank_score <- calculate_integrated_score(rank_scores, method, num_tau_levels)
  end_time <- Sys.time()
  cat("\nTotal vQTL score runtime:", round(difftime(end_time, start_time, units = "secs"), 2), "seconds\n")
      
  return(int_rank_score)
}