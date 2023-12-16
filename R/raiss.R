#' Robust and accurate imputation from summary statistics
#'
#' This function is a part of the statistical library for SNP imputation from:
#' https://gitlab.pasteur.fr/statistical-genetics/raiss/-/blob/master/raiss/stat_models.py
#' It is R implementation of the imputation model described in the paper by Bogdan Pasaniuc,
#' Noah Zaitlen, et al., titled "Fast and accurate imputation of summary
#' statistics enhances evidence of functional enrichment", published in
#' Bioinformatics in 2014.
#'
#' @param zt Vector of known Z scores.
#' @param sig_t Matrix of known linkage disequilibrium (LD) correlation.
#' @param sig_i_t Correlation matrix with rows corresponding to unknown SNPs (to impute)
#'               and columns to known SNPs.
#' @param lamb Regularization term added to the diagonal of the sig_t matrix.
#' @param rcond Threshold for filtering eigenvalues in the pseudo-inverse computation.
#' @param batch Boolean indicating whether batch processing is used.
#'
#' @return A list containing the variance 'var', estimation 'mu', LD score 'ld_score',
#'         condition number 'condition_number', and correctness of inversion
#'         'correct_inversion'.
#' @export
#'
#' @examples
#' # Example usage
#' # Define zt, sig_t, and sig_i_t with appropriate values
#' # result <- raiss_model(zt, sig_t, sig_i_t)
raiss_model <- function(zt, sig_t, sig_i_t, lamb=0.01, rcond=0.01, batch=TRUE, report_condition_number=FALSE) {
  # Translated content from Python function
  sig_t_inv <- invert_sig_t(sig_t, lamb, rcond)
  if (is.null(sig_t_inv)) {
    return(NULL)
  }

  if (batch) {
    
    condition_number <- if(report_condition_number) rep(kappa(sig_t, exact=T, norm="2"), nrow(sig_i_t)) else NA
    correct_inversion <- rep(check_inversion(sig_t, sig_t_inv), nrow(sig_i_t))
  } else {
    condition_number <- if(report_condition_number) kappa(sig_t, exact=T, norm="2") else NA
    correct_inversion <- check_inversion(sig_t, sig_t_inv)
  }

  var_ld_score <- compute_var(sig_i_t, sig_t_inv, lamb, batch)
  var <- var_ld_score$var
  ld_score <- var_ld_score$ld_score

  mu <- compute_mu(sig_i_t, sig_t_inv, zt)
  var_norm <- var_in_boundaries(var, lamb)

  R2 <- ((1 + lamb) - var_norm)
  mu <- mu / sqrt(R2)

  return(list(var=var_norm, mu=mu, ld_score=ld_score, condition_number=condition_number, correct_inversion=correct_inversion))
}

compute_mu <- function(sig_i_t, sig_t_inv, zt) {
  return(sig_i_t %*% (sig_t_inv %*% zt))
}

compute_var <- function(sig_i_t, sig_t_inv, lamb, batch=TRUE) {
  if (batch) {
    var <- (1 + lamb) - rowSums((sig_i_t %*% sig_t_inv) * sig_i_t)
    ld_score <- rowSums(sig_i_t^2)
  } else {
    var <- (1 + lamb) - (sig_i_t %*% (sig_t_inv %*% t(sig_i_t)))
    ld_score <- sum(sig_i_t^2)
  }
  return(list(var=var, ld_score=ld_score))
}

check_inversion <- function(sig_t, sig_t_inv) {
  return(all.equal(sig_t, sig_t %*% (sig_t_inv %*% sig_t), tolerance=1e-5))
}

var_in_boundaries <- function(var, lamb) {
  var[var < 0] <- 0
  var[var > (0.99999 + lamb)] <- 1
  return(var)
}

invert_sig_t <- function(sig_t, lamb, rcond) {
  diag(sig_t) <- 1 + lamb
  sig_t_inv <- MASS::ginv(sig_t, tol=rcond)
  return(sig_t_inv)
}