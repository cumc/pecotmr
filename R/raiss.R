#' Robust and accurate imputation from summary statistics
#'
#' This function is a part of the statistical library for SNP imputation from:
#' https://gitlab.pasteur.fr/statistical-genetics/raiss/-/blob/master/raiss/stat_models.py
#' It is R implementation of the imputation model described in the paper by Bogdan Pasaniuc,
#' Noah Zaitlen, et al., titled "Fast and accurate imputation of summary
#' statistics enhances evidence of functional enrichment", published in
#' Bioinformatics in 2014.
#' @param ref_panel A data frame containing 'chr', 'pos', 'variant_id', 'A1', and 'A2'.
#' @param known_zscores A data frame containing 'chr', 'pos', 'variant_id', 'A1', 'A2', and 'Z' values.
#' @param LD_matrix A square matrix of dimension equal to the number of rows in ref_panel.
#' @param lamb Regularization term added to the diagonal of the LD_matrix in the RAImputation model.
#' @param rcond Threshold for filtering eigenvalues in the pseudo-inverse computation in the RAImputation model.
#' @param R2_threshold R square threshold below which SNPs are filtered from the output.
#' @param minimum_ld Minimum LD score threshold for SNP filtering.
#'
#' @return A data frame that is the result of merging the imputed SNP data with known z-scores.
#' @importFrom dplyr arrange
#' @export
#'
#' @examples
#' # Example usage (assuming appropriate data is available):
#' # result <- raiss(ref_panel, known_zscores, LD_matrix, lamb = 0.01, rcond = 0.01, R2_threshold = 0.6, minimum_ld = 5)
raiss <- function(ref_panel, known_zscores, LD_matrix, lamb = 0.01, rcond = 0.01, R2_threshold = 0.6, minimum_ld = 5) {
  # Check that ref_panel and known_zscores are both increasing in terms of pos
  if (is.unsorted(ref_panel$pos)|| is.unsorted(known_zscores$pos)) {
    stop("ref_panel and known_zscores must be in increasing order of pos.")
  }

  # Define knowns and unknowns
    knowns_id = intersect(known_zscores$variant_id, ref_panel$variant_id)
    knowns =which(ref_panel$variant_id %in% knowns_id)
    unknowns = which(!ref_panel$variant_id %in% knowns_id)
    if(is.data.frame(LD_matrix)){
        LD_matrix = as.matrix(LD_matrix)
    }
  # Extract zt, sig_t, and sig_i_t
  zt <- known_zscores$Z
  sig_t <- LD_matrix[knowns, knowns, drop = FALSE]
  sig_i_t <- LD_matrix[unknowns, knowns, drop = FALSE]

  # Call raiss_model
  results <- raiss_model(zt, sig_t, sig_i_t, lamb, rcond)

  # Format the results
  results <- format_raiss_df(results, ref_panel, unknowns)

  # Filter output
  results <- filter_raiss_output(results, R2_threshold, minimum_ld)

  # Merge with known z-scores
  results <- merge_raiss_df(results, known_zscores)

  return(results)
}

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
raiss_model <- function(zt, sig_t, sig_i_t, lamb=0.01, rcond=0.01, batch=TRUE, report_condition_number=FALSE) {
  sig_t_inv <- invert_mat_recursive(sig_t, lamb, rcond)
  if (!is.numeric(zt) || !is.numeric(sig_t) || !is.numeric(sig_i_t)) {
    stop("zt, sig_t, and sig_i_t must be numeric.")
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

#' @param imp is the output of raiss_model()
#' @param ref_panel is a data frame with columns 'chr', 'pos', 'variant_id', 'ref', and 'alt'.
format_raiss_df <- function(imp, ref_panel, unknowns) {
  result_df <- data.frame(
    chr = ref_panel[unknowns, 'chr'],
    pos = ref_panel[unknowns, 'pos'],
    variant_id = ref_panel[unknowns, 'variant_id'],
    A1 = ref_panel[unknowns, 'A1'],
    A2 = ref_panel[unknowns, 'A2'],
    Z = imp$mu,
    Var = imp$var,
    ld_score = imp$ld_score,
    condition_number = imp$condition_number,
    correct_inversion = imp$correct_inversion
  )

  # Specify the column order
  column_order <- c('chr', 'pos', 'variant_id', "A1", "A2", 'Z', 'Var', 'ld_score', 'condition_number', 
                    'correct_inversion')

  # Reorder the columns
  result_df <- result_df[, column_order]
  return(result_df)
}

merge_raiss_df <- function(raiss_df, known_zscores) {
  # Merge the data frames
  merged_df <- merge(raiss_df, known_zscores, by = c("chr", "pos", "variant_id", "A1", "A2"), all = TRUE)

  # Identify rows that came from known_zscores
  from_known <- !is.na(merged_df$Z.y) & is.na(merged_df$Z.x)

  # Set Var to -1 and ld_score to Inf for these rows
  merged_df$Var[from_known] <- -1
  merged_df$ld_score[from_known] <- Inf

  # If there are overlapping columns (e.g., Z.x and Z.y), resolve them
  # For example, use Z from known_zscores where available, otherwise use Z from raiss_df
  merged_df$Z <- ifelse(from_known, merged_df$Z.y, merged_df$Z.x)

  # Remove the extra columns resulted from the merge (e.g., Z.x, Z.y)
  merged_df <- merged_df[, !colnames(merged_df) %in% c("Z.x", "Z.y")]
  merged_df = arrange(merged_df, pos)
  return(merged_df)
}

filter_raiss_output <- function(zscores, R2_threshold = 0.6, minimum_ld = 5) {
  # Reset the index and subset the data frame
  zscores <- zscores[, c('chr', 'pos', 'variant_id', 'A1', 'A2', 'Z', 'Var', 'ld_score')]
  zscores$imputation_R2 <- 1 - zscores$Var

  # Count statistics before filtering
  NSNPs_bf_filt <- nrow(zscores)
  NSNPs_initial <- sum(zscores$imputation_R2 == 2.0)
  NSNPs_imputed <- sum(zscores$imputation_R2 != 2.0)
  NSNPs_ld_filt <- sum(zscores$ld_score < minimum_ld)
  NSNPs_R2_filt <- sum(zscores$imputation_R2 < R2_threshold)

  # Apply filters
  zscores <- zscores[zscores$imputation_R2 > R2_threshold & zscores$ld_score >= minimum_ld, ]
  NSNPs_af_filt <- nrow(zscores)

  # Print report
  cat("IMPUTATION REPORT\n")
  cat("Number of SNPs:\n")
  cat("before filter:", NSNPs_bf_filt, "\n")
  cat("not imputed:", NSNPs_initial, "\n")
  cat("imputed:", NSNPs_imputed, "\n")
  cat("filtered because of ld:", NSNPs_ld_filt, "\n")
  cat("filtered because of R2:", NSNPs_R2_filt, "\n")
  cat("after filter:", NSNPs_af_filt, "\n")
  return(zscores)
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

invert_mat <- function(mat, lamb, rcond) {
  tryCatch({
    # Modify the diagonal elements of mat
    diag(mat) <- 1 + lamb
    # Compute the pseudo-inverse
    mat_inv <- MASS::ginv(mat, tol = rcond)
    return(mat_inv)
  }, error = function(e) {
    # Second attempt with updated lamb and rcond in case of an error
    diag(mat) <- 1 + lamb * 1.1
    mat_inv <- MASS::ginv(mat, tol = rcond * 1.1)
    return(mat_inv)
  })
}

invert_mat_recursive <- function(mat, lamb, rcond) {
  tryCatch({
    # Modify the diagonal elements of mat
    diag(mat) <- 1 + lamb
    # Compute the pseudo-inverse
    mat_inv <- MASS::ginv(mat, tol = rcond)
    return(mat_inv)
  }, error = function(e) {
    # Recursive call with updated lamb and rcond in case of an error
    invert_mat(mat, lamb * 1.1, rcond * 1.1)
  })
}

invert_mat_eigen <- function(mat, tol = 1e-3) {
    
    eigen_mat <- eigen(mat)
    L <- which(cumsum(eigen_mat$values) / sum(eigen_mat$values) > 1-tol)[1]
    if (is.na(L)) {
      # all eigen values are extremely small
      stop("Cannot invert the input matrix because all its eigen values are negative or close to zero")
    }
    mat_inv <- eigen_mat$vectors[,1:L] %*% 
        diag(1/eigen_mat$values[1:L]) %*% 
        t(eigen_mat$vectors[,1:L])
    
    return(mat_inv)
}