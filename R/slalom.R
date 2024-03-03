#' Slalom Function for Summary Statistics QC for Fine-Mapping Analysis
#'
#' Performs Approximate Bayesian Factor (ABF) analysis, identifies credible sets,
#' and annotates lead variants and DENTIST-S QC based on fine-mapping results.
#'
#' @param dat A data frame containing the fine-mapping results. Must include columns
#'   for variant identifiers, z-scores, and, optionally, standard errors (se) and p-values.
#' @param LD A matrix or data frame representing linkage disequilibrium (LD) information
#'   between variants. The structure and indexing should match the 'variant' column in `dat`.
#' @param abf_prior_variance Numeric, the prior effect size variance for ABF calculations.
#'   Default is 0.4, aligning with the default in the Python implementation.
#' @param nlog10p_dentist_s_threshold Numeric, the -log10 DENTIST-S P value threshold
#'   for identifying outlier variants for prediction. Default is 4.
#' @param r2_threshold Numeric, the r2 threshold for DENTIST-S outlier variants
#'   for prediction. Default is 0.6.
#' @details The function calculates the approximate Bayesian factor (ABF) for each variant,
#'   identifies credible sets at 95% and 99% coverage, and annotates the lead variant based
#'   on specified criteria. It also calculates and annotates DENTIST-S statistics to mark zOutliers. 
#' @return A list containing two elements: `dat`, the input dataframe annotated with ABF
#'   results, credible sets, lead variant, and DENTIST-S statistics; and `dat_summary`, a
#'   summary dataframe with aggregate statistics.
#' @examples
#' # Assuming `dat` is your dataframe with columns 'variant', 'z', 'se', and 'pvalue',
#' # and `LD` is your linkage disequilibrium matrix/dataframe:
#' results <- slalom(dat, LD)
#' @export
#'
slalom <- function(dat, LD, abf_prior_variance = 0.4, nlog10p_dentist_s_threshold = 4, r2_threshold = 0.6) {
  logSumExp <- function(x) {
    max_x <- max(x, na.rm = TRUE)
    sum_exp <- sum(exp(x - max_x), na.rm = TRUE)
    return(max_x + log(sum_exp))
  }
  
  abf <- function(z, se, W = 0.04) {
    V <- se^2
    r <- W / (W + V)
    lbf <- 0.5 * (log(1 - r) + (r * z^2))
    denom <- logSumExp(lbf)
    prob <- exp(lbf - denom)
    return(list(lbf = lbf, prob = prob))
  }
  
  get_cs <- function(variant, prob, coverage = 0.95) {
    ordering <- order(prob, decreasing = TRUE)
    cumprob <- cumsum(prob[ordering])
    idx <- which(cumprob > coverage)[1]
    cs <- variant[ordering][1:idx]
    return(cs)
  }

  # Calculate ABF and probabilities
  abf_results <- abf(dat$z, ifelse("se" %in% names(dat), dat$se, 1), W = abf_prior_variance)
  dat$lbf <- abf_results$lbf
  dat$prob <- abf_results$prob
  
  # Get credible sets
  cs <- get_cs(dat$variant, dat$prob, coverage = 0.95)
  cs_99 <- get_cs(dat$variant, dat$prob, coverage = 0.99)
  dat$cs <- dat$variant %in% cs
  dat$cs_99 <- dat$variant %in% cs_99
  
  # Identify lead variant
  lead_idx_snp <- if ("pvalue" %in% names(dat)) {
    which.min(dat$pvalue)
  } else {
    which.max(dat$prob)
  }
  lead_variant <- dat$variant[lead_idx_snp]
  dat$lead_variant <- FALSE
  dat$lead_variant[lead_idx_snp] <- TRUE
  
  # Annotate LD
  dat$r <- sapply(dat$variant, function(x) LD[lead_variant, x])
  
  lead_z <- dat$z[lead_idx_snp]
  dat$t_dentist_s <- (dat$z - dat$r * lead_z)^2 / (1 - dat$r^2)
  dat$t_dentist_s[dat$t_dentist_s < 0] <- Inf
  dat$t_dentist_s[lead_idx_snp] <- NA
  dat$nlog10p_dentist_s <- -log10(pchisq(dat$t_dentist_s, dat = 1, lower.tail = FALSE))
  dat$r2 <- dat$r^2
  dat$outliers <- (dat$r2 > r2_threshold) & (dat$nlog10p_dentist_s > nlog10p_dentist_s_threshold)
  
  # Summary
  n_r2 <- sum(dat$r2 > r2_threshold)
  n_dentist_s_outlier <- sum(dat$outliers)
  max_pip_idx <- which.max(dat$prob)
  
  dat_summary <- data.frame(
    lead_pip_variant = dat$variant[max_pip_idx],
    n_total = nrow(dat),
    n_r2 = n_r2,
    n_dentist_s_outlier = n_dentist_s_outlier,
    fraction = ifelse(n_r2 > 0, n_dentist_s_outlier / n_r2, 0),
    max_pip = max(dat$prob)
  )
  
  return(list(dat = dat, dat_summary = dat_summary))
}