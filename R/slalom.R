#' Slalom Function for Summary Statistics QC for Fine-Mapping Analysis
#'
#' Performs Approximate Bayesian Factor (ABF) analysis, identifies credible sets,
#' and annotates lead variants based on fine-mapping results. It computes p-values
#' from z-scores assuming a two-sided standard normal distribution.
#'
#' @param zScore Numeric vector of z-scores corresponding to each variant.
#' @param LDmat Square matrix representing linkage disequilibrium (LD) information
#'   between variants. Must have dimensions matching the length of `zScore`.
#' @param standard_error Optional numeric vector of standard errors corresponding
#'   to each z-score. If not provided, a default value of 1 is assumed for all variants.
#' @param abf_prior_variance Numeric, the prior effect size variance for ABF calculations.
#'   Default is 0.04.
#' @param nlog10p_dentist_s_threshold Numeric, the -log10 DENTIST-S P value threshold
#'   for identifying outlier variants for prediction. Default is 4.0.
#' @param r2_threshold Numeric, the r2 threshold for DENTIST-S outlier variants
#'   for prediction. Default is 0.6.
#' @param lead_variant_choice Character, method to choose the lead variant, either
#'   "pvalue" or "abf", with default "pvalue".
#' @return A list containing the annotated LD matrix with ABF results, credible sets,
#'   lead variant, and DENTIST-S statistics; and a summary dataframe with aggregate statistics.
#' @examples
#' # Assuming `zScore` is your vector of z-scores, `LDmat` is your LD matrix,
#' # and optionally `standard_error` is your vector of standard errors:
#' results <- slalom(zScore, LDmat, standard_error)
#' @export
#'
slalom <- function(zScore, LDmat, standard_error = rep(1, length(zScore)), abf_prior_variance = 0.04,
                   nlog10p_dentist_s_threshold = 4.0, r2_threshold = 0.6, lead_variant_choice = "pvalue") {
  if (!is.matrix(LDmat) || nrow(LDmat) != ncol(LDmat) || nrow(LDmat) != length(zScore)) {
    stop("LDmat must be a square matrix matching the length of zScore.")
  }

  pvalue <- 2 * pnorm(abs(zScore), lower.tail = FALSE)

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

  abf_results <- abf(zScore, standard_error, W = abf_prior_variance)
  lbf <- abf_results$lbf
  prob <- abf_results$prob

  get_cs <- function(prob, coverage = 0.95) {
    ordering <- order(prob, decreasing = TRUE)
    cumprob <- cumsum(prob[ordering])
    idx <- which(cumprob > coverage)[1]
    cs <- ordering[1:idx]
    return(cs)
  }

  cs <- get_cs(prob, coverage = 0.95)
  cs_99 <- get_cs(prob, coverage = 0.99)

  lead_idx <- if (lead_variant_choice == "pvalue") {
    which.min(pvalue)
  } else {
    which.max(prob)
  }

  r2 <- LDmat^2
  t_dentist_s <- (zScore - LDmat[, lead_idx] * zScore[lead_idx])^2 / (1 - r2[, lead_idx])
  t_dentist_s[t_dentist_s < 0] <- Inf
  nlog10p_dentist_s <- -log10(1 - pchisq(t_dentist_s, df = 1))
  outliers <- (r2[, lead_idx] > r2_threshold) & (nlog10p_dentist_s > nlog10p_dentist_s_threshold)

  n_r2 <- sum(r2[, lead_idx] > r2_threshold)
  n_dentist_s_outlier <- sum(outliers, na.rm = TRUE)
  max_pip <- max(prob)

  summary <- list(
    lead_pip_variant = lead_idx,
    n_total = length(zScore),
    n_r2 = n_r2,
    n_dentist_s_outlier = n_dentist_s_outlier,
    fraction = ifelse(n_r2 > 0, n_dentist_s_outlier / n_r2, 0),
    max_pip = max_pip,
    cs_95 = cs,
    cs_99 = cs_99
  )
  result <- as.data.frame(list(original_z = zScore, prob = prob, pvalue = pvalue, outliers = outliers, nlog10p_dentist_s = nlog10p_dentist_s))

  return(list(data = result, summary = summary))
}
