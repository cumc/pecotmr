#' @title Implementation of enrichment analysis described in https://doi.org/10.1371/journal.pgen.1006646
#'
#' @description Largely does the same thing as fastenloc https://github.com/xqwen/fastenloc 
#' but uses `susieR` objects as input and outputs parameters to use as prior with `coloc` package.
#'
#' @details Uses output of \code{\link[susieR]{susie}} from the
#'   \code{susieR} package.
#'
#' @param susie_gwas 
#' @param susie_qtl_regions 
#' @param pi_gwas
#' @param pi_qtl
#' @param ImpN
#' @param shrinkage_prior_variance
#' @param num_threads
#' @return A list of enrichment parameter estimates 
#'
#' @examples
#' 
#'# Simulate fake data for gwas_pip
#'n_gwas_pip <- 1000
#'gwas_pip <- runif(n_gwas_pip)
#'names(gwas_pip) <- paste0("snp", 1:n_gwas_pip)
#'gwas_fit <- list(pip=gwas_pip)
#'# Simulate fake data for a single SuSiEFit object
#'simulate_susiefit <- function(n, p) {
#'  pip <- runif(n)
#'  names(pip) <- paste0("snp", 1:n)
#'  alpha <- t(matrix(runif(n * p), nrow = n))
#'  alpha <- t(apply(alpha, 1, function(row) row / sum(row)))
#'  list(
#'    pip = pip,
#'    alpha = alpha,
#'    prior_variance = runif(p)
#'  )
#'}
#'# Simulate multiple SuSiEFit objects
#'n_susie_fits <- 2 
#'susie_fits <- replicate(n_susie_fits, simulate_susiefit(n_gwas_pip, 10), simplify = FALSE)
#'# Add these fits to a list, providing names to each element
#'names(susie_fits) <- paste0("fit", 1:length(susie_fits))
#'# Set other parameters
#'pi_gwas <- 0.01
#'pi_qtl <- 0.01
#'ImpN <- 10
#'shrinkage_prior_variance <- 1
#'num_threads <- 1
#'en <- compute_qtl_enrichment(gwas_fit, susie_fits, pi_gwas, pi_qtl, ImpN, shrinkage_prior_variance)
#' 
#' @seealso \code{\link[susieR]{susie}}
#' @useDynLib qtl_enrichment
#' 
#' @export
#' 
compute_qtl_enrichment <- function(susie_gwas, susie_qtl_regions, 
                           pi_gwas = NULL, pi_qtl = NULL, 
                           ImpN = 10, shrinkage_prior_variance = 0,
                           num_threads = 1) {
# FIXME: revisit this setting
if (is.null(pi_gwas)) {
  pi_gwas = sum(susie_gwas$pip) / length(susie_gwas$pip)
  cat(paste("Estimated pi_gwas is", pi_gwas))
}
if (is.null(pi_qtl)) {
  num_signal = 0
  num_test = 0
  for (d in susie_qtl_regions) {
    num_signal = num_signal + sum(d$pip)
    num_test = num_test + length(d$pip)
  }
  pi_qtl = num_signal/num_test
  cat(paste("Estimated pi_qtl is", pi_qtl))
}

en <- qtl_enrichment_rcpp(r_gwas_pip = susie_gwas$pip, 
                         r_qtl_susie_fit = susie_qtl_regions,
                         pi_gwas = pi_gwas,
                         pi_qtl = pi_qtl,
                         ImpN = ImpN,
                         shrinkage_prior = shrinkage_prior_variance,
                         num_threads = num_threads)
return(en)
}