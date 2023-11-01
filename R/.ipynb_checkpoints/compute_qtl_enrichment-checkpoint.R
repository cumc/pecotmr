#' @title Implementation of enrichment analysis described in https://doi.org/10.1371/journal.pgen.1006646
#'
#' @description Largely follows from fastenloc https://github.com/xqwen/fastenloc 
#' but uses `susieR` fitted objects as input to estimate prior for use with `coloc` package (coloc v5, aka SuSiE-coloc).
#' The main differences are 1) now enrichment is based on all QTL variants whether or not they are inside signal clusters;
#' 2) Causal QTL are sampled from SuSiE single effects, not signal clusters;
#' 3) Allow a variant to be QTL for not only multiple conditions (eg cell types) but also multiple regions (eg genes).
#' Other minor improvements include 1) Make GSL RNG thread-safe; 2) Release memory from QTL binary annotation samples immediately after they are used. 
#' @details Uses output of \code{\link[susieR]{susie}} from the
#'   \code{susieR} package.
#'
#' @param gwas_pip This is a vector of GWAS PIP, genome-wide. 
#' @param susie_qtl_regions This is a list of SuSiE fitted objects per QTL unit analyzed 
#' @param pi_gwas This parameter is highly important if GWAS input does not contain all SNPs interrogated (e.g., in some cases, only fine-mapped geomic regions are included). 
#' Then users must pick a value of total_variants and estimate pi_gwas beforehand by: sum(gwas_pip$pip)/total_variants
#' @param pi_qtl This parameter can be safely left to default if your input QTL data has enough regions to estimate it.
#' @param lambda Similar to the shrinkage parameter used in ridge regression. It takes any non-negative value and shrinks the enrichment estimate towards 0. 
#' When it is set to 0, no shrinkage will be applied. A large value indicates strong shrinkage. The default value is set to 1.0.
#' @param ImpN Rounds of multiple imputation to draw QTL from, default is 25.
#' @param num_threads Number of Simultaneous running CPU threads for multiple imputation, default is 1.
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
#'ImpN <- 10
#'lambda <- 1
#'num_threads <- 1
#'library(pecotmr)
#'en <- compute_qtl_enrichment(gwas_fit, susie_fits, lambda = lambda, ImpN = ImpN, num_threads = num_threads)
#' 
#' @seealso \code{\link[susieR]{susie}}
#' @useDynLib pecotmr
#' @importFrom Rcpp evalCpp 
#' @export
#' 
compute_qtl_enrichment <- function(gwas_pip, susie_qtl_regions, 
                           pi_gwas = NULL, pi_qtl = NULL, 
                           lambda = 1.0, ImpN = 25,
                           num_threads = 1) {
if (is.null(pi_gwas)) {
  warning("Using data to estimate pi_gwas. This will be problematic if your input gwas_pip does not contain genome-wide variants.")
  pi_gwas = sum(gwas_pip$pip) / length(gwas_pip$pip)
  cat(paste("Estimated pi_gwas is:", round(pi_gwas, 5), "\n"))
}
if (is.null(pi_qtl)) {
  warning("Using data to estimate pi_qtl. This will be problematic if either 1) your input susie_qtl_regions is not genome-wide, or 2) your single effects only includes variables inside of credible sets or signal clusters.")
  num_signal = 0
  num_test = 0
  for (d in susie_qtl_regions) {
    num_signal = num_signal + sum(d$pip)
    num_test = num_test + length(d$pip)
  }
  pi_qtl = num_signal/num_test
  cat(paste("Estimated pi_qtl is:", round(pi_qtl, 5), "\n"))
}
if (pi_gwas == 0) stop("Cannot perform enrichment analysis because there is no association signal from GWAS")
if (pi_qtl == 0) stop("Cannot perform enrichment analysis because there is no QTL associated with molecular phenotype")

# FIXME: need to check SNP name overlapping here? Need enough SNP names to overlap

en <- qtl_enrichment_rcpp(r_gwas_pip = gwas_pip$pip, 
                         r_qtl_susie_fit = susie_qtl_regions,
                         pi_gwas = pi_gwas,
                         pi_qtl = pi_qtl,
                         ImpN = ImpN,
                         shrinkage_lambda = lambda,
                         num_threads = num_threads)
return(en)
}
