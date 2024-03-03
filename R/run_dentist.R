#' @title run_dentist: Detecting Errors iN analyses of summary staTISTics
#'
#' @description DENTIST (Detecting Errors iN analyses of summary staTISTics) is a quality control
#' tool for GWAS summary data. It uses linkage disequilibrium (LD) information from a reference
#' panel to identify and correct problematic variants by comparing observed GWAS statistics to
#' predicted values. It can detect errors in genotyping/imputation, allelic errors, and
#' heterogeneity between GWAS and LD reference samples.
#'
#' @param LDmat A square matrix of linkage disequilibrium (LD) values where the dimensions equal the length of the \code{zScore} vector.
#' @param nSample The number of samples used in the GWAS whose summary statistics are being analyzed.
#' @param zScore A numeric vector of Z-scores from the GWAS summary statistics.
#' @param pValueThreshold A numeric threshold for the p-value, below which variants are considered significant for quality control. Default is 5e-8.
#' @param propSVD A numeric value specifying the proportion of SVD components to retain in the analysis. Default is 0.5.
#' @param gcControl Logical; if \code{TRUE}, applies genomic control corrections. Default is FALSE.
#' @param nIter An integer specifying the number of iterations for the DENTIST algorithm. Default is 8.
#' @param gPvalueThreshold A numeric threshold for p-value for grouping variants into significant and null. Default is 0.05.
#' @param ncpus An integer specifying the number of CPU cores to use for parallel computation. Default is 1.
#' @param seed An integer seed for random number generation to ensure reproducible results. Default is 123.
#'
#' @return A list containing the following components:
#' \itemize{
#'   \item{\code{imputedZ}}{A numeric vector of imputed Z-scores for each marker.}
#'   \item{\code{rsq}}{A numeric vector of R-squared values for each marker, indicating the goodness of fit.}
#'   \item{\code{zScore_e}}{A numeric vector of adjusted Z-scores after error detection.}
#'   \item{\code{iterID}}{An integer vector indicating the iteration in which each marker passed the quality control.}
#'   \item{\code{groupingGWAS}}{A binary vector indicating whether each marker is considered problematic (1) or not (0).}
#' }
#'
#' @examples
#' \dontrun{
#' # Simulate some data for demonstration purposes
#' nMarkers <- 1000
#' LDmat <- matrix(runif(nMarkers^2), nrow = nMarkers)
#' diag(LDmat) <- 1 # Ensure the diagonal is 1 for LD matrix
#' zScore <- rnorm(nMarkers)
#'
#' results <- run_dentist(LDmat = LDmat, nSample = 5000, zScore = zScore)
#' }
#'
#' @export
run_dentist <- function(LDmat, nSample, zScore,
                        pValueThreshold = 5e-8, propSVD = 0.4, gcControl = FALSE,
                        nIter = 8, gPvalueThreshold = 0.05, ncpus = 1, seed = 999) {
  # Check that LDmat dimensions match the length of zScore
  if (!is.matrix(LDmat) || nrow(LDmat) != ncol(LDmat) || nrow(LDmat) != length(zScore)) {
    stop("LDmat must be a square matrix with dimensions equal to the length of zScore.")
  }

  results <- dentist_rcpp(LDmat, nSample, zScore,
                   pValueThreshold, propSVD, gcControl, nIter,
                   gPvalueThreshold, ncpus, seed)

  return(results)
}
