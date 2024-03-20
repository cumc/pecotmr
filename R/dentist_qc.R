#' @title Detecting Errors iN analyses of summary staTISTics
#'
#' @description DENTIST (Detecting Errors iN analyses of summary staTISTics) is a quality control
#' tool for GWAS summary data. It uses linkage disequilibrium (LD) information from a reference
#' panel to identify and correct problematic variants by comparing observed GWAS statistics to
#' predicted values. It can detect errors in genotyping/imputation, allelic errors, and
#' heterogeneity between GWAS and LD reference samples.
#'
#' @param zScore A numeric vector of Z-scores from the GWAS summary statistics.
#' @param LDmat A square matrix of linkage disequilibrium (LD) values where the dimensions equal the length of the \code{zScore} vector.
#' @param nSample The number of samples used in the GWAS whose summary statistics are being analyzed.
#' @param pValueThreshold A numeric threshold for the p-value, below which variants are considered significant for quality control. Default is 5e-8.
#' @param propSVD A numeric value specifying the proportion of SVD components to retain in the analysis. Default is 0.5.
#' @param gcControl Logical; if \code{TRUE}, applies genomic control corrections. Default is FALSE.
#' @param nIter An integer specifying the number of iterations for the DENTIST algorithm. Default is 10.
#' @param gPvalueThreshold A numeric threshold for p-value for grouping variants into significant and null. Default is 0.05.
#' @param ncpus An integer specifying the number of CPU cores to use for parallel computation. Default is 1.
#' @param seed An integer seed for random number generation to ensure reproducible results. Default is 123.
#'
#' @return A data frame containing the following columns
#' \itemize{
#'   \item{\code{imputed_z}}{A numeric vector of imputed Z-scores for each marker.}
#'   \item{\code{rsq}}{A numeric vector of R-squared values for each marker, indicating the goodness of fit.}
#'   \item{\code{corrected_z}}{A numeric vector of adjusted Z-scores after error detection.}
#'   \item{\code{iter_to_correct}}{An integer vector indicating the iteration in which each marker passed the quality control.}
#'   \item{\code{is_problematic}}{A binary vector indicating whether each marker is considered problematic (1) or not (0).}
#' }
#'
#' @examples
#' \dontrun{
#' # Simulate some data for demonstration purposes (FIXME: need more serious simulation this does not work)
#' nMarkers <- 1000
#' LDmat <- matrix(runif(nMarkers^2), nrow = nMarkers)
#' LDmat <- (LDmat + t(LDmat)) / 2 # Making the matrix symmetric
#' diag(LDmat) <- 1 # Ensure the diagonal is 1 for LD matrix
#' zScore <- rnorm(nMarkers)
#' zhat <- zScore %*% LDmat
#' results <- dentist(LDmat = LDmat, nSample = 5000, zScore = zhat)
#' }
#'
#' @export
dentist <- function(zScore, LDmat, nSample,
                    pValueThreshold = 5e-8, propSVD = 0.4, gcControl = FALSE,
                    nIter = 10, gPvalueThreshold = 0.05, ncpus = 1, seed = 999, 
                    correct_chen_et_al_bug = TRUE) {
  # Check that LDmat dimensions match the length of zScore
  if (!is.matrix(LDmat) || nrow(LDmat) != ncol(LDmat) || nrow(LDmat) != length(zScore)) {
    stop("LDmat must be a square matrix with dimensions equal to the length of zScore.")
  }

  # Define a custom condition to capture warnings
  warning_handler <- function(w) {
    # Check if the warning message matches the specified pattern
    if (grepl("Adjusted rsq_eigen value exceeding 1", w$message)) {
      # Convert the warning to an error
      stop(w$message)
    }
    # Otherwise, invoke the default warning handler
    invokeRestart("muffleWarning")
  }

  res <- tryCatch(
    {
      dentist_iterative_impute(
        LDmat, nSample, zScore,
        pValueThreshold, propSVD, gcControl, nIter,
        gPvalueThreshold, ncpus, seed, correct_chen_et_al_bug
      )
    },
    warning = warning_handler
  )
  res <- as.data.frame(res)
  res$original_z <- zScore
  #outlier_test <- function(imp_z, orignal_z) {
    # FIXME: add the final QC step here
    #return (res)
  #}
  #res <- res %>% mutate(outlier=outlier_test(imp_z, orignal_z,...whatever it is that you define as outlier_test function)
  return(res)
}
