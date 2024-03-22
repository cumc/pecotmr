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
dentist_impute_single_window <- function(zScore, LDmat, nSample,
                                         pValueThreshold = 5e-8, propSVD = 0.4, gcControl = FALSE,
                                         nIter = 10, gPvalueThreshold = 0.05, ncpus = 1, seed = 999, correct_chen_et_al_bug = TRUE) {
  # Check that number of variants cannot be below 2000
  if (length(zScore) < 2000) {
    warning("The number of variants is below 2000. The function will not be executed.")
    return(NULL) # Return NULL to indicate no results
  }
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
  return(res)
}

minusLogPvalueChisq <- function(stat) {
  p <- pchisq(stat, df = 1, lower.tail = FALSE)
  return(-log10(p))
}

calculate_stat <- function(impOp_zScores, impOp_imputed, impOp_rsq) {
  (impOp_zScores - impOp_imputed)^2 / (1 - impOp_rsq)
}

outlier_test <- function(stat, lambda) {
  ifelse(minusLogPvalueChisq(stat / lambda) > -log10(5e-8), TRUE, FALSE)
}

divide_into_windows <- function(pos, window_size, correct_chen_et_al_bug) {
  windowStartIdx <- c()
  windowEndIdx <- c()
  fillStartIdx <- c()
  fillEndIdx <- c()
  input_pos_start <- min(pos)
  input_pos_end <- max(pos)
  pos_range <- input_pos_end - input_pos_start
  if (window_size <= 0 | (window_size >= pos_range & correct_chen_et_al_bug == TRUE)) {
    windowStartIdx <- 1
    windowEndIdx <- length(pos)
    fillStartIdx <- 1
    fillEndIdx <- length(pos)
  } else {
    for (i in 1:(2 * ceiling(pos_range / window_size))) {
      start_pointer <- ((i - 1) * window_size) * 0.5 + 1
      end_pointer <- start_pointer + window_size
      start_idx <- which.min(abs(pos - start_pointer - input_pos_start))
      windowStartIdx <- c(windowStartIdx, start_idx)
      end_idx <- which.min(abs(pos - end_pointer - input_pos_start))
      windowEndIdx <- c(windowEndIdx, end_idx)
    }
    # avoid the situation where the last window is too small
    cutoff <- which.min(abs(pos - input_pos_end + window_size))
    closest_idx <- which.min(abs(windowStartIdx - cutoff))
    windowStartIdx <- windowStartIdx[windowStartIdx <= cutoff | seq_along(windowStartIdx) == closest_idx]
    windowEndIdx <- windowEndIdx[1:length(windowStartIdx)]

    if (length(windowStartIdx) == 1 & (windowEndIdx[1] < 4122 | correct_chen_et_al_bug == FALSE)) {
      # to avoid only 1 window
      windowEndIdx <- c(windowEndIdx, length(pos))
      start_pointer <- (pos[windowEndIdx] - pos[windowStartIdx]) * 0.25 + pos[windowStartIdx] # position
      start_idx <- which.min(abs(pos - start_pointer))
      windowStartIdx <- c(windowStartIdx, start_idx)
    }
    # decide the fillStartIdx and fillEndIdx
    if (length(windowStartIdx) == 1) { # for single window (correct_chen_et_al_bug must be TRUE)
      fillStartIdx <- 1
      fillEndIdx <- length(pos)
    } else { # for multiple windows
      for (i in 1:length(windowStartIdx)) {
        if (i == 1 & i != length(windowStartIdx)) {
          # for the first window and not the last window, fill Start from 1 and fill End is 0.75 quantile
          fill_start_idx <- 1
          fillStartIdx <- c(fillStartIdx, fill_start_idx)
          fill_end_pointer <- 0.75 * window_size + input_pos_start
          fill_end_idx <- which.max(pos[pos <= fill_end_pointer])
          fillEndIdx <- c(fillEndIdx, fill_end_idx)
        } else if (i != length(windowStartIdx)) {
          # for the windows in the middle, the start is the next variant of the previous fillendIdx, and the end is the 0.75 quantile
          fillStartIdx <- c(fillStartIdx, fillEndIdx[i - 1] + 1)
          fill_end_pointer <- pos[windowStartIdx[i]] + 0.75 * window_size
          fill_end_idx <- which.max(pos[pos <= fill_end_pointer])
          fillEndIdx <- c(fillEndIdx, fill_end_idx)
        }
      }
      # for the last window, the start is the next variant of the previous fillEndIdx, and the end is the last variant
      fillStartIdx <- c(fillStartIdx, fillEndIdx[length(fillEndIdx)] + 1)
      fillEndIdx <- c(fillEndIdx, length(pos))
    }
  }
  # combine window information into a data frame
  window_index <- 1:length(windowStartIdx)
  window_divided_res <- data.frame(
    windowIdx = window_index,
    windowStartIdx = windowStartIdx,
    windowEndIdx = windowEndIdx,
    fillStartIdx = fillStartIdx,
    fillEndIdx = fillEndIdx
  )
  # check if the divided windows are valid
  invalid_ends <- window_divided_res$fillStartIdx[1] != 1 | window_divided_res$fillEndIdx[nrow(window_divided_res)] != length(pos)
  if (nrow(window_divided_res) == 1 & invalid_ends) {
    stop("Invalid window divided!")
  } else if (nrow(window_divided_res) > 1) {
    invalid_middle <- 0
    for (i in 1:(nrow(window_divided_res) - 1)) {
      if (fillStartIdx[i + 1] - fillEndIdx[i] != 1) {
        invalid_middle <- invalid_middle + 1
      }
    }
    if (invalid_ends | invalid_middle != 0) {
      stop("Invalid window divided!")
    }
  }
  return(window_divided_res)
}


merge_windows <- function(imputed_result_by_window, window_divided_res) {
  if (length(imputed_result_by_window) != nrow(window_divided_res)) {
    stop("Different number of windows and imputed results!")
  }
  merged_results <- c()
  for (k in 1:nrow(window_divided_res)) {
    imputed_k <- imputed_result_by_window[[k]]
    imputed_k$index_within_window <- seq(1:nrow(imputed_k))
    imputed_k <- imputed_k %>%
      mutate(index_global = index_within_window + window_divided_res$windowStartIdx[k] - 1)
    extracted_results <- imputed_k %>%
      filter(index_global >= window_divided_res$fillStartIdx[k] & index_global <= window_divided_res$fillEndIdx[k])
    merged_results <- rbind(merged_results, extracted_results)
  }
  return(merged_results)
}

dentist_detect_outliers <- function(sum_stat, LDmat, nSample,
                                    window_size = 2000000, pValueThreshold = 5e-8, propSVD = 0.4, gcControl = FALSE,
                                    nIter = 10, gPvalueThreshold = 0.05, ncpus = 1, seed = 999, correct_chen_et_al_bug = TRUE) {
  # detect for column names and order by pos
  if (!any(tolower(c("pos", "position")) %in% tolower(colnames(sum_stat))) ||
    !any(tolower(c("z", "zscore")) %in% tolower(colnames(sum_stat)))) {
    stop("Input sum_stat is missing either 'pos'/'position' or 'z'/'zscore' column.")
  }
  sum_stat <- sum_stat %>% arrange(pos)
  ### FIXME: formalize this logic somewhere
  ### if correct_chen_et_al_bug = FALSE, then there can never be a single window
  # 1. if window_size is invalid (<=0 or >=range), and
  # 1.1. if correct_chen_et_al_bug==TRUE, run dentist single window, so just window (0-4122)
  # 1.2. if correct_chen_et_al_bug==FALSE, run in the dentist original way, so window (0-4122) and (742-4122)
  # 2. if window_size is valid, and
  # 2.1 if correct_chen_et_al_bug==TRUE, divide the window but don't add additional window if the first window covers all, so just (0-4122)
  # 2.2 if correct_chen_et_al_bug==FALSE, divide the window in the dentist original way, so window (0-4122) and (742-4122)
  if (window_size <= 0 | ((window_size >= max(sum_stat$pos) - min(sum_stat$pos) | is.na(window_size)) & (correct_chen_et_al_bug == TRUE))) {
    print("single window!")
    imputed_result <- dentist_impute_single_window(
      sum_stat$z, LDmat, nSample,
      pValueThreshold, propSVD, gcControl,
      nIter, gPvalueThreshold, ncpus, seed, correct_chen_et_al_bug
    )
  } else {
    # divide windows
    window_divided_res <- divide_into_windows(sum_stat$pos, window_size = window_size, correct_chen_et_al_bug = TRUE)
    print(window_divided_res)
    # compute dentist result for each window
    imputed_result_by_window <- list()
    for (k in 1:nrow(window_divided_res)) {
      zScore_k <- sum_stat$z[window_divided_res$windowStartIdx[k]:window_divided_res$windowEndIdx[k]]
      LDmat_k <- LDmat[
        window_divided_res$windowStartIdx[k]:window_divided_res$windowEndIdx[k],
        window_divided_res$windowStartIdx[k]:window_divided_res$windowEndIdx[k]
      ]
      imputed_result_by_window[[k]] <- dentist_impute_single_window(
        zScore_k, LDmat_k, nSample,
        pValueThreshold, propSVD, gcControl,
        nIter, gPvalueThreshold, ncpus, seed, correct_chen_et_al_bug
      )
    }
    # merge imputed result and generate a final imputed_result (similar to imputed_result above)
    imputed_result <- merge_windows(imputed_result_by_window, window_divided_res)
  }
  # detect outlier
  lambda_original <- 1
  imputed_result <- imputed_result %>%
    mutate(
      stat = calculate_stat(original_z, imputed_z, rsq),
      outlier = outlier_test(stat, lambda_original)
    ) %>%
    filter(!(imputed_z == 0 & rsq == 0))
  return(imputed_result)
}
