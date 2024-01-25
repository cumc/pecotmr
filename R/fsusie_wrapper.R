#' @title  Calculate Purity Measures for Credible Sets
#'
#' As an extension of the internal cal_purity function. This function computes purity metrics (minimum, mean, and median absolute correlations)
#' for each credible set in a list of credible set indices, based on the provided X matrix.
#' The output Purity depends on the method specified: for the 'min' method,
#' it returns a single value for single-element sets or the minimum absolute correlation for others.
#' For other methods, it returns a vector of three values (min, mean, median) for each set.
#'
#' @param l_cs A list of credible set indices, where each element is a vector of indices
#'             corresponding to variables in a credible set.
#' @param X The data matrix used to compute correlations between variables in each credible set.
#' @param method A character string specifying the method to use for calculating purity.
#'               Defaults to 'min'. Other methods return a vector of min, mean, and median
#'               absolute correlations for each credible set.
#' @return A list where each element corresponds to a credible set and contains either a single
#'         purity value (for 'min' method and single-element sets) or a vector of purity metrics
#'         (for other methods and multi-element sets).
#' @export


cal_purity <- function(l_cs, X, method = "min") {
  tt <- list()
  
  for (k in 1:length(l_cs)) {
      cs_indices <- unlist(l_cs[[k]])  # Extract indices for the current credible set
       # Calculate purity based on the specified method
      if (method == "min") { 
          if (length(cs_indices) == 1 ) {
              tt[[k]] <- 1  # Set purity to 1 for non-"min" methods
            } else {
                   x <- abs(cor(X[, cs_indices]))  # Compute the absolute correlation matrix
              x[col(x) == row(x)] <- NA  # Set diagonal elements to NA to exclude them
            tt[[k]] <- min(x, na.rm = TRUE)  # Calculate minimum off-diagonal correlation for "min" method
    # Check if the credible set has only one element and the method is not "min"
    }} else {
           
          if (length(cs_indices) == 1 ) {
      tt[[k]] <- c(1,1,1)  # Set purity to 1 for non-"min" methods
    } else {
      x <- abs(cor(X[, cs_indices]))  # Compute the absolute correlation matrix
      x[col(x) == row(x)] <- NA  # Set diagonal elements to NA to exclude them
        # Calculate min, mean, and median of off-diagonal correlations for other methods
        tt[[k]] <- c( min(x , na.rm = TRUE), 
                     mean(x, na.rm = TRUE), 
                      median(x, na.rm = TRUE))
      }
    }
  }

  return(tt)
}


#'  @title Create Sets Similar to SuSiE Output from susiF Object
#'
#' This function constructs a list that mimics the structure of SuSiE output sets
#' from a susiF object. It includes credible sets (cs) with their names, a purity
#' dataframe, coverage information, and the requested coverage level.
#'
#' @param susiF.obj A susiF object containing the results from a susiF analysis.
#' expected to at least have 'cs' and 'alpha' components.
#' @param requested_coverage A numeric value specifying the desired coverage level for the
#'  credible sets. This is purely for record purpose so should be
#'  manually ensured that it correctly reflect the actual coverage used. Defaults to 0.95.
#' @return A list containing named credible sets (cs), a dataframe of purity metrics
#'         (min.abs.corr, mean.abs.corr, median.abs.corr), an index of credible sets (cs_index),
#'         coverage values for each set, and the requested coverage level. Similar to the SuSiE set output
#' @export


create_set_susiF = function(susiF.obj, X ,requested_coverage = 0.95) {
  # Create 'cs' set with names
  cs_named <- setNames(object = susiF.obj$cs, nm = paste0("L", seq_along(susiF.obj$cs)))
  
  # Create 'purity' data frame without dplyr
  purity_df <- do.call(rbind, lapply(cal_purity(susiF.obj$cs,X = X,method = "susie"), function(x) as.data.frame(t(x))))
  rownames(purity_df) <- names(cs_named)
  colnames(purity_df) <- c("min.abs.corr","mean.abs.corr",	"median.abs.corr")

  # Create 'coverage' without purrr
  coverage_vector <- numeric(length(susiF.obj$alpha))
  for (i in seq_along(susiF.obj$alpha)) {
    alpha_i <- susiF.obj$alpha[[i]]
    cs_i <- susiF.obj$cs[[i]]
    coverage_vector[i] <- sum(alpha_i[cs_i])
  }
  
  # Combine all elements into a list
  sets <- list(
    cs = cs_named,
    purity = purity_df,
    cs_index = 1:length(susiF.obj$cs),
    coverage = coverage_vector,
    requested_coverage = requested_coverage
  )
  
  return(sets)
}



#' @title Wrapper for fsusie Function with Automatic Post-Processing
#'
#' This function serves as a wrapper for the fsusie function, facilitating
#' automatic post-processing such as removing dummy credible sets (cs) that don't meet
#' the minimum purity threshold and calculating correlations for the remaining cs.
#' The function parameters are identical to those of the susiF function.
#'
#' @param X Residual genotype matrix.
#' @param Y Response phenotype matrix.
#' @param pos Genomics position of phenotypes, used for specifying the wavelet model.
#' @param L The maximum number of the credible set.
#' @param prior method to generate the prior.
#' @param max_SNP_EM maximum number of SNP used for learning the prior.
#' @param cov_lev Coverage level for the credible sets.
#' @param max_scale numeric, define the maximum of wavelet coefficients used in the analysis (2^max_scale).
#'        Set 10 true by default.
#' @param min.purity Minimum purity threshold for credible sets to be retained.
#' @param ... Additional arguments passed to the fsusie function.
#' @return A modified fsusie object with the susie sets list, correlations for cs,
#'         and without the dummy cs that do not meet the minimum purity requirement.
#' @importFrom susiF.alpha cal_cor_cs susiF
#' @export
fsusie_wrapper <- function(X, Y, pos, L, prior, max_SNP_EM, cov_lev, min.purity,max_scale, ...) {
  # Run fsusie
  susif.obj <- susiF(X = X, Y = Y, pos = pos, L = L, prior = prior,
                       max_SNP_EM = max_SNP_EM, cov_lev = cov_lev,
                       min.purity = min.purity,max_scale = max_scale, ...)

  # Remove dummy cs based on purity threshold
  if (all(abs(as.numeric(susiF.obj$purity)) < min.purity)) {
    susiF.obj$cs <- list(NULL)
    susiF.obj$sets <- list(cs = list(NULL), requested_coverage = cov_lev)
    susiF.obj$cs_corr <- NULL  # Set cs correlations to NULL if no credible sets meet purity criteria
  } else {
    # Create sets and add correlation for CS if purity criteria are met
    susiF.obj$sets <- create_set_susiF(susiF.obj,X, requested_coverage = cov_lev)
    susiF.obj$cs_corr <- cal_cor_cs(susiF.obj, X)
  }

  return(susiF.obj)
}
