#' @title Fast Centering and Scaling of Matrices
#'
#' @description Does the same thing as \code{scale(x)}, but faster.
#'
#' @details Uses \code{\link[matrixStats]{colSds}} from the
#'   \code{matrixStats} package.
#'
#' @param x A numeric matrix; \code{is.matrix(x)} should be \code{TRUE}.
#'
#' @return A matrix in which the columns are centered to have zero
#'   mean, and they are also scaled to have standard deviation of 1.
#'
#' @examples
#' X <- matrix(1:24,4,6)
#' Y <- scale_faster(X)
#' apply(Y,2,mean) # Should be all zeros
#' apply(Y,2,sd)   # Should be all ones
#'
#' @seealso \code{\link{scale}}
#'
#' @importFrom matrixStats colSds
#' @importFrom Rcpp evalCpp
#'
#' @useDynLib fasterscale
#' 
#' @export
#' 
scale_faster <- function (x) {
  a <- colMeans(x)
  b <- colSds(x)
  
  # The scale_rcpp call should do the same as:
  #
  #   x <- t(t(x) - a)
  #   x <- t(t(x) / b)
  #
  scale_rcpp(x,a,b)
  return(x)
}
