#' Bayesian Multiple Regression with Mixture-of-Normals Prior
#'
#' This function performs Bayesian multiple regression with a mixture-of-normals prior using the `rcpp_mr_ash_rss` function from the C++ implementation.
#'
#' @param bhat Numeric vector of observed effect sizes (standardized).
#' @param shat Numeric vector of standard errors of effect sizes.
#' @param z Numeric vector of Z-scores.
#' @param R Numeric matrix of the correlation matrix.
#' @param var_y Numeric value of the variance of the outcome.
#' @param n Integer value of the sample size.
#' @param sigma2_e Numeric value of the error variance.
#' @param s0 Numeric vector of prior variances for the mixture components.
#' @param w0 Numeric vector of prior weights for the mixture components.
#' @param mu1_init Numeric vector of initial values for the posterior mean of the coefficients.
#' @param tol Numeric value of the convergence tolerance. Default is 1e-8.
#' @param max_iter Integer value of the maximum number of iterations. Default is 1e5.
#' @param update_w0 Logical value indicating whether to update the mixture weights. Default is TRUE.
#' @param update_sigma Logical value indicating whether to update the error variance. Default is TRUE.
#' @param compute_ELBO Logical value indicating whether to compute the Evidence Lower Bound (ELBO). Default is TRUE.
#' @param standardize Logical value indicating whether to standardize the input data. Default is FALSE.
#'
#' @return A list containing the following components:
#' \describe{
#'   \item{mu1}{Numeric vector of the posterior mean of the coefficients.}
#'   \item{sigma2_1}{Numeric vector of the posterior variance of the coefficients.}
#'   \item{w1}{Numeric matrix of the posterior assignment probabilities.}
#'   \item{sigma2_e}{Numeric value of the error variance.}
#'   \item{w0}{Numeric vector of the mixture weights.}
#'   \item{ELBO}{Numeric value of the Evidence Lower Bound (if `compute_ELBO = TRUE`).}
#' }
#'
#' @examples
#' # Generate example data
#' set.seed(985115)
#' n <- 350
#' p <- 16
#' sigmasq_error <- 0.5
#' zeroes <- rbinom(p, 1, 0.6)
#' beta.true <- rnorm(p, 1, sd = 4)
#' beta.true[zeroes] <- 0
#'
#' X <- cbind(matrix(rnorm(n * p), nrow = n))
#' X <- scale(X, center = TRUE, scale = FALSE)
#' y <- X %*% matrix(beta.true, ncol = 1) + rnorm(n, 0, sqrt(sigmasq_error))
#' y <- scale(y, center = TRUE, scale = FALSE)
#'
#' # Calculate sufficient statistics
#' XtX <- t(X) %*% X
#' Xty <- t(X) %*% y
#' yty <- t(y) %*% y
#'
#' # Set the prior
#' K <- 9
#' sigma0 <- c(0.001, .1, .5, 1, 5, 10, 20, 30, .005)
#' omega0 <- rep(1/K, K)
#'
#' # Calculate summary statistics
#' b.hat <- sapply(1:p, function(j) { summary(lm(y ~ X[, j]))$coefficients[-1, 1] })
#' s.hat <- sapply(1:p, function(j) { summary(lm(y ~ X[, j]))$coefficients[-1, 2] })
#' R.hat <- cor(X)
#' var_y <- var(y)
#' sigmasq_init <- 1.5
#'
#' # Run mr_ash_rss
#' out <- mr_ash_rss(b.hat, s.hat, R = R.hat, var_y = var_y, n = n,
#'                   sigma2_e = sigmasq_init, s0 = sigma0, w0 = omega0,
#'                   mu1_init = rep(0, ncol(X)), tol = 1e-8, max_iter = 1e5,
#'                   update_w0 = TRUE, update_sigma = TRUE, compute_ELBO = TRUE,
#'                   standardize = FALSE)
#' # In sample prediction correlations
#' cor(X%*%out1$mu1, y) # 0.9984064
#' @export
mr_ash_rss <- function(bhat, shat, z = numeric(0), R, var_y, n,
                       sigma2_e, s0, w0, mu1_init = numeric(0),
                       tol = 1e-8, max_iter = 1e5,
                       update_w0 = TRUE, update_sigma = TRUE,
                       compute_ELBO = TRUE, standardize = FALSE) {
  if (is.null(var_y)) var_y = Inf
  result <- rcpp_mr_ash_rss(bhat = bhat, shat = shat, z = z, R = R,
                            var_y = var_y, n = n, sigma2_e = sigma2_e,
                            s0 = s0, w0 = w0, mu1_init = mu1_init,
                            tol = tol, max_iter = max_iter,
                            update_w0 = update_w0, update_sigma = update_sigma,
                            compute_ELBO = compute_ELBO, standardize = standardize)
  
  return(result)
}