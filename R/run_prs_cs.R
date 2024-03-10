#' Run PRS-CS: a polygenic prediction method that infers posterior SNP effect sizes under continuous shrinkage (CS) priors
#'
#' This function is a wrapper for the PRS-CS method implemented in C++. It takes marginal effect size estimates from regression and an external LD reference panel
#' and infers posterior SNP effect sizes using Bayesian regression with continuous shrinkage priors.
#'
#' @param bhat A vector of marginal effect sizes.
#' @param LD A list of LD blocks, where each element is a matrix representing an LD block.
#' @param n Sample size of the GWAS.
#' @param a Shape parameter for the prior distribution of psi. Default is 1.
#' @param b Scale parameter for the prior distribution of psi. Default is 0.5.
#' @param phi Global shrinkage parameter. If NULL, it will be estimated automatically. Default is NULL.
#' @param n_iter Number of MCMC iterations. Default is 1000.
#' @param n_burnin Number of burn-in iterations. Default is 500.
#' @param thin Thinning factor for MCMC. Default is 5.
#' @param maf A vector of minor allele frequencies, if available, will standardize the effect sizes by MAF. Default is NULL.
#' @param verbose Whether to print verbose output. Default is FALSE.
#' @param seed Random seed for reproducibility. Default is NULL.
#'
#' @return A list containing the posterior estimates:
#'   - beta_est: Posterior estimates of SNP effect sizes.
#'   - psi_est: Posterior estimates of psi (shrinkage parameters).
#'   - sigma_est: Posterior estimate of the residual variance.
#'   - phi_est: Posterior estimate of the global shrinkage parameter.
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
#' omega0 <- rep(1 / K, K)
#'
#' # Calculate summary statistics
#' b.hat <- sapply(1:p, function(j) {
#'   summary(lm(y ~ X[, j]))$coefficients[-1, 1]
#' })
#' s.hat <- sapply(1:p, function(j) {
#'   summary(lm(y ~ X[, j]))$coefficients[-1, 2]
#' })
#' R.hat <- cor(X)
#' var_y <- var(y)
#' sigmasq_init <- 1.5
#'
#' # Run PRS CS
#' maf <- rep(0.5, length(b.hat)) # fake MAF
#' LD <- list(blk1 = R.hat)
#' out <- prs_cs(b.hat, LD, n, maf = maf)
#' # In sample prediction correlations
#' cor(X %*% out$beta_est, y) # 0.9944553
#' @export
prs_cs <- function(bhat, LD, n,
                   a = 1, b = 0.5, phi = NULL,
                   maf = NULL, n_iter = 1000, n_burnin = 500,
                   thin = 5, verbose = FALSE, seed = NULL) {
  # Check input parameters
  if (missing(LD) || !is.list(LD)) {
    stop("Please provide a valid list of LD blocks using 'LD'.")
  }
  if (missing(n) || n <= 0) {
    stop("Please provide a valid sample size using 'n'.")
  }

  # Check if maf is provided and its length matches that of bhat
  if (!is.null(maf) && length(bhat) != length(maf)) {
    stop("The length of 'bhat' must be the same as 'maf'.")
  }

  # Check if the length of bhat matches the sum of the nrow of all elements in the LD list
  total_rows_in_LD <- sum(sapply(LD, nrow))
  if (length(bhat) != total_rows_in_LD) {
    stop("The length of 'bhat' must be the same as the sum of the number of rows of all elements in the 'LD' list.")
  }

  # Run PRS-CS
  result <- prs_cs_rcpp(
    a = a, b = b, phi = phi, bhat, maf,
    n = n, ld_blk = LD,
    n_iter = n_iter, n_burnin = n_burnin, thin = thin,
    verbose = verbose, seed = seed
  )

  # Return the result as a list
  list(
    beta_est = result$beta_est,
    psi_est = result$psi_est,
    sigma_est = result$sigma_est,
    phi_est = result$phi_est
  )
}
