context("regularized_regression")
library(tidyverse)
library(snpStats)

generate_X_Y <- function(seed=1, num_samples=10, num_features=10, X_rownames=TRUE, y_rownames=TRUE) {
  set.seed(seed)
  X <- scale(
    matrix(rnorm(num_samples * num_features), nrow = num_samples),
    center = TRUE, scale = TRUE)
  
  if (X_rownames) {
    rownames(X) <- paste0("sample", 1:num_samples)
  } else {
    rownames(X) <- NULL
  }
  
  beta = rep(0, num_features)
  beta[1:4] = 1
  y <- X %*% beta + rnorm(num_samples)
  y <- matrix(y, nrow = num_samples, ncol = 1)
  if (y_rownames) {
    rownames(y) <- paste0("sample", 1:num_samples)
  } else {
    rownames(y) <- NULL
  }
  colnames(y) <- c("Outcome")
  
  return(list(X=X, Y=y))
}

generate_susie_obj <- function(X, y) {
    return(susie(X,y, L=10,
        max_iter=500,
        estimate_residual_variance=TRUE,
        estimate_prior_variance=TRUE,
        refine=TRUE,
        compute_univariate_zscore=FALSE,
        min_abs_corr=0.5,
        coverage=0.95))
}

test_that("Check susie_weights susie_fit works as expected", {
    local_mocked_bindings(
        susie_wrapper = function(X, y, ...) list(pip = rep(0.5, 10)),
        coef.susie = function(X, y, ...) rep(TRUE, 10)
    )
    expect_equal(
        length(susie_weights()), 10)
    expect_equal(
        length(susie_weights(susie_fit=list(alpha=T, mu=T, X_column_scale_factors=T))), 9)
})

test_that("Check susie_weights susie_fit runs with sim data no susie_fit", {
    sim <- generate_X_Y(seed = 1)
    X <- sim$X
    y <- sim$Y
    expect_equal(
        length(susie_weights(X=X, y=y)), 10)
})

test_that("Check susie_weights susie_fit runs with sim data with susie_fit", {
    sim <- generate_X_Y(seed = 1)
    X <- sim$X
    y <- sim$Y
    susie_fit <- generate_susie_obj(X, y)
    expect_equal(
        length(susie_weights(X=X, y=y, susie_fit = susie_fit)), 10)
})

test_that("Check susie_weights produces equal output w/ and w/o susie_fit", {
    sim <- generate_X_Y(seed = 1)
    X <- sim$X
    y <- sim$Y
    susie_fit <- generate_susie_obj(X, y)
    res <- susie_weights(X=X, y=y)
    res_susie_fit <- susie_weights(X=X, y=y, susie_fit = susie_fit)
    expect_equal(res, res_susie_fit)
})

test_that("Check glmnet_weights runs", {
    sim <- generate_X_Y(seed = 1, num_samples=30)
    X <- sim$X
    y <- sim$Y
    res <- glmnet_weights(X, y, 0.5) 
    expect_equal(length(res), 10)
})

test_that("Check enet_weights works", {
    sim <- generate_X_Y(seed = 1, num_samples=30)
    X <- sim$X
    y <- sim$Y
    set.seed(1)
    res <- glmnet_weights(X, y, 0.5) 
    set.seed(1)
    enet <- enet_weights(X, y) 
    expect_equal(res, enet)
})

test_that("Check lasso_weights works", {
    sim <- generate_X_Y(seed = 1, num_samples=30)
    X <- sim$X
    y <- sim$Y
    set.seed(1)
    res <- glmnet_weights(X, y, 1) 
    set.seed(1)
    enet <- lasso_weights(X, y) 
    expect_equal(res, enet)
})

generate_mr_ash_inputs <- function(seed=1, n = 360, p = 16, multiple_y=FALSE, k=1) {
  set.seed(seed)
  sigmasq_error <- 0.5
  zeroes <- rbinom(p, 1, 0.6)

  X <- cbind(matrix(rnorm(n * p), nrow = n))
  X <- scale(X, center = TRUE, scale = FALSE)
  if (multiple_y) {
    beta.true <- matrix(rnorm(p * k, 1, sd = 4), nrow = p, ncol = k)
    beta.true[zeroes, ] <- 0
    y <- X %*% beta.true + matrix(rnorm(n * k, 0, sqrt(sigmasq_error)), nrow = n, ncol = k)
  } else {
    beta.true <- rnorm(p, 1, sd = 4)
    beta.true[zeroes] <- 0
    y <- X %*% matrix(beta.true, ncol = 1) + rnorm(n, 0, sqrt(sigmasq_error))
  }
  y <- scale(y, center = TRUE, scale = FALSE)

  # Calculate sufficient statistics
  XtX <- t(X) %*% X
  Xty <- t(X) %*% y
  yty <- t(y) %*% y

  # Set the prior
  K <- 9
  sigma0 <- c(0.001, .1, .5, 1, 5, 10, 20, 30, .005)
  omega0 <- rep(1 / K, K)

  # Calculate summary statistics
  b.hat <- sapply(1:p, function(j) {
    summary(lm(y ~ X[, j]))$coefficients[-1, 1]
  })
  s.hat <- sapply(1:p, function(j) {
    summary(lm(y ~ X[, j]))$coefficients[-1, 2]
  })
  R.hat <- cor(X)
  var_y <- var(y)
  sigmasq_init <- 1.5

  return(list(
    bhat = b.hat, shat = s.hat, R = R.hat, var_y=var_y, n=n,
    sigma2_e=sigmasq_init, s0 = sigma0, w0 = omega0,
    mu1_init = rep(0, ncol(X)), X = X, y=y
  ))
}

test_that("Check mr_ash_rss works", {
  data <- generate_mr_ash_inputs()
  res <- mr_ash_rss(data$bhat, data$shat, data$R, data$var_y, data$n,
    data$sigma2_e, data$s0, data$w0, mu1_init = numeric(0))
  expect_true(all(names(res) %in% c("mu1", "sigma2_1", "w1", "sigma2_e", "w0", "ELBO")))
})

test_that("Check mr_ash_rss error on ncpu", {
  data <- generate_mr_ash_inputs()
  expect_error(mr_ash_rss(data$bhat, data$shat, data$R, data$var_y, data$n,
    data$sigma2_e, data$s0, data$w0, mu1_init = numeric(0), ncpu=-1))
  })

test_that("Check mr_ash_rss works null var_y", {
  data <- generate_mr_ash_inputs()
  res <- mr_ash_rss(data$bhat, data$shat, data$R, NULL, data$n,
    data$sigma2_e, data$s0, data$w0, mu1_init = numeric(0)) 
  expect_true(all(names(res) %in% c("mu1", "sigma2_1", "w1", "sigma2_e", "w0", "ELBO")))
  })

test_that("Check mr_ash_rss_weights works", {
  data <- generate_mr_ash_inputs()
  input <- list(b = data$bhat, seb=data$shat, n=rep(data$n, ncol(data$X)))
  res <- mr_ash_rss_weights(input, data$R, data$var_y,
    data$sigma2_e, data$s0, data$w0, mu1_init = numeric(0)) 
  expect_true(length(res) == ncol(data$R))
})

test_that("Check prs_cs works", {
  data <- generate_mr_ash_inputs()
  maf <- rep(0.5, length(data$bhat))
  LD <- list(blk1 = data$R)
  res <- prs_cs(data$bhat, LD, data$n, maf = maf)
  expect_true(all(names(res) %in% c("beta_est", "psi_est", "sigma_est", "phi_est")))
})

test_that("Check prs_cs missing LD", {
  data <- generate_mr_ash_inputs()
  maf <- rep(0.5, length(data$bhat))
  LD <- NULL
  expect_error(prs_cs(data$bhat, LD, data$n, maf = maf))
})

test_that("Check prs_cs missing n", {
  data <- generate_mr_ash_inputs()
  maf <- rep(0.5, length(data$bhat))
  LD <- list(blk1 = data$R)
  expect_error(prs_cs(data$bhat, LD, NULL, maf = maf))
})

test_that("Check prs_cs bhat and maf mismatch", {
  data <- generate_mr_ash_inputs()
  maf <- rep(0.5, length(data$bhat) + 2)
  LD <- list(blk1 = data$R)
  expect_error(prs_cs(data$bhat, LD, n=data$n, maf = maf))
})

test_that("Check prs_cs bhat and LD mismatch", {
  data <- generate_mr_ash_inputs()
  alt_data <- generate_mr_ash_inputs(seed=2, n = 400, p = 20)
  maf <- rep(0.5, length(alt_data$bhat))
  LD <- list(blk1 = data$R)
  expect_error(prs_cs(alt_data$bhat, LD, n=alt_data$n, maf = maf))
})

test_that("Check prs_cs_weights works", {
  data <- generate_mr_ash_inputs()
  input <- list(b = data$bhat, seb=data$shat, n=data$n)
  maf <- rep(0.5, length(data$bhat))
  res <- prs_cs_weights(input, data$R)
  expect_true(length(res) == ncol(data$R))
})

test_that("Check sdpr works", {
  data <- generate_mr_ash_inputs()
  maf <- rep(0.5, length(data$bhat))
  LD <- list(blk1 = data$R)
  res <- sdpr(data$bhat, LD, data$n)
  expect_true(all(names(res) %in% c("beta_est", "h2")))
})

test_that("Check sdpr bhat and LD mismatch", {
  data <- generate_mr_ash_inputs()
  alt_data <- generate_mr_ash_inputs(seed=2, n = 400, p = 24)
  maf <- rep(0.5, length(data$bhat))
  LD <- list(blk1 = data$R)
  expect_error(sdpr(alt_data$bhat, LD, alt_data$n))
})

test_that("Check sdpr missing n", {
  data <- generate_mr_ash_inputs()
  maf <- rep(0.5, length(data$bhat))
  LD <- list(blk1 = data$R)
  expect_error(sdpr(data$bhat, LD, -1))
})

test_that("Check sdpr negative per_variant_sample_size", {
  data <- generate_mr_ash_inputs()
  maf <- rep(0.5, length(data$bhat))
  LD <- list(blk1 = data$R)
  expect_error(sdpr(data$bhat, LD, data$n, per_variant_sample_size=rep(-1, data$n)))
})

test_that("Check sdpr array violation", {
  data <- generate_mr_ash_inputs()
  maf <- rep(0.5, length(data$bhat))
  LD <- list(blk1 = data$R)
  expect_error(sdpr(data$bhat, LD, data$n, array=sample(1:10, data$n, replace=TRUE)))
})

test_that("Check sdpr_weights works", {
  data <- generate_mr_ash_inputs()
  input <- list(b = data$bhat, seb=data$shat, n=data$n)
  maf <- rep(0.5, length(data$bhat))
  res <- sdpr_weights(input, data$R)
  expect_true(length(res) == ncol(data$R))
})

test_that("Check mrmash_weights mrmash_fit, X, Y, missing", {
  expect_error(mrmash_weights())
})

test_that("Check mrmash_weights mrmash_fit missing", {
  data <- generate_mr_ash_inputs(multiple_y=TRUE, k=4)
  #res <- mrmash_weights(X=data$X, Y=data$y)
  #expect_true(length(res) == ncol(data$X))
})

test_that("Check mvsusie_weights mvsusie_fit, X, Y, missing", {
  expect_error(mvsusie_weights())
})

test_that("Check mvsusie_weights mvsusie_fit missing", {
  data <- generate_mr_ash_inputs(multiple_y=TRUE, k=4)
  res <- mvsusie_weights(X=data$X, Y=data$y)
  expect_true(nrow(res) == ncol(data$X))
})

bayes_input <- function(seed=1, n=100000, p=1000) {
  X <- matrix(rnorm(n), nrow = p)
  Z <- matrix(round(runif(p*3, 0, 0.8), 0), nrow = p)
  set1 <- sample(1:ncol(X), 5)
  set2 <- sample(1:ncol(X), 5)
  sets <- list(set1, set2)
  g <- rowSums(X[, c(set1, set2)])
  e <- rnorm(nrow(X), mean = 0, sd = 1)
  y <- g + e
  return(list(X=X, y=y, Z=Z, sets=sets))
}
test_that("Check bayes_alphabet_weights works", {
  data <- bayes_input()
  res1 <- bayes_alphabet_weights(data$X, data$y, 'bayesN')
  res2 <- bayes_alphabet_weights(data$X, data$y, 'bayesL')
  res3 <- bayes_alphabet_weights(data$X, data$y, 'bayesA')
  res4 <- bayes_alphabet_weights(data$X, data$y, 'bayesC')
  res5 <- bayes_alphabet_weights(data$X, data$y, 'bayesR')
  expect_true(length(res1) == ncol(data$X))
  expect_true(length(res2) == ncol(data$X))
  expect_true(length(res3) == ncol(data$X))
  expect_true(length(res4) == ncol(data$X))
  expect_true(length(res5) == ncol(data$X))
  res1 <- bayes_alphabet_weights(data$X, data$y, 'bayesN', Z = data$Z)
  res2 <- bayes_alphabet_weights(data$X, data$y, 'bayesL', Z = data$Z)
  res3 <- bayes_alphabet_weights(data$X, data$y, 'bayesA', Z = data$Z)
  res4 <- bayes_alphabet_weights(data$X, data$y, 'bayesC', Z = data$Z)
  res5 <- bayes_alphabet_weights(data$X, data$y, 'bayesR', Z = data$Z)
  expect_true(length(res1) == ncol(data$X))
  expect_true(length(res2) == ncol(data$X))
  expect_true(length(res3) == ncol(data$X))
  expect_true(length(res4) == ncol(data$X))
  expect_true(length(res5) == ncol(data$X))
})

test_that("Check bayes_alphabet_weights mismatch X, y", {
  data <- bayes_input()
  alt_data <- bayes_input(p=2000)
  expect_error(bayes_alphabet_weights(data$X, alt_data$y, 'bayesN'))
})

test_that("Check bayes_alphabet_weights mismatch X, Z", {
  data <- bayes_input()
  alt_data <- bayes_input(n=10000, p=2000)
  expect_error(bayes_alphabet_weights(data$X, NULL, 'bayesN', Z = alt_data$Z))
})

test_that("Check bayes_n_weights works", {
  data <- bayes_input()
  res <- bayes_n_weights(data$X, data$y, Z = data$Z)
  expect_true(length(res) == ncol(data$X))
})

test_that("Check bayes_l_weights works", {
  data <- bayes_input()
  res <- bayes_l_weights(data$X, data$y, Z = data$Z)
  expect_true(length(res) == ncol(data$X))
})

test_that("Check bayes_a_weights works", {
  data <- bayes_input()
  res <- bayes_a_weights(data$X, data$y, Z = data$Z)
  expect_true(length(res) == ncol(data$X))
})

test_that("Check bayes_c_weights works", {
  data <- bayes_input()
  res <- bayes_c_weights(data$X, data$y, Z = data$Z)
  expect_true(length(res) == ncol(data$X))
})

test_that("Check bayes_r_weights works", {
  data <- bayes_input()
  res <- bayes_r_weights(data$X, data$y, Z = data$Z)
  expect_true(length(res) == ncol(data$X))
})
