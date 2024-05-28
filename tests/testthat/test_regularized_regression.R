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