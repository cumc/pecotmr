#' @importFrom susieR coef.susie
#' @export
susie_weights <- function(X = NULL, y = NULL, susie_fit = NULL, ...) {
  if (is.null(susie_fit)) {
    # get susie_fit object
    susie_fit <- susie_wrapper(X, y, ...)
  }
  if ("alpha" %in% names(susie_fit) && "mu" %in% names(susie_fit) && "X_column_scale_factors" %in% names(susie_fit)) {
    # This is designed to cope with output from pecotmr::susie_post_processor()
    # We set intercept to 0 and later trim it off anyways
    susie_fit$intercept <- 0
    return(coef.susie(susie_fit)[-1])
  } else {
    return(rep(0, length(susie_fit$pip)))
  }
}

#' @importFrom mr.mash.alpha coef.mr.mash
#' @export
mrmash_weights <- function(...) {
  res <- mrmash_wrapper(...)
  return(coef.mr.mash(res)[-1, ])
}

#' @importFrom mvsusieR mvsusie coef.mvsusie create_mixture_prior
#' @export
mvsusie_weights <- function(mvsusie_fit = NULL, X = NULL, Y = NULL, prior_variance = NULL, residual_variance = NULL, L = 30, mvsusie_max_iter = 200, ...) {
  if (is.null(mvsusie_fit)) {
    message("Did not provide mvsusie_fit, fitting mvSuSiE now")
    if (is.null(X) || is.null(Y)) {
      stop("Both X and Y must be provided if mvsusie_fit is NULL.")
    }
    if (is.null(prior_variance)) prior_variance <- create_mixture_prior(R = ncol(Y))
    if (is.null(residual_variance)) residual_variance <- mr.mash.alpha:::compute_cov_flash(Y)

    mvsusie_fit <- mvsusie(
      X = X, Y = Y, L = L, prior_variance = prior_variance,
      residual_variance = residual_variance, precompute_covariances = F,
      compute_objective = T, estimate_residual_variance = F, estimate_prior_variance = T,
      estimate_prior_method = "EM", max_iter = mvsusie_max_iter,
      n_thread = 1, approximate = F
    )
  }
  return(coef.mvsusie(mvsusie_fit)[-1, ])
}

# Get a reasonable setting for the standard deviations of the mixture
# components in the mixture-of-normals prior based on the data (X, y).
# Input se is an estimate of the residual *variance*, and n is the
# number of standard deviations to return. This code is adapted from
# the autoselect.mixsd function in the ashr package.
#' @importFrom susieR univariate_regression
init_prior_sd <- function(X, y, n = 30) {
  res <- univariate_regression(X, y)
  smax <- 3 * max(res$betahat)
  seq(0, smax, length.out = n)
}

#' @importFrom glmnet cv.glmnet
#' @importFrom stats coef
glmnet_weights <- function(X, y, alpha) {
  eff.wgt <- matrix(0, ncol = 1, nrow = ncol(X))
  sds <- apply(X, 2, sd)
  keep <- sds != 0 & !is.na(sds)
  enet <- cv.glmnet(x = X[, keep], y = y, alpha = alpha, nfold = 5, intercept = T, standardize = F)
  eff.wgt[keep] <- coef(enet, s = "lambda.min")[2:(sum(keep) + 1)]
  return(eff.wgt)
}

#' @export
enet_weights <- function(X, y) glmnet_weights(X, y, 0.5)

#' @export
lasso_weights <- function(X, y) glmnet_weights(X, y, 1)

#' @examples
#' wgt.mr.ash <- mrash_weights(eqtl$X, eqtl$y_res, beta.init = lasso_weights(X, y))
#' @importFrom mr.ash.alpha mr.ash
#' @importFrom stats predict
#' @export
mrash_weights <- function(X, y, init_prior_sd = TRUE, ...) {
  args_list <- list(...)
  if (!"beta.init" %in% names(args_list)) {
    args_list$beta.init <- lasso_weights(X, y)
  }
  fit.mr.ash <- do.call("mr.ash", c(list(X = X, y = y, sa2 = if (init_prior_sd) init_prior_sd(X, y)^2 else NULL), args_list))
  predict(fit.mr.ash, type = "coefficients")[-1]
}
