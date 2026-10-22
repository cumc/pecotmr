# Shared fixtures for the TWAS weight-learner wrappers. Used by both
# test_regularizedRegressionWrappers.R (lasso/enet/scad/... + RSS solvers) and
# test_fineMappingWrappers.R (the SuSiE/mvSuSiE/fSuSiE weight extractors), so
# they live in a helper that testthat auto-loads for every test file.

# Individual-level (X, y) fixture: n samples x p variants with a sparse signal.
.rrwXy <- function(n = 50, p = 6, seed = 1) {
    set.seed(seed)
    X <- matrix(
        rnorm(n * p),
        n,
        p,
        dimnames = list(NULL, sprintf("chr1:%d:A:G", 100L * seq_len(p)))
    )
    b <- rnorm(p)
    b[-(1:2)] <- 0
    y <- as.numeric(X %*% b + rnorm(n))
    list(X = X, y = y, n = n, p = p)
}

# Single-context summary statistics + LD for the RSS weights contract:
# `stat` names its variants, and `ldReference` holds the same LD as an
# LdMatrix for the weight functions to match them to.
.rrwStatLd <- function(n = 50, p = 6, seed = 1) {
    d <- .rrwXy(n, p, seed)
    bhat <- map_dbl(
        seq_len(p),
        function(j) summary(lm(d$y ~ d$X[, j]))$coefficients[2, 1]
    )
    sehat <- map_dbl(
        seq_len(p),
        function(j) summary(lm(d$y ~ d$X[, j]))$coefficients[2, 2]
    )
    zhat <- bhat / sehat
    LD <- cor(d$X)
    stat <- list(
        b = bhat,
        seb = sehat,
        z = zhat,
        cor = bhat,
        n = rep(n, p),
        var_y = 1,
        variantNames = colnames(LD)
    )
    list(stat = stat, LD = LD, ldReference = LdMatrix(LD), p = p, n = n)
}

# Multi-context (variants x conditions) fixture for mr.mash / mvSuSiE wrappers.
.rrwMulti <- function(n = 60, p = 6, K = 3, seed = 2) {
    set.seed(seed)
    X <- matrix(
        rnorm(n * p),
        n,
        p,
        dimnames = list(NULL, sprintf("chr1:%d:A:G", 100L * seq_len(p)))
    )
    B <- matrix(0, p, K)
    B[1, ] <- rnorm(K, sd = 2)
    B[2, ] <- rnorm(K, sd = 2)
    Y <- X %*% B + matrix(rnorm(n * K), n, K)
    colnames(Y) <- paste0("ctx", seq_len(K))
    Z <- exec(
        cbind,
        !!!map(seq_len(K), function(k) {
            map_dbl(seq_len(p), function(j) {
                summary(lm(Y[, k] ~ X[, j]))$coefficients[2, 3]
            })
        })
    )
    colnames(Z) <- colnames(Y)
    list(
        X = X,
        Y = Y,
        LD = cor(X),
        ldReference = LdMatrix(cor(X)),
        stat = list(z = Z, n = n, variantNames = colnames(X)),
        p = p,
        K = K
    )
}

# Variant ids chr1:100:A:G, chr1:200:A:G, ... for p variants, and an LdMatrix
# over `R` named by them, for tests that build their own statistics and LD.
.rrwIds <- function(p) {
    sprintf("chr1:%d:A:G", 100L * seq_len(p))
}

.rrwLdMatrix <- function(R) {
    ids <- .rrwIds(nrow(R))
    LdMatrix(`dimnames<-`(R, list(ids, ids)))
}
