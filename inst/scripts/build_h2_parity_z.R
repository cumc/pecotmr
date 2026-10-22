# =============================================================================
# Rebuild tests/testthat/test_data/h2ParityZ.rds -- the z-scores the tests in
# test_h2EstimationWrappers.R compare against the original ldsc, HDL, LDER
# and gLDSC implementations.
# -----------------------------------------------------------------------------
# These are stored rather than re-simulated inside the test because the draw
# uses eigen(), and eigen() fixes eigenvectors only up to SIGN. The sign a
# given LAPACK build returns is arbitrary, so regenerating from a seed
# produced different data on linux-64 than on osx-arm64 / linux-aarch64 --
# and the expected numbers came from running the original implementations
# (ldsc's ldscore.regressions.Hsq, HDL, LDER, gldsc) on one platform's
# realization. That made three of those tests fail on linux-64 only.
#
# The draw is over the LD in tests/testthat/test_data/h2ParityLd.rds: the 14
# chr22 blocks the bundled ldScoreExample held until 2026-10-20, frozen there
# so the bundled examples can be rebuilt without touching the parity tests.
#
# Rerunning this script invalidates the expected values in
# test_h2EstimationWrappers.R. Re-derive them on the NEW vectors at the same
# time: the data and the expected numbers are a matched pair.
# =============================================================================

library(purrr)

fixture <- readRDS("tests/testthat/test_data/h2ParityLd.rds")

# Each block's correlation matrix, rebuilt from its stored strict upper
# triangle.
blockR <- function(b) {
    R <- diag(length(b$snpIdx))
    R[upper.tri(R)] <- b$upper
    R[lower.tri(R)] <- t(R)[lower.tri(R)]
    R
}

# EXACT copy of the generator the pinned values were produced from.
simZ <- function(fixture, h2, n, seed) {
    M <- nrow(fixture$snpInfo)
    perSnpVar <- rep(h2 / M, M)
    set.seed(seed)
    z <- numeric(M)
    for (b in fixture$blocks) {
        R <- blockR(b)
        idx <- b$snpIdx
        p <- length(idx)
        sigma <- n * (R %*% diag(perSnpVar[idx], p) %*% R) + R
        e <- eigen((sigma + t(sigma)) / 2, symmetric = TRUE)
        z[idx] <- as.vector(e$vectors %*% (sqrt(pmax(e$values, 0)) * rnorm(p)))
    }
    z
}

# Every (h2, n, seed) the test file asks for. Keep in step with the call
# sites of .h2ParityZ() in test_h2EstimationWrappers.R -- a missing entry is
# a hard error there, not a silent re-simulation.
cases <- c(
    list(
        c(0.4, 1e4, 42),
        c(0.15, 4e3, 7),
        c(0.6, 2.5e4, 99),
        c(0.3, 1e4, 5),
        c(0.6, 5e4, 9),
        c(0.4, 5e4, 3)
    ),
    # the HDL recovery test sweeps seeds 1:5 at (h2 = 0.4, n = 1e5)
    map(1:5, function(s) c(0.4, 1e5, s))
)
out <- list()
for (cs in cases) {
    key <- paste(cs[[1]], cs[[2]], cs[[3]], sep = "_")
    out[[key]] <- simZ(fixture, cs[[1]], cs[[2]], cs[[3]])
}
saveRDS(out, "tests/testthat/test_data/h2ParityZ.rds", compress = "xz")
cat("keys:", paste(names(out), collapse = " | "), "\n")
cat("length each:", unique(lengths(out)), "\n")
cat(
    "file size:",
    file.size("tests/testthat/test_data/h2ParityZ.rds"),
    "bytes\n"
)
