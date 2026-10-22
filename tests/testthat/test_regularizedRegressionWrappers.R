context("SS-TWAS: weights, pipeline, and omnibus combination")

# Previous TwasWeights-class tests used the legacy constructor
# `TwasWeights(weights = list(...), variantIds = ..., standardized = ...)`.
# The new `TwasWeights` is a DFrame collection class with (study,
# context, trait, method, entry) columns where each entry is a
# `TwasWeightsRow` S4 object carrying weights / fits / cvResult.
# Class-shape tests for the new collection should live alongside the
# pipeline tests and assert via accessors (`weights`, `studyName`,
# `cvResult`, etc.) — not against legacy slot shapes.
#
# `twasAnalysis()` was collapsed into the unified `twasZ()` dispatcher
# (task #37); its tests are removed here.
# `twasWeightsSumstatPipeline()` was removed without replacement in the
# S4 refactor (twasWeightsPipeline now dispatches directly on
# `QtlSumStats` / `QtlDataset` / `MultiStudyQtlDataset`).
#
# What remains: tests of the internal SuSiE-RSS weight extractors that
# are still present in `R/twasWeights.R` (`.susieRssExtractWeights`,
# `susieRssWeights`, `susieInfRssWeights`).

# =============================================================================
# SuSiE-RSS weight extraction
# =============================================================================

test_that("mrmashWeights fitRetention: slim omits the full fit, full keeps it", {
    skip_if_not_installed("mr.mashr")
    fakeFit <- list(w0 = c(a_1 = 0.5, a_2 = 0.5), V = diag(2))
    ddpm <- list(U = list(a = diag(2)))
    # Mock coef extraction so we exercise only the retain payload logic, not a
    # real mr.mash fit. coef.mr.mash(fit)[-1, ] -> drop the intercept row.
    local_mocked_bindings(
        coef.mr.mash = function(object, ...) rbind(c(0, 0), c(0.1, 0.2)),
        .package = "mr.mashr"
    )
    fitSlim <- attr(
        mrmashWeights(
            mrmashFit = fakeFit,
            fitRetention = "slim",
            dataDrivenPriorMatrices = ddpm
        ),
        "fit"
    )
    expect_setequal(names(fitSlim), c("dataDrivenPriorMatrices", "w0", "V"))
    expect_null(fitSlim$fit) # slim: no full fit
    expect_identical(fitSlim$dataDrivenPriorMatrices, ddpm)
    expect_identical(fitSlim$w0, fakeFit$w0)

    fitFull <- attr(
        mrmashWeights(
            mrmashFit = fakeFit,
            fitRetention = "full",
            dataDrivenPriorMatrices = ddpm
        ),
        "fit"
    )
    expect_true("fit" %in% names(fitFull)) # full: the whole fit retained
    expect_identical(fitFull$fit, fakeFit)
    expect_identical(fitFull$w0, fakeFit$w0) # slim fields still present
})


# =============================================================================
# Two-stage SuSiE-RSS fitting
# =============================================================================

# === Tests migrated from test_mrmashWrapper.R (mr.mash + glasso/glmnet coef helpers) ===

test_that("computeW0 returns uniform weights when ncomps == 1", {
    Bhat <- matrix(c(1, 0, 0, 2, 0, 0), nrow = 3, ncol = 2)
    result <- pecotmr:::computeW0(Bhat, ncomps = 1)
    expect_equal(result, 1)
})


test_that("computeW0 handles all-zero Bhat by returning uniform weights", {
    # When Bhat is all zero, prop_nonzero = 0
    # w0 = c(1, 0, ..., 0) => sum(w0 != 0) < 2 => fallback to uniform
    Bhat <- matrix(0, nrow = 5, ncol = 3)
    result <- pecotmr:::computeW0(Bhat, ncomps = 4)
    expect_equal(result, rep(1 / 4, 4))
    expect_equal(sum(result), 1)
})


test_that("computeW0 distributes weight based on nonzero rows when ncomps > 1", {
    # 2 out of 4 rows have nonzero entries
    Bhat <- matrix(0, nrow = 4, ncol = 2)
    Bhat[1, 1] <- 1
    Bhat[3, 2] <- 2
    result <- pecotmr:::computeW0(Bhat, ncomps = 3)
    expect_equal(length(result), 3)
    expect_equal(sum(result), 1, tolerance = 1e-10)
    # First element should be (1 - prop_nonzero) = 0.5
    expect_equal(result[1], 0.5)
})

# =========================================================================
# mrmashWrapper.R: rescaleCovW0 (lines 300-329)
# =========================================================================

test_that("rescaleCovW0 removes null component and renormalizes", {
    w0 <- c(
        null = 0.3,
        XtX_1 = 0.2,
        XtX_2 = 0.1,
        FLASH_1 = 0.15,
        FLASH_2 = 0.25
    )
    result <- pecotmr:::rescaleCovW0(w0)
    expect_false("null" %in% names(result))
    expect_equal(sum(result), 1, tolerance = 1e-10)
})


test_that("rescaleCovW0 handles all-zero non-null weights", {
    w0 <- c(null = 1.0, XtX_1 = 0, XtX_2 = 0, FLASH_1 = 0)
    result <- pecotmr:::rescaleCovW0(w0)
    # All non-null weights are zero -> equal weights
    expect_equal(sum(result), 1, tolerance = 1e-10)
    expect_true(all(result == result[1])) # all equal
})


test_that("rescaleCovW0 groups correctly by prior group prefix", {
    w0 <- c(
        null = 0.5,
        PCA_1 = 0.1,
        PCA_2 = 0.2,
        tFLASH_1 = 0.1,
        tFLASH_2 = 0.1
    )
    result <- pecotmr:::rescaleCovW0(w0)
    expect_true("PCA" %in% names(result))
    expect_true("tFLASH" %in% names(result))
    expect_equal(sum(result), 1, tolerance = 1e-10)
})

# =========================================================================
# mrmashWrapper.R: computeGrid (mr.mash's own exported grid builder)
# =========================================================================

test_that("computeGrid is mr.mash's grid, not a reimplementation of it", {
    skip_if_not_installed("mr.mashr")
    # gridMin/gridMax/autoselectMixsd used to be ported copies here, and
    # gridMin had drifted: upstream grid_min() is min(Shat)/10, the port
    # dropped the /10, so every grid started ten times too high.
    set.seed(3)
    bhat <- matrix(rnorm(20, sd = 2), nrow = 10, ncol = 2)
    sbhat <- matrix(abs(rnorm(20, mean = 0.5, sd = 0.1)), nrow = 10, ncol = 2)
    expect_equal(
        pecotmr:::computeGrid(bhat, sbhat),
        mr.mashr::autoselect.mixsd(
            list(Bhat = bhat, Shat = sbhat),
            mult = sqrt(2)
        )^2
    )
    # The ported floor would have been 10^2 times this one, since the grid
    # scales variances.
    ported <- min(sbhat)
    expect_lt(min(pecotmr:::computeGrid(bhat, sbhat)), ported^2)
})


test_that("compute_grid produces a valid grid from summary statistics", {
    set.seed(42)
    bhat <- matrix(rnorm(20, sd = 2), nrow = 10, ncol = 2)
    sbhat <- matrix(abs(rnorm(20, mean = 0.5, sd = 0.1)), nrow = 10, ncol = 2)
    result <- pecotmr:::computeGrid(bhat, sbhat)
    expect_true(is.numeric(result))
    expect_true(length(result) > 0)
    expect_true(all(result > 0))
})


test_that("compute_grid handles NA and zero sbhat values", {
    bhat <- c(1, 2, NA, 4, 5)
    sbhat <- c(0.5, 0, NA, 0.3, 0.8)
    result <- pecotmr:::computeGrid(bhat, sbhat)
    expect_true(is.numeric(result))
    expect_true(length(result) > 0)
})

# =========================================================================
# mrmashWrapper.R: mrmashWrapper input validation (lines 99-130)
# =========================================================================

# Note: Cannot mock requireNamespace via local_mocked_bindings because it is a
# base R function, not in pecotmr's namespace. We skip these tests and instead
# test the downstream validation that we CAN exercise.

test_that("mrmashWrapper errors when X and Y are not matrices", {
    skip_if_not_installed("glmnet")
    skip_if_not_installed("mr.mashr")
    expect_error(
        mrmashWrapper(data.frame(x = 1:3), matrix(1:6, nrow = 3, ncol = 2)),
        "matrices"
    )
})


test_that("mrmashWrapper errors when X and Y row counts differ", {
    skip_if_not_installed("glmnet")
    skip_if_not_installed("mr.mashr")
    expect_error(
        mrmashWrapper(
            matrix(1:6, nrow = 3, ncol = 2),
            matrix(1:8, nrow = 4, ncol = 2)
        ),
        "Assertion on 'Y'.*Must have exactly 3 rows"
    )
})


test_that("MrmashPriorParam refuses a priorGrid that is not a vector", {
    # The typed slot catches this when the record is built, so the wrapper's
    # own is.vector() guard became unreachable and was deleted.
    expect_error(
        MrmashPriorParam(priorGrid = matrix(1:4, nrow = 2)),
        'invalid object for slot "priorGrid"'
    )
    expect_equal(
        MrmashPriorParam(priorGrid = c(0.1, 0.2))$priorGrid,
        c(0.1, 0.2)
    )
})


test_that("mrmashWrapper errors when no prior matrices and canonical_priorMatrices is FALSE", {
    skip_if_not_installed("glmnet")
    skip_if_not_installed("mr.mashr")
    X <- matrix(rnorm(12), nrow = 3, ncol = 4)
    Y <- matrix(rnorm(6), nrow = 3, ncol = 2)
    expect_error(
        mrmashWrapper(
            X,
            Y,
            dataDrivenPriorMatrices = NULL,
            prior = MrmashPriorParam(canonicalPriorMatrices = FALSE)
        ),
        "dataDrivenPriorMatrices"
    )
})


test_that("mrmashWrapper warns when Y has missing and B_init_method is glasso", {
    skip_if_not_installed("glmnet")
    skip_if_not_installed("mr.mashr")
    set.seed(42)
    n <- 20
    p <- 5
    r <- 2
    X <- matrix(rnorm(n * p), nrow = n, ncol = p)
    Y <- matrix(rnorm(n * r), nrow = n, ncol = r)
    Y[1, 1] <- NA # introduce missing values
    colnames(Y) <- c("cond1", "cond2")

    # Should produce warning about glasso and NAs, then likely fail on the
    # downstream mr.mashr call, but the warning is what we test
    expect_warning(
        tryCatch(
            mrmashWrapper(
                X,
                Y,
                bInitMethod = "glasso",
                dataDrivenPriorMatrices = list(U = list(matrix(1, 2, 2))),
                prior = MrmashPriorParam(canonicalPriorMatrices = FALSE)
            ),
            error = function(e) NULL
        ),
        "glasso"
    )
})

# =========================================================================
# mrmashWrapper.R: computeCoefficientsGlasso (lines 211-240)
# =========================================================================

test_that("computeCoefficientsGlasso runs without Xnew", {
    skip_if_not_installed("glmnet")
    set.seed(42)
    n <- 50
    p <- 5
    r <- 3
    X <- matrix(rnorm(n * p), nrow = n, ncol = p)
    Y <- matrix(rnorm(n * r), nrow = n, ncol = r)
    colnames(Y) <- paste0("cond", seq_len(r))
    result <- pecotmr:::computeCoefficientsGlasso(
        X,
        Y,
        standardize = FALSE,
        numThreads = 1,
        Xnew = NULL
    )
    expect_true("Bhat" %in% names(result))
    expect_true("Ytrain" %in% names(result))
    expect_equal(nrow(result$Bhat), p)
    expect_equal(ncol(result$Bhat), r)
    expect_null(result$Yhat_new)
})


test_that("computeCoefficientsGlasso runs with Xnew", {
    skip_if_not_installed("glmnet")
    set.seed(42)
    n <- 50
    p <- 5
    r <- 3
    X <- matrix(rnorm(n * p), nrow = n, ncol = p)
    Y <- matrix(rnorm(n * r), nrow = n, ncol = r)
    colnames(Y) <- paste0("cond", seq_len(r))
    Xnew <- matrix(rnorm(10 * p), nrow = 10, ncol = p)
    result <- pecotmr:::computeCoefficientsGlasso(
        X,
        Y,
        standardize = FALSE,
        numThreads = 1,
        Xnew = Xnew
    )
    expect_true("Yhat_new" %in% names(result))
    expect_equal(nrow(result$Yhat_new), 10)
    expect_equal(ncol(result$Yhat_new), r)
    expect_equal(colnames(result$Yhat_new), colnames(Y))
})

# =========================================================================
# mrmashWrapper.R: computeCoefficientsUnivGlmnet (lines 243-281)
# =========================================================================

test_that("computeCoefficientsUnivGlmnet runs without Xnew", {
    skip_if_not_installed("glmnet")
    set.seed(42)
    n <- 60
    p <- 5
    r <- 2
    X <- matrix(rnorm(n * p), nrow = n, ncol = p)
    Y <- matrix(rnorm(n * r), nrow = n, ncol = r)
    colnames(Y) <- paste0("cond", seq_len(r))
    result <- pecotmr:::computeCoefficientsUnivGlmnet(
        X,
        Y,
        alpha = 0.5,
        standardize = FALSE,
        Xnew = NULL
    )
    expect_true("Bhat" %in% names(result))
    expect_true("intercept" %in% names(result))
    expect_equal(nrow(result$Bhat), p)
    expect_equal(ncol(result$Bhat), r)
    expect_null(result$Yhat_new)
})


test_that("computeCoefficientsUnivGlmnet runs with Xnew", {
    skip_if_not_installed("glmnet")
    set.seed(42)
    n <- 60
    p <- 5
    r <- 2
    X <- matrix(rnorm(n * p), nrow = n, ncol = p)
    Y <- matrix(rnorm(n * r), nrow = n, ncol = r)
    colnames(Y) <- paste0("cond", seq_len(r))
    Xnew <- matrix(rnorm(8 * p), nrow = 8, ncol = p)
    result <- pecotmr:::computeCoefficientsUnivGlmnet(
        X,
        Y,
        alpha = 0.5,
        standardize = FALSE,
        Xnew = Xnew
    )
    expect_true("Yhat_new" %in% names(result))
    expect_equal(nrow(result$Yhat_new), 8)
    expect_equal(ncol(result$Yhat_new), r)
    expect_equal(colnames(result$Yhat_new), colnames(Y))
})


test_that("computeCoefficientsUnivGlmnet handles NA in Y", {
    skip_if_not_installed("glmnet")
    set.seed(42)
    n <- 60
    p <- 5
    r <- 2
    X <- matrix(rnorm(n * p), nrow = n, ncol = p)
    Y <- matrix(rnorm(n * r), nrow = n, ncol = r)
    colnames(Y) <- paste0("cond", seq_len(r))
    Y[1:5, 1] <- NA # introduce missing values in one condition
    result <- pecotmr:::computeCoefficientsUnivGlmnet(
        X,
        Y,
        alpha = 0.5,
        standardize = FALSE,
        Xnew = NULL
    )
    expect_true("Bhat" %in% names(result))
    expect_equal(nrow(result$Bhat), p)
})

# =========================================================================
# mrmashWrapper.R: mrmashWrapper seed warning (line 107-108)
# =========================================================================

# Note: Cannot mock base::exists() via local_mocked_bindings.
# The seed-check message on line 107-108 would require removing .Random.seed
# from the global environment, which is not safe to do in tests.

# =============================================================================
# Real-fit coverage for the solver wrappers in R/regularizedRegressionWrappers.R
# -----------------------------------------------------------------------------
# The fine-mapping / TWAS pipelines MOCK these wrappers, so their bodies are
# otherwise untested. Here we drive each wrapper on a SMALL real fixture and
# assert the return shape (weight length == #variants; matrix for multivariate;
# attr(.,"fit") unless fitRetention is "none"). MCMC iterations kept tiny.
# fsusieWeights and the mock-based mrmash/mvsusie payload tests live in
# test_rrMrmashMvsusie.R and are not duplicated here.
# =============================================================================

# (Shared .rrwXy / .rrwStatLd / .rrwMulti fixtures live in helper-rrwFixtures.R
#  so test_fineMappingWrappers.R can reuse them for the SuSiE weight extractors.)

# -------------------------------- individual --------------------------------

# Capture what reaches cv.glmnet, which is all fitGlmnetCv() is responsible
# for -- the engine itself is glmnet's. The fake fit is not a usable
# cv.glmnet object, so an extractor downstream of the call errors; `try`
# swallows that, since only the captured arguments are under test.
.rrw_captureGlmnet <- function(expr) {
    captured <- NULL
    local_mocked_bindings(
        cv.glmnet = function(...) {
            captured <<- list(...)
            structure(list(), class = "cv.glmnet")
        },
        .package = "glmnet"
    )
    try(force(expr), silent = TRUE)
    captured
}

test_that("fitGlmnetCv: NULL formals leave glmnet's own defaults", {
    skip_if_not_installed("glmnet")
    f <- .rrwXy()
    args <- .rrw_captureGlmnet(fitGlmnetCv(f$X, f$y, alpha = 0.5))
    expect_equal(args$alpha, 0.5)
    # An explicit NULL would reach cv.glmnet and override its default with
    # nothing, so an unset formal must be ABSENT from the call.
    for (nm in c("family", "standardize", "intercept", "nfolds")) {
        expect_false(is_in(nm, names(args)))
    }
})

test_that("fitGlmnetCv: set formals use glmnet's own argument names", {
    skip_if_not_installed("glmnet")
    f <- .rrwXy()
    args <- .rrw_captureGlmnet(fitGlmnetCv(
        f$X,
        f$y,
        alpha = 1,
        family = "mgaussian",
        standardize = TRUE,
        intercept = FALSE,
        nfolds = 3
    ))
    expect_equal(args$family, "mgaussian")
    expect_true(args$standardize)
    expect_false(args$intercept)
    expect_equal(args$nfolds, 3)
    expect_false(is_in("nfold", names(args)))
})

test_that("fitGlmnetCv: methodArgs win over the caller's defaults", {
    skip_if_not_installed("glmnet")
    f <- .rrwXy()
    args <- .rrw_captureGlmnet(fitGlmnetCv(
        f$X,
        f$y,
        alpha = 0.5,
        nfolds = 5,
        methodArgs = GlmnetOptions(nfolds = 3)
    ))
    expect_equal(args$nfolds, 3)
    expect_error(
        fitGlmnetCv(f$X, f$y, alpha = 0.5, methodArgs = list(nfolds = 3)),
        "GlmnetOptions"
    )
})

test_that("GlmnetOptions: nfold is refused and names nfolds instead", {
    expect_error(GlmnetOptions(nfold = 5), "fold count is `nfolds`")
})

test_that("GlmnetOptions: parallel is refused outright", {
    # glmnet's `parallel` needs a foreach backend pecotmr does not declare,
    # and weight fitting already runs inside BiocParallel workers over the
    # same CV folds -- so there is no value a caller could usefully set.
    # cv.glmnet takes `...`, so nothing else would report it.
    expect_error(GlmnetOptions(parallel = TRUE), "cannot be set here")
    expect_error(GlmnetOptions(parallel = FALSE), "foreach backend")
    # A refusal must not swallow the names that ARE settable.
    expect_equal(GlmnetOptions(nfolds = 3)$nfolds, 3)
})

test_that("fitGlmnetCv: the engine is never asked to parallelise", {
    skip_if_not_installed("glmnet")
    f <- .rrwXy()
    args <- .rrw_captureGlmnet(fitGlmnetCv(f$X, f$y, alpha = 0.5))
    expect_false(is_in("parallel", names(args)))
})

test_that("glmnetWeights: the fold count reaches the engine as nfolds", {
    skip_if_not_installed("glmnet")
    f <- .rrwXy()
    args <- .rrw_captureGlmnet(glmnetWeights(f$X, f$y, alpha = 0.5))
    expect_equal(args$nfolds, 5)
    expect_false(args$standardize)
    expect_true(args$intercept)
    # Regression: pecotmr passed `nfold = 5`, which reached cv.glmnet only by
    # partial matching, so a user who spelled the engine's own `nfolds`
    # supplied BOTH and glmnet failed with a message about `family`.
    w <- glmnetWeights(f$X, f$y, 0.5, methodArgs = GlmnetOptions(nfolds = 3))
    expect_length(as.numeric(w), f$p)
})

test_that("mr.mash glmnet initialisation routes through fitGlmnetCv", {
    skip_if_not_installed("glmnet")
    f <- .rrwXy()
    Y <- cbind(f$y, f$y + rnorm(length(f$y)))
    gl <- .rrw_captureGlmnet(computeCoefficientsGlasso(
        f$X,
        Y,
        standardize = TRUE,
        numThreads = 1L
    ))
    expect_equal(gl$family, "mgaussian")
    expect_equal(gl$alpha, 1)
    expect_true(gl$standardize)
    # `parallel = FALSE` was a restatement of cv.glmnet's own default.
    expect_false(is_in("parallel", names(gl)))
    uv <- .rrw_captureGlmnet(pecotmr:::.linreg(
        1L,
        f$X,
        Y,
        alpha = 0.5,
        standardize = FALSE,
        Xnew = NULL
    ))
    expect_equal(uv$family, "gaussian")
    expect_equal(uv$alpha, 0.5)
    expect_false(uv$standardize)
    expect_false(is_in("parallel", names(uv)))
})

test_that("lassoWeights / enetWeights (glmnet) return length-p weights", {
    skip_if_not_installed("glmnet")
    f <- .rrwXy()
    expect_length(as.numeric(lassoWeights(f$X, f$y)), f$p)
    expect_length(as.numeric(enetWeights(f$X, f$y)), f$p)
})

test_that("scadWeights / mcpWeights (ncvreg) return length-p weights", {
    skip_if_not_installed("ncvreg")
    f <- .rrwXy()
    expect_length(as.numeric(scadWeights(f$X, f$y)), f$p)
    expect_length(as.numeric(mcpWeights(f$X, f$y)), f$p)
})

test_that("l0learnWeights returns length-p weights", {
    skip_if_not_installed("L0Learn")
    f <- .rrwXy()
    expect_length(as.numeric(l0learnWeights(f$X, f$y)), f$p)
})

test_that("mrashWeights returns length-p weights and can retain the fit", {
    skip_if_not_installed("susieR")
    skip_if_not_installed("glmnet")
    f <- .rrwXy()
    # mr.ash's own default grid is too narrow for this fixture's two strong
    # effects, and mr.ash says so.
    expect_warning(
        w <- mrashWeights(f$X, f$y, fitRetention = "slim"),
        "larger setting of the prior variance"
    )
    expect_length(w, f$p)
    expect_false(is.null(attr(w, "fit")))
})

test_that("qgg Bayes-alphabet weights (N/L/A/C/R) return length-p weights", {
    skip_if_not_installed("qgg")
    f <- .rrwXy()
    mc <- list(nit = 200, nburn = 20, nthin = 1)
    expect_length(exec(bayesNWeights, !!!c(list(f$X, f$y), mc)), f$p)
    expect_length(exec(bayesLWeights, !!!c(list(f$X, f$y), mc)), f$p)
    expect_length(exec(bayesAWeights, !!!c(list(f$X, f$y), mc)), f$p)
    expect_length(exec(bayesCWeights, !!!c(list(f$X, f$y), mc)), f$p)
    expect_length(exec(bayesRWeights, !!!c(list(f$X, f$y), mc)), f$p)
})

test_that("buildMrmashPriorMatrices exposes expand_covs' zeromat", {
    skip_if_not_installed("mr.mashr")
    # compute_canonical_covs() was already fully exposed (singletons,
    # hetgrid) but expand_covs()'s `zeromat` was hardcoded TRUE, so the
    # null component could not be dropped.
    seen <- NULL
    # The default must match the real expand_covs(): pecotmr no longer
    # passes zeromat unless the caller set it, so an omitted argument has to
    # fall through to the engine's own default here too.
    local_mocked_bindings(
        expand_covs = function(mats, grid, zeromat = TRUE) {
            seen <<- zeromat
            mats
        },
        .package = "mr.mashr"
    )
    set.seed(5)
    Bhat <- matrix(rnorm(12), 6, 2)
    Shat <- matrix(abs(rnorm(12, 0.5, 0.1)), 6, 2)
    invisible(buildMrmashPriorMatrices(
        Bhat,
        Shat,
        expandCovs = MrmashExpandCovsOptions(zeromat = FALSE)
    ))
    expect_false(seen)
    invisible(buildMrmashPriorMatrices(Bhat, Shat))
    expect_true(seen)
})

test_that("bayesAlphabetWeights forwards every MCMC control to gbayes", {
    # `nthin` was declared and documented (default 5) but never put in
    # callArgs, so every fit silently ran at gbayes's own nthin = 1 while
    # its siblings nit and nburn were forwarded.
    seen <- NULL
    local_mocked_bindings(
        gbayes = function(...) {
            seen <<- list(...)
            list(bm = rep(0, 5))
        },
        .package = "qgg"
    )
    f <- .rrwXy(p = 5)
    invisible(bayesAlphabetWeights(
        f$X,
        f$y,
        method = "bayesN",
        nit = 300,
        nburn = 30,
        nthin = 7
    ))
    expect_equal(seen$nit, 300)
    expect_equal(seen$nburn, 30)
    expect_equal(seen$nthin, 7)
})

test_that("bayesAlphabetWeights validates matching row counts before fitting", {
    skip_if_not_installed("qgg")
    f <- .rrwXy()
    expect_error(
        bayesAlphabetWeights(f$X, f$y[-1], method = "bayesN"),
        "y.*Must have length 50"
    )
    expect_error(
        bayesAlphabetWeights(
            f$X,
            f$y,
            method = "bayesN",
            Z = matrix(1, f$n - 1, 1)
        ),
        "Z.*Must have exactly 50 rows"
    )
})

test_that("bayesBWeights / bLassoWeights (BGLR) return length-p weights", {
    skip_if_not_installed("BGLR")
    f <- .rrwXy()
    expect_length(
        bayesBWeights(f$X, f$y, nIter = 200, burnIn = 20, thin = 1),
        f$p
    )
    expect_length(
        bLassoWeights(f$X, f$y, nIter = 200, burnIn = 20, thin = 1),
        f$p
    )
})

test_that("dprVbWeights returns length-p weights and retains the fit", {
    skip_if_not_installed("RcppDPR")
    f <- .rrwXy()
    w <- dprVbWeights(f$X, f$y, fitRetention = "slim")
    expect_length(w, f$p)
    expect_false(is.null(attr(w, "fit")))
})

test_that("dprGibbsWeights returns length-p weights", {
    skip_if_not_installed("RcppDPR")
    f <- .rrwXy()
    invisible(capture.output(w <- dprGibbsWeights(f$X, f$y, sStep = 200)))
    expect_length(w, f$p)
})

test_that("dprAdaptiveGibbsWeights returns length-p weights", {
    skip_if_not_installed("RcppDPR")
    f <- .rrwXy()
    invisible(capture.output(
        w <- dprAdaptiveGibbsWeights(
            f$X,
            f$y,
            methodArgs = DprOptions(s_step = 100)
        )
    ))
    expect_length(w, f$p)
})

test_that("mrmashWeights fits from (X, Y) and returns p x K weights", {
    skip_if_not_installed("mr.mashr")
    skip_if_not_installed("glmnet")
    set.seed(3)
    m <- .rrwMulti(n = 60, p = 6, K = 3)
    w <- suppressMessages(mrmashWeights(
        X = m$X,
        Y = m$Y,
        methodArgs = MrmashOptions(canonicalPriorMatrices = TRUE)
    ))
    expect_equal(dim(w), c(m$p, m$K))
    expect_true(all(is.finite(w)))
})


# ----------------------------- RSS solvers (C++) ----------------------------

test_that("lassosumRss returns a p x nlambda beta matrix", {
    f <- .rrwStatLd()
    out <- lassosumRss(f$stat$b, f$LD, f$n)
    expect_equal(nrow(out$beta), f$p)
    expect_equal(ncol(out$beta), length(out$lambda))
    expect_length(out$conv, length(out$lambda))
})

test_that("penalizedRss traces a solution path for MCP / SCAD / L0", {
    f <- .rrwStatLd()
    for (pen in c("MCP", "SCAD")) {
        out <- penalizedRss(f$stat$b, f$LD, f$n, penalty = pen)
        expect_equal(nrow(out$beta), f$p)
    }
    outL0 <- penalizedRss(
        f$stat$b,
        f$LD,
        f$n,
        penalty = "L0",
        lambda0 = 0.01,
        lambda = c(0)
    )
    expect_equal(nrow(outL0$beta), f$p)
})

test_that("prsCs returns posterior betaEst of length p", {
    f <- .rrwStatLd()
    out <- prsCs(f$stat$b, f$LD, f$n, nIter = 100, nBurnin = 20, thin = 1)
    expect_length(out$betaEst, f$p)
    expect_true(all(is.finite(out$betaEst)))
})

test_that("sdpr returns betaEst of length p", {
    f <- .rrwStatLd()
    out <- sdpr(
        f$stat$b,
        f$LD,
        f$n,
        iter = 100,
        burn = 20,
        thin = 1,
        verbose = FALSE
    )
    expect_length(out$betaEst, f$p)
})

test_that("RSS solvers validate their R-matrix / sample-size / length inputs", {
    f <- .rrwStatLd()
    # A list (the old block-list contract) is now rejected: R must be a matrix.
    expect_error(prsCs(f$stat$b, list(blk1 = f$LD), f$n), "as a matrix")
    expect_error(lassosumRss(f$stat$b, list(blk1 = f$LD), f$n), "as a matrix")
    expect_error(penalizedRss(f$stat$b, list(blk1 = f$LD), f$n), "as a matrix")
    expect_error(prsCs(f$stat$b, f$LD, -1), "sample size")
    expect_error(prsCs(f$stat$b[-1], f$LD, f$n), "bhat.*Must have length 6")
    expect_error(sdpr(f$stat$b[-1], f$LD, f$n), "bhat.*Must have length 6")
    expect_error(sdpr(f$stat$b, f$LD, f$n, M = 2), "at least 4")
})

# --------------------------- RSS weight wrappers ----------------------------

test_that("lassosumRssWeights returns length-p weights and records the selection", {
    f <- .rrwStatLd()
    w <- lassosumRssWeights(f$stat, f$ldReference)
    expect_length(w, f$p)
    expect_equal(unname(attr(w, "lassosum_selection")["mode"]), "ldQuadratic")
    expect_length(
        lassosumRssWeights(f$stat, f$ldReference, selection = "minFbeta"),
        f$p
    )
})

test_that("scadRssWeights / mcpRssWeights / l0learnRssWeights return length-p weights", {
    f <- .rrwStatLd()
    expect_length(scadRssWeights(f$stat, f$ldReference), f$p)
    expect_length(mcpRssWeights(f$stat, f$ldReference), f$p)
    expect_length(l0learnRssWeights(f$stat, f$ldReference), f$p)
})

test_that("mrashRssWeights attaches the mr.ash.rss fit on request", {
    f <- .rrwStatLd()
    none <- mrashRssWeights(f$stat, f$ldReference)
    expect_null(attr(none, "fit"))
    slim <- mrashRssWeights(f$stat, f$ldReference, fitRetention = "slim")
    expect_equal(as.numeric(slim), as.numeric(none))
    expect_equal(as.numeric(attr(slim, "fit")$mu1), as.numeric(none))
})

test_that("the raw RSS solvers are internal behind their weight adapters", {
    # The adapters are the public route; the Options constructors still read
    # the now-internal solvers' formals, so they keep rejecting a stray name.
    exported <- getNamespaceExports("pecotmr")
    solvers <- c("prsCs", "sdpr", "lassosumRss", "penalizedRss")
    adapters <- c(
        "prsCsWeights",
        "sdprWeights",
        "lassosumRssWeights",
        "scadRssWeights",
        "mcpRssWeights",
        "l0learnRssWeights"
    )
    expect_false(any(is_in(solvers, exported)))
    expect_true(all(is_in(adapters, exported)))
    expect_error(PrsCsOptions(notAnOption = 1), "unknown argument")
    expect_error(SdprOptions(notAnOption = 1), "unknown argument")
    expect_error(LassosumOptions(notAnOption = 1), "unknown argument")
    expect_error(PenalizedRssOptions(notAnOption = 1), "unknown argument")
})

test_that("prsCsWeights and sdprWeights attach the sampler's estimates", {
    # What the samplers estimate besides the weights reaches the caller as
    # the "fit" attribute; asking for it leaves the weights unchanged.
    f <- .rrwStatLd()
    prsOpts <- PrsCsOptions(nIter = 100, nBurnin = 20, thin = 1, seed = 42L)
    prsNone <- prsCsWeights(f$stat, f$ldReference, methodArgs = prsOpts)
    expect_null(attr(prsNone, "fit"))
    prsSlim <- prsCsWeights(
        f$stat,
        f$ldReference,
        fitRetention = "slim",
        methodArgs = prsOpts
    )
    prsFit <- attr(prsSlim, "fit")
    expect_named(prsFit, c("betaEst", "psiEst", "sigmaEst", "phiEst"))
    expect_equal(as.numeric(prsSlim), as.numeric(prsNone))
    expect_equal(prsFit$betaEst, as.numeric(prsNone))
    expect_length(prsFit$psiEst, f$p)
    sdprOpts <- SdprOptions(
        iter = 100,
        burn = 20,
        thin = 1,
        verbose = FALSE,
        seed = 42L
    )
    sdprNone <- sdprWeights(f$stat, f$ldReference, methodArgs = sdprOpts)
    expect_null(attr(sdprNone, "fit"))
    sdprFull <- sdprWeights(
        f$stat,
        f$ldReference,
        fitRetention = "full",
        methodArgs = sdprOpts
    )
    expect_named(attr(sdprFull, "fit"), c("betaEst", "h2"))
    expect_length(attr(sdprFull, "fit")$h2, 1L)
    expect_equal(as.numeric(sdprFull), as.numeric(sdprNone))
    expect_error(
        prsCsWeights(f$stat, f$ldReference, fitRetention = "all"),
        "fitRetention"
    )
})

test_that("lassosumRssWeights keeps its search: path at slim, beta at full", {
    f <- .rrwStatLd()
    s <- c(0.5, 1)
    none <- lassosumRssWeights(f$stat, f$ldReference, s = s)
    expect_null(attr(none, "fit"))
    slim <- lassosumRssWeights(
        f$stat,
        f$ldReference,
        s = s,
        fitRetention = "slim"
    )
    expect_named(attr(slim, "fit"), "path")
    path <- attr(slim, "fit")$path
    expect_named(
        path,
        c("s", "lambda", "fbeta", "loss", "conv", "nparams", "selected")
    )
    # One row per (s, lambda): lassosumRss() traces 20 lambdas by default.
    expect_equal(nrow(path), length(s) * 20L)
    expect_equal(sum(path$selected), 1L)
    expect_equal(as.numeric(slim), as.numeric(none))
    # The selected row is the one the selection attribute names.
    sel <- attr(slim, "lassosum_selection")
    expect_equal(path$s[path$selected], as.numeric(sel[["s"]]))
    expect_equal(path$lambda[path$selected], as.numeric(sel[["lambda"]]))
    full <- lassosumRssWeights(
        f$stat,
        f$ldReference,
        s = s,
        fitRetention = "full"
    )
    beta <- attr(full, "fit")$beta
    expect_equal(dim(beta), c(f$p, nrow(path)))
    expect_equal(unname(beta[, path$selected]), as.numeric(full))
    expect_equal(colSums(beta != 0), path$nparams)
})

test_that("scad, mcp and l0learn RSS weights pass fitRetention through", {
    f <- .rrwStatLd()
    cols <- c("s", "lambda", "fbeta", "loss", "conv", "nparams", "selected")
    scad <- scadRssWeights(f$stat, f$ldReference, s = 1, fitRetention = "slim")
    mcp <- mcpRssWeights(f$stat, f$ldReference, s = 1, fitRetention = "slim")
    expect_named(attr(scad, "fit")$path, cols)
    expect_named(attr(mcp, "fit")$path, cols)
    l0 <- l0learnRssWeights(
        f$stat,
        f$ldReference,
        s = 1,
        lambda0 = c(0.01, 0.1),
        fitRetention = "full"
    )
    path <- attr(l0, "fit")$path
    expect_named(path, c("s", "lambda0", cols[-1]))
    # Plain L0 traces the single lambda = 0, so one row per lambda0.
    expect_equal(path$lambda0, c(0.01, 0.1))
    expect_equal(ncol(attr(l0, "fit")$beta), 2L)
})

test_that("the RSS solvers run on a factor of the LD as on the LD itself", {
    # The factor route keeps F beta rather than R_s beta, so R_s is never
    # formed; the coordinate updates, the L0 swaps and the path are the same.
    # The dense route applies the shrinkage itself too, as fitting the
    # pre-shrunk matrix would.
    set.seed(8)
    X <- matrix(rbinom(80 * 120, 2, 0.3), 80, 120)
    F <- pecotmr:::.ldFactorFromDosage(X)
    R <- crossprod(F) / (nrow(F) - 1)
    beta <- replace(numeric(120), c(5, 40, 90), c(0.3, -0.2, 0.25))
    z <- as.vector(R %*% beta) + rnorm(120, sd = 0.05)
    lambda <- exp(seq(log(0.1), log(0.001), length.out = 8))
    for (s in c(0.2, 0.9)) {
        Rs <- (1 - s) * R + s * diag(120)
        dense <- pecotmr:::lassosumRssCpp(z, list(R), s, lambda, 1e-6, 1000L)
        expect_equal(
            dense,
            pecotmr:::lassosumRssCpp(z, list(Rs), 0, lambda, 1e-6, 1000L),
            tolerance = 1e-10
        )
        expect_equal(
            pecotmr:::lassosumRssFactorCpp(z, F, s, lambda, 1e-6, 1000L),
            dense,
            tolerance = 1e-8
        )
        for (pen in c("lasso", "MCP", "SCAD", "L0")) {
            args <- list(
                lambdaR = lambda,
                penaltyStr = pen,
                gamma = 3.7,
                alpha = 1,
                lambda0 = 0.001,
                lambda2 = 0,
                thr = 1e-6,
                maxiter = 1000L,
                maxSwaps = 10L
            )
            dense <- exec(pecotmr:::penalizedRssCpp, z, list(R), s, !!!args)
            expect_equal(
                dense,
                exec(pecotmr:::penalizedRssCpp, z, list(Rs), 0, !!!args),
                tolerance = 1e-10,
                label = pen
            )
            expect_equal(
                exec(pecotmr:::penalizedRssFactorCpp, z, F, s, !!!args),
                dense,
                tolerance = 1e-8,
                label = pen
            )
        }
    }
})

test_that("SDPR on a factor of the LD draws as it does on the matrix", {
    # The factor route holds SDPR's matrices in low-rank form, B = W W', and
    # runs SDPR's own sampler on them, so a seeded chain matches the dense one
    # to rounding. optLlk = 2 reaches B through the Woodbury identity while
    # its low-rank part is narrower than p (one sample size, three sizes, or
    # arrays), and otherwise forms it as the dense route does (three sizes
    # with arrays).
    set.seed(7)
    X <- matrix(rbinom(40 * 120, 2, 0.3), 40, 120)
    F <- pecotmr:::.ldFactorFromDosage(X)
    R <- crossprod(F) / (nrow(F) - 1)
    p <- ncol(F)
    effects <- replace(numeric(p), c(5, 40, 90), c(0.3, -0.25, 0.2))
    bhat <- as.vector(R %*% effects) + rnorm(p, sd = 0.02)
    sizes <- sample(c(400, 500, 600), p, replace = TRUE)
    arrays <- sample(0:2, p, replace = TRUE)
    base <- list(
        bhatR = bhat,
        n = 500L,
        perVariantSampleSize = NULL,
        array = NULL,
        a = 0.1,
        c = 1,
        M = 50L,
        a0k = 0.5,
        b0k = 0.5,
        iter = 200L,
        burn = 50L,
        thin = 5L,
        numThreads = 1L,
        optLlk = 1L,
        verbose = FALSE,
        seed = 11L
    )
    cases <- list(
        llk1 = list(),
        llk2 = list(optLlk = 2L),
        sizes = list(optLlk = 2L, perVariantSampleSize = sizes),
        arrays = list(optLlk = 2L, array = arrays),
        sizesArrays = list(
            optLlk = 2L,
            perVariantSampleSize = sizes,
            array = arrays
        )
    )
    for (nm in names(cases)) {
        args <- list_modify(base, !!!cases[[nm]])
        expect_equal(
            exec(pecotmr:::sdprFactorCpp, factorR = F, !!!args),
            exec(pecotmr:::sdprCpp, LD = list(R), !!!args),
            tolerance = 1e-8,
            label = nm
        )
    }
    expect_error(
        exec(pecotmr:::sdprFactorCpp, factorR = F[, -1], !!!base),
        "The LD factor has 119 columns but there are 120 effects"
    )
})

test_that("prsCsWeights and sdprWeights take (stat, ldReference)", {
    f <- .rrwStatLd()
    expect_length(
        prsCsWeights(
            f$stat,
            f$ldReference,
            methodArgs = PrsCsOptions(nIter = 100, nBurnin = 20, thin = 1)
        ),
        f$p
    )
    expect_length(
        sdprWeights(
            f$stat,
            f$ldReference,
            methodArgs = SdprOptions(
                iter = 100,
                burn = 20,
                thin = 1,
                verbose = FALSE
            )
        ),
        f$p
    )
})

test_that("prsCsWeights and sdprWeights hand the engine z / sqrt(n)", {
    # Both upstream tools fit the standardized marginal effects, so a z-only
    # stat (the QtlSumStats shape) and b with its standard error both reach
    # the engine on that scale.
    f <- .rrwStatLd()
    want <- as.numeric(f$stat$z) / sqrt(f$n)
    local_mocked_bindings(
        .prsCsRun = function(bhat, ld, n, ...) list(betaEst = bhat),
        .sdprRun = function(bhat, ld, n, ...) list(betaEst = bhat)
    )
    ids <- f$stat$variantNames
    zOnly <- list(z = f$stat$z, n = f$n, variantNames = ids)
    withSe <- list(
        b = 3 * f$stat$z,
        seb = rep(3, f$p),
        n = rep(f$n, f$p),
        variantNames = ids
    )
    for (stat in list(zOnly, withSe)) {
        expect_equal(prsCsWeights(stat, f$ldReference), want)
        expect_equal(sdprWeights(stat, f$ldReference), want)
    }
})

# Statistics over the first p variants of the test genotype panel, with the
# same LD as an LdMatrix: a weight method must read the panel's LD (through
# its standardized genotypes) as it reads the matrix.
.rrwPanelStat <- function(p = 40) {
    panel <- readGenotypes(test_path("test_data", "test_variants.bed"))
    ids <- names(SummarizedExperiment::rowRanges(panel))[seq_len(p)]
    R <- pecotmr:::.ldFromRef(panel, ids, "test")
    set.seed(3)
    effects <- replace(numeric(p), c(3, 20), c(4, -3))
    z <- as.vector(R %*% effects) + rnorm(p)
    list(
        panel = panel,
        ldMatrix = LdMatrix(R, nSamples = ncol(panel)),
        stat = list(z = z, n = rep(500, p), variantNames = ids)
    )
}

test_that("RSS weights read a genotype panel's LD as they read the matrix", {
    # The panel route gives the solvers the panel's factor and never forms
    # the LD matrix; it must land on the matrix route's weights.
    d <- .rrwPanelStat()
    fns <- list(
        lassosum = lassosumRssWeights,
        scad = scadRssWeights,
        mcp = mcpRssWeights,
        l0learn = l0learnRssWeights
    )
    for (nm in names(fns)) {
        expect_equal(
            fns[[nm]](d$stat, d$panel),
            fns[[nm]](d$stat, d$ldMatrix),
            tolerance = 1e-6,
            label = nm
        )
    }
    # mr.ash.rss takes only a matrix, formed from the panel.
    expect_equal(
        mrashRssWeights(d$stat, d$panel),
        mrashRssWeights(d$stat, d$ldMatrix),
        tolerance = 1e-8
    )
    # SDPR runs its own sampler on the low-rank form of its matrices, so a
    # seeded chain lands on the matrix route's.
    for (optLlk in 1:2) {
        opts <- SdprOptions(
            iter = 100,
            burn = 20,
            thin = 1,
            optLlk = optLlk,
            verbose = FALSE,
            seed = 1L
        )
        expect_equal(
            sdprWeights(d$stat, d$panel, methodArgs = opts),
            sdprWeights(d$stat, d$ldMatrix, methodArgs = opts),
            tolerance = 1e-6,
            label = str_c("sdpr optLlk ", optLlk)
        )
    }
})

test_that("an RSS weight method matches its LD reference by variant id", {
    f <- .rrwStatLd()
    unnamed <- f$stat
    unnamed$variantNames <- NULL
    expect_error(
        lassosumRssWeights(unnamed, f$ldReference),
        "must carry `variantNames`"
    )
    expect_error(
        mrashRssWeights(unnamed, f$ldReference),
        "must carry `variantNames`"
    )
    # The extractors need the ids only to check against a reference.
    expect_equal(pecotmr:::.rssVariantCount(unnamed, NULL, "t"), f$p)
    expect_equal(pecotmr:::.rssVariantCount(f$stat, f$ldReference, "t"), f$p)
    expect_true(is.na(pecotmr:::.rssVariantCount(list(), NULL, "t")))
})

test_that(".formalDefaults evaluates a function's own defaults", {
    expect_equal(
        pecotmr:::.formalDefaults(prsCs, c("a", "phi", "nIter")),
        list(a = 1, phi = NULL, nIter = 1000)
    )
})

test_that("mrashRssWeights defaults to mr.ash's own prior grid", {
    # On the standardized scale mr.ash.rss fits, its sufficient statistics
    # are those of standardized genotypes, whose grid mr.ash builds itself.
    set.seed(1)
    X <- scale(matrix(rnorm(200 * 10), 200))
    fit <- susieR::mr.ash(X, rnorm(200), verbose = FALSE, max.iter = 5)
    expect_equal(
        pecotmr:::.mrashRssDefaultGrid(200, rnorm(10), numeric(0), NULL),
        fit$data$sa2
    )
    f <- .rrwStatLd()
    expect_length(mrashRssWeights(f$stat, f$ldReference), f$p)
    # With b, seb and varY the unstandardized statistics set the scale.
    z <- f$stat$z
    shat <- f$stat$seb
    w <- 2 * (f$n - 1) / (z^2 + f$n - 2) / shat^2
    expect_equal(
        pecotmr:::.mrashRssDefaultGrid(f$n, z, shat, 2),
        (2^((0:24) / 25) - 1)^2 / median(w) * f$n
    )
})

test_that("mrashRssWeights returns posterior-mean weights of length p", {
    skip_if_not_installed("susieR")
    f <- .rrwStatLd()
    w <- mrashRssWeights(
        f$stat,
        f$ldReference,
        varY = 1,
        sigma2E = 1,
        s0 = c(0, 0.01, 0.1, 0.5, 1),
        w0 = rep(1 / 5, 5)
    )
    expect_length(w, f$p)
    expect_true(all(is.finite(w)))
})


test_that("mrmashRssWeights fits mr.mash.rss and returns p x K weights", {
    skip_if_not_installed("mr.mashr")
    m <- .rrwMulti(n = 60, p = 6, K = 3)
    w <- mrmashRssWeights(m$stat, m$ldReference)
    expect_equal(dim(w), c(m$p, m$K))
    expect_true(all(is.finite(w)))
})

test_that("mrmashRssWeights errors on single-context stat$z", {
    skip_if_not_installed("mr.mashr")
    f <- .rrwStatLd()
    oneCol <- list(
        z = matrix(f$stat$z, ncol = 1),
        n = f$n,
        variantNames = f$stat$variantNames
    )
    expect_error(mrmashRssWeights(oneCol, f$ldReference), ">= 2 columns")
})


# ------------------------------- pure helpers -------------------------------

test_that(".rssMarginalCor reads cor / z / b and validates length", {
    f <- .rrwStatLd()
    expect_length(
        pecotmr:::.rssMarginalCor(
            list(cor = f$stat$cor),
            n = f$n,
            p = f$p
        ),
        f$p
    )
    expect_equal(
        pecotmr:::.rssMarginalCor(list(z = f$stat$z), n = f$n, p = f$p),
        as.numeric(f$stat$z) / sqrt(f$n)
    )
    expect_equal(
        pecotmr:::.rssMarginalCor(list(b = f$stat$b), n = f$n, p = f$p),
        as.numeric(f$stat$b)
    )
    # With its standard error, b is a z-score on any scale.
    expect_equal(
        pecotmr:::.rssMarginalCor(
            list(b = 2 * f$stat$z, seb = rep(2, f$p)),
            n = f$n,
            p = f$p
        ),
        as.numeric(f$stat$z) / sqrt(f$n)
    )
    expect_error(
        pecotmr:::.rssMarginalCor(list(), n = f$n, p = f$p),
        "one of"
    )
    expect_error(
        pecotmr:::.rssMarginalCor(
            list(z = f$stat$z[-1]),
            n = f$n,
            p = f$p
        ),
        "must equal"
    )
})

test_that(".lassosumClampCor scales values with |cor| >= 1 below 1", {
    expect_equal(pecotmr:::.lassosumClampCor(c(0.1, 0.5)), c(0.1, 0.5))
    expect_lt(max(abs(pecotmr:::.lassosumClampCor(c(0.5, 1.5, -2)))), 1)
})

test_that(".lassosumFirstMax returns the first index of the maximum", {
    expect_equal(pecotmr:::.lassosumFirstMax(c(1, 3, 3, 2)), 2L)
    expect_equal(pecotmr:::.lassosumFirstMax(c(5, 1, 2)), 1L)
})

test_that(".lassosumSelectMinFbeta picks the minimum-fbeta candidate", {
    set.seed(7)
    cb <- matrix(rnorm(6 * 4), 6, 4)
    r <- pecotmr:::.lassosumSelectMinFbeta(
        cb,
        data.frame(fbeta = c(3, 1, 2, 4))
    )
    expect_equal(r$index, 2)
    expect_equal(r$mode, "minFbeta")
    expect_equal(r$beta, cb[, 2])
})

test_that(".lassosumSelectLdQuadratic scores candidates by c'b / sqrt(b'Rb)", {
    f <- .rrwStatLd()
    set.seed(8)
    cb <- matrix(rnorm(f$p * 4), f$p, 4)
    r <- pecotmr:::.lassosumSelectLdQuadratic(cb, f$stat$b, list(R = f$LD))
    expect_equal(r$mode, "ldQuadratic")
    expect_true(r$index %in% seq_len(4))
    expect_equal(r$beta, cb[, r$index])
})

test_that("computeCovDiag returns a diagonal condition covariance", {
    m <- .rrwMulti(n = 60, p = 6, K = 3)
    cv <- computeCovDiag(m$Y)
    expect_equal(dim(cv), c(m$K, m$K))
    expect_equal(cv[upper.tri(cv)], rep(0, sum(upper.tri(cv))))
    expect_equal(unname(diag(cv)), unname(apply(m$Y, 2, var)))
})

test_that("buildMrmashPriorMatrices builds an expanded S0 list and a prior grid", {
    skip_if_not_installed("mr.mashr")
    set.seed(9)
    res <- buildMrmashPriorMatrices(
        Bhat = matrix(rnorm(18), 6, 3),
        Shat = matrix(0.2, 6, 3),
        K = 3
    )
    expect_true(is.list(res$S0))
    expect_gt(length(res$S0), 1)
    expect_true(is.numeric(res$priorGrid))
    expect_true(all(map_lgl(
        res$S0,
        function(s) all(dim(s) == c(3, 3))
    )))
})

test_that("buildMrmashPriorMatrices errors without canonical or data-driven priors", {
    skip_if_not_installed("mr.mashr")
    expect_error(
        buildMrmashPriorMatrices(
            Bhat = matrix(rnorm(6), 3, 2),
            Shat = matrix(0.2, 3, 2),
            K = 2,
            canonicalPriorMatrices = FALSE
        ),
        "dataDrivenPriorMatrices"
    )
})


# ===========================================================================
# mr.mash weight tests (relocated from test_rrMrmashMvsusie.R)
# ===========================================================================

# ---- mrmashWeights ----
test_that("mrmashWeights errors when mr.mashr package is not available", {
    skip_if(
        requireNamespace("mr.mashr", quietly = TRUE),
        "mr.mashr is installed; skipping missing-package test"
    )

    expect_error(
        mrmashWeights(
            mrmashFit = NULL,
            X = matrix(1, 10, 5),
            Y = matrix(1, 10, 3)
        ),
        "mr\\.mash\\.alpha"
    )
})

test_that("mrmashWeights errors when X and Y are NULL and fit is NULL", {
    skip_if_not(
        requireNamespace("mr.mashr", quietly = TRUE),
        "mr.mashr not installed"
    )
    expect_error(
        mrmashWeights(mrmashFit = NULL, X = NULL, Y = NULL),
        "Both X and Y must be provided"
    )
})

test_that("mrmashWeights retaining attaches {dataDrivenPriorMatrices, w0, V}", {
    skip_if_not(
        requireNamespace("mr.mashr", quietly = TRUE),
        "mr.mashr not installed"
    )
    # These are exactly the parts fineMappingPipeline needs to rebuild the
    # mvSuSiE reweighted mixture prior (w0 -> rescaleCovW0, original $U) and the
    # residual variance (V); the heavy mu1 coefficient matrix is not retained.
    ddpm <- list(U = list(comp = diag(2)))
    fakeFit <- structure(
        list(w0 = c(null = 0.4, comp_grid1 = 0.6), V = diag(2) * 2),
        class = "mr.mash"
    )
    fakeCoef <- matrix(0.1, nrow = 5, ncol = 2)
    local_mocked_bindings(
        coef.mr.mash = function(object, ...) fakeCoef,
        .package = "mr.mashr"
    )
    w <- mrmashWeights(
        mrmashFit = fakeFit,
        dataDrivenPriorMatrices = ddpm,
        fitRetention = "slim"
    )
    fit <- attr(w, "fit")
    expect_true(is.list(fit))
    expect_identical(fit$dataDrivenPriorMatrices, ddpm)
    expect_identical(fit$w0, fakeFit$w0)
    expect_identical(fit$V, fakeFit$V)
    # The default ("none") leaves the weights free of the fit attribute.
    expect_null(attr(
        mrmashWeights(mrmashFit = fakeFit, dataDrivenPriorMatrices = ddpm),
        "fit"
    ))
})

test_that("MCMC / optimiser control arguments are guarded", {
    expect_error(
        bayesBWeights(NULL, NULL, nIter = 0),
        "nIter.*Must be >= 1"
    )
    expect_error(
        bayesBWeights(NULL, NULL, thin = 0),
        "thin.*Must be >= 1"
    )
    expect_error(
        bayesBWeights(NULL, NULL, probIn = 2),
        "probIn.*is not <= 1"
    )
    expect_error(
        bLassoWeights(NULL, NULL, burnIn = -1),
        "burnIn.*Must be >= 0"
    )
    expect_error(
        dprGibbsWeights(NULL, NULL, sStep = 0),
        "sStep.*Must be >= 1"
    )
    expect_error(
        dprAdaptiveGibbsWeights(NULL, NULL, fitRetention = "sometimes"),
        "must be one of"
    )
    expect_error(
        ncvregWeights(NULL, NULL, penalty = "SCAD", nfolds = 0),
        "nfolds.*Must be >= 1"
    )
    expect_error(
        computeCoefficientsGlasso(NULL, NULL, standardize = NA, numThreads = 1),
        "standardize.*May not be NA"
    )
    expect_error(
        computeCoefficientsGlasso(
            NULL,
            NULL,
            standardize = TRUE,
            numThreads = 1.5
        ),
        "numThreads.*integerish"
    )
})

test_that("RSS solver control arguments are guarded", {
    expect_error(
        lassosumRss(NULL, NULL, NULL, thr = -1),
        "thr.*is not >= 0"
    )
    expect_error(
        lassosumRss(NULL, NULL, NULL, maxiter = 0),
        "maxiter.*Must be >= 1"
    )
    expect_error(
        scadRssWeights(stat = "nope", ldReference = NULL),
        "stat.*Must be of type 'list'"
    )
    expect_error(
        scadRssWeights(stat = list(), ldReference = NULL, s = -1),
        "s.*is not >= 0"
    )
    expect_error(
        mcpRssWeights(stat = "nope", ldReference = NULL),
        "stat.*Must be of type 'list'"
    )
    expect_error(
        mrashRssWeights(stat = "nope", ldReference = NULL),
        "stat.*Must be of type 'list'"
    )
    expect_error(
        mrmashRssWeights(stat = "nope", ldReference = NULL),
        "stat.*Must be of type 'list'"
    )
})

test_that("MrmashPriorParam carries prior options and rejects a typo", {
    pa <- MrmashPriorParam(canonicalPriorMatrices = TRUE)
    expect_s4_class(pa, "MrmashPriorParam")
    expect_true(pa$canonicalPriorMatrices)
    # pecotmr's own builder, so no `...`: R rejects an unknown name. Note
    # `canonicalPrior` would be ACCEPTED, by R's partial matching.
    expect_error(
        MrmashPriorParam(canonicalPriorMatrix = TRUE),
        "unused argument"
    )
    expect_true(MrmashPriorParam(canonicalPrior = TRUE)$canonicalPriorMatrices)
})

test_that("MrmashPriorParam nests each engine function's own arguments", {
    skip_if_not_installed("mr.mashr")
    # hetgrid and singletons belong to compute_canonical_covs, zeromat to
    # expand_covs. Nesting them says where each goes AND gets each checked
    # against that function's live formals rather than a transcribed list.
    pa <- MrmashPriorParam(
        canonicalCovs = MrmashCanonicalCovsOptions(
            hetgrid = c(0, 0.5),
            singletons = FALSE
        ),
        expandCovs = MrmashExpandCovsOptions(zeromat = FALSE)
    )
    expect_equal(pa$canonicalCovs$hetgrid, c(0, 0.5))
    expect_false(pa$canonicalCovs$singletons)
    expect_false(pa$expandCovs$zeromat)
    # Each bundle rejects the other's name, which a single flat bag could not.
    expect_error(
        MrmashCanonicalCovsOptions(zeromat = FALSE),
        "unknown argument"
    )
    expect_error(
        MrmashExpandCovsOptions(singletons = FALSE),
        "unknown argument"
    )
    # Neither carries a default: pecotmr used to copy the engine's own.
    expect_length(MrmashCanonicalCovsOptions(), 0L)
    expect_length(MrmashExpandCovsOptions(), 0L)
})

test_that("mrmashWrapper forwards methodArgs under mr.mash's own names", {
    skip_if_not_installed("mr.mashr")
    ma <- MrmashOptions(max_iter = 10, tol = 0.5)
    expect_equal(ma$max_iter, 10)
    # The renamed pecotmr spellings are gone; mr.mash's names are the contract.
    expect_error(MrmashOptions(maxIter = 10), "unknown argument")
})

test_that("mrmashWrapper refuses to let methodArgs override derived values", {
    skip_if_not_installed("mr.mashr")
    skip_if_not_installed("glmnet")
    set.seed(1)
    X <- matrix(rnorm(60 * 4), 60, 4)
    Y <- matrix(rnorm(60 * 2), 60, 2)
    expect_error(
        mrmashWrapper(
            X,
            Y,
            dataDrivenPriorMatrices = list(U = list(diag(2))),
            methodArgs = MrmashOptions(S0 = list(diag(2)))
        ),
        "supplied by pecotmr"
    )
})

test_that("a flat prior name is routed into the prior group", {
    # The pipeline spells per-token kwargs flat, so mrmashWeights routes
    # prior-construction names into mrmashWrapper's `prior` record.
    split <- pecotmr:::.mrmashSplitPriorArgs(
        list(canonicalPriorMatrices = TRUE, max_iter = 5)
    )
    expect_equal(names(split$prior), "canonicalPriorMatrices")
    expect_equal(names(split$rest), "max_iter")
})
