# Tests migrated from test_h2ClassesSumstats.R

# === Tests migrated from test_h2ClassesSumstats.R (H2Estimate) ===

test_that("H2Estimate constructs with all slots", {
    obj <- new(
        "H2Estimate",
        estimate = 0.3,
        estimateSe = 0.05,
        intercept = 1.01,
        interceptSe = 0.02,
        localBlocks = NULL,
        enrichment = NULL,
        annotationJackknifeCoefs = NULL,
        scoreStats = NULL,
        method = "lder",
        nSnps = 10000L,
        traitName = "height"
    )
    expect_s4_class(obj, "H2Estimate")
    expect_equal(heritabilityEstimate(obj), 0.3)
    expect_equal(methodNames(obj), "lder")
})

# show() smoke test, moved here from test_showMethods.R so the test
# tree mirrors R/.
test_that("show(H2Estimate) does not error", {
    h2 <- new(
        "H2Estimate",
        estimate = 0.3,
        estimateSe = 0.05,
        intercept = 1.0,
        interceptSe = 0.01,
        localBlocks = NULL,
        enrichment = NULL,
        annotationJackknifeCoefs = NULL,
        scoreStats = NULL,
        method = "lder",
        nSnps = 1000L,
        traitName = "test"
    )
    expect_output(show(h2), "H2Estimate")
})

test_that("H2Estimate types its optional slots", {
    # Each optional slot is a type-or-NULL union, so a tibble fills the table
    # slots while a value of the wrong type is refused when the object is
    # built.
    h2 <- function(...) {
        new(
            "H2Estimate",
            estimate = 0.3,
            estimateSe = 0.05,
            intercept = 1.0,
            interceptSe = 0.01,
            method = "lder",
            nSnps = 1000L,
            traitName = "test",
            ...
        )
    }
    local <- tibble::tibble(blockId = 1L, h2Local = 0.1, h2LocalSe = 0.01)
    expect_identical(localH2Blocks(h2(localBlocks = local)), local)
    expect_null(localH2Blocks(h2()))
    expect_error(
        h2(localBlocks = list(blockId = 1L)),
        'invalid object for slot "localBlocks"'
    )
    expect_error(
        h2(enrichment = "none"),
        'invalid object for slot "enrichment"'
    )
    expect_error(
        h2(annotationJackknifeCoefs = c(0.1, 0.2)),
        'invalid object for slot "annotationJackknifeCoefs"'
    )
    expect_error(
        h2(scoreStats = c(z = 1.5)),
        'invalid object for slot "scoreStats"'
    )
})
