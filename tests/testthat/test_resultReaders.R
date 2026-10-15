# ---------------------------------------------------------------------------
# Published result readers: readFineMappingResult() / readTwasWeights()
#
# The reader and the writer are each other's only real test. Each direction
# had only ever been checked against a stored fixture before, which is how the
# two came to disagree about a reserved key while both "passed" -- so most of
# what is here is a round trip, and the comparisons are against the ORIGINAL
# object rather than against a file.
# ---------------------------------------------------------------------------

# @noRd
.rdr_tmp <- function(name) file.path(tempdir(), name)

# The ways a round trip legitimately differs from its input, so a test can say
# "everything else is identical" and mean it:
#   dataSource / standardized  ADDED on the way back (one is derived by the
#     writer, the other supplied to the reader)
#   cs_<n>  the published form is an integer index, and "susie_0" -- the
#     not-captured sentinel -- publishes as NA, which is the canonical form
#     every reader of the column already treats as "in no credible set"
#   z_supplied  resolved from NA to FALSE: a file publishing ES/SE and no EZ
#     says the z was NOT supplied by the source, which is what the fixture's
#     NA meant before the flag existed
# @noRd
.RDR_NORMALIZED <- c(
    "dataSource",
    "standardized",
    "cs_95",
    "cs_70",
    "cs_50",
    "z_supplied"
)

# Which columns of two per-variant tables disagree, aligned on variant id.
#
# Aligned rather than compared in place because a published file is sorted
# genomically and the stored table is not, and restricted to the variants
# `b` carries because an export is topLoci-FILTERED by default -- the
# question is whether what was published came back intact, not whether a
# filtered file holds every variant of the fit.
# @noRd
.rdr_diff <- function(a, b, tolerance = 1e-6) {
    a <- a[is_in(a$variant_id, b$variant_id), ]
    b <- b[match(a$variant_id, b$variant_id), ]
    shared <- intersect(names(a), names(b))
    bad <- character(0)
    for (nm in shared) {
        va <- a[[nm]]
        vb <- b[[nm]]
        ok <- if (is.numeric(va) && is.numeric(vb)) {
            isTRUE(all.equal(va, vb, tolerance = tolerance))
        } else {
            identical(as.character(va), as.character(vb))
        }
        if (!ok) {
            bad <- c(bad, nm)
        }
    }
    setdiff(bad, .RDR_NORMALIZED)
}

# A minimal published table written by hand, as a user converting another
# tool's output would produce it.
# @noRd
.rdr_writeTable <- function(df, name) {
    path <- .rdr_tmp(name)
    utils::write.table(
        df,
        path,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        na = "NA"
    )
    path
}

# ---- round trips ----------------------------------------------------------

test_that("a fine-mapping result round-trips through a published BED", {
    data(qtlFineMappingExample, envir = environment())
    path <- .rdr_tmp("rt_fm.bed")
    exportToBed(qtlFineMappingExample, path, columns = "all")
    back <- readFineMappingResult(path, weightStandardized = TRUE)
    expect_s4_class(back, "QtlFineMappingResult")
    expect_equal(nrow(back), nrow(qtlFineMappingExample))
    expect_equal(
        .rdr_diff(
            as.data.frame(topLoci(qtlFineMappingExample)),
            as.data.frame(topLoci(back))
        ),
        character(0)
    )
    # The marginal view is the other half of the same table and has its own
    # published fields (ES/SE/EZ/LP), so it is checked separately.
    expect_equal(
        .rdr_diff(
            as.data.frame(marginalEffects(qtlFineMappingExample)),
            as.data.frame(marginalEffects(back))
        ),
        character(0)
    )
    # The analysis WINDOW, without which genomicRegion() silently degrades to
    # the span of the variants that survived QC.
    expect_equal(
        genomicRegion(back),
        genomicRegion(qtlFineMappingExample)
    )
})

test_that("a fine-mapping result round-trips through a published VCF", {
    data(qtlFineMappingExample, envir = environment())
    path <- .rdr_tmp("rt_fm.vcf")
    exportToVcf(qtlFineMappingExample, path, context = "context_1")
    back <- readFineMappingResult(path, weightStandardized = TRUE)
    expect_s4_class(back, "QtlFineMappingResult")
    expect_equal(
        .rdr_diff(
            as.data.frame(topLoci(qtlFineMappingExample)),
            as.data.frame(topLoci(back)),
            tolerance = 1e-5
        ),
        character(0)
    )
    expect_equal(genomicRegion(back), genomicRegion(qtlFineMappingExample))
})

test_that("the two formats agree with each other", {
    # The point of one shared vocabulary: whichever format a result was
    # published in, it reads back the same. A `region_start` that is 0-based
    # in the table and 1-based in the VCF header was caught by exactly this.
    data(qtlFineMappingExample, envir = environment())
    bed <- .rdr_tmp("agree.bed")
    vcf <- .rdr_tmp("agree.vcf")
    exportToBed(qtlFineMappingExample, bed, columns = "all")
    exportToVcf(qtlFineMappingExample, vcf, context = "context_1")
    fromBed <- readFineMappingResult(bed, weightStandardized = TRUE)
    fromVcf <- readFineMappingResult(vcf, weightStandardized = TRUE)
    expect_equal(genomicRegion(fromBed), genomicRegion(fromVcf))
    expect_equal(
        .rdr_diff(
            as.data.frame(topLoci(fromBed)),
            as.data.frame(topLoci(fromVcf)),
            tolerance = 1e-5
        ),
        character(0)
    )
})

test_that("a published VCF restores the per-effect fit the header carried", {
    # `V` and `sets$cs_index` are per EFFECT, so no per-variant column can
    # hold them -- the `##SAMPLE` block is why they survive at all, and
    # colocalization's credible-set filter reads them off the fit.
    data(qtlFineMappingExample, envir = environment())
    path <- .rdr_tmp("fit.vcf")
    exportToVcf(qtlFineMappingExample, path, context = "context_1")
    fit <- susieFit(readFineMappingResult(path, weightStandardized = TRUE))
    original <- susieFit(qtlFineMappingExample)
    expect_equal(as.numeric(fit$V), as.numeric(original$V), tolerance = 1e-8)
    expect_equal(
        as.integer(fit$sets$cs_index),
        as.integer(original$sets$cs_index)
    )
    # A BED has no header to carry them, so there the fit is honestly absent
    # rather than half-invented.
    bed <- .rdr_tmp("fit.bed")
    exportToBed(qtlFineMappingExample, bed)
    expect_null(susieFit(readFineMappingResult(bed, weightStandardized = TRUE)))
})

test_that("TWAS weights round-trip through a published BED", {
    data(twasWeightsExample, envir = environment())
    path <- .rdr_tmp("rt_w.bed")
    expect_message(
        exportToBed(twasWeightsExample, path),
        "fine-mapping-derived"
    )
    back <- readTwasWeights(path)
    # Fine-mapping-derived rows are skipped on export, so fewer rows come
    # back than went out -- that is the writer's rule, not a reader loss.
    expect_lt(nrow(back), nrow(twasWeightsExample))
    md <- S4Vectors::mcols(back)
    for (i in seq_len(nrow(back))) {
        key <- list(
            studyName = as.character(md$study[[i]]),
            context = as.character(md$context[[i]]),
            trait = as.character(md$trait[[i]]),
            method = as.character(md$method[[i]])
        )
        was <- exec(resolveWeights, twasWeightsExample, !!!key)
        now <- exec(resolveWeights, back, !!!key)
        expect_setequal(now$variantIds, was$variantIds)
        expect_equal(
            as.numeric(now$weights)[order(now$variantIds)],
            as.numeric(was$weights)[order(was$variantIds)],
            tolerance = 1e-10
        )
    }
    expect_equal(
        genomicRegion(back),
        genomicRegion(twasWeightsExample)[
            match(
                as.character(md$method),
                as.character(S4Vectors::mcols(twasWeightsExample)$method)
            )
        ]
    )
})

# ---- multi-row tables -----------------------------------------------------

test_that("a merged table reads back as a multi-row collection", {
    data(qtlFineMappingExample, qtlFineMappingLbfExample, envir = environment())
    combined <- combineFineMappingResults(
        qtlFineMappingExample,
        qtlFineMappingLbfExample
    )
    path <- .rdr_tmp("merged.bed")
    exportToBed(combined, path)
    back <- readFineMappingResult(path, weightStandardized = TRUE)
    expect_equal(nrow(back), nrow(combined))
    md <- S4Vectors::mcols(back)
    expect_setequal(
        as.character(md$trait),
        as.character(S4Vectors::mcols(combined)$trait)
    )
})

test_that("a per-result value that disagrees within one tuple is an error", {
    # Two results concatenated under one identity. Picking either window
    # would silently publish one fit's region as the other's.
    df <- data.frame(
        `#chr` = c("chr1", "chr1"),
        start = c(99L, 199L),
        end = c(100L, 200L),
        ID = c("chr1:100:A:G", "chr1:200:C:T"),
        PIP = c(0.9, 0.1),
        study = c("s", "s"),
        context = c("c", "c"),
        trait = c("t", "t"),
        method = c("susie", "susie"),
        region_start = c(0L, 500L),
        region_end = c(1000L, 1500L),
        check.names = FALSE
    )
    path <- .rdr_writeTable(df, "clash.bed")
    expect_error(
        readFineMappingResult(path),
        "`region_start` takes 2 different values"
    )
})

# ---- supplying what the file does not record -------------------------------

test_that("identity arguments fill in a file that did not record them", {
    df <- data.frame(
        `#chr` = "chr1",
        start = 99L,
        end = 100L,
        ID = "chr1:100:A:G",
        PIP = 0.9,
        check.names = FALSE
    )
    path <- .rdr_writeTable(df, "bare.bed")
    back <- readFineMappingResult(
        path,
        studyName = "s",
        context = "c",
        trait = "t",
        method = "susie"
    )
    expect_s4_class(back, "QtlFineMappingResult")
    md <- S4Vectors::mcols(back)
    expect_equal(as.character(md$study), "s")
    expect_equal(as.character(md$method), "susie")
    expect_equal(as.numeric(topLoci(back)$pip), 0.9)
})

test_that("an identity argument that disagrees with the file is an error", {
    df <- data.frame(
        `#chr` = "chr1",
        start = 99L,
        end = 100L,
        ID = "chr1:100:A:G",
        PIP = 0.9,
        study = "recorded",
        context = "c",
        trait = "t",
        method = "susie",
        check.names = FALSE
    )
    path <- .rdr_writeTable(df, "disagree.bed")
    expect_error(
        readFineMappingResult(path, studyName = "supplied"),
        "was given as 'supplied' but the file records 'recorded'"
    )
    # Dropping the argument uses the file's value.
    expect_equal(
        as.character(S4Vectors::mcols(readFineMappingResult(path))$study),
        "recorded"
    )
})

test_that("a result keyed by axes the file lacks names them", {
    df <- data.frame(
        `#chr` = "chr1",
        start = 99L,
        end = 100L,
        ID = "chr1:100:A:G",
        PIP = 0.9,
        context = "c",
        check.names = FALSE
    )
    path <- .rdr_writeTable(df, "noaxes.bed")
    expect_error(
        readFineMappingResult(path),
        "does not record study, trait, method.*`studyName`, `trait`, `method`"
    )
})

# ---- permissiveness -------------------------------------------------------

test_that("a table carrying only PIP is enough for a fine-mapping result", {
    # TWAS needs only the weight, MR beta/se/PIP, colocalization the LBF. A
    # file converted from another tool has some and not others, so the reader
    # requires only what the class cannot exist without.
    df <- data.frame(
        `#chr` = c("chr1", "chr1"),
        start = c(99L, 199L),
        end = c(100L, 200L),
        ID = c("chr1:100:A:G", "chr1:200:C:T"),
        PIP = c(0.9, 0.05),
        study = "s",
        context = "c",
        trait = "t",
        method = "susie",
        check.names = FALSE
    )
    path <- .rdr_writeTable(df, "pip_only.bed")
    back <- readFineMappingResult(path)
    tl <- topLoci(back, signalCutoff = 0)
    expect_equal(length(tl), 2L)
    expect_equal(as.numeric(S4Vectors::mcols(tl)$pip), c(0.9, 0.05))
    # Absent columns are absent, not NA-filled: nothing invents a beta.
    expect_true(all(is.na(marginalEffects(back)$beta)))
})

test_that("a table carrying only a weight is enough for TwasWeights", {
    df <- data.frame(
        `#chr` = c("chr1", "chr1"),
        start = c(99L, 199L),
        end = c(100L, 200L),
        ID = c("chr1:100:A:G", "chr1:200:C:T"),
        WGT = c(0.3, -0.2),
        study = "s",
        context = "c",
        trait = "t",
        method = "lasso",
        weight_standardized = TRUE,
        check.names = FALSE
    )
    path <- .rdr_writeTable(df, "wgt_only.bed")
    back <- readTwasWeights(path)
    w <- resolveWeights(back)
    expect_equal(as.numeric(w$weights), c(0.3, -0.2))
    expect_true(weightStandardized(back))
})

test_that("ids are rebuilt from coordinates when ID does not encode them", {
    df <- data.frame(
        `#chr` = "chr1",
        start = 99L,
        end = 100L,
        ID = "rs123",
        REF = "A",
        ALT = "G",
        PIP = 0.9,
        study = "s",
        context = "c",
        trait = "t",
        method = "susie",
        check.names = FALSE
    )
    path <- .rdr_writeTable(df, "rsid.bed")
    back <- readFineMappingResult(path)
    expect_equal(variantIds(back), "chr1:100:A:G")
})

test_that("a table with no way to identify its variants says so", {
    df <- data.frame(
        ID = "rs123",
        PIP = 0.9,
        check.names = FALSE
    )
    path <- .rdr_writeTable(df, "noid.bed")
    expect_error(
        readFineMappingResult(path, studyName = "s"),
        "must identify its variants"
    )
})

test_that("a weights file read as a fine-mapping result points back", {
    data(twasWeightsExample, envir = environment())
    wt <- .rdr_tmp("x_wt.bed")
    suppressMessages(exportToBed(twasWeightsExample, wt))
    expect_error(
        readFineMappingResult(wt, weightStandardized = TRUE),
        "no `PIP` column.*readTwasWeights"
    )
})

test_that("a fine-mapping table reads as weights, because WGT is the weight", {
    # Not a category error in this direction: a susie posterior mean IS the
    # TWAS coefficient, which is why it publishes as WGT in both schemas.
    data(qtlFineMappingExample, envir = environment())
    fm <- .rdr_tmp("x_fm.bed")
    exportToBed(qtlFineMappingExample, fm)
    asWeights <- readTwasWeights(fm, weightStandardized = TRUE)
    expect_s4_class(asWeights, "TwasWeights")
    expect_equal(
        as.numeric(resolveWeights(asWeights)$weights),
        as.numeric(
            S4Vectors::mcols(
                topLoci(qtlFineMappingExample)
            )$posterior_effect
        )[order(variantIds(
            readFineMappingResult(fm, weightStandardized = TRUE)
        ))],
        tolerance = 1e-8
    )
})

# ---- the weight scale -----------------------------------------------------

test_that("a published weight without its scale has to be told", {
    data(twasWeightsExample, envir = environment())
    path <- .rdr_tmp("scale.bed")
    suppressMessages(exportToBed(twasWeightsExample, path))
    # This fixture records the flag, so it reads without one.
    expect_s4_class(readTwasWeights(path), "TwasWeights")
    df <- utils::read.delim(path, check.names = FALSE)
    df$weight_standardized <- NA
    stripped <- .rdr_writeTable(df, "scale_stripped.bed")
    expect_error(
        readTwasWeights(stripped),
        "publishes WGT but not `weight_standardized`"
    )
    expect_s4_class(
        readTwasWeights(stripped, weightStandardized = FALSE),
        "TwasWeights"
    )
})

# ---- opt-in Wakefield log Bayes factors ------------------------------------

test_that("priorVariance derives log Bayes factors only when asked", {
    df <- data.frame(
        `#chr` = c("chr1", "chr1"),
        start = c(99L, 199L),
        end = c(100L, 200L),
        ID = c("chr1:100:A:G", "chr1:200:C:T"),
        PIP = c(0.9, 0.05),
        ES = c(0.3, -0.1),
        SE = c(0.05, 0.05),
        study = "s",
        context = "c",
        trait = "t",
        method = "susie",
        check.names = FALSE
    )
    path <- .rdr_writeTable(df, "wakefield.bed")
    plain <- readFineMappingResult(path)
    expect_true(all(is.na(
        S4Vectors::mcols(
            topLoci(plain, signalCutoff = 0)
        )$logBF
    )))
    derived <- readFineMappingResult(path, priorVariance = 0.04)
    lbf <- S4Vectors::mcols(topLoci(derived, signalCutoff = 0))$logBF
    expect_equal(
        lbf,
        .approxBayesFactor(df$ES / df$SE, df$SE, W = 0.04)$lbf,
        tolerance = 1e-10
    )
    # Never silently over published values.
    df$LBF <- c(5, 0)
    withLbf <- .rdr_writeTable(df, "wakefield_has_lbf.bed")
    expect_error(
        readFineMappingResult(withLbf, priorVariance = 0.04),
        "already carries LBF"
    )
})

test_that("priorVariance without a z or a standard error says so", {
    df <- data.frame(
        `#chr` = "chr1",
        start = 99L,
        end = 100L,
        ID = "chr1:100:A:G",
        PIP = 0.9,
        study = "s",
        context = "c",
        trait = "t",
        method = "susie",
        check.names = FALSE
    )
    path <- .rdr_writeTable(df, "wakefield_bare.bed")
    expect_error(
        readFineMappingResult(path, priorVariance = 0.04),
        "publishes no z \\(EZ, or ES/SE\\)"
    )
})

# ---- EZ, and what it means --------------------------------------------------

test_that("EZ round-trips as a SUPPLIED z and is not invented", {
    # The spec defines EZ as the z "provided if it was used to derive ES and
    # SE". A z pecotmr computed as ES/SE must not come back claiming to be
    # the source's, or a reader would believe ES/SE were derived.
    df <- data.frame(
        `#chr` = "chr1",
        start = 99L,
        end = 100L,
        ID = "chr1:100:A:G",
        PIP = 0.9,
        ES = 0.3,
        SE = 0.05,
        study = "s",
        context = "c",
        trait = "t",
        method = "susie",
        check.names = FALSE
    )
    noEz <- readFineMappingResult(.rdr_writeTable(df, "ez_no.bed"))
    me <- marginalEffects(noEz)
    expect_equal(as.numeric(me$z), 6, tolerance = 1e-8)
    expect_false(any(as.logical(me$z_supplied)))
    df$EZ <- 5.5
    withEz <- readFineMappingResult(.rdr_writeTable(df, "ez_yes.bed"))
    me <- marginalEffects(withEz)
    expect_equal(as.numeric(me$z), 5.5)
    expect_true(all(as.logical(me$z_supplied)))
})

# ---- class choice, format choice, region ------------------------------------

test_that("a result with no context or trait is a GWAS result", {
    data(gwasFineMappingExample, envir = environment())
    path <- .rdr_tmp("gwas.bed")
    exportToBed(gwasFineMappingExample, path)
    back <- readFineMappingResult(path, weightStandardized = TRUE)
    expect_s4_class(back, "GwasFineMappingResult")
    expect_equal(
        .rdr_diff(
            as.data.frame(topLoci(gwasFineMappingExample)),
            as.data.frame(topLoci(back))
        ),
        character(0)
    )
})

test_that("the format follows the extension unless it is named", {
    expect_equal(pecotmr:::.prdFormat("x/y.bed.gz", "auto"), "bed")
    expect_equal(pecotmr:::.prdFormat("x/y.tsv", "auto"), "bed")
    expect_equal(pecotmr:::.prdFormat("x/y.vcf.gz", "auto"), "vcf")
    expect_equal(pecotmr:::.prdFormat("x/y.bcf", "auto"), "vcf")
    # A file named something else is read as whatever it is said to be.
    expect_equal(pecotmr:::.prdFormat("x/y.txt", "vcf"), "vcf")
})

test_that("a tabix-indexed table can be read by region", {
    data(qtlFineMappingExample, envir = environment())
    path <- .rdr_tmp("region.bed.gz")
    exportToBed(qtlFineMappingExample, path)
    all <- readFineMappingResult(path, weightStandardized = TRUE)
    full <- topLoci(all, signalCutoff = 0)
    window <- GenomicRanges::GRanges(
        as.character(GenomicRanges::seqnames(full))[[1L]],
        IRanges::IRanges(
            start = GenomicRanges::start(full)[[1L]],
            end = GenomicRanges::start(full)[[3L]]
        )
    )
    part <- readFineMappingResult(
        path,
        region = window,
        weightStandardized = TRUE
    )
    expect_lt(length(variantIds(part)), length(variantIds(all)))
    expect_true(all(is_in(variantIds(part), variantIds(all))))
})

test_that("an un-indexed file ignores a region loudly", {
    df <- data.frame(
        `#chr` = "chr1",
        start = 99L,
        end = 100L,
        ID = "chr1:100:A:G",
        PIP = 0.9,
        study = "s",
        context = "c",
        trait = "t",
        method = "susie",
        check.names = FALSE
    )
    path <- .rdr_writeTable(df, "noindex.bed")
    expect_warning(
        readFineMappingResult(path, region = "chr1:1-1000"),
        "not bgzipped \\+ tabix-indexed"
    )
})

test_that("a missing file names itself", {
    expect_error(
        readFineMappingResult(.rdr_tmp("absent_xyz.bed")),
        "does not exist"
    )
})

# ---- the ##SAMPLE key table -------------------------------------------------

test_that("every ##SAMPLE key is alphanumeric", {
    # Rsamtools' header parser splits the `<...>` on `,(?=[[:alnum:]]+=)`, so
    # an underscored key is not recognised as a key at all: it glues onto the
    # previous value and every field after it on the line is lost. Six of
    # these were unreadable for that reason.
    keys <- pecotmr:::.PUB_SAMPLE_KEYS
    expect_true(all(str_detect(keys, "^[A-Za-z0-9]+$")))
    expect_equal(length(unique(keys)), length(keys))
})

test_that("a per-fit field with no ##SAMPLE key stops the writer", {
    expect_error(
        pecotmr:::.pubSampleKeys(c("study", "brand_new_field")),
        "brand_new_field.*no\\s+##SAMPLE key"
    )
})

test_that("the ##SAMPLE block survives a VariantAnnotation read", {
    data(qtlFineMappingExample, envir = environment())
    path <- .rdr_tmp("sample_block.vcf")
    exportToVcf(qtlFineMappingExample, path, context = "context_1")
    meta <- pecotmr:::.prdVcfSampleMeta(VariantAnnotation::readVcf(path))
    # Named in the shared vocabulary, not the header's CamelCase: the reader
    # inverts the key table so everything downstream sees one set of names.
    expect_true(all(is_in(
        c("study", "context", "trait", "method", "region_start", "cs_index"),
        colnames(meta)
    )))
    expect_equal(meta$method, "susie")
})

test_that("a header-only file reads back as an empty collection", {
    # An export of a result whose every entry was screened out writes its
    # schema and no rows. That is a legitimate file, so it reads as an empty
    # collection rather than erroring about identity axes it has no rows to
    # record.
    header <- str_flatten(pecotmr:::.pubSchema(character(0)), "\t")
    path <- .rdr_tmp("empty.bed")
    writeLines(header, path)
    back <- readFineMappingResult(path)
    expect_s4_class(back, "QtlFineMappingResult")
    expect_equal(nrow(back), 0L)
    weights <- readTwasWeights(path)
    expect_s4_class(weights, "TwasWeights")
    expect_equal(nrow(weights), 0L)
})

test_that("an empty result survives an export and a read", {
    data(qtlFineMappingExample, envir = environment())
    path <- .rdr_tmp("empty_rt.bed")
    # A signal cutoff above 1 retains nothing, so the table is vacuous.
    exportToBed(qtlFineMappingExample, path, signalCutoff = 2)
    expect_equal(nrow(readFineMappingResult(path)), 0L)
})
