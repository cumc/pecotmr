# Tests for S4 classes (h2_classes.R), GwasSumStats (h2_sumstats.R),
# and AnnotationMatrix (h2Annotations.R)

# =============================================================================
# Test data helpers
# =============================================================================

makeTestGRanges <- function(n = 10) {
    GenomicRanges::GRanges(
        seqnames = rep("chr1", n),
        ranges = IRanges::IRanges(
            start = seq(1000, by = 100, length.out = n),
            width = 1L
        )
    )
}

makeTestSumstatsDf <- function(n = 50) {
    set.seed(42)
    data.frame(
        SNP = paste0("rs", seq_len(n)),
        CHR = rep("1", n),
        BP = seq(1000, by = 100, length.out = n),
        A1 = rep("A", n),
        A2 = rep("G", n),
        Z = rnorm(n),
        N = rep(10000, n),
        stringsAsFactors = FALSE
    )
}

makeTestLdBlocks <- function() {
    blocks_gr <- GenomicRanges::GRanges(
        seqnames = c("chr1", "chr1"),
        ranges = IRanges::IRanges(start = c(1, 5001), end = c(5000, 10000))
    )
    blocks_gr
}

makeTestSnpInfo <- function(n = 10) {
    data.frame(
        SNP = paste0("rs", seq_len(n)),
        CHR = rep("1", n),
        BP = seq(1000, by = 100, length.out = n),
        A1 = rep("A", n),
        A2 = rep("G", n),
        stringsAsFactors = FALSE
    )
}

makeTestAnnotationMeta <- function() {
    data.frame(
        name = c("base", "enhancer", "promoter"),
        tier = c("baseline", "candidate", "candidate"),
        type = c("binary", "binary", "continuous"),
        stringsAsFactors = FALSE
    )
}

# Bridge helper: turn the legacy per-study data.frame shape into a
# single-row GwasSumStats collection using the new API. Keeps the bulk of
# the per-accessor tests below readable.
.testGenotypeHandle <- function() {
    new(
        "GenotypeHandle",
        path = "/tmp/test",
        format = "plink2",
        snpInfo = data.frame(),
        nSamples = 0L,
        sampleIds = character(),
        pgenPtr = NULL
    )
}

.dfToSumstatsGr <- function(df) {
    chrs <- as.character(df$CHR)
    if (!all(grepl("^chr", chrs))) {
        chrs <- paste0("chr", sub("^chr", "", chrs))
    }
    gr <- GenomicRanges::GRanges(
        seqnames = chrs,
        ranges = IRanges::IRanges(start = as.integer(df$BP), width = 1L)
    )
    mc <- df[, setdiff(colnames(df), c("CHR", "BP")), drop = FALSE]
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(mc)
    gr
}

makeGwasSumStatsFromDf <- function(
    df,
    traitName = "test",
    genome = "hg19",
    varY = NA_real_
) {
    required <- c("SNP", "A1", "A2", "Z", "N")
    missingCols <- setdiff(required, colnames(df))
    if (length(missingCols) > 0L) {
        stop("Missing required columns: ", paste(missingCols, collapse = ", "))
    }
    keep <- stats::complete.cases(df[, required, drop = FALSE])
    if (!all(keep)) {
        message(sprintf(
            "Removed %d SNPs with missing required-column values.",
            sum(!keep)
        ))
    }
    df <- df[keep, , drop = FALSE]
    if (is.null(varY)) {
        varY <- NA_real_
    }
    GwasSumStats(
        studyName = traitName,
        entry = list(.dfToSumstatsGr(df)),
        genome = genome,
        ldReference = .testGenotypeHandle(),
        varY = varY
    )
}

# A handle over a real (fake-path) variant set, so tests can build the
# per-block LD panels a block-parallel pipeline produces: same genotype source,
# each narrowed to its own variants.
.blockGenotypeHandle <- function(bp, path = "/tmp/test.pgen") {
    new(
        "GenotypeHandle",
        path = path,
        format = "plink2",
        snpInfo = data.frame(
            SNP = paste0("chr1:", bp, ":A:G"),
            CHR = rep("1", length(bp)),
            BP = as.integer(bp),
            A1 = rep("G", length(bp)),
            A2 = rep("A", length(bp)),
            fileIdx = seq_along(bp),
            stringsAsFactors = FALSE
        ),
        nSamples = 10L,
        sampleIds = paste0("s", seq_len(10)),
        pgenPtr = NULL
    )
}

# One LD block's worth of GWAS summary statistics: `n` variants starting at
# `start`, carrying a QC record and an LD panel narrowed to just those
# variants (what a per-block pipeline step writes out).
makeGwasBlock <- function(
    blockId,
    start,
    n = 5L,
    genome = "hg19",
    qcOptions = list(mafCutoff = 0.01),
    ldReference = .blockGenotypeHandle(bp)
) {
    bp <- seq(start, by = 100L, length.out = n)
    df <- data.frame(
        SNP = paste0("chr1:", bp, ":A:G"),
        CHR = rep("1", n),
        BP = bp,
        A1 = rep("A", n),
        A2 = rep("G", n),
        Z = seq_len(n) / 10,
        N = rep(1000L, n),
        stringsAsFactors = FALSE
    )
    obj <- GwasSumStats(
        studyName = "trait1",
        entry = list(.dfToSumstatsGr(df)),
        genome = genome,
        ldReference = ldReference,
        blockId = blockId,
        qcInfo = list(
            timestamp = NA_character_,
            options = qcOptions,
            entryAudit = list(list(variantsIn = n, block = blockId))
        )
    )
    obj
}


# =============================================================================
# LdMatrix fixtures for buildLdEigen() / buildLdScore()
# =============================================================================

# One block's variants as LdMatrix() takes them: SNP ids, alleles and an AF
# column, at 100 bp spacing.
.testLdVariants <- function(n, chrom, startBp) {
    bp <- seq(startBp, by = 100L, length.out = n)
    data.frame(
        chrom = chrom,
        pos = as.integer(bp),
        SNP = paste0(chrom, ":", bp, ":C:T"),
        A1 = "T",
        A2 = "C",
        AF = seq(0.2, 0.8, length.out = n)
    )
}

.testArBlock <- function(n, rho) {
    rho^abs(outer(seq_len(n), seq_len(n), "-"))
}

# One LD block: an LdMatrix with AR(1) LD.
makeTestLdMatrix <- function(
    n = 6L,
    chrom = "chr1",
    startBp = 1000L,
    rho = 0.5,
    nSamples = 500L
) {
    LdMatrix(
        .testArBlock(n, rho),
        variants = .testLdVariants(n, chrom, startBp),
        nSamples = nSamples
    )
}

# Several consecutive LD blocks along one chromosome, one LdMatrix each.
makeTestLdMatrixBlocks <- function(
    sizes = c(4L, 3L),
    chrom = "chr1",
    startBp = 1000L,
    rho = 0.5,
    nSamples = 500L
) {
    starts <- startBp + 100L * c(0L, cumsum(head(sizes, -1L)))
    map2(
        sizes,
        starts,
        function(n, s) makeTestLdMatrix(n, chrom, s, rho, nSamples)
    )
}

# A genotype panel of 100 samples cut into three LD blocks of 60, 140 and 100
# variants, with the same blocks' LD as LdMatrix objects: a statistic read
# from the panel should agree with one built from the matrices. The middle
# block has more variants than the panel has samples, so its LD has 40 zero
# eigenvalues that the SVD of the panel's factor does not list.
.h2PanelFixture <- function() {
    panel <- readGenotypes(test_path("test_data", "test_variants.bed"))
    gr <- SummarizedExperiment::rowRanges(panel)
    pos <- GenomicRanges::start(gr)
    cuts <- list(1:60, 61:200, 201:300)
    blocks <- data.frame(
        chrom = "chr21",
        start = map_int(cuts, function(i) pos[[min(i)]]),
        end = map_int(cuts, function(i) pos[[max(i)]])
    )
    ld <- map(cuts, function(i) {
        LdMatrix(
            unname(computeLd(panel, snpIdx = i)),
            variants = gr[i],
            nSamples = ncol(panel)
        )
    })
    list(panel = panel, blocks = blocks, ld = ld)
}

# Two baseline annotations (binary and continuous) and one candidate over the
# variants of `ref`.
.h2PanelAnnotations <- function(ref) {
    set.seed(11)
    p <- length(ref)
    mat <- cbind(
        coding = rbinom(p, 1, 0.3),
        cons = runif(p),
        cand = rbinom(p, 1, 0.2)
    )
    AnnotationMatrix(
        mat,
        GenomicRanges::granges(ref, use.mcols = FALSE),
        data.frame(
            name = colnames(mat),
            tier = c("baseline", "baseline", "candidate"),
            type = c("binary", "continuous", "binary")
        ),
        genome = "hg19"
    )
}

# GWAS z-scores over the fixture's blocks with heritability `h2`:
# z = sqrt(n) R beta + e, e ~ N(0, R), block by block.
.h2PanelZ <- function(ld, n, h2, seed = 3) {
    set.seed(seed)
    M <- sum(map_int(ld, length))
    list_c(map(ld, function(x) {
        R <- ldMatrix(x)
        beta <- rnorm(nrow(R), sd = sqrt(h2 / M))
        e <- eigen(R, symmetric = TRUE)
        noise <- e$vectors %*% (sqrt(pmax(e$values, 0)) * rnorm(nrow(R)))
        as.vector(sqrt(n) * R %*% beta + noise)
    }))
}
