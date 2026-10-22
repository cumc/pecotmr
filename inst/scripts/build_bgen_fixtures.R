# =============================================================================
# Copy the BGEN reader's test fixtures from bgen's own examples.
# -----------------------------------------------------------------------------
# complex.bgen and haplotypes.bgen are the bgen project's example files
# (github.com/gavinband/bgen), distributed under the Boost Software License
# 1.0, whose text is copied beside them.
#
#   complex.bgen     haploid and diploid samples, phased and unphased
#                    records, multi-allelic records and ploidy above 2;
#                    complex.vcf holds its probabilities as text, the truth
#                    the tests compare against.
#   haplotypes.bgen  phased haplotypes, with haplotypes.haps as its truth.
#
# The other BGEN fixtures -- each compression, odd bit depths, missing
# samples, a file without sample ids, layout 1 -- are written by the tests
# themselves (tests/testthat/helper-bgen.R).
#
# Run from the package root, with the bgen repository checked out at
# $BGEN_REPO (default ~/ghq/github.com/gavinband/bgen):
#   Rscript inst/scripts/build_bgen_fixtures.R
# =============================================================================

bgenRepo <- Sys.getenv(
    "BGEN_REPO",
    file.path(Sys.getenv("HOME"), "ghq", "github.com", "gavinband", "bgen")
)
source <- file.path(bgenRepo, "example")
dest <- file.path("tests", "testthat", "test_data", "bgen")
files <- c(
    "complex.bgen",
    "complex.bgen.bgi",
    "complex.sample",
    "complex.vcf",
    "haplotypes.bgen",
    "haplotypes.bgen.bgi",
    "haplotypes.haps"
)
stopifnot(dir.exists(source), all(file.exists(file.path(source, files))))
dir.create(dest, showWarnings = FALSE, recursive = TRUE)
stopifnot(all(file.copy(file.path(source, files), dest, overwrite = TRUE)))
stopifnot(file.copy(
    file.path(bgenRepo, "LICENSE_1_0.txt"),
    dest,
    overwrite = TRUE
))
