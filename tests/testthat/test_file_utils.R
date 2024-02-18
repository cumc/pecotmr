context("file_utils")
library(tidyverse)

test_that("read_pvar dummy data works",{
    dummy_path <- gsub("//", "/", tempfile(pattern = "dummy_pvar", tmpdir = tempdir(), fileext = ".pvar"))
    dummy <- data.frame("#CHROM" = c(1, 2, 3, 4, 5),
        "ID" = c("rs1", "rs2", "rs3", "rs4", "rs5"),
        "POS" = c(100, 200, 300, 400, 500),
        "REF" = c("A", "T", "C", "G", "A"),
        "ALT" = c("T", "C", "G", "A", "T"))
    colnames(dummy) <- c("#CHROM", "ID", "POS", "REF", "ALT")
    cat(c("#DUMMY HEADER 1", "#DUMMY HEADER 2", "#DUMMY HEADER 3"), file = dummy_path, sep = "\n")
    write_delim(
        dummy, dummy_path, delim = "\t", col_names = TRUE, append = TRUE)
    expect_equal(colnames(read_pvar(dummy_path)), c("chrom", "id", "pos", "alt", "ref"))
    file.remove(dummy_path)
})

test_that("read_bim dummy data works",{
    example_path <- "test_data/protocol_example.genotype.bed"
    res <- read_bim(example_path)
    expect_equal(colnames(res), c("chrom", "id", "gpos", "pos", "a1", "a0"))
    expect_equal(nrow(res), 100)
})

test_that("read_psam dummy data works",{
    dummy_path <- gsub("//", "/", tempfile(pattern = "dummy_psam", tmpdir = tempdir(), fileext = ".psam"))
    dummy <- data.frame("#CHROM" = c(1, 2, 3, 4, 5),
        "IID" = c("rs1", "rs2", "rs3", "rs4", "rs5"),
        "SID" = c(100, 200, 300, 400, 500),
        "PAT" = c("A", "T", "C", "G", "A"),
        "MAT" = c("T", "C", "G", "A", "T"),
        "SEX" = c(1, 2, 1, 2, 1))
    write_delim(
        dummy, dummy_path, delim = "\t", col_names = TRUE, append = TRUE)
    res <- read_psam(dummy_path)
    expect_equal(colnames(res), c("FID", "IID", "SID", "PAT", "MAT", "SEX"))
    file.remove(dummy_path)
})

test_that("read_fam dummy data works",{
    example_path <- "test_data/protocol_example.genotype.bed"
    res <- read_fam(example_path)
    expect_equal(nrow(res), 100)
})

test_that("open_pgen dummy data works",{
    example_path <- "test_data/dummy_data.pgen"
    res <- open_pgen(example_path)
    expect_equal(res$class, "pgen")
})

test_that("read_pgen dummy data works",{
    example_path <- "test_data/dummy_data.pgen"
    res <- open_pgen(example_path)
    expect_equal(res$class, "pgen")
})

test_that("open_bed dummy data works",{
    example_path <- "test_data/protocol_example.genotype.bed"
    res <- open_bed(example_path)
    expect_equal(res$class, "pgen")
})

test_that("find_valid_file_path works",{
    ref_path <- "test_data/protocol_example.genotype.bed"
    expect_error(
        find_valid_file_path(paste0(ref_path, "s"), "protocol_example.genotype.bamf"),
        "Both reference and target file paths do not work. Tried paths: 'test_data/protocol_example.genotype.beds' and 'test_data/protocol_example.genotype.bamf'")
    expect_equal(
        find_valid_file_path(ref_path, "abc"),
        ref_path)
    expect_equal(
        find_valid_file_path(ref_path, "protocol_example.genotype.bim"),
        "test_data/protocol_example.genotype.bim")
    expect_equal(
        find_valid_file_path(ref_path, "test_data/protocol_example.genotype.bim"),
        "test_data/protocol_example.genotype.bim")
})
