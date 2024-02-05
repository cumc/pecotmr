context("LD")
library(tidyverse)

generate_dummy_data <- function() {
  region <- data.frame(
    chrom = "chr1",
    start = c(1000),
    end = c(1190)
  )
  meta_df <- data.frame(
    chrom = "chr1",
    start = c(1000, 1200, 1400, 1600, 1800),
    end = c(1200, 1400, 1600, 1800, 2000),
    path = c(
      "./test_data/LD_block_1.chr1_1000_1200.float16.txt.xz,./test_data/LD_block_1.chr1_1000_1200.float16.bim",
      "./test_data/LD_block_2.chr1_1200_1400.float16.txt.xz,./test_data/LD_block_2.chr1_1200_1400.float16.bim",
      "./test_data/LD_block_3.chr1_1400_1600.float16.txt.xz,./test_data/LD_block_3.chr1_1400_1600.float16.bim",
      "./test_data/LD_block_4.chr1_1600_1800.float16.txt.xz,./test_data/LD_block_4.chr1_1600_1800.float16.bim",
      "./test_data/LD_block_5.chr1_1800_2000.float16.txt.xz,./test_data/LD_block_5.chr1_1800_2000.float16.bim"
    ))
  return(list(region = region, meta = meta_df))
}


test_that("Check that we correctly retrieve the names from the matrix",{
  data <- generate_dummy_data()
  region <- data$region
  LD_meta_file_path <- gsub("//", "/", tempfile(pattern = "ld_meta_file", tmpdir = tempdir(), fileext = ".RDS"))
  write_delim(data$meta, LD_meta_file_path, delim = "\t")
  res <- load_LD_matrix(LD_meta_file_path, region)
  variants <- unlist(
    c("1:1000:A:G", "1:1040:A:G", "1:1080:A:G", "1:1120:A:G", "1:1160:A:G"))
  expect_equal(
    unlist(res$combined_LD_variants),
    variants)
  file.remove(LD_meta_file_path)
})

test_that("Check that the LD block has the appropriate rownames and colnames",{
  data <- generate_dummy_data()
  region <- data$region
  LD_meta_file_path <- gsub("//", "/", tempfile(pattern = "ld_meta_file", tmpdir = tempdir(), fileext = ".RDS"))
  write_delim(data$meta, LD_meta_file_path, delim = "\t")
  res <- load_LD_matrix(LD_meta_file_path, region)
  variants <- unlist(
    c("1:1000:A:G", "1:1040:A:G", "1:1080:A:G", "1:1120:A:G", "1:1160:A:G"))
  expect_identical(rownames(res$combined_LD_matrix), variants)
  expect_identical(colnames(res$combined_LD_matrix), variants)
  file.remove(LD_meta_file_path)
})

test_that("Check that the LD block contains the correct information",{
  data <- generate_dummy_data()
  region <- data$region
  LD_meta_file_path <- gsub("//", "/", tempfile(pattern = "ld_meta_file", tmpdir = tempdir(), fileext = ".RDS"))
  write_delim(data$meta, LD_meta_file_path, delim = "\t")
  res <- load_LD_matrix(LD_meta_file_path, region)
  # Variant names
  variants <- unlist(
    c("1:1000:A:G", "1:1040:A:G", "1:1080:A:G", "1:1120:A:G", "1:1160:A:G"))
  # Check LD Block 1
  ld_block_one <- res$combined_LD_matrix
  ld_block_one_original <- as.matrix(
    read_delim(
      "test_data/LD_block_1.chr1_1000_1200.float16.txt.xz",
      delim = " ", col_names = F))
  rownames(ld_block_one_original) <- colnames(ld_block_one_original) <- variants
  expect_equal(ld_block_one, ld_block_one_original)
  file.remove(LD_meta_file_path)
})