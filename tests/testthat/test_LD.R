context("LD")

# Use input data that was generated using ChatGPT 4
# Generated:
#   - LD_meta_file.csv
#   - region.csv
#   - LD_block_#.chr1_#_#.float16.txt.gz
#   - LD_block_#.bim

test_that("Check that we correctly retrieve the names from the matrix",{
  region <- read_delim("test_data/region.csv", delim = ",")
  meta <- read_delim("test_data/LD_meta_file.csv", delim = ",")
  res <- load_LD_matrix(meta, region)
  variants <- unlist(
    c(
      paste0("chr1:100", seq(0,4), ":A:G"),
      paste0("chr1:300", seq(0,4), ":A:G"),
      paste0("chr1:500", seq(0,4), ":A:G")
    ))
  expect_identical(res$variants_df$variants, variants)
})

test_that("Check that the LD block has the appropriate rownames and colnames",{
  region <- read_delim("test_data/region.csv", delim = ",")
  meta <- read_delim("test_data/LD_meta_file.csv", delim = ",")
  res <- load_LD_matrix(meta, region)
  variants <- unlist(
    c(
      paste0("chr1:100", seq(0,4), ":A:G"),
      paste0("chr1:300", seq(0,4), ":A:G"),
      paste0("chr1:500", seq(0,4), ":A:G")
    ))
  expect_identical(rownames(res$LD), variants)
  expect_identical(colnames(res$LD), variants)
})

test_that("Check that the LD block contains the correct information",{
  region <- read_delim("test_data/region.csv", delim = ",")
  meta <- read_delim("test_data/LD_meta_file.csv", delim = ",")
  res <- load_LD_matrix(meta, region)
  # Variant names
  variants <- unlist(
    c(
      paste0("chr1:100", seq(0,4), ":A:G"),
      paste0("chr1:300", seq(0,4), ":A:G"),
      paste0("chr1:500", seq(0,4), ":A:G")
    ))
  # Check LD Block 1
  ld_block_one <- res$LD %>%
    subset(., rownames(.) %in% variants[1:5]) %>%
    subset(select=variants[1:5])
  ld_block_one_original <- as.matrix(
    read_delim(
      "/Users/travyseedwards/Documents/Research/PECOTMR/12-20-2023/tests/testthat/test_data/LD_block_1.chr1_1000_2000.float16.txt.xz",
      delim = "\t", col_names = T)[-1])
  rownames(ld_block_one_original) <- colnames(ld_block_one_original) <- gsub("_", ":", colnames(ld_block_one_original))
  expect_equal(ld_block_one, ld_block_one_original)
  # Check LD Block 2
  ld_block_two <- res$LD %>%
    subset(., rownames(.) %in% variants[6:10]) %>%
    subset(select=variants[6:10])
  ld_block_two_original <- as.matrix(
    read_delim(
      "/Users/travyseedwards/Documents/Research/PECOTMR/12-20-2023/tests/testthat/test_data/LD_block_2.chr1_3000_4000.float16.txt.xz",
      delim = "\t", col_names = T)[-1])
  rownames(ld_block_two_original) <- colnames(ld_block_two_original) <- gsub("_", ":", colnames(ld_block_two_original))
  expect_equal(ld_block_two, ld_block_two_original)
  # Check LD Block 3
  ld_block_three <- res$LD %>%
    subset(., rownames(.) %in% variants[11:15]) %>%
    subset(select=variants[11:15])
  ld_block_three_original <- as.matrix(
    read_delim(
      "/Users/travyseedwards/Documents/Research/PECOTMR/12-20-2023/tests/testthat/test_data/LD_block_3.chr1_5000_6000.float16.txt.xz",
      delim = "\t", col_names = T)[-1])
  rownames(ld_block_three_original) <- colnames(ld_block_three_original) <- gsub("_", ":", colnames(ld_block_three_original))
  expect_equal(ld_block_three, ld_block_three_original)
})