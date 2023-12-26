context("LD")

# Use input data that was generated using ChatGPT 4
# Generated:
#   - LD_meta_file.csv
#   - region.csv
#   - LD_block_#.chr1_#_#.float16.txt.gz
#   - LD_block_#.bim

test_that("Check that we correctly retrieve the names from the matrix",{
  region <- read_delim("./test_data/region.csv", delim = ",")
  meta <- read_delim("./test_data/LD_meta_file.csv", delim = ",")
  res <- load_LD_matrix(meta, region)
  variants <- unlist(
    c(
      c("1:1000:A:G", "1:1040:A:G", "1:1080:A:G", "1:1120:A:G", "1:1160:A:G"),
      c("1:1400:A:G", "1:1440:A:G", "1:1480:A:G", "1:1520:A:G", "1:1560:A:G"),
      c("1:1800:A:G", "1:1840:A:G", "1:1880:A:G", "1:1920:A:G", "1:1960:A:G")
    ))
  expect_identical(
    unlist(res[[1]]$variants_df$variants),
    variants[1:5])
  expect_identical(
    unlist(res[[2]]$variants_df$variants),
    variants[7:10])
  expect_identical(
    unlist(res[[3]]$variants_df$variants),
    variants[12:15])
})

test_that("Check that the LD block has the appropriate rownames and colnames",{
  region <- read_delim("test_data/region.csv", delim = ",")
  meta <- read_delim("test_data/LD_meta_file.csv", delim = ",")
  res <- load_LD_matrix(meta, region)
  variants <- unlist(
    c(
      c("1:1000:A:G", "1:1040:A:G", "1:1080:A:G", "1:1120:A:G", "1:1160:A:G"),
      c("1:1400:A:G", "1:1440:A:G", "1:1480:A:G", "1:1520:A:G", "1:1560:A:G"),
      c("1:1800:A:G", "1:1840:A:G", "1:1880:A:G", "1:1920:A:G", "1:1960:A:G")
    ))
  expect_identical(rownames(res[[1]]$LD), variants[1:5])
  expect_identical(colnames(res[[1]]$LD), variants[1:5])
  expect_identical(rownames(res[[2]]$LD), variants[7:10])
  expect_identical(colnames(res[[2]]$LD), variants[7:10])
  expect_identical(rownames(res[[3]]$LD), variants[12:15])
  expect_identical(colnames(res[[3]]$LD), variants[12:15])
})

test_that("Check that the LD block contains the correct information",{
  region <- read_delim("test_data/region.csv", delim = ",")
  meta <- read_delim("test_data/LD_meta_file.csv", delim = ",")
  res <- load_LD_matrix(meta, region)
  # Variant names
  variants <- unlist(
    c(
      c("1:1000:A:G", "1:1040:A:G", "1:1080:A:G", "1:1120:A:G", "1:1160:A:G"),
      c("1:1400:A:G", "1:1440:A:G", "1:1480:A:G", "1:1520:A:G", "1:1560:A:G"),
      c("1:1800:A:G", "1:1840:A:G", "1:1880:A:G", "1:1920:A:G", "1:1960:A:G")
    ))
  # Check LD Block 1
  ld_block_one <- res[[1]]$LD
  ld_block_one_original <- as.matrix(
    read_delim(
      "test_data/LD_block_1.chr1_1000_1200.float16.txt.xz",
      delim = " ", col_names = F))
  rownames(ld_block_one_original) <- colnames(ld_block_one_original) <- variants[1:5]
  expect_equal(ld_block_one, ld_block_one_original)
  # Check LD Block 2
  ld_block_two <- res[[2]]$LD
  ld_block_two_original <- as.matrix(
    read_delim(
      "test_data/LD_block_3.chr1_1400_1600.float16.txt.xz",
      delim = " ", col_names = F))[2:5,2:5]
  rownames(ld_block_two_original) <- colnames(ld_block_two_original) <- variants[7:10]
  expect_equal(ld_block_two, ld_block_two_original)
  # Check LD Block 3
  ld_block_three <- res[[3]]$LD
  ld_block_three_original <- as.matrix(
    read_delim(
      "test_data/LD_block_5.chr1_1800_2000.float16.txt.xz",
      delim = " ", col_names = F))[2:5,2:5]
  rownames(ld_block_three_original) <- colnames(ld_block_three_original) <- variants[12:15]
  expect_equal(ld_block_three, ld_block_three_original)
})