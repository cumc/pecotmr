context("LD")

# Use LD matrices that were generated
# Use simulated data from PLINK

test_that("Check that we correctly retrieve the names from the matrix",{
  region <- read_delim("test_data/region.csv", delim = ",")
  meta <- read_delim("test_data/LD_meta_file.csv", delim = ",")
  load_LD_matrix(meta, region)
})



#test_that("Check that the block is correctly acquired",{
#  set.seed(1)
#})