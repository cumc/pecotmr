context("mr")
library(tidyverse)

# Define a function to randomly generate ref:alt pairs
library(testthat)
generate_format_mock_data <- function(seed = 1, num_variants = 10, empty_sets = F) {
generate_ref_alt <- function(num_variants){
  ref_alleles <- c("A", "T", "G", "C")
  alt_alleles <- c("A", "T", "G", "C")
  pairs <- paste0(sample(ref_alleles, num_variants, replace = TRUE), ":", sample(alt_alleles, num_variants, replace = TRUE))
  # Ensure that ref is not the same as alt
  while(any(sapply(strsplit(pairs, ":"), function(x) x[1] == x[2]))) {
    pairs <- paste0(sample(ref_alleles, num_variants, replace = TRUE), ":", sample(alt_alleles, num_variants, replace = TRUE))
  }
  return(pairs)
}
top_loci_mock <- data.frame(
  variant_id = paste0("chr1:", 1:num_variants, ":", generate_ref_alt(num_variants)),
  betahat = rnorm(num_variants),
  sebetahat = runif(num_variants, 0.05, 0.1),
  #maf = runif(n_entries, 0.1, 0.5),
  pip = runif(num_variants, 0, 1),
  cs_coverage_0.95 = sample(0:2, num_variants, replace = TRUE),
  cs_coverage_0.7 = sample(0:2, num_variants, replace = TRUE),
  cs_coverage_0.5 = sample(0:2, num_variants, replace = TRUE)
)
susie_result_mock <- list(
  susie_results = list(
    condition1 = list(
      top_loci = top_loci_mock,
      region_info = list(region_name = "Gene1")
    )
  )
) 
gwas_sumstats_db_mock <- data.frame(
  variant_id = paste0("chr1:", 1:num_variants, ":", generate_ref_alt(num_variants)),
  beta = rnorm(num_variants),
  se = runif(num_variants, 0.05, 0.1)
)
return(
        list(
            susie_result = susie_result_mock,
            gwas_sumstats_db = gwas_sumstats_db_mock
          )
    )

}


test_that("calc_I2 works with dummy data",{
    # Test with Q > 1e-3 and positive I2
    expect_equal(calc_I2(10, c(1, 2, 3)), (10 - 3 + 1)/10)

    # Test with Q exactly 1e-3 and I2 should be 0
    expect_equal(calc_I2(1e-3, c(1, 2)), 0)

    # Test with Q < 1e-3 and I2 should be 0
    expect_equal(calc_I2(1e-4, c(1, 2, 3, 4)), 0)

    # Test with negative Q and I2 should be 0
    expect_equal(calc_I2(-10, c(1, 2, 3)), 0)

    # Test with Q leading to negative I2, I2 should be 0
    expect_equal(calc_I2(5, c(1, 2, 3, 4, 5, 6)), 0)

    # Test with Q as 0, I2 should be 0
    expect_equal(calc_I2(0, c(1, 2)), 0)

    # Test with Est as empty vector, I2 should be 0
    expect_equal(calc_I2(10, c()), 1.1)
})

test_that("mr_format functions with normal parameters", {
  input_data <- generate_format_mock_data()  

  condition <- "condition1"
  coverage <- "cs_coverage_0.95"
  allele_qc <- TRUE
  
  res <- mr_format(input_data$susie_result, condition, input_data$gwas_sumstats_db, coverage, allele_qc)
   expect_true(all(c("gene_name", "variant_id", "bhat_x", "sbhat_x", "cs", "pip", "bhat_y", "sbhat_y") %in% names(res)))
 })
test_that("mr_format returns a dataframe with NAs for zero coverage in top_loci", {
  # Assume generate_format_mock_data generates mock data including susie_result with top_loci
  # and gwas_sumstats_db as per the structure required by the mr_format function.
  input_data <- generate_format_mock_data()

  condition <- "condition1"
  coverage <- "cs_coverage_0.95"
  allele_qc <- TRUE

  # Create a mock susie_result with zero coverage in top_loci
  susie_result_mock <- input_data$susie_result
  susie_result_mock[["susie_results"]][[condition]][["top_loci"]][[coverage]] <- rep(0, nrow(susie_result_mock[["susie_results"]][[condition]][["top_loci"]]))
  
  # Run the mr_format function with the mock data where coverage is zero
  result <- mr_format(susie_result_mock, condition, input_data$gwas_sumstats_db, coverage = coverage, allele_qc = allele_qc)
  
  # Check if the result is a dataframe with the expected NAs
  expect_true(is.data.frame(result))
  expect_true(all(is.na(result[,-1])))
})                   
# Test whether function handles non-existent top_loci
test_that("mr_format returns a dataframe with NAs with non-existent top_loci", {
  input_data <- generate_format_mock_data()

  condition <- "condition1"
  coverage <- "cs_coverage_0.95"
  allele_qc <- TRUE

 
  # Create a mock susie_result with zero coverage in top_loci
  susie_result_mock <- input_data$susie_result
  susie_result_mock[["susie_results"]][[condition]][["top_loci"]] <- list()

  
  # Run the mr_format function with the mock data where coverage is zero
  result <- mr_format(susie_result_mock, condition, input_data$gwas_sumstats_db, coverage = coverage, allele_qc = allele_qc)

  # Check if the result is a dataframe with the expected NAs
  expect_true(is.data.frame(result))
  expect_true(all(is.na(result[,-1])))
})

generate_mock_mr_formatted_input <- function(num_variants = NULL, generate_full_dataset = TRUE) {
  generate_ref_alt <- function(num_variants){
  ref_alleles <- c("A", "T", "G", "C")
  alt_alleles <- c("A", "T", "G", "C")
  pairs <- paste0(sample(ref_alleles, num_variants, replace = TRUE), ":", sample(alt_alleles, num_variants, replace = TRUE))
  # Ensure that ref is not the same as alt
  while(any(sapply(strsplit(pairs, ":"), function(x) x[1] == x[2]))) {
    pairs <- paste0(sample(ref_alleles, num_variants, replace = TRUE), ":", sample(alt_alleles, num_variants, replace = TRUE))
  }
  return(pairs)
  }
  if (generate_full_dataset) {
    data.frame(
      gene_name = rep("Gene1", num_variants),
      #cs = sample(1:2, num_variants, replace = TRUE),
      cs = as.integer(rep(1,num_variants)),
      variant_id = paste0("1:", 1:num_variants, ":", generate_ref_alt(num_variants)),
      bhat_x = rnorm(num_variants),
      sbhat_x = runif(num_variants, 0.1, 0.2),
      bhat_y = rnorm(num_variants),
      sbhat_y = runif(num_variants, 0.1, 0.2),
      pip = runif(num_variants, 0, 1),
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
     gene_name = "Gene1",
     variant_id = as.character(NA),
     bhat_x = as.numeric(NA),
     sbhat_x = as.numeric(NA),
     cs = as.numeric(NA),
     pip = as.numeric(NA),
     bhat_y = as.numeric(NA),
     sbhat_y = as.numeric(NA),
     stringsAsFactors = FALSE # Optional, to prevent factors
     )
    }
  }
# Test 1: normal input data
test_that("mr_analysis returns expected output with normal inputs", {
  input_data <- generate_mock_mr_formatted_input(num_variants = 10, generate_full_dataset = TRUE)
  result <- mr_analysis(input_data, cpip_cutoff=0.5)
  expect_true(is.data.frame(result))
  expect_gt(nrow(result), 0)
  expect_true(all(c("gene_name", "num_CS", "num_IV", "cpip", "meta_eff", "se_meta_eff", "meta_pval", "Q", "Q_pval", "I2") %in% names(result)))
})
# Test 2: All NA values except gene_name
test_that("mr_analysis returns null output for all NA input except gene_name", {
  input_data <- generate_mock_mr_formatted_input(generate_full_dataset = FALSE)
  result <- mr_analysis(input_data, cpip_cutoff=0.5)
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 1)
  expect_true(all(is.na(result[,-1])))
})
# Test 3: No significant cpip values
test_that("mr_analysis handles no significant cpip values correctly", {
  input_data <- generate_mock_mr_formatted_input(num_variants = 5, generate_full_dataset = TRUE)
  input_data$pip <- runif(nrow(input_data), 0, 0.1) # Setting low pip values to ensure cpip < cpip_cutoff
  result <- mr_analysis(input_data, cpip_cutoff = 0.5)
  expect_true(is.data.frame(result))
  #Expecting an empty or null output depending on how you handle this case
  expect_equal(nrow(result), 1)
  expect_true(all(is.na(result[,-1])))
})     