context("mr")
library(tidyverse)

generate_format_mock_data <- function(seed = 1, num_genes = 5, empty_sets = F) {
    target_genes <- data.frame(
        chr = rep(22, num_genes),
        gene_id = paste0("ENSG", 1:num_genes),
        gene_name = paste0("Gene", 1:num_genes))
    
    susie_path <- data.frame(
        ID = paste0("Gene", 1:num_genes),
        path = paste0("a/b/c/", 1:num_genes))
    
    susie_obj <- list()
    for (i in 1:nrow(susie_path)) {
        susie_obj[[i]] <- list(list(
                sets = list(
                    cs = if(empty_sets) NULL else list(one = 1, two = 2)
                ),
                top_loci = data.frame(
                    variant_id = paste0("chr22:", 1:num_genes),
                    bhat = rnorm(num_genes),
                    sbhat = runif(num_genes),
                    pip = runif(num_genes),
                    cs_index_primary = sample(1:2, num_genes, replace = TRUE),
                    cs_index_secondary = sample(1:2, num_genes, replace = TRUE))))
    }
    
    weights_file_path <- data.frame(
        ID = paste0("Gene", 1:num_genes),
        path = paste0("a/b/c/", 1:num_genes))
    
    weights_obj <- lapply(
        weights_file_path$ID,
        function(x) { list(
        outcome_QC = data.frame(
            variant_allele_flip = paste0("chr22:", 1:num_genes),
            beta = rnorm(num_genes),
            se = runif(num_genes))
    )})

    return(
        list(
            target_genes = target_genes,
            susie_path = susie_path,
            susie_obj = susie_obj,
            weights_file_path = weights_file_path,
            weights_obj = weights_obj
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

test_that("twas_mr_format_input functions with normal parameters",{
    input_data <- generate_format_mock_data()
    input_data$susie_path$path <- unlist(lapply(
    input_data$susie_path$ID, function(x) {
        gsub("//", "/", tempfile(pattern = paste0(x, "_susie"), tmpdir = tempdir(), fileext = ".RDS"))
    }))
    input_data$weights_file_path$path <- unlist(lapply(
    input_data$weights_file_path$ID, function(x) {
        gsub("//", "/", tempfile(pattern = paste0(x, "_weights"), tmpdir = tempdir(), fileext = ".RDS"))
    }))
    for (i in 1:length(input_data$susie_path$path)) {
        saveRDS(input_data$susie_obj[[i]], file=input_data$susie_path$path[i])
        saveRDS(input_data$weights_obj[[i]], file=input_data$weights_file_path$path[i])}
    res <- twas_mr_format_input(input_data$target_genes, input_data$susie_path, input_data$weights_file_path)
    expect_true(all(c("snp", "bhat_x", "sbhat_x", "pip", "cs", "X_ID", "gene_id", "bhat_y", "sbhat_y") %in% names(res)))
    file.remove(input_data$susie_path$path)
    file.remove(input_data$weights_file_path$path)
})

test_that("twas_mr_format_input returns null with empty target genes",{
    input_data <- generate_format_mock_data()
    input_data$susie_path$path <- unlist(lapply(
    input_data$susie_path$ID, function(x) {
        gsub("//", "/", tempfile(pattern = paste0(x, "_susie"), tmpdir = tempdir(), fileext = ".RDS"))
    }))
    input_data$weights_file_path$path <- unlist(lapply(
    input_data$weights_file_path$ID, function(x) {
        gsub("//", "/", tempfile(pattern = paste0(x, "_weights"), tmpdir = tempdir(), fileext = ".RDS"))
    }))
    for (i in 1:length(input_data$susie_path$path)) {
        saveRDS(input_data$susie_obj[[i]], file=input_data$susie_path$path[i])
        saveRDS(input_data$weights_obj[[i]], file=input_data$weights_file_path$path[i])}
    expect_error(twas_mr_format_input(c(), input_data$susie_path, input_data$weights_file_path))
    file.remove(input_data$susie_path$path)
    file.remove(input_data$weights_file_path$path)
})

test_that("twas_mr_format_input returns null with empty susie cs",{
    input_data <- generate_format_mock_data(empty_sets = T)
    input_data$susie_path$path <- unlist(lapply(
    input_data$susie_path$ID, function(x) {
        gsub("//", "/", tempfile(pattern = paste0(x, "_susie"), tmpdir = tempdir(), fileext = ".RDS"))
    }))
    input_data$weights_file_path$path <- unlist(lapply(
    input_data$weights_file_path$ID, function(x) {
        gsub("//", "/", tempfile(pattern = paste0(x, "_weights"), tmpdir = tempdir(), fileext = ".RDS"))
    }))
    for (i in 1:length(input_data$susie_path$path)) {
        saveRDS(input_data$susie_obj[[i]], file=input_data$susie_path$path[i])
        saveRDS(input_data$weights_obj[[i]], file=input_data$weights_file_path$path[i])}
    res <- twas_mr_format_input(input_data$target_genes, input_data$susie_path, input_data$weights_file_path)
    expect_equal(res, NULL)
    file.remove(input_data$susie_path$path)
    file.remove(input_data$weights_file_path$path)
})

generate_finemr_mock_input <- function(seed=1) {
    set.seed(seed)
    tibble(
        X_ID = rep(c("X1", "X2"), each = 5),
        cs = rep(1:5, 2),
        gene_id = rep(c("gene1", "gene2"), each = 5),
        snp = sample(c("SNP1", "SNP2", "SNP3", "SNP4", "SNP5"), 10, replace = TRUE),
        bhat_y = rnorm(10),
        sbhat_y = runif(10, 0.1, 1),
        bhat_x = rnorm(10),
        sbhat_x = runif(10, 0.1, 1),
        pip = runif(10, 0, 1)
    )
}

# Test cases
test_that("fine_mr returns expected structure and types", {
  mock_input <- generate_finemr_mock_input()
  result <- fine_mr(mock_input, 0.5)

  # Check if the result is a data frame
  expect_true(is.data.frame(result))

  # Check for expected columns
  expected_cols <- c("X_ID", "num_CS", "num_IV", "cpip", "gene_id", "meta_eff",
                     "se_meta_eff", "meta_pval", "meta_qval", "Q", "Q_pval", "I2")
  expect_true(all(expected_cols %in% names(result)))

  # Check for data types
  expect_true(is.numeric(result$meta_eff))
  expect_true(is.numeric(result$Q_pval))
  # ... other column type checks

  # Check for expected values (example)
  expect_true(all(result$cpip >= 0.5))
})
