context("raiss")
library(tidyverse)
library(MASS)

generate_dummy_data <- function(seed=1, ref_panel_ordered=TRUE, known_zscores_ordered=TRUE) {
    set.seed(seed)

    n_variants <- 100
    ref_panel <- data.frame(
        chrom = rep(1, n_variants),
        pos = seq(1, n_variants * 10, 10),
        variant_id = paste0("rs", seq_len(n_variants)),
        A0 = sample(c("A", "T", "G", "C"), n_variants, replace = TRUE),
        A1 = sample(c("A", "T", "G", "C"), n_variants, replace = TRUE)
    )

    n_known <- 50
    known_zscores <- data.frame(
        chrom = rep(1, n_known),
        pos = sample(ref_panel$pos, n_known),
        variant_id = sample(ref_panel$variant_id, n_known),
        A0 = sample(c("A", "T", "G", "C"), n_known, replace = TRUE),
        A1 = sample(c("A", "T", "G", "C"), n_known, replace = TRUE),
        Z = rnorm(n_known)
    )

    LD_matrix <- matrix(rnorm(n_variants^2), nrow = n_variants, ncol = n_variants)
    diag(LD_matrix) <- 1 
    known_zscores <- if (known_zscores_ordered) known_zscores[order(known_zscores$pos),] else known_zscores
    ref_panel <- if (ref_panel_ordered) ref_panel else ref_panel[order(ref_panel$pos, decreasing = TRUE),]
    return(list(ref_panel=ref_panel, known_zscores=known_zscores, LD_matrix=LD_matrix))
}

test_that("Input validation for raiss works correctly", {
    input_data_ref_panel_unordered <- generate_dummy_data(ref_panel_ordered=FALSE)
    input_data_zscores_unordered <- generate_dummy_data(known_zscores_ordered=FALSE)
    expect_error(raiss(input_data_ref_panel_unordered$ref_panel, input_data$known_zscores, input_data$LD_matrix))
    expect_error(raiss(input_data$ref_panel, input_data_zscores_unordered$known_zscores, input_data$LD_matrix))
})

test_that("Default parameters for raiss work correctly", {
    # TODO - ask Gao about merging on removed columns
    input_data <- generate_dummy_data()
    result <- raiss(input_data$ref_panel, input_data$known_zscores, input_data$LD_matrix)
    expect_true(is.data.frame(result))
})

test_that("Test Default Parameters for raiss_model", {
  zt <- c(1.2, 0.5)
  sig_t <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
  sig_i_t <- matrix(c(0.5, 0.2, 0.2, 0.5), nrow = 2)

  result <- raiss_model(zt, sig_t, sig_i_t)

  expect_is(result, "list")
  expect_true(all(names(result) %in% c("var", "mu", "ld_score", "condition_number", "correct_inversion")))
})

test_that("Test with Different lamb Values for raiss_model", {
  zt <- c(1.2, 0.5)
  sig_t <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
  sig_i_t <- matrix(c(0.5, 0.2, 0.2, 0.5), nrow = 2)

  lamb_values <- c(0.01, 0.05, 0.1)
  for (lamb in lamb_values) {
    result <- raiss_model(zt, sig_t, sig_i_t, lamb)
    expect_is(result, "list")
  }
})

test_that("Test Batch processing in raiss_model", {
  # TODO 
})

test_that("Report Condition Number in raiss_model", {
  zt <- c(1.2, 0.5)
  sig_t <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
  sig_i_t <- matrix(c(0.5, 0.2, 0.2, 0.5), nrow = 2)

  result_with_cn <- raiss_model(zt, sig_t, sig_i_t, report_condition_number = TRUE)
  result_without_cn <- raiss_model(zt, sig_t, sig_i_t, report_condition_number = FALSE)

  expect_is(result_with_cn, "list")
  expect_is(result_without_cn, "list")
})

test_that("Input Validation of raiss_model", {

  zt <- c(1.2, 0.5)
  sig_t <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
  sig_i_t <- matrix(c(0.5, 0.2, 0.2, 0.5), nrow = 2)
  zt_invalid <- "not a numeric vector"
  sig_t_invalid <- "not a matrix"
  sig_i_t_invalid <- "not a matrix"

  expect_error(raiss_model(zt_invalid, sig_t, sig_i_t))
  expect_error(raiss_model(zt, sig_t_invalid, sig_i_t))
  expect_error(raiss_model(zt, sig_t, sig_i_t_invalid))
})

test_that("Boundary Conditions of raiss_model", {
  zt <- c(1.2, 0.5)
  sig_t <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
  sig_i_t <- matrix(c(0.5, 0.2, 0.2, 0.5), nrow = 2)

  zt_empty <- numeric(0)
  sig_t_empty <- matrix(numeric(0), nrow = 0)
  sig_i_t_empty <- matrix(numeric(0), nrow = 0)

  expect_error(raiss_model(zt_empty, sig_t, sig_i_t))
  expect_error(raiss_model(zt, sig_t_empty, sig_i_t))
  expect_error(raiss_model(zt, sig_t, sig_i_t_empty))
})

test_that("Test with Different rcond Values for raiss_model", {
  zt <- c(1.2, 0.5)
  sig_t <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
  sig_i_t <- matrix(c(0.5, 0.2, 0.2, 0.5), nrow = 2)

  rcond_values <- c(0.01, 0.05, 0.1)
  for (rcond in rcond_values) {
    result <- raiss_model(zt, sig_t, sig_i_t, lamb = 0.01, rcond = rcond)
    expect_is(result, "list")
    expect_true(all(names(result) %in% c("var", "mu", "ld_score", "condition_number", "correct_inversion")))
  }
})

test_that("format_raiss_df returns correctly formatted data frame", {
  imp <- list(
    mu = rnorm(5),
    var = runif(5),
    ld_score = rnorm(5),
    condition_number = runif(5),
    correct_inversion = sample(c(TRUE, FALSE), 5, replace = TRUE)
  )
  
  ref_panel <- data.frame(
    chrom = sample(1:22, 10, replace = TRUE),
    pos = sample(1:10000, 10),
    variant_id = paste0("rs", 1:10),
    A0 = sample(c("A", "T", "G", "C"), 10, replace = TRUE),
    A1 = sample(c("A", "T", "G", "C"), 10, replace = TRUE)
  )

  unknowns <- sample(1:nrow(ref_panel), 5)

  result <- format_raiss_df(imp, ref_panel, unknowns)

  expect_true(is.data.frame(result))
  expect_equal(ncol(result), 10)
  expect_equal(colnames(result), c('chrom', 'pos', 'variant_id', 'A0', 'A1', 'Z', 'Var', 'ld_score', 'condition_number', 'correct_inversion'))

  for (col in c('chrom', 'pos', 'variant_id', 'A0', 'A1')) {
    expect_equal(setNames(unlist(result[col]), NULL), unlist(ref_panel[unknowns, col, drop = TRUE]))
  }
  for (col in c('Z', 'Var', 'ld_score', 'condition_number', 'correct_inversion')) {
    expected_col <- if (col == "Z") "mu" else if (col == "Var") "var" else col
    expect_equal(setNames(unlist(result[col]), NULL), setNames(unlist(imp[expected_col]), NULL))
  }
})

test_that("Merge operation is correct for merge_raiss_df", {
    raiss_df_example <- data.frame(
        chrom = c("chr21", "chr22"),
        pos = c(123, 456),
        variant_id = c("var1", "var2"),
        A0 = c("A", "T"),
        A1 = c("T", "A"),
        Z = c(0.5, 1.5),
        Var = c(0.2, 0.3),
        ld_score = c(10, 20),
        imputation_R2 = c(0.8, 0.7))

    known_zscores_example <- data.frame(
        chrom = c("chr21", "chr22"),
        pos = c(123, 456),
        variant_id = c("var1", "var2"),
        A0 = c("A", "T"),
        A1 = c("T", "A"),
        Z = c(0.5, 1.5))

    merged_df <- merge_raiss_df(raiss_df_example, known_zscores_example)
    expect_equal(nrow(merged_df), 2)
    expect_true(all(c("chr21", "chr22") %in% merged_df$chrom))
})

generate_fro_test_data <- function(seed=1) {
    set.seed(seed)
    return(data.frame(
        chrom = paste0("chr", rep(22, 10)),
        pos = seq(1, 100, 10),
        variant_id = 1:10,
        A0 = rep("A", 10),
        A1 = rep("T", 10),
        Z = rnorm(10),
        Var = runif(10, 0, 1),
        ld_score = rnorm(10, 5, 2)
    ))
}

test_that("Correct columns are selected in filter_raiss_output", {
    test_data <- generate_fro_test_data()
    output <- filter_raiss_output(test_data)
    expect_true(all(c('variant_id', 'A0', 'A1', 'Z', 'Var', 'ld_score') %in% names(output)))
})

test_that("imputation_R2 is calculated correctly in filter_raiss_output", {
    test_data <- generate_fro_test_data()
    output <- filter_raiss_output(test_data)
    expected_R2 <- 1 - test_data[which(test_data$ld_score >= 5),]$Var
    expect_equal(output$imputation_R2, expected_R2[which(expected_R2 > 0.6)])
})

test_that("Filtering is applied correctly in filter_raiss_output", {
    test_data <- generate_fro_test_data()
    R2_threshold <- 0.6
    minimum_ld <- 5
    output <- filter_raiss_output(test_data, R2_threshold, minimum_ld)

    expect_true(all(output$imputation_R2 > R2_threshold))
    expect_true(all(output$ld_score >= minimum_ld))
})

test_that("Function returns the correct subset in filter_raiss_output", {
    test_data <- generate_fro_test_data()
    test_data$imputation_R2 <- 1 - test_data$Var
    output <- filter_raiss_output(test_data)

    manual_filter <- test_data[test_data$imputation_R2 > 0.6 & test_data$ld_score >= 5, ]

    expect_equal(nrow(output), nrow(manual_filter))
    expect_equal(sum(output$variant_id != manual_filter$variant_id), 0)
})

test_that("compute_mu basic functionality", {
    sig_i_t <- matrix(c(1, 2, 3, 4), nrow = 2)
    sig_t_inv <- matrix(c(5, 6, 7, 8), nrow = 2)
    zt <- matrix(c(9, 10, 11, 12), nrow = 2)

    expected_result <- matrix(c(517, 766, 625, 926), nrow = 2)
    result <- compute_mu(sig_i_t, sig_t_inv, zt)
    expect_equal(result, expected_result)
})

generate_mock_data_for_compute_var <- function(seed=1) {
    return(
        list(
            sig_i_t_1 = matrix(c(1, 2, 3, 4), nrow = 2),
            sig_t_inv_1 = matrix(c(5, 6, 7, 8), nrow = 2),
            lamb_1 = 0.5))
}

test_that("compute_var returns correct output for batch = TRUE", {
    input_data <- generate_mock_data_for_compute_var()
    result <- compute_var(input_data$sig_i_t_1, input_data$sig_t_inv_1, input_data$lamb_1, batch = TRUE)
    expect_is(result, "list")
    expect_length(result, 2)
    expect_true("var" %in% names(result))
    expect_true("ld_score" %in% names(result))
})

test_that("compute_var returns correct output for batch = FALSE", {
    input_data <- generate_mock_data_for_compute_var()
    result <- compute_var(input_data$sig_i_t_1, input_data$sig_t_inv_1, input_data$lamb_1, batch = FALSE)
    expect_is(result, "list")
    expect_length(result, 2)
    expect_true("var" %in% names(result))
    expect_true("ld_score" %in% names(result))
})

test_that("check_inversion correctly identifies inverse matrices in", {
  sig_t <- matrix(c(1, 2, 3, 4), nrow=2, ncol=2)
  sig_t_inv <- solve(sig_t)  
  expect_true(check_inversion(sig_t, sig_t_inv))
})

test_that("var_in_boundaries sets boundaries correctly", {
  lamb_test <- 0.05
  var <- c(-1, 0, 0.5, 1.04, 1.05)  

  result <- var_in_boundaries(var, lamb_test)

  expect_equal(result[1], 0)                   # Value less than 0 should be set to 0
  expect_equal(result[2], 0)                   # Value within lower boundary should remain unchanged
  expect_equal(result[3], 0.5)                 # Value within boundaries should remain unchanged
  expect_equal(result[4], 1.04)                   # Value greater than 0.99999 + lamb should be set to 1
  expect_equal(result[5], 1)                   # Value greater than 0.99999 + lamb should be set to 1
})

test_that("invert_mat computes correct pseudo-inverse", {
  mat <- matrix(c(1, 2, 3, 4), nrow = 2)
  lamb <- 0.5
  rcond <- 1e-7
  result <- invert_mat(mat, lamb, rcond)
  expect_true(is.matrix(result))
})

test_that("invert_mat handles errors and retries", {
  mat <- matrix(c(0, 0, 0, 0), nrow = 2) 
  lamb <- 0.1
  rcond <- 1e-7
  result <- invert_mat(mat, lamb, rcond)
  expect_true(is.matrix(result))
})

test_that("invert_mat_recursive correctly inverts a valid square matrix", {
  mat <- matrix(c(2, -1, -1, 2), nrow = 2)
  lamb <- 0.5
  rcond <- 0.01
  result <- invert_mat_recursive(mat, lamb, rcond)
  expect_true(is.matrix(result))
  expect_equal(dim(result), dim(mat))
})

test_that("invert_mat_recursive handles non-square matrices appropriately", {
  mat <- matrix(1:6, nrow = 2)
  lamb <- 0.5
  rcond <- 0.01
  expect_silent(invert_mat_recursive(mat, lamb, rcond))
})

test_that("invert_mat_recursive handles errors and performs recursive call correctly", {
  mat <- "not a matrix"
  lamb <- 0.5
  rcond <- 0.01
  expect_error(invert_mat_recursive(mat, lamb, rcond))
})

# Test with Different Tolerance Levels
test_that("invert_mat_eigen behaves differently with varying tolerance levels", {
  mat <- matrix(c(1, 0, 0, 1e-4), nrow = 2)
  tol_high <- 1e-2
  tol_low <- 1e-6
  result_high_tol <- invert_mat_eigen(mat, tol_high)
  result_low_tol <- invert_mat_eigen(mat, tol_low)
  expect_true(!is.logical(all.equal(result_high_tol, result_low_tol)))
})

test_that("invert_mat_eigen handles non-square matrices", {
  mat <- matrix(1:6, nrow = 2)
  expect_error(invert_mat_eigen(mat))
})

test_that("invert_mat_eigen returns the same matrix for an identity matrix", {
    mat <- diag(2)
    expected <- mat
    actual <- invert_mat_eigen(mat)
    expect_equal(actual, expected)
})

test_that("invert_mat_eigen returns a zero matrix for a zero matrix input", {
    mat <- matrix(0, nrow = 2, ncol = 2)
    expected <- mat
    actual <- invert_mat_eigen(mat)
    expect_equal(actual, expected)
})

test_that("invert_mat_eigen handles matrices with negative eigenvalues", {
    mat <- matrix(c(-2, 0, 0, -3), nrow = 2)
    expect_silent(invert_mat_eigen(mat))
})