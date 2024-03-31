context("misc")
library(tidyverse)
library(data.table)

test_that("Test compute_maf freq 0.5",{
    expect_equal(compute_maf(rep(1, 20)), 0.5)
})

test_that("Test compute_maf freq 0.6",{
    expect_equal(compute_maf(rep(1.2, 20)), 0.4)
})

test_that("Test compute_maf freq 0.3",{
    expect_equal(compute_maf(rep(0.6, 20)), 0.3)
})

test_that("Test compute_maf with NA",{
    set.seed(1)
    generate_small_dataset <- function(sample_size = 20) {
        vals <- c(1.2, NA)
        return(sample(vals, sample_size, replace = TRUE))
    }
    expect_equal(compute_maf(generate_small_dataset()), 0.4)
})

test_that("test compute_missing",{
    small_dataset <- c(rep(NA, 20), rep(1, 80))
    expect_equal(compute_missing(small_dataset), 0.2)
})

test_that("Test compute_non_missing_y",{
    small_dataset <- c(rep(NA, 20), rep(1, 80))
    expect_equal(compute_non_missing_y(small_dataset), 80)
})

test_that("Test compute_all_missing_y",{
    small_dataset <- c(rep(NA, 20), rep(1, 80))
    expect_equal(compute_all_missing_y(small_dataset), F)
})

test_that("Test mean_impute",{
    dummy_data <- matrix(c(1,2,NA,1,2,3), nrow=3, ncol=2)
    expect_equal(mean_impute(dummy_data)[3,1], 1.5)
})

test_that("Test is_zero_variance",{
    dummy_data <- matrix(c(1,2,3,1,1,1), nrow=3, ncol=2)
    col <- which(apply(dummy_data, 2, is_zero_variance))
    expect_equal(col, 2)
})

test_that("Test filter_X",{
    dummy_data <- matrix(
        c(1,NA,NA,NA, 0,0,1,1, 2,2,2,2, 1,1,1,2, 2,2,0,1, 0,1,1,2),
        # Missing Rate, MAF thresh, Zero Var, Var Thresh, Regular values
        nrow=4, ncol=6)
    var_thres <- 0.3
    expect_equal(filter_X(dummy_data, 0.70, 0.3, var_thres = 0.3), matrix(c(2,2,0,1, 0,1,1,2), nrow=4, ncol=2))
})

test_that("Test filter_Y non-matrix",{
    dummy_data <- matrix(c(1,NA,NA,NA, 1,1,2,NA), nrow=4, ncol=2)
    res <- filter_Y(as.data.frame(dummy_data), 3)
    expect_equal(length(res$Y), 3)
    expect_equal(res$rm_rows, NULL)
})

test_that("Test filter_Y is-matrix",{
    dummy_data <- matrix(c(1,NA,NA,NA, 1,1,2,NA, 2,1,2,NA), nrow=4, ncol=3)
    expect_equal(nrow(filter_Y(dummy_data, 3)$Y), 3)
    expect_equal(ncol(filter_Y(dummy_data, 3)$Y), 2)
    expect_equal(length(filter_Y(dummy_data, 3)$rm_rows), 1)
})

test_that("Test format_variant_id",{
    expect_equal(format_variant_id(c("chr1_123_G_C", "chr1_132_A_T")), c("chr1:123:G:C", "chr1:132:A:T"))
})