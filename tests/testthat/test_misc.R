context("misc")

dummy_geno_data <- function(number_of_samples = 100, number_of_snps = 1000, with_indels = FALSE, missing_rate = 0, maf = 0, var_thresh = 0) {
    set.seed(1)
    X <- matrix(
        sample(c(0,1,2), number_of_samples*number_of_snps, replace = TRUE),
        nrow=number_of_snps, ncol=number_of_samples)
    colnames(X) <- paste0("chr1:", seq(1000,1000+number_of_snps), "_G_C")
    rownames(X) <- paste0("sample_", seq(1, number_of_samples))

    if (missing_rate > 0) {
        X[sample(length(X), missing_rate*length(X))] <- NA
    }

    return(X)
}

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
    expect_equal(compute_all_missing_y(small_dataset), 20)
})

test_that("Test mean_impute",{
    dummy_data <- matrix(c(1,2,NA,1,2,3), nrow=3, ncol=2)
    expect_equal(mean_impute(dummy_data)[3,1], 1.5)
})

test_that("Test is_zero_variance",{
    dummy_data <- matrix(c(1,2,3,1,1,1), nrow=3, ncol=2)
    col <- which(apply(X, 2, is_zero_variance))
    expect_equal(col, 2)
})

test_that("Test filter_X",{
    dummy_data <- matrix(
        c(1,NA,NA,NA, 0,0,1,1, 2,2,2,2, 1,1,1,2, 2,2,0,1, 0,1,1,2),
        # Missing Rate, MAF thresh, Zero Var, Var Thresh, Regular values
        nrow=4, ncol=6)
    var_thres <- 0.3
    expect_equal(filter_X(dummy_data, 0.70, 0.5, 0.3), matrix(c(2,2,0,1, 0,1,1,2), nrow=4, ncol=2))
})

test_that("Test filter_Y non-matrix",{
    dummy_data <- matrix(c(1,NA,NA,NA, 1,1,2,NA), nrow=4, ncol=2)
    expect_equal(nrow(filter_Y(dummy_data, 3)$Y), 3)
    expect_equal(ncol(filter_Y(dummy_data, 3)$Y), 1)
    expect_equal(length(filter_Y(dummy_data, 3)$rm_rows), 1)
})

test_that("Test filter_Y is-matrix",{
    dummy_data <- matrix(c(1,NA,NA,NA, 1,1,2,NA, 0,1,2,NA), nrow=4, ncol=2)
    expect_equal(nrow(filter_Y(dummy_data, 3)$Y), 3)
    expect_equal(ncol(filter_Y(dummy_data, 3)$Y), 2)
    expect_equal(length(filter_Y(dummy_data, 3)$rm_rows), 1)
})

test_that("Test format_variant_id",{
    expect_equal(format_variant_id(c("chr1_123_G_C", "chr1_132_A_T")), c("chr1:123:G:C", "chr1:132:A:T"))
})

test_that("Test load_genotype_data",{
  set.seed(1)
})

test_that("Test load_covariate_data",{
  set.seed(1)
})

test_that("Test load_phenotype_data",{
  set.seed(1)
})

test_that("Test filter_by_common_samples",{
    common_samples <- c("Sample_1", "Sample_2", "Sample_3")
    dat <- as.data.frame(matrix(c(1,2,3,4,5,6,7,8), nrow=4, ncol=2))
    rownames(dat) <- c("Sample_1", "Sample_2", "Sample_3", "Sample_4")
    colnames(dat) <- c("chr1:122:G:C", "chr1:123:G:C")
    expect_equal(nrow(filter_by_common_samples(dat, common_samples)), 3)
    expect_equal(rownames(filter_by_common_samples(dat, common_samples)), common_samples)
})

test_that("Test prepare_data_list",{
    # Create dummy data
    dummy_geno_data <- as.data.frame(matrix(
        c(1,NA,NA,NA, 0,0,1,1, 2,2,2,2, 1,1,1,2, 2,2,0,1, 0,1,1,2),
        # Missing Rate, MAF thresh, Zero Var, Var Thresh, Regular values
        nrow=4, ncol=6))
    rownames(dummy_geno_data) <- c("Sample_1", "Sample_2", "Sample_3", "Sample_14")
    colnames(dummy_geno_data) <- c("chr1:122:G:C", "chr1:123:G:C", "chr1:124:G:C", "chr1:125:G:C", "chr1:126:G:C", "chr1:127:G:C")
    dummy_pheno_data <- as.data.frame(matrix(c(1,NA,NA,NA, 1,1,2,NA, 0,1,2,NA), nrow=4, ncol=2))
    rownames(dummy_pheno_data) <- c("Sample_10", "Sample_1", "Sample_2", "Sample_3")
    colnames(dummy_pheno_data) <- c("pheno_1", "pheno_2")
    dummy_covar_data <- as.data.frame(matrix(c(70,71,72,73, 28,30,15,20, 1,2,3,4), nrow=4, ncol=3))
    rownames(dummy_covar_data) <- c("Sample_1", "Sample_2", "Sample_3", "Sample_4")
    colnames(dummy_covar_data) <- c("covar_1", "covar_2", "covar_3")
    # Set parameters
    imiss_cutoff <- 0.70
    maf_cutoff <- 0.5
    mac_cutoff <- 2.4
    xvar_cutoff <- 0.3
    keep_samples <- c("Sample_1", "Sample_2", "Sample_3")
    res <- prepare_data_list(dummy_geno_data, dummy_pheno_data, dummy_covar_data, imiss_cutoff, maf_cutoff, mac_cutoff, xvar_cutoff, keep_samples=keep_samples)
    # Check that Covar, X, and Y have the same number of rows
    expect_equal(nrow(res$covar), 3)
    expect_equal(nrow(res$X), 3)
    expect_equal(nrow(res$Y), 3)
    # Check that filter_X occured
    expect_equal(ncol(res$X), 2)
    # Check that Covar, X, and Y have the same samples
    expect_equal(rownames(res$covar), rownames(res$X))
    expect_equal(rownames(res$covar), rownames(res$Y))
    expect_equal(rownames(res$X), rownames(res$Y))
})

test_that("Test prepare_X_matrix",{
    dummy_geno_data <- as.data.frame(matrix(
        c(1,NA,NA,NA,2, 0,0,1,1,0, 2,2,2,2,2, 1,1,1,2,2, 2,2,0,1,2, 0,1,1,2,2),
        # Missing Rate, MAF thresh, Zero Var, Var Thresh, Regular values
        nrow=5, ncol=6))
    rownames(dummy_geno_data) <- c("Sample_1", "Sample_2", "Sample_3", "Sample_4", "Sample_5")
    colnames(dummy_geno_data) <- c("chr1:122:G:C", "chr1:123:G:C", "chr1:124:G:C", "chr1:125:G:C", "chr1:126:G:C", "chr1:127:G:C")
    dummy_covar_data <- as.data.frame(
        matrix(
            c(70,71,72,73, 28,30,15,20, 1,2,3,4),
            nrow=4, ncol=3))
    rownames(dummy_covar_data) <- c("Sample_1", "Sample_2", "Sample_3", "Sample_4")
    colnames(dummy_covar_data) <- c("covar_1", "covar_2", "covar_3")
    dummy_data_list <- data.frame(
        covar = dummy_covar_data)
    # Set parameters
    imiss_cutoff <- 0.70
    maf_cutoff <- 0.5
    mac_cutoff <- 2.4
    xvar_cutoff <- 0.3
    expect_equal(
        prepare_X_matrix(dummy_geno_data, dummy_data_list, imiss_cutoff, maf_cutoff, mac_cutoff, xvar_cutoff),
        matrix(c(2,2,0,1, 0,1,1,2), nrow=4, ncol=2))
})

test_that("Test add_X_residuals",{
    dummy_geno_data <- matrix(
        c(2,2,0,1, 0,1,1,2),
        nrow=4, ncol=2)
    dummy_covar_data <- matrix(
        c(70,71,72,73, 28,30,15,20, 1,2,3,4),
        nrow=4, ncol=3)
    dummy_data_list %>%
        mutate(
            lm_res_X = map2(X, covar, ~ .lm.fit(x = cbind(1, .y), y = .x)$residuals %>% as.matrix()))
  set.seed(1)
})

test_that("Test add_Y_residuals",{
  set.seed(1)
})

test_that("Test load_regional_association_data",{
  set.seed(1)
})

test_that("Test load_regional_univariate_data",{
  set.seed(1)
})

test_that("Test load_regional_regression_data",{
  set.seed(1)
})

test_that("Test load_regional_multivariate_data",{
  set.seed(1)
})