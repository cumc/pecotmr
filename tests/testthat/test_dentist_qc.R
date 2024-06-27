context("dentist_qc")
library(MASS)
library(corpcor)

generate_dentist_data <- function(seed=42, n_snps = 100, sample_size = 100, n_outliers = 5, start_pos = 1000000, end_pos = 4000000) {
    set.seed(seed)
    cor_matrix <- matrix(0, nrow = n_snps, ncol = n_snps)
    for (i in 1:(n_snps - 1)) {
        for (j in (i + 1):n_snps) {
            cor_matrix[i, j] <- runif(1, 0.2, 0.8)
            cor_matrix[j, i] <- cor_matrix[i, j]
        }
    }
    diag(cor_matrix) <- 1
    ld_matrix <- cov2cor(make.positive.definite(cor_matrix))
    z_scores <- mvrnorm(n = 1, mu = rep(0, n_snps), Sigma = ld_matrix)
    outlier_indices <- sample(1:n_snps, n_outliers)
    z_scores[outlier_indices] <- rnorm(n_outliers, mean = 0, sd = 5)
    sumstat <- data.frame(
        position = unlist(lapply(seq(start_pos,end_pos,length.out = n_snps), round)),
        z = z_scores
    )
    return(list(sumstat = sumstat, LD_mat = ld_matrix, nSample = sample_size))
}

test_that("Test dentist works", {
    data <- generate_dentist_data()
    expect_warning(expect_equal(length(dentist(data$sumstat, data$LD_mat, data$nSample)$imputed_z), 100))
})

test_that("Test dentist works correct_chen_et_al_bug = F", {
    data <- generate_dentist_data()
    expect_warning(expect_equal(length(dentist(data$sumstat, data$LD_mat, data$nSample, correct_chen_et_al_bug = F)$imputed_z), 100))
})

test_that("Test dentist stops when missing position", {
    data <- generate_dentist_data()
    colnames(data$sumstat) <- c("something", "z")
    expect_error(dentist(data$sumstat, data$LD_mat, data$nSample, correct_chen_et_al_bug = F))
})

test_that("Test dentist stops when missing zscore", {
    data <- generate_dentist_data()
    colnames(data$sumstat) <- c("position", "something")
    expect_error(dentist(data$sumstat, data$LD_mat, data$nSample, correct_chen_et_al_bug = F))
})

generate_dentist_single_window_data <- function(seed=42, n_snps = 100, sample_size = 100, n_outliers = 5) {
    set.seed(seed)
    cor_matrix <- matrix(0, nrow = n_snps, ncol = n_snps)
    for (i in 1:(n_snps - 1)) {
        for (j in (i + 1):n_snps) {
            cor_matrix[i, j] <- runif(1, 0.2, 0.8)
            cor_matrix[j, i] <- cor_matrix[i, j]
        }
    }
    diag(cor_matrix) <- 1
    ld_matrix <- cov2cor(make.positive.definite(cor_matrix))
    z_scores <- mvrnorm(n = 1, mu = rep(0, n_snps), Sigma = ld_matrix)
    outlier_indices <- sample(1:n_snps, n_outliers)
    z_scores[outlier_indices] <- rnorm(n_outliers, mean = 0, sd = 5)
    return(list(z_scores = z_scores, LD_mat = ld_matrix, nSample = sample_size))
}

test_that("Test dentist_single_window works", {
    data <- generate_dentist_single_window_data()
    expect_warning(expect_equal(length(dentist_single_window(data$z_scores, data$LD_mat, data$nSample)$original_z), 100))
})

test_that("Test dentist_single_window warning when < 2000 variants", {
    data <- generate_dentist_single_window_data()
    expect_warning(dentist_single_window(data$z_scores, data$LD_mat, data$nSample))
})

test_that("Test dentist_single_window works < 1.0 duprThreshold", {
    data <- generate_dentist_single_window_data()
    expect_warning(expect_equal(length(dentist_single_window(data$z_scores, data$LD_mat, data$nSample, duprThreshold = 0.99)$original_z), 100))
})

test_that("Test dentist_single_window stops with zscore/LD matrix dimension mismatch", {
    data <- generate_dentist_single_window_data()
    expect_warning(expect_error(dentist_single_window(generate_dentist_single_window_data()$z_scores, generate_dentist_single_window_data(n_snps = 80)$LD_mat, data$nSample)))
})

#add_dups_back_dentist <- function(zScore, dentist_output, find_dup_output) {
generate_add_dups_back_dentist_data <- function(seed=42, n_snps = 1000, sample_size = 1000, n_corr = 20, n_outliers = 5) {
    seed <- 42
    n_snps <- 100
    sample_size <- 100
    n_corr <- 20
    n_outliers <- 5
    rThreshold <- 0.2

    set.seed(seed)
    cor_matrix <- matrix(0, nrow = n_snps, ncol = n_snps)
    for (i in 1:(n_snps - 1)) {
        for (j in (i + 1):n_snps) {
            cor_matrix[i, j] <- runif(1, 0.2, 0.8)
            cor_matrix[j, i] <- cor_matrix[i, j]
        }
    }

    diag(cor_matrix) <- 1

    ld_matrix <- cov2cor(make.positive.definite(cor_matrix))

    z_scores <- mvrnorm(n = 1, mu = rep(0, n_snps), Sigma = ld_matrix)
    outlier_indices <- sample(1:n_snps, n_outliers)
    z_scores[outlier_indices] <- rnorm(n_outliers, mean = 0, sd = 5)

    find_dup_output <- find_duplicate_variants(z_scores, ld_matrix, 0.5)
    org_z_scores <- z_scores
    z_scores <- find_dup_output$filteredZ
    ld_matrix <- find_dup_output$filteredLD
    print(find_dup_output$dupBearer)
    dentist_output <- dentist_single_window(z_scores, ld_matrix, sample_size)

    return(list(z_scores = org_z_scores, dentist_output = dentist_output, find_dup_output = find_dup_output))
}

test_that("Test add_dups_back_dentist works", {
    zScore <- c(1.2, 2.3, 2.4, 1.4, 5.6)
    dentist_output <- data.frame(
        original_z = c(1.2, 2.3, 5.6),
        imputed_z = c(1.1, 2.1, 5.1),
        iter_to_correct = c(1, 1, 3),
        rsq = c(0.9, 0.8,  0.5),
        z_diff = c(0.1, 0.2, 0.5)
    )
    find_dup_output <- data.frame(
        dupBearer = c(-1, -1, 2, 1, -1),
        sign = c(1, -1, 1, -1, 1)
    )

    expected_output <- data.frame(
        original_z = zScore,
        imputed_z = c(1.1, 2.1, 2.1, -1.1, 5.1),
        iter_to_correct = c(1, 1, 1, 1, 3),
        rsq = c(0.9, 0.8, 0.8, 0.9, 0.5),
        z_diff = c(0.1, 0.2, 0.2, 0.1, 0.5),
        is_duplicate = c(FALSE, FALSE, TRUE, TRUE, FALSE)
    )
    res <- add_dups_back_dentist(zScore, dentist_output, find_dup_output)
    expect_equal(expected_output, res)
})

test_that("Test add_dups_back_dentist stops when nrow mismatch ", {
    z_scores <- rep(0, 5)
    dentist_output <- list(
        original_z = c(1, 2, 3, 4, 5),
        imputed_z = c(1, 2, 3, 4, 5),
        iter_to_correct = c(1, 2, 3, 4, 5),
        rsq = c(1, 2, 3, 4, 5),
        z_diff = c(1, 2, 3, 4, 5)
    )
    find_dup_output <- list(
        dupBearer = c(-1, -1, -1, -1, -1, -1),
        sign = c(1, 2, 3, 4, 5, 6)
    )
    expect_error(add_dups_back_dentist(z_scores, dentist_output, find_dup_output))
})

test_that("Test add_dups_back_dentist stops when zScore and find_dup_output nrow mismatch ", {
    z_scores <- rep(0, 6)
    dentist_output <- list(
        original_z = c(1, 2, 3, 4, 5),
        imputed_z = c(1, 2, 3, 4, 5),
        iter_to_correct = c(1, 2, 3, 4, 5),
        rsq = c(1, 2, 3, 4, 5),
        z_diff = c(1, 2, 3, 4, 5)
    )
    find_dup_output <- list(
        dupBearer = c(-1, -1, -1, -1, -1),
        sign = c(1, 2, 3, 4, 5)
    )
    expect_error(add_dups_back_dentist(z_scores, dentist_output, find_dup_output))
})

test_that("Test divide_into_windows works", {
    res <- divide_into_windows(seq(2000000,5000000,100000), 2000000, F)
    expect_equal(nrow(res), 2)
})

test_that("Test divide_into_windows works window_size <= 0", {
    res <- divide_into_windows(seq(2000000,5000000,100000), 0, F)
    expect_equal(nrow(res), 1)
})

test_that("Test divide_into_windows works correct_chen_et_al_bug", {
    res <- divide_into_windows(seq(2000000,5000000,100000), 2000000, T)
    expect_equal(nrow(res), 2)
})

#merge_windows <- function(dentist_result_by_window, window_divided_res) {
test_that("Test merge_windows works", {
    data <- generate_dentist_data(n_snps = 1000, sample_size = 1000, start_pos = 0, end_pos = 2000)
    window_divided_res <- divide_into_windows(data$sumstat$position, window_size = 1000, correct_chen_et_al_bug = TRUE)
    pValueThreshold <- 5.0369e-8
    propSVD <- 0.4
    gcControl <- FALSE
    nIter <- 10
    gPvalueThreshold <- 0.05
    duprThreshold <- 0.99
    ncpus <- 1
    seed <- 999
    correct_chen_et_al_bug <- TRUE
    dentist_result_by_window <- list()
    suppressWarnings({
        for (k in 1:nrow(window_divided_res)) {
        zScore_k <- data$sumstat$z[window_divided_res$windowStartIdx[k]:window_divided_res$windowEndIdx[k]]
        LD_mat_k <- data$LD_mat[
            window_divided_res$windowStartIdx[k]:window_divided_res$windowEndIdx[k],
            window_divided_res$windowStartIdx[k]:window_divided_res$windowEndIdx[k]
        ]
        dentist_result_by_window[[k]] <- dentist_single_window(
            zScore_k, LD_mat_k, 100,
            pValueThreshold, propSVD, gcControl,
            nIter, gPvalueThreshold, duprThreshold,
            ncpus, seed, correct_chen_et_al_bug
        )
        }
    })
    res <- merge_windows(dentist_result_by_window, window_divided_res)
    expect_equal(nrow(res), 1000)
})

test_that("Test merge_windows stops with window and imputed mismatch", {
    expect_error(merge_windows(rep(0, 5), data.frame(windowStartIdx = rep(0,2), windowEndIdx = rep(0, 2))))
})