context("misc")
library(tidyverse)
library(data.table)

dummy_geno_data <- function(
    number_of_samples = 10, number_of_snps = 10, sample_start_id = 1,
    number_missing = 10, number_low_maf = 10, number_zero_var = 10, number_var_thresh = 10) {
    set.seed(1)
    # Create portion of Matrix with satisfactory values
    X <- matrix(
        sample(c(0,1,2), number_of_samples*number_of_snps, replace = TRUE),
        nrow=number_of_samples, ncol=number_of_snps)
    # Create portion of Matrix that should get pruned
    ## Missing Rate
    if (number_missing > 0) {
        X_missing <- rbind(
            matrix(
                sample(c(0,1,2), (number_of_samples-3)*number_of_snps, replace = TRUE),
                nrow=number_of_samples-3, ncol=number_of_snps),
            matrix(
                rep(NA, 3*number_of_snps), nrow=3, ncol=number_of_snps))
        X <- cbind(X, X_missing)
    }
    ## MAF
    if (number_low_maf > 0) {
        X_maf <- matrix(
            rep(0.1, number_of_samples*number_of_snps), nrow=number_of_samples, ncol=number_of_snps)
        X <- cbind(X, X_maf)
    }
    ## Zero Variance
    if (number_zero_var > 0) {
        X_zerovar <- matrix(
            rep(1, number_of_samples*number_of_snps), nrow=number_of_samples, ncol=number_of_snps)
        X <- cbind(X, X_zerovar)
    }
    ## Variance Threshold, just one row
    if (number_var_thresh > 0) {
        X_varthresh <- matrix(
            c(rep(1, (number_of_samples - 1)), 2), nrow=number_of_samples, ncol=1)
        X <- cbind(X, X_varthresh)
    }
    colnames(X) <- paste0(
        "chr1:",
        seq(1000,1000+number_of_snps+number_missing+number_low_maf+number_zero_var+number_var_thresh-1),
        "_G_C")
    rownames(X) <- paste0("Sample_", seq(sample_start_id, number_of_samples + sample_start_id - 1))
    return(X)
}

dummy_pheno_data <- function(number_of_samples = 10, number_of_phenotypes = 10, randomize = FALSE, y_as_matrix = F, sample_start_id = 1) {
    # Create dummy phenotype bed file
    # columns: Chrom, Start, End, Sample_1, Sample_2, ..., Sample_N
    start_matrix <- matrix(
        c(
            rep("chr1", number_of_phenotypes),
            seq(100, 100+number_of_phenotypes-1),
            seq(101, 101+number_of_phenotypes-1)
        ),
        nrow=number_of_phenotypes, ncol=3)
    end_matrix <- matrix(
        rnorm(number_of_samples*number_of_phenotypes), nrow=number_of_phenotypes, ncol=number_of_samples)
    pheno_data <- cbind(start_matrix, end_matrix)
    sample_ids <- paste0("Sample_", seq(sample_start_id, number_of_samples + sample_start_id - 1))
    colnames(pheno_data) <- c("#chr", "start", "end", sample_ids)
    colnames(end_matrix) <- sample_ids
    if (randomize) {
        end_matrix <- end_matrix[sample(nrow(end_matrix)),]
    }
    pheno_data <- t(pheno_data)
    pheno_data <- if (y_as_matrix) pheno_data else lapply(seq_len(ncol(pheno_data)), function(i) pheno_data[,i,drop=FALSE])
    return(pheno_data)
}

dummy_covar_data <- function(number_of_samples = 10, number_of_covars = 10, row_na = FALSE, randomize = FALSE, sample_start_id = 1) {
    covar <- matrix(
        sample(1:20, number_of_samples*number_of_covars, replace = TRUE),
        nrow=number_of_samples, ncol=number_of_covars)
    colnames(covar) <- paste0("Covar_", seq(1, number_of_covars))
    rownames(covar) <- paste0("Sample_", seq(sample_start_id, number_of_samples + sample_start_id - 1))
    if (randomize) {
        covar <- covar[sample(nrow(covar)),]
    }
    if (row_na) {
        covar[sample(length(covar),1), 1:number_of_covars] <- NA
    }
    return(covar)
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

test_that("Test load_genotype_data",{
  res <- load_genotype_data("test_data/protocol_example.genotype")
  sample_ids <- read_delim(
    "test_data/protocol_example.genotype.fam", delim = "\t", col_names = F
  ) %>% pull(X1)
  expect_equal(nrow(res), length(sample_ids))
  expect_equal(rownames(res), sample_ids)
})

test_that("Test load_genotype_region",{
  res <- load_genotype_region(
    "test_data/protocol_example.genotype")
  sample_ids <- read_delim(
    "test_data/protocol_example.genotype.fam", delim = "\t", col_names = F
  ) %>% pull(X1)
  expect_equal(nrow(res), length(sample_ids))
  expect_equal(rownames(res), sample_ids)
})

test_that("Test load_genotype_region no indels",{
  res <- load_genotype_region(
    "test_data/protocol_example.genotype", keep_indel = F)
  bim_file <- read_delim(
    "test_data/protocol_example.genotype.bim", delim = "\t", col_names = F
  )
  sample_ids <- read_delim(
    "test_data/protocol_example.genotype.fam", delim = "\t", col_names = F
  ) %>% pull(X1)
  expect_equal(nrow(res), length(sample_ids))
  expect_equal(rownames(res), sample_ids)
  indels <- with(bim_file, grepl("[^ATCG]", X5) | grepl("[^ATCG]", X6) | nchar(X5) > 1 | nchar(X6) > 1)
  expect_equal(
    nrow(bim_file[!indels, ]),
    ncol(res)
  )
})

test_that("Test load_genotype_region with region",{
  res <- load_genotype_region(
    "test_data/protocol_example.genotype",
    region = "chr22:20689453-20845958")
  sample_ids <- read_delim(
    "test_data/protocol_example.genotype.fam", delim = "\t", col_names = F
  ) %>% pull(X1)
  snp_ids <- read_delim(
    "test_data/protocol_example.genotype.bim", delim = "\t", col_names = F
  ) %>% pull(X2)
  expect_equal(nrow(res), length(sample_ids))
  expect_equal(rownames(res), sample_ids)
  expect_equal(ncol(res), 8)
  expect_equal(colnames(res), snp_ids[1:8])
})

test_that("Test load_genotype_region with region and no indels",{
  res <- load_genotype_region(
    "test_data/protocol_example.genotype",
    region = "chr22:20689453-20845958", keep_indel = F)
  bim_file <- read_delim(
    "test_data/protocol_example.genotype.bim", delim = "\t", col_names = F
  )[1:8, ]
  sample_ids <- read_delim(
    "test_data/protocol_example.genotype.fam", delim = "\t", col_names = F
  ) %>% pull(X1)
  expect_equal(nrow(res), length(sample_ids))
  expect_equal(rownames(res), sample_ids)
  indels <- with(bim_file, grepl("[^ATCG]", X5) | grepl("[^ATCG]", X6) | nchar(X5) > 1 | nchar(X6) > 1)
  expect_equal(
    nrow(bim_file[!indels, ]),
    ncol(res))
  expect_equal(colnames(res), bim_file[!indels, ]$X2)
})

test_that("Test load_genotype_region equals load_genotype_data",{
  res <- load_genotype_data("test_data/protocol_example.genotype")
  res_alt <- load_genotype_region(
    "test_data/protocol_example.genotype")
  expect_equal(dim(res), dim(res_alt))
  expect_equal(rownames(res), rownames(res_alt))
  expect_equal(colnames(res), colnames(res_alt))
  expect_equal(res_alt, res)
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

test_that("Test prepare_data_list multiple pheno",{
    # Create dummy data
    ## Prepare Genotype Data
    dummy_geno_data <- matrix(
        c(1,NA,NA,NA, 0,0,1,1, 2,2,2,2, 1,1,1,2, 2,2,0,1, 0,1,1,2),
        # Missing Rate, MAF thresh, Zero Var, Var Thresh, Regular values
        nrow=4, ncol=6)
    rownames(dummy_geno_data) <- c("Sample_1", "Sample_2", "Sample_3", "Sample_14")
    colnames(dummy_geno_data) <- c("chr1:122:G:C", "chr1:123:G:C", "chr1:124:G:C", "chr1:125:G:C", "chr1:126:G:C", "chr1:127:G:C")
    ## Prepare Phenotype Data
    dummy_pheno_data_one <- matrix(c("chr1", "222", "223", "1","1","2",NA), nrow=7, ncol=1)
    rownames(dummy_pheno_data_one) <- c("#chr", "start", "end", "Sample_3", "Sample_1", "Sample_2", "Sample_10")
    dummy_pheno_data_two <- matrix(c("chr1", "222", "223", "2","1","2",NA), nrow=7, ncol=1)
    rownames(dummy_pheno_data_two) <- c("#chr", "start", "end", "Sample_3", "Sample_1", "Sample_2", "Sample_10")
    ## Prepare Covariate Data
    dummy_covar_data <- matrix(c(70,71,72,73, 28,30,15,20, 1,2,3,4), nrow=4, ncol=3)
    rownames(dummy_covar_data) <- c("Sample_1", "Sample_2", "Sample_3", "Sample_4")
    colnames(dummy_covar_data) <- c("covar_1", "covar_2", "covar_3")
    # Set parameters
    imiss_cutoff <- 0.70
    maf_cutoff <- 0.1
    mac_cutoff <- 1.8
    xvar_cutoff <- 0.3
    keep_samples <- c("Sample_1", "Sample_2", "Sample_3")
    res <- prepare_data_list(
        dummy_geno_data, list(dummy_pheno_data_one, dummy_pheno_data_two), list(dummy_covar_data, dummy_covar_data),
        imiss_cutoff, maf_cutoff, mac_cutoff, xvar_cutoff, phenotype_header = 3, keep_samples=keep_samples)
    # Check that Covar, X, and Y have the same number of rows
    expect_equal(nrow(res$covar[[1]]), 3)
    expect_equal(nrow(res$X[[1]]), 3)
    expect_equal(length(res$Y[[1]]), 3)
    # Check that filter_X occured
    expect_equal(ncol(res$X[[1]]), 2)
    # Check that Covar, X, and Y have the same samples
    expect_equal(rownames(res$covar[[1]]), rownames(res$X[[1]]))
    expect_equal(rownames(res$covar[[1]]), names(res$Y[[1]]))
    expect_equal(rownames(res$X[[1]]), names(res$Y[[1]]))
})

test_that("Test prepare_data_list",{
    # Create dummy data
    ## Prepare Genotype Data
    dummy_geno_data <- matrix(
        c(1,NA,NA,NA, 0,0,1,1, 2,2,2,2, 1,1,1,2, 2,2,0,1, 0,1,1,2),
        # Missing Rate, MAF thresh, Zero Var, Var Thresh, Regular values
        nrow=4, ncol=6)
    rownames(dummy_geno_data) <- c("Sample_1", "Sample_2", "Sample_3", "Sample_14")
    colnames(dummy_geno_data) <- c("chr1:122:G:C", "chr1:123:G:C", "chr1:124:G:C", "chr1:125:G:C", "chr1:126:G:C", "chr1:127:G:C")
    ## Prepare Phenotype Data
    dummy_pheno_data <- matrix(
        c(
            rep("chr1", 4),
            rep(10, 4),
            rep(11, 4),
            1, NA, NA, NA,
            1, 1, 2, NA,
            2, 1, 2, NA
        ), ncol = 6, nrow = 4
    )
    rownames(dummy_pheno_data) <- c("Pheno_1", "Pheno_2", "Pheno_3", "Pheno_4")
    colnames(dummy_pheno_data) <- c("chrom", "start", "end", "Sample_1", "Sample_2", "Sample_3")
    dummy_pheno_data <- t(dummy_pheno_data)
    ## Prepare Covariate Data
    dummy_covar_data <- matrix(c(70,71,72,73, 28,30,15,20, 1,2,3,4), nrow=4, ncol=3)
    rownames(dummy_covar_data) <- c("Sample_1", "Sample_2", "Sample_3", "Sample_4")
    colnames(dummy_covar_data) <- c("covar_1", "covar_2", "covar_3")
    # Set parameters
    imiss_cutoff <- 0.70
    maf_cutoff <- 0.1
    mac_cutoff <- 1.8
    xvar_cutoff <- 0.3
    keep_samples <- c("Sample_1", "Sample_2", "Sample_3")
    res <- prepare_data_list(
        dummy_geno_data, list(dummy_pheno_data), list(dummy_covar_data), imiss_cutoff, maf_cutoff, mac_cutoff, xvar_cutoff,
        phenotype_header = 3, keep_samples=keep_samples)
    # Check that Covar, X, and Y have the same number of rows
    expect_equal(nrow(res$covar[[1]]), 3)
    expect_equal(nrow(res$X[[1]]), 3)
    expect_equal(nrow(res$Y[[1]]), 3)
    # Check that filter_X occured
    expect_equal(ncol(res$X[[1]]), 2)
    # Check that Covar, X, and Y have the same samples
    expect_equal(rownames(res$covar[[1]]), rownames(res$X[[1]]))
    expect_equal(rownames(res$covar[[1]]), rownames(res$Y[[1]]))
    expect_equal(rownames(res$X[[1]]), rownames(res$Y[[1]]))
})

test_that("Test prepare_X_matrix",{
    dummy_geno_data <- matrix(
        c(1,NA,NA,NA,2, 0,0,1,1,0, 2,2,2,2,2, 1,1,1,2,1, 2,2,0,1,2, 0,1,1,2,2),
        # Missing Rate, MAF thresh, Zero Var, Var Thresh, Regular values
        nrow=5, ncol=6)
    rownames(dummy_geno_data) <- c("Sample_1", "Sample_2", "Sample_3", "Sample_4", "Sample_5")
    colnames(dummy_geno_data) <- c("chr1:122:G:C", "chr1:123:G:C", "chr1:124:G:C", "chr1:125:G:C", "chr1:126:G:C", "chr1:127:G:C")
    dummy_covar_data <- matrix(
        c(70,71,72,73,74, 28,30,15,20,22, 1,2,3,4,5),
        nrow=5, ncol=3)
    rownames(dummy_covar_data) <- c("Sample_1", "Sample_2", "Sample_3", "Sample_4", "Sample_5")
    colnames(dummy_covar_data) <- c("covar_1", "covar_2", "covar_3")
    dummy_data_list <- tibble(
        covar = list(dummy_covar_data))
    # Set parameters
    imiss_cutoff <- 0.70
    maf_cutoff <- 0.3
    mac_cutoff <- 1.8
    xvar_cutoff <- 0.3
    res <- prepare_X_matrix(dummy_geno_data, dummy_data_list, imiss_cutoff, maf_cutoff, mac_cutoff, xvar_cutoff)
    target <- matrix(c(2,2,0,1,2, 0,1,1,2,2), nrow=5, ncol=2)
    rownames(target) <- c("Sample_1", "Sample_2", "Sample_3", "Sample_4", "Sample_5")
    colnames(target) <- c("chr1:126:G:C", "chr1:127:G:C")
    expect_equal(res, target)
})

test_that("Test add_X_residuals",{
    dummy_geno_data <- matrix(
        c(2,2,0,1, 0,1,1,2),
        nrow=4, ncol=2)
    dummy_covar_data <- matrix(
        c(70,71,72,73, 28,30,15,20, 1,2,3,4),
        nrow=4, ncol=3)
    dummy_data_list <- tibble(
        X = list(dummy_geno_data),
        covar = list(dummy_covar_data))
    res <- add_X_residuals(dummy_data_list)
    res_X <- .lm.fit(x = cbind(1, dummy_covar_data), y = dummy_geno_data)$residuals %>% as.matrix()
    res_X_mean <- apply(res_X, 2, mean)
    res_X_sd <- apply(res_X, 2, sd)
    expect_equal(res$lm_res_X[[1]], res_X)
    expect_equal(res$X_resid_mean[[1]], res_X_mean)
    expect_equal(res$X_resid_sd[[1]], res_X_sd)
})

# test_that("Test add_Y_residuals y_as_matrix=T",{
#     dummy_pheno_data <- matrix(
#         rnorm(8),
#         nrow=4, ncol=4)
#     dummy_covar_data <- matrix(
#         c(70,71,72,73, 28,30,15,20, 1,2,3,4),
#         nrow=4, ncol=3)
#     dummy_data_list <- tibble(
#         Y = list(dummy_pheno_data),
#         covar = list(dummy_covar_data))
#     conditions <- c("cond_1", "cond_2", "cond_3", "cond_4")
#     res_Y <- .lm.fit(x = cbind(1, dummy_covar_data), y = dummy_pheno_data)$residuals %>% as.matrix()
#     res_Y_mean <- apply(res_Y, 2, mean)
#     res_Y_sd <- apply(res_Y, 2, sd)
#     res <- add_Y_residuals(dummy_data_list, conditions, y_as_matrix=T)
#     expect_equal(res$lm_res[[1]], res_Y)
#     expect_equal(res$Y_resid_mean[[1]], res_Y_mean)
#     expect_equal(res$Y_resid_sd[[1]], res_Y_sd)
# })

test_that("Test add_Y_residuals y_as_matrix=F",{
    dummy_pheno_data <- rnorm(4)
    dummy_covar_data <- matrix(
        c(70,71,72,73, 28,30,15,20, 1,2,3,4),
        nrow=4, ncol=3)
    dummy_data_list <- tibble(
        Y = list(dummy_pheno_data),
        covar = list(dummy_covar_data))
    conditions <- c("cond_1")
    res_Y <- .lm.fit(x = cbind(1, dummy_covar_data), y = dummy_pheno_data)$residuals %>% as.matrix()
    res_Y_mean <- apply(res_Y, 2, mean)
    res_Y_sd <- apply(res_Y, 2, sd)
    res <- add_Y_residuals(dummy_data_list, conditions, y_as_matrix=F)
    expect_equal(res$lm_res[[1]], res_Y)
    expect_equal(res$Y_resid_mean[[1]], res_Y_mean)
    expect_equal(res$Y_resid_sd[[1]], res_Y_sd)
})

# For load_regional_association_data tests
## We mock the load data functionality for ease of use
test_that("Test load_regional_association_data complete overlap",{
    geno_data <- dummy_geno_data(
                number_of_samples = 10, number_of_snps = 10, sample_start_id = 1,
                number_missing = 0, number_low_maf = 0, number_zero_var = 0, number_var_thresh = 0)
    covar_data <- dummy_covar_data(
            number_of_samples = 10, number_of_covars = 5, row_na = FALSE, randomize = FALSE, sample_start_id = 1)
    pheno_data <- dummy_pheno_data(
            number_of_samples = 10, number_of_phenotypes = 1, randomize = FALSE, sample_start_id = 1)
    local_mocked_bindings(
        load_genotype_region = function(...) geno_data,
        load_covariate_data = function(...) list(covar_data),
        load_phenotype_data = function(...) pheno_data
    )
    res <- load_regional_association_data(
        "dummy_geno.bed.gz",
        "dummy_pheno.bed.gz",
        "dummy_covar.txt.gz",
        "chr1:1000-2000",
        "cond_1",
        imiss_cutoff = 0.70,
        maf_cutoff = 0.1,
        mac_cutoff = (0.1*10*2),
        xvar_cutoff = 0.2,
        y_as_matrix = F,
        phenotype_header = 3,
        keep_samples = NULL)
    expect_equal(nrow(res$X), 10)
    expect_equal(ncol(res$X), 10)
    colnames(geno_data) <- gsub("_", ":", colnames(geno_data))
    expect_equal(res$X[order(as.numeric(gsub("Sample_", "", rownames(res$X)))), , drop = FALSE], geno_data)
    expect_equal(length(res$Y[[1]]), 10)
    expect_equal(
        as.vector(res$Y[[1]][order(as.numeric(gsub("Sample_", "", names(res$Y[[1]]))))]),
        as.numeric(as.vector(asplit(pheno_data[[1]], 2)[[1]])[4:13]))
    expect_equal(nrow(res$covar[[1]]), 10)
    expect_equal(ncol(res$covar[[1]]), 5)
    expect_equal(res$covar[[1]][order(as.numeric(gsub("Sample_", "", rownames(res$covar[[1]])))), , drop = FALSE], covar_data)
})

# test_that("Test load_regional_association_data complete overlap y_matrix=T",{
#     geno_data <- dummy_geno_data(
#                 number_of_samples = 10, number_of_snps = 10, sample_start_id = 1,
#                 number_missing = 0, number_low_maf = 0, number_zero_var = 0, number_var_thresh = 0)
#     covar_data <- dummy_covar_data(
#             number_of_samples = 10, number_of_covars = 5, row_na = FALSE, randomize = FALSE, sample_start_id = 1)
#     pheno_data <- dummy_pheno_data(
#             number_of_samples = 10, number_of_phenotypes = 2, randomize = FALSE, y_as_matrix = T, sample_start_id = 1)
#     local_mocked_bindings(
#         load_genotype_data = function(...) geno_data,
#         load_covariate_data = function(...) list(covar_data),
#         load_phenotype_data = function(...) list(pheno_data)
#     )
#     res <- load_regional_association_data(
#         "dummy_geno.bed.gz",
#         "dummy_pheno.bed.gz",
#         "dummy_covar.txt.gz",
#         "chr1:1000-2000",
#         c("cond_1", "cond_2"),
#         imiss_cutoff = 0.70,
#         maf_cutoff = 0.1,
#         mac_cutoff = (0.1*10*2),
#         xvar_cutoff = 0.2,
#         y_as_matrix = T,
#         keep_samples = NULL)
#     expect_equal(nrow(res$X), 10)
#     expect_equal(ncol(res$X), 10)
#     expect_equal(res$X[order(as.numeric(gsub("Sample_", "", rownames(res$X)))), , drop = FALSE], geno_data)
#     expect_equal(nrow(res$Y[[1]]), 10)
#     expect_equal(ncol(res$Y[[1]]), 2)
#     expect_equal(
#         res$Y[[1]][order(as.numeric(gsub("Sample_", "", rownames(res$Y[[1]])))), , drop = FALSE],
#         pheno_data[4:13,])
#     expect_equal(nrow(res$covar[[1]]), 10)
#     expect_equal(ncol(res$covar[[1]]), 5)
#     expect_equal(res$covar[[1]][order(as.numeric(gsub("Sample_", "", rownames(res$covar[[1]])))), , drop = FALSE], covar_data)
# })

test_that("Test load_regional_association_data fewer covar samples",{
    geno_data <- dummy_geno_data(
                number_of_samples = 10, number_of_snps = 10, sample_start_id = 1,
                number_missing = 0, number_low_maf = 0, number_zero_var = 0, number_var_thresh = 0)
    covar_data <- dummy_covar_data(
            number_of_samples = 10, number_of_covars = 5, row_na = FALSE, randomize = FALSE, sample_start_id = 3)
    pheno_data <- dummy_pheno_data(
            number_of_samples = 10, number_of_phenotypes = 1, randomize = FALSE, y_as_matrix = F, sample_start_id = 1)
    local_mocked_bindings(
        load_genotype_region = function(...) geno_data,
        load_covariate_data = function(...) list(covar_data),
        load_phenotype_data = function(...) pheno_data
    )
    res <- load_regional_association_data(
        "dummy_geno.bed.gz",
        "dummy_pheno.bed.gz",
        "dummy_covar.txt.gz",
        "chr1:1000-2000",
        "cond_1",
        imiss_cutoff = 0.70,
        maf_cutoff = 0.1,
        mac_cutoff = (0.1*10*2),
        xvar_cutoff = 0.2,
        y_as_matrix = F,
        phenotype_header = 3,
        keep_samples = NULL)
    expect_equal(nrow(res$X), 8)
    expect_equal(ncol(res$X), 9)
    colnames(geno_data) <- gsub("_", ":", colnames(geno_data))
    expect_equal(
        res$X[order(as.numeric(gsub("Sample_", "", rownames(res$X)))), , drop = FALSE],
        geno_data[3:10,-6])
    expect_equal(length(res$Y[[1]]), 8)
    expect_equal(
        res$Y[[1]][order(as.numeric(gsub("Sample_", "", names(res$Y[[1]]))))],
        setNames(
            as.numeric(pheno_data[[1]][6:13,]),
            names(pheno_data[[1]][6:13,])))
    expect_equal(nrow(res$covar[[1]]), 8)
    expect_equal(ncol(res$covar[[1]]), 5)
    expect_equal(
        res$covar[[1]][order(as.numeric(gsub("Sample_", "", rownames(res$covar[[1]])))), , drop = FALSE],
        covar_data[1:8,])
})

test_that("Test load_regional_association_data slight overlap across geno, pheno, covar",{
    geno_data <- dummy_geno_data(
                number_of_samples = 10, number_of_snps = 10, sample_start_id = 1,
                number_missing = 0, number_low_maf = 0, number_zero_var = 0, number_var_thresh = 0)
    covar_data <- dummy_covar_data(
            number_of_samples = 10, number_of_covars = 5, row_na = FALSE, randomize = FALSE, sample_start_id = 3)
    pheno_data <- dummy_pheno_data(
            number_of_samples = 10, number_of_phenotypes = 1, randomize = FALSE, y_as_matrix = F, sample_start_id = 7)
    local_mocked_bindings(
        load_genotype_region = function(...) geno_data,
        load_covariate_data = function(...) list(covar_data),
        load_phenotype_data = function(...) pheno_data
    )
    res <- load_regional_association_data(
        "dummy_geno.bed.gz",
        "dummy_pheno.bed.gz",
        "dummy_covar.txt.gz",
        "chr1:1000-2000",
        "cond_1",
        imiss_cutoff = 0.70,
        maf_cutoff = 0.1,
        mac_cutoff = (0.1*10*2),
        xvar_cutoff = 0.2,
        y_as_matrix = F,
        phenotype_header = 3,
        keep_samples = NULL)
    expect_equal(nrow(res$X), 4)
    expect_equal(ncol(res$X), 3)
    colnames(geno_data) <- gsub("_", ":", colnames(geno_data))
    expect_equal(
        res$X[order(as.numeric(gsub("Sample_", "", rownames(res$X)))), , drop = FALSE],
        geno_data[7:10,c(2,4,7)])
    expect_equal(length(res$Y[[1]]), 4)
    expect_equal(
        res$Y[[1]][order(as.numeric(gsub("Sample_", "", names(res$Y[[1]]))))],
        setNames(
            as.numeric(pheno_data[[1]][4:7,]),
            names(pheno_data[[1]][4:7,])))
    expect_equal(nrow(res$covar[[1]]), 4)
    expect_equal(ncol(res$covar[[1]]), 5)
    expect_equal(
        res$covar[[1]][order(as.numeric(gsub("Sample_", "", rownames(res$covar[[1]])))), , drop = FALSE],
        covar_data[5:8,])
})

test_that("Test load_regional_association_data no overlap",{
    geno_data <- dummy_geno_data(
                number_of_samples = 10, number_of_snps = 10, sample_start_id = 1,
                number_missing = 0, number_low_maf = 0, number_zero_var = 0, number_var_thresh = 0)
    covar_data <- dummy_covar_data(
            number_of_samples = 10, number_of_covars = 5, row_na = FALSE, randomize = FALSE, sample_start_id = 11)
    pheno_data <- dummy_pheno_data(
            number_of_samples = 10, number_of_phenotypes = 1, randomize = FALSE, y_as_matrix = F, sample_start_id = 21)
    local_mocked_bindings(
        load_genotype_region = function(...) geno_data,
        load_covariate_data = function(...) list(covar_data),
        load_phenotype_data = function(...) pheno_data
    )
    expect_error(
        load_regional_association_data(
            "dummy_geno.bed.gz",
            "dummy_pheno.bed.gz",
            "dummy_covar.txt.gz",
            "chr1:1000-2000",
            "cond_1",
            imiss_cutoff = 0.70,
            maf_cutoff = 0.1,
            mac_cutoff = (0.1*10*2),
            xvar_cutoff = 0.2,
            y_as_matrix = F,
            phenotype_header = 3,
            keep_samples = NULL),
        "No common complete samples between genotype and phenotype/covariate data")
})

test_that("Test load_regional_association_data unordered samples",{
    geno_data <- dummy_geno_data(
                number_of_samples = 10, number_of_snps = 10, sample_start_id = 1,
                number_missing = 0, number_low_maf = 0, number_zero_var = 0, number_var_thresh = 0)
    covar_data <- dummy_covar_data(
            number_of_samples = 10, number_of_covars = 5, row_na = FALSE, randomize = TRUE, sample_start_id = 1)
    pheno_data <- dummy_pheno_data(
            number_of_samples = 10, number_of_phenotypes = 1, randomize = TRUE, y_as_matrix = F, sample_start_id = 1)
    local_mocked_bindings(
        load_genotype_region = function(...) geno_data,
        load_covariate_data = function(...) list(covar_data),
        load_phenotype_data = function(...) pheno_data
    )
    res <- load_regional_association_data(
        "dummy_geno.bed.gz",
        "dummy_pheno.bed.gz",
        "dummy_covar.txt.gz",
        "chr1:1000-2000",
        c("cond_1"),
        imiss_cutoff = 0.70,
        maf_cutoff = 0.1,
        mac_cutoff = (0.1*10*2),
        xvar_cutoff = 0.2,
        y_as_matrix = F,
        phenotype_header = 3,
        keep_samples = NULL)
    expect_equal(nrow(res$X), 10)
    expect_equal(ncol(res$X), 10)
    colnames(geno_data) <- gsub("_", ":", colnames(geno_data))
    expect_equal(res$X[order(as.numeric(gsub("Sample_", "", rownames(res$X)))), , drop = FALSE], geno_data)
    expect_equal(length(res$Y[[1]]), 10)
    expect_equal(
        res$Y[[1]][order(as.numeric(gsub("Sample_", "", names(res$Y[[1]]))))],
        setNames(
            as.numeric(pheno_data[[1]][4:13,]),
            names(pheno_data[[1]][4:13,])))
    expect_equal(nrow(res$covar[[1]]), 10)
    expect_equal(ncol(res$covar[[1]]), 5)
    expect_equal(
        res$covar[[1]][order(as.numeric(gsub("Sample_", "", rownames(res$covar[[1]])))), , drop = FALSE],
        covar_data[order(as.numeric(gsub("Sample_", "", rownames(covar_data)))), , drop = FALSE])
})

test_that("Test load_regional_univariate_data",{
    geno_data <- dummy_geno_data(
                number_of_samples = 10, number_of_snps = 10, sample_start_id = 1,
                number_missing = 0, number_low_maf = 0, number_zero_var = 0, number_var_thresh = 0)
    covar_data <- dummy_covar_data(
            number_of_samples = 10, number_of_covars = 5, row_na = FALSE, randomize = FALSE, sample_start_id = 1)
    pheno_data <- dummy_pheno_data(
            number_of_samples = 10, number_of_phenotypes = 1, randomize = FALSE, sample_start_id = 1)
    local_mocked_bindings(
        load_genotype_region = function(...) geno_data,
        load_covariate_data = function(...) list(covar_data),
        load_phenotype_data = function(...) pheno_data
    )
    res <- load_regional_univariate_data(
        "dummy_geno.bed.gz",
        "dummy_pheno.bed.gz",
        "dummy_covar.txt.gz",
        "chr1:1000-2000",
        "cond_1",
        imiss_cutoff = 0.70,
        maf_cutoff = 0.1,
        mac_cutoff = (0.1*10*2),
        xvar_cutoff = 0.2,
        phenotype_header = 3,
        keep_samples = NULL)
    expect_equal(nrow(res$X), 10)
    expect_equal(ncol(res$X), 10)
    colnames(geno_data) <- gsub("_", ":", colnames(geno_data))
    expect_equal(res$X[order(as.numeric(gsub("Sample_", "", rownames(res$X)))), , drop = FALSE], geno_data)
    expect_equal(length(res$residual_Y$cond_1), 10)
})

test_that("Test load_regional_regression_data",{
    geno_data <- dummy_geno_data(
                number_of_samples = 10, number_of_snps = 10, sample_start_id = 1,
                number_missing = 0, number_low_maf = 0, number_zero_var = 0, number_var_thresh = 0)
    covar_data <- dummy_covar_data(
            number_of_samples = 10, number_of_covars = 5, row_na = FALSE, randomize = FALSE, sample_start_id = 1)
    pheno_data <- dummy_pheno_data(
            number_of_samples = 10, number_of_phenotypes = 1, randomize = FALSE, y_as_matrix = F, sample_start_id = 1)
    local_mocked_bindings(
        load_genotype_region = function(...) geno_data,
        load_covariate_data = function(...) list(covar_data),
        load_phenotype_data = function(...) pheno_data
    )
    res <- load_regional_regression_data(
        "dummy_geno.bed.gz",
        "dummy_pheno.bed.gz",
        "dummy_covar.txt.gz",
        "chr1:1000-2000",
        "cond_1",
        imiss_cutoff = 0.70,
        maf_cutoff = 0.1,
        mac_cutoff = (0.1*10*2),
        xvar_cutoff = 0.2,
        y_as_matrix = F,
        phenotype_header = 3,
        keep_samples = NULL)
    expect_equal(nrow(res$X_data[[1]]), 10)
    expect_equal(ncol(res$X_data[[1]]), 10)
    colnames(geno_data) <- gsub("_", ":", colnames(geno_data))
    expect_equal(res$X_data[[1]][order(as.numeric(gsub("Sample_", "", rownames(res$X_data[[1]])))), , drop = FALSE], geno_data)
    expect_equal(length(res$Y[[1]]), 10)
    expect_equal(
        res$Y[[1]][order(as.numeric(gsub("Sample_", "", names(res$Y[[1]]))))],
        setNames(
            as.numeric(pheno_data[[1]][4:13,]),
            names(pheno_data[[1]][4:13,])))
    expect_equal(nrow(res$covar[[1]]), 10)
    expect_equal(ncol(res$covar[[1]]), 5)
    expect_equal(res$covar[[1]][order(as.numeric(gsub("Sample_", "", rownames(res$covar[[1]])))), , drop = FALSE], covar_data)
})

#test_that("Test load_regional_multivariate_data",{
#  set.seed(1)
#})