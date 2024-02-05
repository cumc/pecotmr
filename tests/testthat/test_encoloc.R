context("encoloc")
library(tidyverse)

generate_mock_data_for_enrichment <- function(seed=1, num_files=2) {
    gwas_finemapped_data <- as.vector(paste0("gwas_file", 1:num_files))
    xqtl_finemapped_data <- "xqtl_file.rds"
    return(list(gwas_finemapped_data = gwas_finemapped_data,
                xqtl_finemapped_data = xqtl_finemapped_data))
}

generate_mock_susie_fit <- function(seed=1, num_samples = 10, num_features=10) {
    set.seed(seed)
    start_index <- sample(1:100000, 1)
    susie_fit <- list(
        pip = setNames(runif(num_features), paste0("rs", start_index:(num_features+start_index-1))),
        variant_names = paste0("chr22:", start_index:(num_features + start_index -1), ":A:C"),
        lbf_variable = matrix(runif(num_features * 10), nrow = 10, ncol=num_features),
        alpha = t(matrix(runif(num_samples * num_features), nrow = num_samples)),
        V = runif(10), 
        prior_variance = runif(num_features))
    return(susie_fit)
}

test_that("xqtl_enrichment_wrapper works with dummy input single threaded",{
    local_mocked_bindings(
        qtl_enrichment_rcpp = function(...) TRUE)
    input_data <- generate_mock_data_for_enrichment()
    input_data$gwas_finemapped_data <- unlist(lapply(
        input_data$gwas_finemapped_data, function(x) {
            gsub("//", "/", tempfile(pattern = x, tmpdir = tempdir(), fileext = ".rds"))
    }))
    input_data$xqtl_finemapped_data <- gsub("//", "/", tempfile(pattern = "xqtl_file", tmpdir = tempdir(), fileext = ".rds"))
    saveRDS(list(gene=list(susie_fit = generate_mock_susie_fit(seed=1))), input_data$xqtl_finemapped_data)
    for (i in 1:length(input_data$gwas_finemapped_data)) {
        saveRDS(list(susie_fit = generate_mock_susie_fit(seed=i)), input_data$gwas_finemapped_data[i])
    }
    res <- xqtl_enrichment_wrapper(
        input_data$xqtl_finemapped_data,input_data$gwas_finemapped_data, 
        gwas_finemapping_obj = "susie_fit",
        xqtl_finemapping_obj = "susie_fit",
        pi_gwas = 0.5, pi_qtl = 0.5, 
        lambda = 1.0, ImpN = 25,
        num_threads = 1)
    expect_true(res)
    file.remove(input_data$gwas_finemapped_data)
    file.remove(input_data$xqtl_finemapped_data)
})

test_that("xqtl_enrichment_wrapper works with dummy input multi threaded",{
    local_mocked_bindings(
        qtl_enrichment_rcpp = function(...) TRUE)
    input_data <- generate_mock_data_for_enrichment()
    input_data$gwas_finemapped_data <- unlist(lapply(
        input_data$gwas_finemapped_data, function(x) {
            gsub("//", "/", tempfile(pattern = x, tmpdir = tempdir(), fileext = ".rds"))
    }))
    input_data$xqtl_finemapped_data <- gsub("//", "/", tempfile(pattern = "xqtl_file", tmpdir = tempdir(), fileext = ".rds"))
    saveRDS(list(gene=list(susie_fit = generate_mock_susie_fit(seed=1))), input_data$xqtl_finemapped_data)
    for (i in 1:length(input_data$gwas_finemapped_data)) {
        saveRDS(list(susie_fit = generate_mock_susie_fit(seed=i)), input_data$gwas_finemapped_data[i])
    }
    res <- xqtl_enrichment_wrapper(
        input_data$xqtl_finemapped_data,input_data$gwas_finemapped_data, 
        gwas_finemapping_obj = "susie_fit",
        xqtl_finemapping_obj = "susie_fit",
        pi_gwas = 0.5, pi_qtl = 0.5, 
        lambda = 1.0, ImpN = 25,
        num_threads = 2)
    expect_true(res)
    file.remove(input_data$gwas_finemapped_data)
    file.remove(input_data$xqtl_finemapped_data)
})

test_that("xqtl_enrichment_wrapper works with dummy input single and multi threaded",{
    local_mocked_bindings(
        qtl_enrichment_rcpp = function(...) TRUE)
    input_data <- generate_mock_data_for_enrichment()
    input_data$gwas_finemapped_data <- unlist(lapply(
        input_data$gwas_finemapped_data, function(x) {
            gsub("//", "/", tempfile(pattern = x, tmpdir = tempdir(), fileext = ".rds"))
    }))
    input_data$xqtl_finemapped_data <- gsub("//", "/", tempfile(pattern = "xqtl_file", tmpdir = tempdir(), fileext = ".rds"))
    saveRDS(list(gene=list(susie_fit = generate_mock_susie_fit(seed=1))), input_data$xqtl_finemapped_data)
    for (i in 1:length(input_data$gwas_finemapped_data)) {
        saveRDS(list(susie_fit = generate_mock_susie_fit(seed=i)), input_data$gwas_finemapped_data[i])
    }
    res_single <- xqtl_enrichment_wrapper(
        input_data$xqtl_finemapped_data,input_data$gwas_finemapped_data, 
        gwas_finemapping_obj = "susie_fit",
        xqtl_finemapping_obj = "susie_fit",
        pi_gwas = 0.5, pi_qtl = 0.5, 
        lambda = 1.0, ImpN = 25,
        num_threads = 1)
    res_multi <- xqtl_enrichment_wrapper(
        input_data$xqtl_finemapped_data,input_data$gwas_finemapped_data, 
        gwas_finemapping_obj = "susie_fit",
        xqtl_finemapping_obj = "susie_fit",
        pi_gwas = 0.5, pi_qtl = 0.5, 
        lambda = 1.0, ImpN = 25,
        num_threads = 2)
    expect_equal(res_single, res_multi)
    file.remove(input_data$gwas_finemapped_data)
    file.remove(input_data$xqtl_finemapped_data)
})

test_that("coloc_wrapper works with dummy input",{
    input_data <- generate_mock_data_for_enrichment()
    input_data$gwas_finemapped_data <- unlist(lapply(
        input_data$gwas_finemapped_data, function(x) {
            gsub("//", "/", tempfile(pattern = x, tmpdir = tempdir(), fileext = ".rds"))
    }))
    input_data$xqtl_finemapped_data <- gsub("//", "/", tempfile(pattern = "xqtl_file", tmpdir = tempdir(), fileext = ".rds"))
    saveRDS(list(gene=list(susie_fit = generate_mock_susie_fit(seed=1))), input_data$xqtl_finemapped_data)
    for (i in 1:length(input_data$gwas_finemapped_data)) {
        saveRDS(list(susie_fit = generate_mock_susie_fit(seed=i)), input_data$gwas_finemapped_data[i])
    }
    res <- coloc_wrapper(input_data$xqtl_finemapped_data, input_data$gwas_finemapped_data, 
                     xqtl_finemapping_obj =  "susie_fit", gwas_finemapping_obj=  "susie_fit", 
                     xqtl_varname_obj= c("susie_fit","variant_names"), gwas_varname_obj = c("susie_fit","variant_names"))
    expect_true(all(names(res) %in% c("summary","results","priors","analysis_region")))
    file.remove(input_data$gwas_finemapped_data)
    file.remove(input_data$xqtl_finemapped_data)
})

test_that("coloc_wrapper error with non-unique",{ })
