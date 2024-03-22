context("compute_qtl_enrichment")

generate_mock_data <- function(seed=1, num_pips = 1000, num_susie_fits = 2) {
  # Simulate fake data for gwas_pip
  n_gwas_pip <- num_pips
  gwas_pip <- runif(n_gwas_pip)
  names(gwas_pip) <- paste0("snp", 1:n_gwas_pip)
  gwas_fit <- list(pip=gwas_pip)

  # Simulate fake data for a single SuSiEFit object
  simulate_susiefit <- function(n, p) {
    pip <- runif(n)
    names(pip) <- paste0("snp", 1:n)
    alpha <- t(matrix(runif(n * p), nrow = n))
    alpha <- t(apply(alpha, 1, function(row) row / sum(row)))
    list(
      pip = pip,
      alpha = alpha,
      prior_variance = runif(p)
    )
  }
  
  # Simulate multiple SuSiEFit objects
  n_susie_fits <- num_susie_fits
  susie_fits <- replicate(n_susie_fits, simulate_susiefit(n_gwas_pip, 10), simplify = FALSE)
  # Add these fits to a list, providing names to each element
  names(susie_fits) <- paste0("fit", 1:length(susie_fits))
  return(list(gwas_fit=gwas_fit, susie_fits=susie_fits))
}

test_that("compute_qtl_enrichment dummy data single-threaded works",{
  local_mocked_bindings(
      qtl_enrichment_rcpp = function(...) TRUE)
  input_data <- generate_mock_data(seed=1, num_pips=10)
  expect_warning(
    compute_qtl_enrichment(input_data$gwas_fit$pip, input_data$susie_fits, lambda = 1, ImpN = 10, num_threads = 1),
    "pi_gwas is not provided. Estimating pi_gwas from the data. Note that this estimate may be biased if the input gwas_pip does not contain genome-wide variants.")
  expect_warning(
    compute_qtl_enrichment(input_data$gwas_fit$pip, input_data$susie_fits, lambda = 1, ImpN = 10, num_threads = 1),
    "pi_qtl is not provided. Estimating pi_qtl from the data. Note that this estimate may be biased if either 1) the input susie_qtl_regions does not have enough data, or 2) the single effects only include variables inside of credible sets or signal clusters.")
  res <- compute_qtl_enrichment(input_data$gwas_fit$pip, input_data$susie_fits, pi_gwas=0.5141, pi_qtl=0.49819, lambda = 1, ImpN = 10, num_threads = 1)
  expect_true(length(res) > 0)
})

test_that("compute_qtl_enrichment dummy data single thread and multi-threaded are equivalent",{
  local_mocked_bindings(
      qtl_enrichment_rcpp = function(...) TRUE)
  input_data <- generate_mock_data(seed=1, num_pips=10)
  res_single <- compute_qtl_enrichment(input_data$gwas_fit$pip, input_data$susie_fits, pi_gwas=0.5141, pi_qtl=0.49819, lambda = 1, ImpN = 10, num_threads = 1)
  res_multi <- compute_qtl_enrichment(input_data$gwas_fit$pip, input_data$susie_fits, pi_gwas=0.5141, pi_qtl=0.49819, lambda = 1, ImpN = 10, num_threads = 2)
  expect_equal(res_single, res_multi)
})