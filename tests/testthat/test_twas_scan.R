context("twas_scan")

generate_mock_data <- function(seed=1, num_snps=100, empty_sets = F, gwas_mismatch = F, LD_mismatch = F) {
  set.seed(seed)

  random_weights <- function(num_true_effects = 10, num_snps = num_snps, seed = 1) {
    # Generate true effect sizes for a subset of variants
    true_effects = rnorm(num_true_effects, mean = 0, sd = 1)  # assuming a normal distribution

    # Set zero effects for the remaining variants
    zero_effects = rep(0, num_snps - num_true_effects)

    # Combine and shuffle the effect sizes
    set.seed(seed)
    effect_sizes = sample(c(true_effects, zero_effects))
    noise = rnorm(num_snps, mean = 0, sd = 0.1)
    effect_sizes = effect_sizes + noise
    return(effect_sizes)
  }

  weights_all_matrix <- matrix(
    c(random_weights(num_true_effects = 10, num_snps = num_snps, seed = 1),
    random_weights(num_true_effects = 10, num_snps = num_snps, seed = 2),
    random_weights(num_true_effects = 10, num_snps = num_snps, seed = 3),
    random_weights(num_true_effects = 10, num_snps = num_snps, seed = 4)),
    ncol = 4, nrow=100)
  
  sumstat_A1 <- sample(c("A", "T", "G", "C"), num_snps, replace = TRUE)
  sumstat_A2 <- unlist(lapply(sumstat_A1, function(x) {
    if (x == "A") {
      return(sample(c("G", "C"), 1))
    } else if (x == "T") {
      return(sample(c("G", "C"), 1))
    } else if (x == "G") {
      return(sample(c("A", "T"), 1))
    } else if (x == "C") {
      return(sample(c("A", "T"), 1))
    }
  }))
  gwas_sumstats_db <- data.frame(
    chr = rep(22, num_snps),
    pos = sample(1e6:5e6, num_snps, replace = TRUE),
    A1.sumstats = sumstat_A1,
    A2.sumstats = sumstat_A2,
    beta = rnorm(num_snps),
    se = runif(num_snps, 0.01, 0.1),
    z = rnorm(num_snps)
  ) %>% arrange(chr, pos)
  
  LD_matrix <- matrix(runif(num_snps^2), nrow = num_snps, ncol = num_snps)
  colnames(LD_matrix) <- rownames(LD_matrix) <- variants_id_all <- paste(
    gwas_sumstats_db$chr, gwas_sumstats_db$pos, gwas_sumstats_db$A1.sumstats, gwas_sumstats_db$A2.sumstats, sep = ":")

  extract_variants_objs <- variants_id_all

  gwas_sumstats_db$variant_allele_flip <- if (!gwas_mismatch) variants_id_all else paste(
    12, gwas_sumstats_db$pos, gwas_sumstats_db$A1.sumstats, gwas_sumstats_db$A2.sumstats, sep = ":")

  colnames(LD_matrix) <- rownames(LD_matrix) <- if (!LD_mismatch) colnames(LD_matrix) else paste(
    3, gwas_sumstats_db$pos, gwas_sumstats_db$A1.sumstats, gwas_sumstats_db$A2.sumstats, sep = ":")

  return(list(
    weights_all_matrix = weights_all_matrix,
    gwas_sumstats_db = gwas_sumstats_db,
    LD_matrix = LD_matrix,
    extract_variants_objs = extract_variants_objs))
}

test_that("Confirm twas_scan works with simulated data",{
  data <- generate_mock_data()
  res <- twas_analysis(data$weights_all_matrix, data$gwas_sumstats_db, data$LD_matrix, data$extract_variants_objs)
  expect_equal(length(res), 4)
  expect_true(all(unlist(lapply(res, function(x) {"z" %in% names(x)}))))
  expect_true(all(unlist(lapply(res, function(x) {"pval" %in% names(x)}))))
})

test_that("twas_analysis raises error if empty gwas",{
  data <- generate_mock_data(gwas_mismatch = T, LD_mismatch = F)
  expect_error(
    twas_analysis(data$weights_all_matrix, data$gwas_sumstats_db, data$LD_matrix, data$extract_variants_objs),
    "No GWAS summary statistics found for the specified variants.")
})

test_that("twas_analysis raises error if empty LD_matrix",{
  data <- generate_mock_data(gwas_mismatch = F, LD_mismatch = T)
  expect_error(
    twas_analysis(data$weights_all_matrix, data$gwas_sumstats_db, data$LD_matrix, data$extract_variants_objs),
    "LD matrix subset extraction failed.")
})