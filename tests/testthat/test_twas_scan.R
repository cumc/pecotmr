context("twas_scan")

generate_mock_data <- function(seed=1, num_snps=100, empty_sets = F) {
  set.seed(seed)
  weights_path <- data.frame(ID = "random_Gene", path = "a/b/c/")

  # Mock data for 'region'
  chroms <- rep(22, num_snps)
  starts <- sample(1e6:5e6, num_snps, replace = TRUE)
  ends <- starts + sample(1e4:1e5, num_snps, replace = TRUE)
  region <- data.frame(chrom = chroms, start = starts, end = ends, ID = "random_Gene")

  # Mock data for 'GWAS_data'
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
  GWAS_data <- data.frame(
    chr = rep(22, num_snps),
    pos = sample(1e6:5e6, num_snps, replace = TRUE),
    A1.sumstats = sumstat_A1,
    A2.sumstats = sumstat_A2,
    beta = rnorm(num_snps),
    se = runif(num_snps, 0.01, 0.1),
    z = rnorm(num_snps)
  ) %>% arrange(chr, pos)

  # Mock data for 'LD_meta_file'
  LD_matrix <- matrix(runif(num_snps^2), nrow = num_snps, ncol = num_snps)
  colnames(LD_matrix) <- rownames(LD_matrix) <- variants_id_all <-paste(
    GWAS_data$chr, GWAS_data$pos, GWAS_data$A1.sumstats, GWAS_data$A2.sumstats, sep = ":")
  ld_res <- list(
    LD = LD_matrix,
    variants_id_all = variants_id_all)

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

  # Mock data for QTL weights
  qtl_weights <- list(list(
      sets = if (empty_sets) NULL else list(
        set1 = "test",
        set2 = "test",
        set3 = "test"
      ),
      variant_names = variants_id_all,
      susie_weights = random_weights(num_true_effects = 10, num_snps = num_snps, seed = 1),
      lasso_weights = random_weights(num_true_effects = 10, num_snps = num_snps, seed = 1),
      enet_weights = random_weights(num_true_effects = 10, num_snps = num_snps, seed = 1),
      mrash_weights = random_weights(num_true_effects = 10, num_snps = num_snps, seed = 1)
  ))
  return(list(gwas = GWAS_data, LD = ld_res, region = region, qtl_weights = qtl_weights, weights_path = weights_path))
}

test_that("Confirm twas_scan works with simulated data",{
  data <- generate_mock_data()
  weight_path <- gsub("//", "/", tempfile(pattern = "qtl_weights", tmpdir = tempdir(), fileext = ".RDS"))
  saveRDS(data$qtl_weights, file=weight_path)
  data$weights_path$path <- weight_path
  local_mocked_bindings(
      load_LD_matrix = function(...) data$LD,
      allele_qc = function(...) data$gwas,
  )
  res <- twas_scan(data$weights_path, data$region, data$gwas, data$LD)
  expect_equal(nrow(res$gene_weights_pq), nrow(data$region))
  expect_equal(
    res$outcome_QC,
    data$gwas %>% 
      mutate(variant_allele_flip = paste(chr, pos, A1.sumstats, A2.sumstats, sep = ":"))
  )
  expect_equal(nrow(res$twas_z_format), nrow(data$region))
  file.remove(weight_path)
})


test_that("Confirm twas_scan works with null qtl_weights",{
  data <- generate_mock_data(empty_sets = T)
  weight_path <- gsub("//", "/", tempfile(pattern = "qtl_weights", tmpdir = tempdir(), fileext = ".RDS"))
  saveRDS(data$qtl_weights, file=weight_path)
  data$weights_path$path <- weight_path
  local_mocked_bindings(
      load_LD_matrix = function(...) data$LD,
      allele_qc = function(...) data$gwas,
  )
  expect_output(
    twas_scan(data$weights_path, data$region, data$gwas, data$LD),
    "The 'qtl_weights' is NULL, so no output is generated.",
    fixed = TRUE)
  file.remove(weight_path)
})