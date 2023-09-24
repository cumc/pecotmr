library(foreach)
library(doParallel)

# Function to generate fake input data
generate_fake_data <- function() {
  gwas_pip <- seq(1, 100)
  eqtl_pip <- list()
  for(i in 1:20) {
    eqtl_pip[[i]] <- list(gwas_pip = sample(gwas_pip, 10))
  }
  return(list(eqtl_pip = eqtl_pip, gwas_pip = gwas_pip))
}

# Fake function to simulate impute_qtn
impute_qtn <- function(eqtl) {
  return(sample(0:(length(eqtl$gwas_pip) - 1), 1))
}

# Fake function to simulate run_EM
run_EM <- function(eqtl_sample) {
  return(runif(4))
}

# Function for enrichment estimation
# gwas_pip is a vector of GWAS PIP from SuSiE analysis for as many SNPs as available.
# eqtl_pip is a list of genes level SuSiE, each has the alpha matrix, pi0 estimates and pip vector.
# Here, the 
# NOTE: difference from previous implementation is that here we estimate number of signals not using CS filter but using sigma_0 filter for effect sizes
enrich_est <- function(eqtl_pip, gwas_pip, pi_gwas=NULL, pi_eqtl=NULL, ImpN, ncores=8, prior_variance=1) {
  
  registerDoParallel(cores = ncores)
  snp_index <- 1:length(gwas_pip)
  all_eqtl_pip <- do.call(c, lapply(eqtl_pip, function(x) x$pip))

  if (is.null(pi_gwas)) {
    pi_gwas <- sum(gwas_pip) / length(gwas_pip)
  }
  if (is.null(pi_eqtl)) {
    pi_eqtl <- lapply(1:length(eqtl_pip), function(i) sum(i)/length(i))
  }
  pi_eqtl <- 

  # FIXME: quit on error if pi_gwas or pi_eqtl is less than 1

  results <- foreach(k = 1:ImpN, .combine = 'rbind') %dopar% {
    eqtl_sample <- integer(length(gwas_pip))
    for (i in 1:length(eqtl_pip)) {
      rst <- impute_qtn(eqtl_pip[[i]])
      if (rst >= 0) {
        snp <- eqtl_pip[[i]]$gwas_pip[rst + 1]
        eqtl_sample[snp_index[snp]] <- 1
      }
    }
    rst <- run_EM(eqtl_sample, gwas_pip)
    return(rst)
  }

  a0_vec <- results[, 1]
  a1_vec <- results[, 2]
  v0_vec <- results[, 3]
  v1_vec <- results[, 4]

  a0_est <- mean(a0_vec)
  a1_est <- mean(a1_vec)
  var0 <- mean(v0_vec)
  var1 <- mean(v1_vec)

  bv0 <- var(a0_vec) * ImpN / (ImpN - 1)
  bv1 <- var(a1_vec) * ImpN / (ImpN - 1)

  sd0 <- sqrt(var0 + bv0)
  sd1 <- sqrt(var1 + bv1)

  a1_est_ns <- a1_est
  sd1_ns <- sd1
  
  if (prior_variance > 0) {
    post_var <- 1 / (1 / prior_variance + 1 / (sd1 * sd1))
    a1_est <- (a1_est_ns * prior_variance) / (prior_variance + sd1_ns * sd1_ns)
    sd1 <- sqrt(post_var)
  }

  a0_est <- log(P_gwas / (1 + P_eqtl * exp(a1_est) - P_eqtl - P_gwas))

  p1 <- (1-P_eqtl) * exp(a0_est) / (1 + exp(a0_est))
  p2 <- P_eqtl / (1 + exp(a0_est + a1_est))
  p12 <- P_eqtl * exp(a0_est + a1_est) / (1 + exp(a0_est + a1_est))

  return(list(a0_est = a0_est, a1_est = a1_est, a1_est_ns = a1_est_ns, sd1_ns = sd1_ns, sd1 = sd1,
              p1 = p1, p2 = p2, p12 = p12, bv = c(bv0, bv1)))
}

# Generate fake data
fake_data <- generate_fake_data()
eqtl_pip <- fake_data$eqtl_pip
gwas_pip <- fake_data$gwas_pip
ImpN <- 10  # Number of imputations

# Run the enrichment estimation function
result <- enrich_est(eqtl_pip, gwas_pip, ImpN)
print(result)


gwas_pip = vector(1:5)
names(gwas_pip) = c("a","b","c","d","e")
eqtl_pip[[1]]$pip = vector(1:5)
names(eqtl_pip[[1]]$pip) = c("a","b","c","d","e")
eqtl_pip[[1]]$alpha = matrix(1:50,5,10)