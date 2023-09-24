# Generate some simulated data
set.seed(123)


# Simulate data for gwas_pip
n_gwas_pip <- 1000
gwas_pip <- runif(n_gwas_pip)
names(gwas_pip) <- paste0("snp", 1:n_gwas_pip)

# Simulate data for a single SuSiEFit object
simulate_susiefit <- function(n, p) {
  pip = runif(n)
  names(pip) = paste0("snp", 1:n)
  alpha = t(matrix(runif(n * p), nrow = n))
  alpha <- t(apply(alpha, 1, function(row) row / sum(row)))

  list(
    pip = pip,
    alpha = alpha,
    prior_variance = runif(p)
  )
}

# Simulate multiple SuSiEFit objects
n_susie_fits <- 10
susie_fits <- replicate(n_susie_fits, simulate_susiefit(n_gwas_pip, 2), simplify = FALSE)

# Add these fits to a list, providing names to each element
names(susie_fits) <- paste0("fit", 1:length(susie_fits))


# Set other parameters
pi_gwas <- 0.01
pi_qtl <- 0.01
ImpN <- 10
prior_variance <- 1
num_threads <- 1


# Now susie_fits is ready to be passed to C++

# Load the Rcpp library
library(Rcpp)
# Compile the C++ code into an R-accessible shared object
sourceCpp("src/qtl_enrichment.cpp")

# Call the function
output <- qtl_enrichment(r_gwas_pip = gwas_pip, 
                         r_qtl_susie_fit = susie_fits,
                         pi_gwas = pi_gwas,
                         pi_qtl = pi_qtl,
                         ImpN = ImpN,
                         prior_variance = prior_variance,
                         num_threads = num_threads)
# Inspect the output
print(output)