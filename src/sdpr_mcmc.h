#include <iostream>
#include <vector>
#include <armadillo.hpp>
#include <cmath>
#include <random>
#include <unordered_map>

class beta_distribution {
public:
// Constructor that takes alpha and beta parameters
beta_distribution(double alpha, double beta)
	: a_distribution(alpha, 1.0), b_distribution(beta, 1.0) {
}

// Generates a beta-distributed random number
template <typename RNG>
double operator()(RNG& rng) {
	double x = a_distribution(rng);
	double y = b_distribution(rng);
	return x / (x + y);
}

private:
std::gamma_distribution<double> a_distribution;
std::gamma_distribution<double> b_distribution;
};

/**
 * @class mcmc_data
 * @brief The per-SNP data the MCMC algorithm reads, and where each LD block's
 * SNPs sit in it.
 *
 * - beta_mrg: the marginal beta values for each SNP, deflated by `c`.
 * - boundary: the [start, end) positions of each LD block's SNPs.
 * - sz: the sample size of each SNP.
 * - array: the genotyping array of each SNP (0 = every cohort, 1 or 2 = the
 *   cohorts genotyped on that array). Pairs of SNPs on arrays 1 and 2 share
 *   no samples, so the optLlk = 2 likelihood zeroes their LD.
 *
 * Example usage:
 * @code
 * std::vector<double> beta_mrg = {0.1, 0.2, 0.3, 0.4, 0.5};
 * std::vector<size_t> block_sizes = {2, 3};
 * std::vector<double> sz = {1000, 1000, 1500, 1500, 1500};
 * std::vector<int> array = {1, 1, 2, 2, 2};
 *
 * mcmc_data data(beta_mrg, block_sizes, sz, array);
 * @endcode
 */
class mcmc_data {
public:
std::vector<double> beta_mrg;
std::vector<std::pair<size_t, size_t> > boundary;
std::vector<double> sz;
std::vector<int> array;

mcmc_data(const std::vector<double>& beta_mrg,
          const std::vector<size_t>& block_sizes,
          const std::vector<double>& sz,
          const std::vector<int>&    array)
	: beta_mrg(beta_mrg), sz(sz), array(array) {
	size_t start = 0;
	for (size_t size : block_sizes) {
		boundary.emplace_back(start, start + size);
		start += size;
	}
}
};

// What the sampler reads from one LD block, R. B is the matrix of the
// likelihood: n R (R + aI)^{-1} R for optLlk = 1, R B2^{-1} R for optLlk = 2
// (B2 the scaled R of solve_ldmat), and calc_b_tmp the matching
// A' beta_mrg. The sampler asks B for its diagonal, a product B beta, and the
// submatrix over the causal SNPs, and R for beta'R beta (the heritability).
// Each block type below answers those from one way of holding the matrices,
// so the sampler runs unchanged on either.

// B and R held whole, each symmetrized once.
struct SdprDenseBlock {
	arma::mat B;
	arma::mat R;
	arma::vec diag_B;
	arma::vec calc_b_tmp;

	arma::vec B_times(const arma::vec& beta) const { return B * beta; }
	arma::mat B_sub(const arma::uvec& idx) const { return B.submat(idx, idx); }
	double R_quad(const arma::vec& beta) const {
		return arma::dot(R * beta, beta);
	}
};

// R = Z Z' with Z (p x r) from a factor of R, and B = W W'. B has rank at
// most r under either likelihood, since it is R times a matrix times R, so a
// product costs O(pr) and neither p x p matrix is formed.
struct SdprFactorBlock {
	arma::mat W;
	arma::mat Z;
	arma::vec diag_B;
	arma::vec calc_b_tmp;

	arma::vec B_times(const arma::vec& beta) const {
		return W * (W.t() * beta);
	}
	arma::mat B_sub(const arma::uvec& idx) const {
		arma::mat W_c = W.rows(idx);
		return W_c * W_c.t();
	}
	double R_quad(const arma::vec& beta) const {
		arma::vec t = Z.t() * beta;
		return arma::dot(t, t);
	}
};

class MCMC_state {
public:
double alpha;
double eta;
double h2;
double N;
arma::vec beta;
arma::vec b;
std::vector<int> cls_assgn;
std::vector<double> V;
std::vector<double> p;
std::vector<double> log_p;
std::vector<double> cluster_var;
std::vector<unsigned> suff_stats;
std::vector<double> sumsq;
// Each block's terms of the eta update, from sample_beta().
std::vector<double> num;
std::vector<double> denom;
MCMC_state(size_t num_snp, size_t max_cluster, size_t num_blocks, \
           double a0, double b0, unsigned int seed) {
	a0k = a0; b0k = b0;
	// Changed May 20 2021
	// N (sz) is absorbed into A, B; so set to 1.
	N = 1.0;

	n_snp = num_snp;
	M = max_cluster;
	alpha = 1;
	eta = 1;
	beta = arma::vec(num_snp, arma::fill::zeros);
	b = arma::vec(num_snp, arma::fill::zeros);
	p.assign(max_cluster, 1.0/max_cluster);
	log_p.assign(max_cluster, 0);
	for (size_t i=0; i<max_cluster; i++) {
		log_p[i] = logf(p[i] + 1e-40);
	}

	cluster_var.assign(max_cluster, 0.0);
	suff_stats.assign(max_cluster, 0);
	sumsq.assign(max_cluster, 0.0);
	V.assign(max_cluster, 0.0);
	num.assign(num_blocks, 0.0);
	denom.assign(num_blocks, 0.0);
	// Initialize all SNPs to the null cluster (k=0). The original SDPR
	// used random initialization (uniform over 0..M-1), but this causes
	// the first sample_beta() call to allocate an enormous dense matrix
	// (nearly all SNPs are "causal"), crashing with "Mat::init() too large".
	// Starting from null is standard MCMC practice and lets the sampler
	// discover causal assignments organically.
	cls_assgn.assign(num_snp, 0);
	r.seed(seed);
}

void sample_sigma2();
template <class Block>
void calc_b(size_t j, const mcmc_data &dat, const std::vector<Block> &blocks);
template <class Block>
void sample_assignment(size_t j, const mcmc_data &dat, \
                       const std::vector<Block> &blocks);
void update_suffstats();
void sample_V();
void update_p();
void sample_alpha();
template <class Block>
void sample_beta(size_t j, const mcmc_data &dat, \
                 const std::vector<Block> &blocks);
template <class Block>
void compute_h2(const mcmc_data &dat, const std::vector<Block> &blocks);
void sample_eta();

private:
double a0k;
double b0k;
size_t M, n_snp;
std::mt19937 r;
};

class MCMC_samples {
public:
arma::vec beta;
double h2;
MCMC_samples(size_t num_snps) {
	beta = arma::vec(num_snps, arma::fill::zeros);
	h2 = 0;
}
};

/**
 * @brief The dense LD blocks of the SDPR likelihood (solve_ldmat).
 *
 * @param dat The per-SNP data, with beta_mrg already deflated by `c`.
 * @param ref_ld_mat The LD matrix of each block, taken over.
 * @param a Factor to shrink the reference LD matrix (optLlk = 1).
 * @param sz The sample size of the GWAS (optLlk = 1).
 * @param opt_llk 1 for equation 6 (slightly shrink the correlation of SNPs)
 *   and 2 for equation 5 (SNPs genotyped on different arrays in a separate
 *   cohort).
 */
std::vector<SdprDenseBlock> sdpr_dense_blocks(const mcmc_data &dat,
                                              std::vector<arma::mat> ref_ld_mat,
                                              double a, unsigned sz,
                                              int opt_llk);

/**
 * @brief The one LD block of the SDPR likelihood from a factor of its LD,
 * R = F'F / (m - 1), never forming R or B.
 *
 * The arguments are sdpr_dense_blocks()'s with the factor F (m x p) in place
 * of the LD matrices. B and A' beta_mrg agree with the dense block's to
 * rounding: optLlk = 1 has them in closed form from the eigenvalues of R;
 * optLlk = 2 applies B2^{-1} to r + 1 vectors, through the Woodbury identity
 * while B2's low-rank part is narrower than p, and otherwise by forming B2.
 */
SdprFactorBlock sdpr_factor_block(const mcmc_data &dat, const arma::mat &F,
                                  double a, unsigned sz, int opt_llk);

/**
 * @brief Perform Markov Chain Monte Carlo (MCMC) for estimating effect sizes.
 *
 * @param data The `mcmc_data` object, with beta_mrg deflated by `c`.
 * @param blocks The LD blocks, one per `data.boundary` entry
 *   (sdpr_dense_blocks() or sdpr_factor_block()).
 * @param M Max number of variance components. Default is 1000.
 * @param a0k Hyperparameter for inverse gamma distribution. Default is 0.5.
 * @param b0k Hyperparameter for inverse gamma distribution. Default is 0.5.
 * @param iter Number of iterations for MCMC. Default is 1000.
 * @param burn Number of burn-in iterations for MCMC. Default is 200.
 * @param thin Thinning interval for MCMC. Default is 5.
 * @param n_threads Number of threads to use. Default is 1.
 * @param verbose Whether to print verbose output. Default is true.
 * @param seed Seed of the sampler's random number generator.
 *
 * @return An `std::unordered_map` with "beta", the estimated effect sizes,
 *   and "h2", the estimated heritability.
 */
template <class Block>
std::unordered_map<std::string, arma::vec> mcmc(
	const mcmc_data&          data,
	const std::vector<Block>& blocks,
	size_t                    M,
	double                    a0k,
	double                    b0k,
	int                       iter,
	int                       burn,
	int                       thin,
	unsigned                  n_threads,
	bool                      verbose,
	unsigned int              seed
	);
