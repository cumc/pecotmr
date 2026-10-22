#include <cpp11.hpp>
#include <cpp11armadillo.hpp>
#include <random>
#include <unordered_map>
#include "sdpr_mcmc.h"

using namespace cpp11;
using namespace arma;

namespace {

// The per-SNP data as R passes it: the marginal betas deflated by `c`, a NULL
// perVariantSampleSize is n for every SNP, and a NULL array is array 1 for
// every SNP.
mcmc_data sdpr_data(const doubles& bhatR, const std::vector<size_t>& block_sizes,
                    int n, sexp perVariantSampleSize, sexp array, double c) {
	std::vector<double> bhat = cpp11::as_cpp<std::vector<double>>(bhatR);
	// Deflation correction
	for (size_t i = 0; i < bhat.size(); i++) {
		bhat[i] /= c;
	}
	std::vector<double> sz;
	std::vector<int> arr;
	if (perVariantSampleSize != R_NilValue) {
		sz = cpp11::as_cpp<std::vector<double>>(perVariantSampleSize);
	} else {
		sz = std::vector<double>(bhat.size(), n);
	}
	if (array != R_NilValue) {
		arr = cpp11::as_cpp<std::vector<int>>(array);
	} else {
		arr = std::vector<int>(bhat.size(), 1);
	}
	return mcmc_data(bhat, block_sizes, sz, arr);
}

// The seed R passes, or a drawn one for NULL.
unsigned int sdpr_seed(sexp seed) {
	if (seed != R_NilValue) {
		return cpp11::as_cpp<unsigned int>(seed);
	}
	return std::random_device{}();
}

cpp11::writable::list sdpr_output(std::unordered_map<std::string, vec>& results) {
	using namespace cpp11::literals;
	writable::list output({
		"betaEst"_nm = as_doubles(results["beta"]),
		"h2"_nm = as_doubles(results["h2"])
	});
	return output;
}

}  // namespace

// cpp11 interface function
[[cpp11::register]]
cpp11::writable::list sdprCpp(
	const doubles& bhatR,
	const list& LD,
	int n,
	sexp perVariantSampleSize = R_NilValue,
	sexp array = R_NilValue,
	double a = 0.1,
	double c = 1.0,
	int M = 1000,
	double a0k = 0.5,
	double b0k = 0.5,
	int iter = 1000,
	int burn = 200,
	int thin = 5,
	int numThreads = 1,
	int optLlk = 1,
	bool verbose = true,
	sexp seed = R_NilValue
	) {
	// Convert list to std::vector<arma::mat>
	std::vector<mat> ref_ld_mat;
	std::vector<size_t> block_sizes;
	for (int i = 0; i < LD.size(); i++) {
		ref_ld_mat.push_back(as_Mat(doubles_matrix<>(LD[i])));
		block_sizes.push_back(ref_ld_mat.back().n_rows);
	}

	mcmc_data data = sdpr_data(bhatR, block_sizes, n, perVariantSampleSize, array, c);
	std::vector<SdprDenseBlock> blocks =
		sdpr_dense_blocks(data, std::move(ref_ld_mat), a, n, optLlk);

	std::unordered_map<std::string, vec> results = mcmc(
		data, blocks, M, a0k, b0k, iter, burn, thin, numThreads, verbose,
		sdpr_seed(seed)
		);
	return sdpr_output(results);
}

/**
 * @brief SDPR from a factor F of the LD, R = F'F / (m - 1), never forming R.
 *
 * The arguments are sdprCpp()'s with `factorR` (m x p) in place of the LD
 * blocks; the window is one block. The sampler is SDPR's own, run on the
 * low-rank form of its matrices (sdpr_factor_block()), so it draws as
 * sdprCpp() does on the matrix the factor gives.
 */
[[cpp11::register]]
cpp11::writable::list sdprFactorCpp(
	const doubles& bhatR,
	const doubles_matrix<>& factorR,
	int n,
	sexp perVariantSampleSize = R_NilValue,
	sexp array = R_NilValue,
	double a = 0.1,
	double c = 1.0,
	int M = 1000,
	double a0k = 0.5,
	double b0k = 0.5,
	int iter = 1000,
	int burn = 200,
	int thin = 5,
	int numThreads = 1,
	int optLlk = 1,
	bool verbose = true,
	sexp seed = R_NilValue
	) {
	mat F = as_Mat(factorR);
	if (F.n_cols != static_cast<uword>(bhatR.size()))
		cpp11::stop("The LD factor has %d columns but there are %d effects.",
		            (int)F.n_cols, (int)bhatR.size());

	mcmc_data data = sdpr_data(bhatR, {F.n_cols}, n, perVariantSampleSize, array, c);
	std::vector<SdprFactorBlock> blocks = {
		sdpr_factor_block(data, F, a, n, optLlk)
	};

	std::unordered_map<std::string, vec> results = mcmc(
		data, blocks, M, a0k, b0k, iter, burn, thin, numThreads, verbose,
		sdpr_seed(seed)
		);
	return sdpr_output(results);
}
