/**
 * @file prscs_mcmc.cpp
 * @brief cpp11 entry points for the PRS-CS sampler (prscs_mcmc.h).
 */

#include <cpp11.hpp>
#include <cpp11armadillo.hpp>
#include <cmath>
#include <utility>
#include "prscs_mcmc.h"

using namespace cpp11;
using namespace arma;

namespace {

// Run the sampler on `ld_blk` with the arguments as R passes them: a NULL
// `maf` is a vector of zeros, a NULL `phi` is learned from the data, and a
// NULL `seed` draws one.
template <class LdBlock>
cpp11::writable::list prs_cs_run(double a, double b, sexp phi, doubles bhat, sexp maf, int n,
                                 const std::vector<LdBlock>& ld_blk, int nIter, int nBurnin,
                                 int thin, bool verbose, sexp seed) {
	std::vector<double> bhat_vec = cpp11::as_cpp<std::vector<double>>(bhat);
	std::vector<double> maf_vec;
	if (maf != R_NilValue) {
		maf_vec = cpp11::as_cpp<std::vector<double>>(maf);
	} else {
		maf_vec = std::vector<double>(bhat_vec.size(), 0.0); // Populate with zeros if maf is NULL
	}

	// Use stack variable to avoid heap allocation and memory leak risk.
	double phi_val = 0.0;
	double* phi_ptr = nullptr;
	if (phi != R_NilValue) {
		phi_val = cpp11::as_cpp<double>(phi);
		phi_ptr = &phi_val;
	}

	unsigned int seed_val = 0;
	if (seed != R_NilValue) {
		seed_val = cpp11::as_cpp<unsigned int>(seed);
	} else {
		seed_val = std::random_device{}();
	}

	std::map<std::string, vec> output = prs_cs_mcmc(a, b, phi_ptr, bhat_vec, maf_vec, n, ld_blk,
	                                                nIter, nBurnin, thin, verbose, seed_val);

	using namespace cpp11::literals;
	writable::list result({
		"betaEst"_nm = as_doubles(output["beta_est"]),
		"psiEst"_nm = as_doubles(output["psi_est"]),
		"sigmaEst"_nm = cpp11::as_sexp(output["sigma_est"](0)),
		"phiEst"_nm = cpp11::as_sexp(output["phi_est"](0))
	});
	return result;
}

}  // namespace

/**
 * @brief PRS-CS from LD matrices.
 *
 * @param a Shape parameter for the prior distribution of psi.
 * @param b Scale parameter for the prior distribution of psi.
 * @param phi Global shrinkage parameter. If NULL, it will be estimated automatically.
 * @param bhat Vector of effect sizes.
 * @param maf Vector of minor allele frequencies. If NULL, it is assumed to be a vector of zeros.
 * @param n Sample size.
 * @param ldBlk List of LD blocks.
 * @param nIter Number of MCMC iterations.
 * @param nBurnin Number of burn-in iterations.
 * @param thin Thinning interval.
 * @param verbose Whether to print verbose output.
 * @param seed Random seed. If NULL, no seed is set.
 * @return A list containing the posterior estimates.
 */
[[cpp11::register]]
cpp11::writable::list prsCsCpp(double a, double b, sexp phi,
                       doubles bhat, sexp maf,
                       int n, list ldBlk,
                       int nIter, int nBurnin, int thin,
                       bool verbose, sexp seed) {
	std::vector<DenseLdBlock> ld_blk_vec;
	for (int i = 0; i < ldBlk.size(); ++i) {
		ld_blk_vec.emplace_back(as_Mat(doubles_matrix<>(ldBlk[i])));
	}
	return prs_cs_run(a, b, phi, bhat, maf, n, ld_blk_vec, nIter, nBurnin, thin, verbose, seed);
}

/**
 * @brief PRS-CS from a factor F of one LD block, R = F'F / (m - 1), never forming R.
 *
 * The arguments are prsCsCpp()'s with `factorR` in place of the LD blocks. The chain
 * targets the same posterior as prsCsCpp() on R but draws differently (FactorLdBlock).
 */
[[cpp11::register]]
cpp11::writable::list prsCsFactorCpp(double a, double b, sexp phi,
                             doubles bhat, sexp maf,
                             int n, const doubles_matrix<>& factorR,
                             int nIter, int nBurnin, int thin,
                             bool verbose, sexp seed) {
	mat F = as_Mat(factorR);
	if (F.n_cols != static_cast<uword>(bhat.size()))
		cpp11::stop("The LD factor has %d columns but there are %d effects.",
		            (int)F.n_cols, (int)bhat.size());
	std::vector<FactorLdBlock> ld_blk_vec;
	ld_blk_vec.emplace_back(F / std::sqrt(static_cast<double>(F.n_rows) - 1.0));
	return prs_cs_run(a, b, phi, bhat, maf, n, ld_blk_vec, nIter, nBurnin, thin, verbose, seed);
}
