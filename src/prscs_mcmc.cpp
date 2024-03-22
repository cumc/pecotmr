/**
 * @file prs_cs_wrapper.cpp
 * @brief Rcpp wrapper for the prs_cs function.
 */

#include <RcppArmadillo.h>
#include "prscs_mcmc.h"

// [[Rcpp::depends(RcppArmadillo)]]
/**
 * @brief Rcpp wrapper for the prs_cs function.
 *
 * @param a Shape parameter for the prior distribution of psi.
 * @param b Scale parameter for the prior distribution of psi.
 * @param phi Global shrinkage parameter. If nullptr, it will be estimated automatically.
 * @param bhat Vector of effect sizes.
 * @param maf Vector of minor allele frequencies. If nullptr, it is assumed to be a vector of zeros.
 * @param n Sample size.
 * @param ld_blk List of LD blocks.
 * @param n_iter Number of MCMC iterations.
 * @param n_burnin Number of burn-in iterations.
 * @param thin Thinning interval.
 * @param verbose Whether to print verbose output.
 * @param seed Random seed. If nullptr, no seed is set.
 * @return A list containing the posterior estimates.
 */
// [[Rcpp::export]]
Rcpp::List prs_cs_rcpp(double a, double b, Rcpp::Nullable<double> phi,
                       Rcpp::NumericVector bhat, Rcpp::Nullable<Rcpp::NumericVector> maf,
                       int n, Rcpp::List ld_blk,
                       int n_iter, int n_burnin, int thin,
                       bool verbose, Rcpp::Nullable<unsigned int> seed) {
	// Convert Rcpp types to C++ types
	std::vector<double> bhat_vec = Rcpp::as<std::vector<double> >(bhat);
	std::vector<double> maf_vec;
	if (maf.isNotNull()) {
		maf_vec = Rcpp::as<std::vector<double> >(maf.get());
	} else {
		maf_vec = std::vector<double>(bhat_vec.size(), 0.0); // Populate with zeros if maf is NULL
	}

	std::vector<arma::mat> ld_blk_vec;
	for (int i = 0; i < ld_blk.size(); ++i) {
		ld_blk_vec.push_back(Rcpp::as<arma::mat>(ld_blk[i]));
	}

	double* phi_ptr = nullptr;
	if (phi.isNotNull()) {
		phi_ptr = new double(Rcpp::as<double>(phi));
	}

	unsigned int seed_val = 0;
	if (seed.isNotNull()) {
		seed_val = Rcpp::as<unsigned int>(seed);
	} else {
		seed_val = std::random_device{}();
	}

	std::map<std::string, arma::vec> output = prs_cs_mcmc(a, b, phi_ptr, bhat_vec, maf_vec, n, ld_blk_vec,
	                                                      n_iter, n_burnin, thin, verbose, seed_val);

	// Convert the output to an Rcpp::List
	Rcpp::List result;
	result["beta_est"] = output["beta_est"];
	result["psi_est"] = output["psi_est"];
	result["sigma_est"] = output["sigma_est"](0);
	result["phi_est"] = output["phi_est"](0);

	// Clean up dynamically allocated memory
	delete phi_ptr;
	return result;
}
