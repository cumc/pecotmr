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
 * @param sumstats Dictionary containing summary statistics.
 * @param n Sample size.
 * @param ld_blk List of LD blocks.
 * @param n_iter Number of MCMC iterations.
 * @param n_burnin Number of burn-in iterations.
 * @param thin Thinning interval.
 * @param beta_std Whether to standardize the effect sizes.
 * @param verbose Whether to print verbose output.
 * @param seed Random seed. If nullptr, no seed is set.
 * @return A list containing the posterior estimates.
 */
// [[Rcpp::export]]
Rcpp::List prs_cs_rcpp(double a, double b, Rcpp::Nullable<double> phi, Rcpp::List sumstats,
                       int n, Rcpp::List ld_blk, 
                       int n_iter, int n_burnin, int thin, 
                       bool beta_std, bool verbose, Rcpp::Nullable<int> seed) {
    // Convert Rcpp types to C++ types
    std::vector<std::vector<double>> sumstats_vec;
    sumstats_vec.push_back(Rcpp::as<std::vector<double>>(sumstats["BETA"]));
    sumstats_vec.push_back(Rcpp::as<std::vector<double>>(sumstats["MAF"]));

    std::vector<arma::mat> ld_blk_vec;
    for (int i = 0; i < ld_blk.size(); ++i) {
        ld_blk_vec.push_back(Rcpp::as<arma::mat>(ld_blk[i]));
    }


    double* phi_ptr = nullptr;
    if (phi.isNotNull()) {
        phi_ptr = new double(Rcpp::as<double>(phi));
    }

    int* seed_ptr = nullptr;
    if (seed.isNotNull()) {
        seed_ptr = new int(Rcpp::as<int>(seed));
    }

    // Call the prs_cs function
    std::map<std::string, arma::vec> output = prs_cs_mcmc(a, b, phi_ptr, sumstats_vec, n, ld_blk_vec,
                                                     n_iter, n_burnin, thin, beta_std, verbose, seed_ptr);

    // Convert the output to an Rcpp::List
    Rcpp::List result;
    result["beta_est"] = output["beta_est"];
    result["psi_est"] = output["psi_est"];
    result["sigma_est"] = output["sigma_est"](0);
    result["phi_est"] = output["phi_est"](0);

    // Clean up dynamically allocated memory
    delete phi_ptr;
    delete seed_ptr;

    return result;
}