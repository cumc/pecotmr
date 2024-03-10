/**
 * @file mcmc.hpp
 * @brief Markov Chain Monte Carlo (MCMC) sampler for polygenic prediction with continuous shrinkage (CS) priors.
 */

#ifndef MCMC_HPP
#define MCMC_HPP

#include <armadillo>
#include <vector>
#include <string>
#include <cmath>
#include <map>
#include <iomanip>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

/**
 * @brief Evaluate the function psi(x, alpha, lambda).
 *
 * @param x Input value.
 * @param alpha Shape parameter.
 * @param lambda Scale parameter.
 * @return Value of psi(x, alpha, lambda).
 */
double fpsi(double x, double alpha, double lambda) {
	double f = -alpha * (std::cosh(x) - 1.0) - lambda * (std::exp(x) - x - 1.0);
	return f;
}

/**
 * @brief Evaluate the derivative of fpsi(x, alpha, lambda).
 *
 * @param x Input value.
 * @param alpha Shape parameter.
 * @param lambda Scale parameter.
 * @return Value of fdpsi(x, alpha, lambda).
 */
double fdpsi(double x, double alpha, double lambda) {
	double f = -alpha * std::sinh(x) - lambda * (std::exp(x) - 1.0);
	return f;
}

/**
 * @brief Evaluate the function g(x, sd, td, f1, f2).
 *
 * @param x Input value.
 * @param sd Parameter sd.
 * @param td Parameter td.
 * @param f1 Parameter f1.
 * @param f2 Parameter f2.
 * @return Value of g(x, sd, td, f1, f2).
 */
double g(double x, double sd, double td, double f1, double f2) {
	if (x >= -sd && x <= td) {
		return 1.0;
	} else if (x > td) {
		return f1;
	} else {
		return f2;
	}
}

/**
 * @brief Generate random variates from the generalized inverse Gaussian distribution.
 *
 * @param p Shape parameter.
 * @param a Scale parameter.
 * @param b Scale parameter.
 * @return Random variate from the generalized inverse Gaussian distribution.
 */
double gigrnd(double p, double a, double b) {
	double lambda = p;
	double omega = std::sqrt(a * b);
	bool swap = false;
	if (lambda < 0) {
		lambda = -lambda;
		swap = true;
	}

	double alpha = std::sqrt(std::pow(omega, 2) + std::pow(lambda, 2)) - lambda;

// find t
	double x = -fpsi(1.0, alpha, lambda);
	double t;
	if (x >= 0.5 && x <= 2.0) {
		t = 1.0;
	} else if (x > 2.0) {
		if (alpha == 0 && lambda == 0) {
			t = 1.0;
		} else {
			t = std::sqrt(2.0 / (alpha + lambda));
		}
	} else if (x < 0.5) {
		if (alpha == 0 && lambda == 0) {
			t = 1.0;
		} else {
			t = std::log(4.0 / (alpha + 2.0 * lambda));
		}
	}

// find s
	x = -fpsi(-1.0, alpha, lambda);
	double s;
	if (x >= 0.5 && x <= 2.0) {
		s = 1.0;
	} else if (x > 2.0) {
		if (alpha == 0 && lambda == 0) {
			s = 1.0;
		} else {
			s = std::sqrt(4.0 / (alpha * std::cosh(1.0) + lambda));
		}
	} else if (x < 0.5) {
		if (alpha == 0 && lambda == 0) {
			s = 1.0;
		} else if (alpha == 0) {
			s = 1.0 / lambda;
		} else if (lambda == 0) {
			s = std::log(1.0 + 1.0 / alpha + std::sqrt(1.0 / std::pow(alpha, 2) + 2.0 / alpha));
		} else {
			s = std::min(1.0 / lambda, std::log(1.0 + 1.0 / alpha + std::sqrt(1.0 / std::pow(alpha, 2) + 2.0 / alpha)));
		}
	}



	double eta = -fpsi(t, alpha, lambda);
	double zeta = -fdpsi(t, alpha, lambda);
	double theta = -fpsi(-s, alpha, lambda);
	double xi = fdpsi(-s, alpha, lambda);

	double p_r = 1.0 / xi;
	double r = 1.0 / zeta;

	double td = t - r * eta;
	double sd = s - p_r * theta;
	double q = td + sd;

	gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
	double rnd = 0.0;
	while (true) {
		double U = gsl_rng_uniform(rng);
		double V = gsl_rng_uniform(rng);
		double W = gsl_rng_uniform(rng);

		if (U < q / (p_r + q + r)) {
			rnd = -sd + q * V;
		} else if (U < (q + r) / (p_r + q + r)) {
			rnd = td - r * std::log(V);
		} else {
			rnd = -sd + p_r * std::log(V);
		}

		double f1 = std::exp(-eta - zeta * (rnd - t));
		double f2 = std::exp(-theta + xi * (rnd + s));

		if (W * g(rnd, sd, td, f1, f2) <= std::exp(fpsi(rnd, alpha, lambda))) {
			break;
		}
	}

	gsl_rng_free(rng);

	rnd = std::exp(rnd) * (lambda / omega + std::sqrt(1.0 + std::pow(lambda / omega, 2)));

	if (swap) {
		rnd = 1.0 / rnd;
	}

	double result = rnd / std::sqrt(a / b);

	// Check if the result is zero and replace it with a small value
	if (result == 0.0) {
		result = std::numeric_limits<double>::min();
	}
	if (result > 1.0) {
		result = 1.0;
	}

	return result;
}

/**
 * @brief Markov Chain Monte Carlo (MCMC) sampler for polygenic prediction with continuous shrinkage (CS) priors.
 *
 * @param a Shape parameter for the prior distribution of psi.
 * @param b Scale parameter for the prior distribution of psi.
 * @param phi Global shrinkage parameter. If nullptr, it will be estimated automatically.
 * @param sumstats Two vectors containing BETA and MAF summary statistics.
 * @param n Sample size.
 * @param ld_blk List of LD blocks.
 * @param n_iter Number of MCMC iterations.
 * @param n_burnin Number of burn-in iterations.
 * @param thin Thinning interval.
 * @param beta_std Whether to standardize the effect sizes.
 * @param verbose Whether to print verbose output.
 * @param seed Random seed. If nullptr, no seed is set.
 * @return A map containing the posterior estimates.
 */
std::map<std::string, arma::vec> prs_cs_mcmc(double a, double b, double* phi, const std::vector<std::vector<double> >& sumstats,
                                             int n, const std::vector<arma::mat>& ld_blk,
                                             int n_iter, int n_burnin, int thin,
                                             bool beta_std, bool verbose, int* seed) {
	if (verbose) {
		std::cout << "Running Markov Chain Monte Carlo (MCMC) sampler..." << std::endl;
	}

	// Seed the random number generator
	gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
	if (seed != nullptr) {
		gsl_rng_set(rng, *seed);
	}

	// Derived statistics
	arma::vec beta_mrg(sumstats[0]);
	arma::vec maf(sumstats[1]);
	int n_pst = (n_iter - n_burnin) / thin;
	int p = beta_mrg.n_elem;

	// Initialization
	arma::vec beta(p, arma::fill::zeros);
	arma::vec psi(p, arma::fill::ones);
	double sigma = 1.0;
	bool phi_updt = (phi == nullptr);
	if (phi_updt) {
		phi = new double(1.0);
	}

	arma::vec beta_est(p, arma::fill::zeros);
	arma::vec psi_est(p, arma::fill::zeros);
	double sigma_est = 0.0;
	double phi_est = 0.0;

	// MCMC
	for (int itr = 1; itr <= n_iter; ++itr) {
		if (verbose && itr % 100 == 0) {
			std::cout << "Iteration " << std::setw(4) << itr << " of " << n_iter << std::endl;
		}

		int mm = 0;
		double quad = 0.0;
		for (int kk = 0; kk < ld_blk.size(); ++kk) {
			if (ld_blk[kk].n_rows == 0) {
				continue;
			}

			arma::uvec idx_blk = arma::regspace<arma::uvec>(mm, mm + ld_blk[kk].n_rows - 1);
			mm += ld_blk[kk].n_rows;

			//std::cout << psi << std::endl;
			//std::cout << idx_blk << std::endl;
			//std::cout << ld_blk[kk] << std::endl;

			arma::mat dinvt = ld_blk[kk] + arma::diagmat(1.0 / psi.elem(idx_blk));
			//std::cout << "dinvt " << dinvt << std::endl;

			arma::mat dinvt_chol = arma::chol(dinvt);

			//std::cout << "dinvt chol " << dinvt_chol << std::endl;

			arma::vec beta_tmp = arma::solve(arma::trimatl(dinvt_chol.t()), beta_mrg(idx_blk)) +
			                     arma::vec(ld_blk[kk].n_rows).transform([&](double) {
				return gsl_ran_gaussian(rng, 1.0);
			}) * std::sqrt(sigma / n);
			std::cout << "beta_tmp " << beta_tmp << std::endl;
			beta(idx_blk) = arma::solve(arma::trimatu(dinvt_chol), beta_tmp);

			quad += arma::as_scalar(beta(idx_blk).t() * dinvt * beta(idx_blk));
			std::cout << "quad " << quad << std::endl;
		}
		std::cout << "beta " << beta << std::endl;

		double err = std::max(n / 2.0 * (1.0 - 2.0 * arma::dot(beta, beta_mrg) + quad),
		                      n / 2.0 * arma::sum(arma::pow(beta, 2) / psi));

		sigma = 1.0 / gsl_ran_gamma(rng, (n + p) / 2.0, 1.0 / err);

		arma::vec delta = arma::vec(p);
		for (int jj = 0; jj < p; ++jj) {
			delta(jj) = gsl_ran_gamma(rng, a + b, 1.0 / (psi(jj) + *phi));
		}

		// FIXME: psi can be so small (close to zero) as in the @example data demo which causes inf in dinvt and solve() to break
		// This happens when beta are all very small and delta all very big
		// std::cout << "sigma " << sigma << std::endl;
		// std::cout << "delta " << delta << std::endl;
		for (int jj = 0; jj < p; ++jj) {
			psi(jj) = gigrnd(a - 0.5, 2.0 * delta(jj), n * std::pow(beta(jj), 2) / sigma);
		}

		if (phi_updt) {
			double w = gsl_ran_gamma(rng, 1.0, 1.0 / (*phi + 1.0));
			*phi = gsl_ran_gamma(rng, p * b + 0.5, 1.0 / (arma::sum(delta) + w));
		}

		// Posterior
		if (itr > n_burnin && (itr % thin == 0)) {
			beta_est += beta / n_pst;
			psi_est += psi / n_pst;
			sigma_est += sigma / n_pst;
			phi_est += *phi / n_pst;
		}
	}

	gsl_rng_free(rng);

	if (phi_updt) {
		delete phi;
	}

	// Convert standardized beta to per-allele beta
	if (!beta_std) {
		beta_est /= arma::sqrt(2.0 * maf % (1.0 - maf));
	}

	// Prepare the output map
	std::map<std::string, arma::vec> output;
	output["beta_est"] = beta_est;
	output["psi_est"] = psi_est;
	output["sigma_est"] = arma::vec(1).fill(sigma_est);
	output["phi_est"] = arma::vec(1).fill(phi_est);

	// Print estimated phi
	if (verbose && phi_updt) {
		std::cout << "Estimated global shrinkage parameter: " << phi_est << std::endl;
	}

	if (verbose) {
		std::cout << "MCMC sampling completed." << std::endl;
	}

	return output;
}
#endif // MCMC_HPP