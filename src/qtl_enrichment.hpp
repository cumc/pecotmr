#ifndef QTL_ENRICHMENT_HPP
#define QTL_ENRICHMENT_HPP
#include <RcppArmadillo.h> // need to include this before RcppGSL otherwise it complains about conflicts
#include <RcppGSL.h>
#include <vector>
#include <string>
#include <map>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <omp.h>
#include <cmath>
#include <cstdio>

// Enable C++11
// [[Rcpp::plugins(cpp11)]]
// Enable openmp
// [[Rcpp::plugins(openmp)]]
// Import Armadillo
// [[Rcpp::depends(RcppArmadillo)]]
// Import GSL
// [[Rcpp::depends(RcppGSL)]]

class SuSiEFit {
public:
std::vector<double> pip;
std::vector<std::string> variable_names;
arma::mat alpha;
std::vector<double> prior_variance;

SuSiEFit(SEXP r_susie_fit) {
	Rcpp::List susie_fit(r_susie_fit);

	Rcpp::NumericVector pip_vec = Rcpp::as<Rcpp::NumericVector>(susie_fit["pip"]);
	pip = Rcpp::as<std::vector<double> >(pip_vec);
	variable_names = Rcpp::as<std::vector<std::string> >(pip_vec.names());
	alpha = Rcpp::as<arma::mat>(susie_fit["alpha"]);
	prior_variance = Rcpp::as<std::vector<double> >(susie_fit["prior_variance"]);

	if (alpha.n_rows != prior_variance.size()) {
		Rcpp::stop("The number of rows in alpha must match the length of prior_variance.");
	}

	// Check if all elements in prior_variance are not greater than 0
	if (std::all_of(prior_variance.begin(), prior_variance.end(), [](double x) {
			return x <= 0;
		})) {
		Rcpp::stop("At least one element in prior_variance must be greater than 0.");
	}

	// Filter out rows with prior_variance = 0
	std::vector<arma::uword> valid_rows;
	for (size_t i = 0; i < prior_variance.size(); ++i) {
		if (prior_variance[i] > 0) {
			valid_rows.push_back(i);
		}
	}
	alpha = alpha.rows(arma::uvec(valid_rows));
	prior_variance.erase(std::remove(prior_variance.begin(), prior_variance.end(), 0), prior_variance.end());

	// Add a check to make sure each row of alpha sums to 1
	for (arma::uword i = 0; i < alpha.n_rows; ++i) {
		double row_sum = arma::sum(alpha.row(i));
		if (std::abs(row_sum - 1.0) > 1e-6) {
			Rcpp::stop("Row " + std::to_string(i + 1) + " of single effect PIP matrix (alpha) does not sum to 1. It is: " + std::to_string(row_sum));
		}
	}
}

std::vector<std::string> impute_qtn(const gsl_rng *r) const {
	std::vector<std::string> qtn_names;

	for (arma::uword i = 0; i < alpha.n_rows; ++i) {
		std::vector<double> alpha_row(alpha.colptr(i), alpha.colptr(i) + alpha.n_cols);
		std::unique_ptr<unsigned int[]> sample(new unsigned int[alpha_row.size()]);
		gsl_ran_multinomial(r, alpha_row.size(), 1, alpha_row.data(), sample.get());

		for (size_t j = 0; j < alpha_row.size(); ++j) {
			if (sample[j] == 1) {
				qtn_names.push_back(variable_names[j]);
				break;
			}
		}
	}

	return qtn_names;
}
};

std::vector<double> run_EM(
	const std::vector<double> &gwas_pip,
	const std::vector<int> &   qtl_sample,
	double                     pi_qtl,
	double                     pi_gwas,
	int                        max_iter = 1000,
	double                     a1_tol = 0.01)
{
	double a0 = log(pi_gwas / (1 - pi_gwas));
	double a1 = 0;
	double var0 = 0;
	double var1 = 0;
	double r1 = exp(a0 + a1);
	double r0 = exp(a0);
	double r_null = pi_gwas / (1 - pi_gwas);
	double total_snp = gwas_pip.size();
	int iter = 0;

	while (true) {
		iter++;
		double pseudo_count = 1.0;
		double e0g0 = pseudo_count * (1 - pi_gwas) * (1 - pi_qtl);
		double e0g1 = pseudo_count * (1 - pi_qtl) * pi_gwas;
		double e1g0 = pseudo_count * (1 - pi_gwas) * pi_qtl;
		double e1g1 = pseudo_count * pi_gwas * pi_qtl;

		for (size_t i = 0; i < gwas_pip.size(); i++) {
			double val = gwas_pip[i];
			if (val == 1)
				val = 1 - 1e-8;
			val = val / (1 - val);

			if (qtl_sample[i] == 0) {
				val = r0 * (val / r_null);
				val = val / (1 + val);
				e0g1 += val;
				e0g0 += 1 - val;
			}

			if (qtl_sample[i] == 1) {
				val = r1 * (val / r_null);
				val = val / (1 + val);
				e1g1 += val;
				e1g0 += 1 - val;
			}
		}

		e0g0 += total_snp - (e0g0 + e0g1 + e1g0 + e1g1);

		double a1_new = log(e1g1 * e0g0 / (e1g0 * e0g1));
		a1 = a1_new;
		a0 = log(e0g1 / e0g0);
		r0 = exp(a0);
		r1 = exp(a0 + a1);
		var1 = (1.0 / e0g0 + 1.0 / e1g0 + 1.0 / e1g1 + 1.0 / e0g1);
		var0 = (1.0 / e0g1 + 1.0 / e0g0);

		if (fabs(a1_new - a1) < a1_tol || iter >= max_iter) {
			break;
		}
		if (iter % 100 == 0) {
			Rcpp::Rcout << "EM Iteration " << iter << ": a0 = " << a0 << ", a1 = " << a1 << std::endl;
		}
	}
	if (iter == max_iter) {
		Rcpp::Rcout << "WARNING: EM algorithm did not converge after " << iter << "iterations!" << std::endl;
	}

	std::vector<double> av;
	av.push_back(a0);
	av.push_back(a1);
	av.push_back(var0);
	av.push_back(var1);

	return av;
}

std::map<std::string, double> qtl_enrichment_workhorse(
	const std::vector<SuSiEFit> &   qtl_susie_fits,
	const std::vector<double> &     gwas_pip,
	const std::vector<std::string> &gwas_variable_names,
	double                          pi_gwas,
	double                          pi_qtl,
	int                             ImpN,
	double                          shrinkage_lambda,
	int                             num_threads = 4)
{

	std::vector<double> a0_vec(ImpN, 0.0);
	std::vector<double> v0_vec(ImpN, 0.0);
	std::vector<double> a1_vec(ImpN, 0.0);
	std::vector<double> v1_vec(ImpN, 0.0);

	std::map<std::string, int> gwas_variant_index;

	for (size_t i = 0; i < gwas_variable_names.size(); ++i) {
		gwas_variant_index[gwas_variable_names[i]] = i;
	}

	Rcpp::Rcout << "Data loaded successfully!" << std::endl;

    #pragma omp parallel for num_threads(num_threads)
	for (int k = 0; k < ImpN; k++) {
		// Initialize the GSL RNG for this thread
		const gsl_rng_type *T;
		gsl_rng *r;
		gsl_rng_env_setup();
		T = gsl_rng_default;
		r = gsl_rng_alloc(T);
		// Set the seed for this thread's RNG
		gsl_rng_set(r, static_cast<unsigned long int>(k));

		std::vector<int> qtl_sample(gwas_pip.size(), 0);

		for (size_t i = 0; i < qtl_susie_fits.size(); i++) {
			std::vector<std::string> variants = qtl_susie_fits[i].impute_qtn(r);
			for (const auto &variant : variants) {
				qtl_sample[gwas_variant_index[variant]] = 1;
			}
		}
		gsl_rng_free(r);

		std::vector<double> rst = run_EM(gwas_pip, qtl_sample, pi_qtl, pi_gwas);

	#pragma omp critical
		{
			a0_vec[k] = rst[0];
			a1_vec[k] = rst[1];
			v0_vec[k] = rst[2];
			v1_vec[k] = rst[3];
		}
	}

	Rcpp::Rcout << "EM algorithm completed!" << std::endl;

	double a0_est = 0;
	double a1_est = 0;
	double var0 = 0;
	double var1 = 0;
	for (int k = 0; k < ImpN; k++) {
		a0_est += a0_vec[k];
		a1_est += a1_vec[k];
		var0 += v0_vec[k];
		var1 += v1_vec[k];
	}
	a0_est /= ImpN;
	a1_est /= ImpN;

	double bv0 = 0;
	double bv1 = 0;
	for (int k = 0; k < ImpN; k++) {
		bv0 += pow(a0_vec[k] - a0_est, 2.0);
		bv1 += pow(a1_vec[k] - a1_est, 2.0);
	}
	bv0 /= (ImpN - 1);
	bv1 /= (ImpN - 1);
	var0 /= ImpN;
	var1 /= ImpN;

	double sd0 = sqrt(var0 + bv0 * (ImpN + 1) / ImpN);
	double sd1 = sqrt(var1 + bv1 * (ImpN + 1) / ImpN);

	double a1_est_ns = a1_est;
	double sd1_ns = sd1;

	// Apply shrinkage
	double pv = (shrinkage_lambda == 0) ? -1 : 1.0/shrinkage_lambda;
	if (pv > 0) {
		double post_var = 1.0 / (1.0 / pv + 1 / (sd1 * sd1));
		a1_est = (a1_est_ns * pv) / (pv + sd1_ns * sd1_ns);
		sd1 = sqrt(post_var);
	}

	a0_est = log(pi_gwas / (1 + pi_qtl * exp(a1_est) - pi_qtl - pi_gwas));

	double p1 = (1 - pi_qtl) * exp(a0_est) / (1 + exp(a0_est));
	double p2 = pi_qtl / (1 + exp(a0_est + a1_est));
	double p12 = pi_qtl * exp(a0_est + a1_est) / (1 + exp(a0_est + a1_est));

	double pi1_e = exp(a0_est + a1_est) / (1 + exp(a0_est + a1_est));
	double pi1_ne = exp(a0_est) / (1 + exp(a0_est));

	// Create the map to store output
	std::map<std::string, double> output_map;
	output_map["Intercept"] = a0_est;
	output_map["sd (intercept)"] = sd0;
	output_map["Enrichment (no shrinkage)"] = a1_est_ns;
	output_map["Enrichment (w/ shrinkage)"] = a1_est;
	output_map["sd (no shrinkage)"] = sd1_ns;
	output_map["sd (w/ shrinkage)"] = sd1;
	output_map["Alternative (coloc) p1"] = p1;
	output_map["Alternative (coloc) p2"] = p2;
	output_map["Alternative (coloc) p12"] = p12;
	output_map["pi1_e"] = pi1_e;
	output_map["pi1_ne"] = pi1_ne;

	return output_map;
}

#endif // QTL_ENRICHMENT_HPP