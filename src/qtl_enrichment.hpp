#ifndef QTL_ENRICHMENT_HPP
#define QTL_ENRICHMENT_HPP

#include <RcppArmadillo.h>
#include <vector>
#include <string>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

// Enable C++11
// [[Rcpp::plugins(cpp11)]]

// Import Armadillo
// [[Rcpp::depends(RcppArmadillo)]]

class SuSiEFit {
public:
    std::vector<double> pip;
    std::vector<std::string> pip_names;
    arma::mat alpha;
    std::vector<double> prior_variance;

    SuSiEFit(SEXP r_susie_fit) {
        Rcpp::List susie_fit(r_susie_fit);

        pip = Rcpp::as<std::vector<double>>(susie_fit["pip"]);
        pip_names = Rcpp::as<std::vector<std::string>>(Rcpp::rownames(susie_fit["pip"]));
        alpha = Rcpp::as<arma::mat>(susie_fit["alpha"]);
        prior_variance = Rcpp::as<std::vector<double>>(susie_fit["prior_variance"]);

        if (alpha.n_rows != prior_variance.size()) {
            Rcpp::stop("The number of rows in alpha must match the length of prior_variance.");
        }
    }

    std::vector<std::string> impute_qtn(const gsl_rng *r) {
        std::vector<std::string> qtn_names;

        for (arma::uword i = 0; i < alpha.n_rows; ++i) {
            std::vector<double> alpha_row(alpha.colptr(i), alpha.colptr(i) + alpha.n_cols);
            unsigned int *sample = new unsigned int[alpha_row.size()];

            gsl_ran_multinomial(r, alpha_row.size(), 1, alpha_row.data(), sample);

            for (size_t j = 0; j < alpha_row.size(); ++j) {
                if (sample[j] == 1) {
                    qtn_names.push_back(pip_names[j]);
                    break;
                }
            }
            delete[] sample;
        }

        return qtn_names;
    }
};

#endif // QTL_ENRICHMENT_HPP