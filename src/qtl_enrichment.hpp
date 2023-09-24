#ifndef QTL_ENRICHMENT_HPP
#define QTL_ENRICHMENT_HPP

#include <RcppArmadillo.h>
#include <vector>
#include <string>

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
};

#endif // QTL_ENRICHMENT_HPP