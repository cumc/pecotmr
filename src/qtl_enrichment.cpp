#include "qtl_enrichment.hpp"

// [[Rcpp::export]]
void qtl_enrichment(SEXP r_gwas_pip, SEXP r_susie_fit_list) {
  std::vector<double> gwas_pip = Rcpp::as<std::vector<double>>(r_gwas_pip);
  std::vector<std::string> gwas_pip_names = Rcpp::as<std::vector<std::string>>(Rcpp::rownames(r_gwas_pip));
  
  Rcpp::List susie_fit_list(r_susie_fit_list);
  std::vector<SuSiEFit> susie_fits;
  
  for (int i = 0; i < susie_fit_list.size(); ++i) {
    SuSiEFit susie_fit(susie_fit_list[i]);
    susie_fits.push_back(susie_fit);
  }

/*
  // For demonstration, print out the first element's pip values and names
  Rcpp::Rcout << "First susie_fit pip values: ";
  for (double value : susie_fits[0].pip) {
    Rcpp::Rcout << value << " ";
  }
  Rcpp::Rcout << "\nFirst susie_fit pip names: ";
  for (std::string name : susie_fits[0].pip_names) {
    Rcpp::Rcout << name << " ";
  }
  Rcpp::Rcout << std::endl;
*/
}