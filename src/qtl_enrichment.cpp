#include "qtl_enrichment.hpp"

// [[Rcpp::export]]
Rcpp::List qtl_enrichment(
    SEXP r_gwas_pip, SEXP r_susie_fit_list,
    double pi_gwas = 0, double pi_eqtl = 0,
    int ImpN = 10, double prior_variance = 1, 
    int num_threads = 4)
{
    // Convert r_gwas_pip to C++ type
    std::vector<double> gwas_pip = Rcpp::as<std::vector<double>>(r_gwas_pip);
    std::vector<std::string> gwas_pip_names = Rcpp::as<std::vector<std::string>>(Rcpp::rownames(r_gwas_pip));

    // Convert r_susie_fit_list to C++ type
    Rcpp::List susie_fit_list(r_susie_fit_list);
    std::vector<SuSiEFit> susie_fits;

    for (int i = 0; i < susie_fit_list.size(); ++i) {
        SuSiEFit susie_fit(susie_fit_list[i]);
        susie_fits.push_back(susie_fit);
    }

    std::map<std::string, double> output = qtl_enrichment_workhorse(susie_fits, gwas_pip, gwas_pip_names, pi_gwas, pi_eqtl, ImpN, prior_variance, num_threads);

    // Convert std::map to Rcpp::List
    Rcpp::List output_list;
    for (auto const& element : output) {
        output_list[element.first] = element.second;
    }

    return output_list;
}