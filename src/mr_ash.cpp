#include "mr_ash.h"
// [[Rcpp::export]]
List rcpp_mr_ash_rss(const NumericVector& bhat, const NumericVector& shat, const NumericVector& z,
                     const NumericMatrix& R, double var_y, int n, double sigma2_e, const NumericVector& s0,
                     const NumericVector& w0, const NumericVector& mu1_init, double tol = 1e-8,
                     int max_iter = 1e5, bool update_w0 = true, bool update_sigma = true,
                     bool compute_ELBO = true, bool standardize = false) {
  // Convert input types
  vec bhat_vec = as<vec>(bhat);
  vec shat_vec = as<vec>(shat);
  vec z_vec = as<vec>(z);
  mat R_mat = as<mat>(R);
  vec s0_vec = as<vec>(s0);
  vec w0_vec = as<vec>(w0);
  vec mu1_init_vec = as<vec>(mu1_init);
  
  // Call the C++ function
  List result = mr_ash_rss(bhat_vec, shat_vec, z_vec, R_mat, var_y, n, sigma2_e, s0_vec, w0_vec,
                           mu1_init_vec, tol, max_iter, update_w0, update_sigma, compute_ELBO,
                           standardize);
  
  // Return the result
  return result;
}