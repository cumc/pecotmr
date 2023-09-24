#include <Rcpp.h>
using namespace Rcpp;

// Center and scale the columns of matrix X. Note that matrix X is
// modified in place. Input argument "a" should be the vector of
// column means, and input "b" should be the standard deviations.
//
// [[Rcpp::export]]
void scale_rcpp (NumericMatrix& X, const NumericVector& a,
		 const NumericVector& b) {
  for (int j = 0; j < X.ncol(); j++)
    for (int i = 0; i < X.nrow(); i++)
      X(i,j) = (X(i,j) - a(j)) / b(j);
}
