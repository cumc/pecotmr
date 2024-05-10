elnet <- function(lambda1, lambda2, diag, X, r, thr, x, yhat, trace, maxiter) {
    .Call('_lassosum_elnet', PACKAGE = 'lassosum', lambda1, lambda2, diag, X, r, thr, x, yhat, trace, maxiter)
}

repelnet <- function(lambda1, lambda2, diag, X, r, thr, x, yhat, trace, maxiter, startvec, endvec) {
    .Call('_lassosum_repelnet', PACKAGE = 'lassosum', lambda1, lambda2, diag, X, r, thr, x, yhat, trace, maxiter, startvec, endvec)
}

List runElnet(arma::vec& lambda, double shrink, const std::string fileName,
              arma::vec& r, int N, int P, 
			  arma::Col<int>& col_skip_pos, arma::Col<int>& col_skip, 
			  arma::Col<int>& keepbytes, arma::Col<int>& keepoffset, 
			  double thr, arma::vec& x, int trace, int maxiter, 
			  arma::Col<int>& startvec, arma::Col<int>& endvec) {
  // a) read bed file
  // b) standardize genotype matrix
  // c) multiply by constatant factor
  // d) perfrom elnet

  // Rcout << "ABC" << std::endl;

  int i,j;
  arma::mat genotypes = genotypeMatrix(fileName, N, P, col_skip_pos, col_skip, keepbytes,
                             keepoffset, 1);
  // Rcout << "DEF" << std::endl;

  if (genotypes.n_cols != r.n_elem) {
    throw std::runtime_error("Number of positions in reference file is not "
                             "equal the number of regression coefficients");
  }

  arma::vec sd = normalize(genotypes);

  genotypes *= sqrt(1.0 - shrink);

  arma::Col<int> conv(lambda.n_elem);
  int len = r.n_elem;

  arma::mat beta(len, lambda.n_elem);
  arma::mat pred(genotypes.n_rows, lambda.n_elem); pred.zeros();
  arma::vec out(lambda.n_elem);
  arma::vec loss(lambda.n_elem);
  arma::vec diag(r.n_elem); diag.fill(1.0 - shrink); 
  // Rcout << "HIJ" << std::endl;

	for(j=0; j < diag.n_elem; j++) {
		if(sd(j) == 0.0) diag(j) = 0.0;
	}
  // Rcout << "LMN" << std::endl;

  arma::vec fbeta(lambda.n_elem);
  arma::vec yhat(genotypes.n_rows);
  // yhat = genotypes * x;

  
  // Rcout << "Starting loop" << std::endl;
  for (i = 0; i < lambda.n_elem; ++i) {
    if (trace > 0)
      Rcout << "lambda: " << lambda(i) << "\n" << std::endl;
    out(i) =
        repelnet(lambda(i), shrink, diag,genotypes, r, thr, x, yhat, trace-1, maxiter, 
                 startvec, endvec);
    beta.col(i) = x;
    for(j=0; j < beta.n_rows; j++) {
      if(sd(j) == 0.0) beta(j,i)=beta(j,i) * shrink;
    }
    if (out(i) != 1) {
      throw std::runtime_error("Not converging.....");
    }
    pred.col(i) = yhat;
    loss(i) = arma::as_scalar(arma::sum(arma::pow(yhat, 2)) -
                              2.0 * arma::sum(x % r));
    fbeta(i) =
        arma::as_scalar(loss(i) + 2.0 * arma::sum(arma::abs(x)) * lambda(i) +
                        arma::sum(arma::pow(x, 2)) * shrink);
  }
  return List::create(Named("lambda") = lambda, 
					  Named("beta") = beta,
                      Named("conv") = out,
					  Named("pred") = pred,
                      Named("loss") = loss, 
					  Named("fbeta") = fbeta, 
					  Named("sd")= sd);