#include <algorithm>
#include <cmath>
#include <thread>
#include <chrono>
#include <fstream>
#include <numeric>
#include <random>
#include <x86intrin.h>
#include "sse_mathfun.h"
#include "function_pool.h"
#include "sdpr_mcmc.h"

using namespace std::chrono;

using std::cout; using std::endl;
using std::thread; using std::ref;
using std::vector; using std::ofstream;
using std::string; using std::min;

#define square(x) ((x)*(x))

void MCMC_state::sample_sigma2() {
	std::gamma_distribution<double> dist;
	for (size_t i=1; i<M; i++) {
		double a = suff_stats[i] / 2.0 + a0k;
		double b = 1.0 / (sumsq[i] / 2.0 + b0k);
		dist = std::gamma_distribution<double>(a, 1.0/b);
		cluster_var[i] = 1.0 / dist(r);
		if (std::isinf(cluster_var[i])) {
			cluster_var[i] = 1e5;
			std::cerr << "Cluster variance is infinite." << std::endl;
		}
		else if (cluster_var[i] == 0) {
			cluster_var[i] = 1e-10;
			std::cerr << "Cluster variance is zero." << std::endl;
		}
	}
}

void MCMC_state::calc_b(size_t j, const mcmc_data &dat, const ldmat_data &ldmat_dat) {
	size_t start_i = dat.boundary[j].first;
	size_t end_i = dat.boundary[j].second;
	arma::vec b_j = b.subvec(start_i, end_i-1);
	arma::vec beta_j = beta.subvec(start_i, end_i-1);

	arma::vec diag = ldmat_dat.B[j].diag();

	// diag(B) * beta
	b_j = beta_j % diag;

	// eta^2 * (diag(B) * beta) - eta^2 * B * beta
	b_j = eta*eta * (b_j - ldmat_dat.B[j] * beta_j);

	// eta^2 * (diag(B) * beta) - eta^2 * B * beta + eta * A^T * beta_mrg
	b_j += eta * ldmat_dat.calc_b_tmp[j];
}

void MCMC_state::sample_assignment(size_t j, const mcmc_data &dat, const ldmat_data &ldmat_dat) {
	size_t start_i = dat.boundary[j].first;
	size_t end_i = dat.boundary[j].second;

	vector<vector<float> > prob(end_i-start_i, vector<float>(M));
	vector<vector<float> > tmp(end_i-start_i, vector<float>(M));
	vector<float> Bjj(end_i-start_i);
	vector<float> bj(end_i-start_i);
	vector<float> rnd(end_i-start_i);

	float max_elem, log_exp_sum = 0;

	std::uniform_real_distribution<float> unif(0.0, 1.0);

	for (size_t i=0; i<end_i-start_i; i++) {
		Bjj[i] = ldmat_dat.B[j](i, i);
		bj[i] = b(start_i+i);

		prob[i][0] = log_p[0];
		rnd[i] = unif(r);
	}

	// N = 1.0 after May 21 2021
	float C = pow(eta, 2.0) * N;

	// auto vectorized
	for (size_t i=0; i<end_i-start_i; i++) {
		for (size_t k=1; k<M; k++) {
			prob[i][k] = C * Bjj[i] * cluster_var[k] + 1;
		}
	}

	// unable to auto vectorize due to log
	// explicitly using SSE
	__m128 _v, _m;
	for (size_t i=0; i<end_i-start_i; i++) {
		size_t k = 1;
		for (; k<M; k+=4) { // require M >= 4
			_v = log_ps(_mm_loadu_ps(&prob[i][k]));
			_mm_storeu_ps(&tmp[i][k], _v);
		}

		for (; k<M; k++) {
			tmp[i][k] = logf(prob[i][k]);
		}
	}

	// auto vectorized
	for (size_t i=0; i<end_i-start_i; i++) {
		for (size_t k=1; k<M; k++) {
			prob[i][k] = -0.5*tmp[i][k] + log_p[k] + square(N*bj[i]) * cluster_var[k] / (2*prob[i][k]);
		}
	}

	for (size_t i=0; i<end_i-start_i; i++) {
		// SSE version to find max
		_v = _mm_loadu_ps(&prob[i][0]);
		size_t k = 4;
		for (; k<M; k+=4) {
			_v = _mm_max_ps(_v, _mm_loadu_ps(&prob[i][k]));
		}

		for (size_t m=0; m<3; m++) {
			_v = _mm_max_ps(_v, _mm_shuffle_ps(_v, _v, 0x93));
		}

		_mm_store_ss(&max_elem, _v);

		for (; k<M; k++) {
			max_elem = (max_elem > prob[i][k]) ? max_elem : prob[i][k];
		}

		// SSE version log exp sum
		_m = _mm_load1_ps(&max_elem);
		_v = exp_ps(_mm_sub_ps(_mm_loadu_ps(&prob[i][0]), _m));

		k = 4;
		for (; k<M; k+=4) {
			_v = _mm_add_ps(_v, exp_ps(_mm_sub_ps(_mm_loadu_ps(&prob[i][k]), _m)));
		}

		_v = _mm_hadd_ps(_v, _v);
		_v = _mm_hadd_ps(_v, _v);
		_mm_store_ss(&log_exp_sum, _v);

		for (; k<M; k++) {
			log_exp_sum += expf(prob[i][k] - max_elem);
		}

		log_exp_sum = max_elem + logf(log_exp_sum);

		cls_assgn[i+start_i] = M-1;
		for (size_t k=0; k<M-1; k++) {
			rnd[i] -= expf(prob[i][k] - log_exp_sum);
			if (rnd[i] < 0) {
				cls_assgn[i+start_i] = k;
				break;
			}
		}
	}
}

void MCMC_state::update_suffstats() {
	std::fill(suff_stats.begin(), suff_stats.end(), 0.0);
	std::fill(sumsq.begin(), sumsq.end(), 0.0);
	for (size_t i=0; i<n_snp; i++) {
		suff_stats[cls_assgn[i]]++;
		double tmp = beta(i);
		sumsq[cls_assgn[i]] += square(tmp);
	}
}

void MCMC_state::sample_V() {
	vector<double> a(M-1);

	a[M-2] = suff_stats[M-1];
	for (int i=M-3; i>=0; i--) {
		a[i] = suff_stats[i+1] + a[i+1];
	}

	for (size_t i=0; i<M-1; i++) {
		beta_distribution dist(1 + suff_stats[i], alpha + a[i]);
		V[i] = dist(r);
	}
	V[M-1] = 1;
}

void MCMC_state::update_p() {
	vector<double> cumprod(M-1);

	cumprod[0] = 1 - V[0];

	for (size_t i=1; i<M-1; i++) {
		cumprod[i] = cumprod[i-1] * (1 - V[i]);

		if (V[i] == 1) {
			std::fill(cumprod.begin()+i+1, cumprod.end(), 0.0);
			break;
		}
	}

	p[0] = V[0];
	for (size_t i=1; i<M-1; i++) {
		p[i] = cumprod[i-1] * V[i];
	}

	double sum = std::accumulate(p.begin(), p.end()-1, 0.0);
	if (1 - sum > 0) {
		p[M-1] = 1 - sum;
	}
	else {
		p[M-1] = 0;
	}

	for (size_t i=0; i<M; i++) {
		log_p[i] = logf(p[i] + 1e-40);
	}
}

void MCMC_state::sample_alpha() {
	double sum = 0, m = 0;
	for (size_t i=0; i<M; i++) {
		if (V[i] != 1) {
			sum += log(1 - V[i]);
			m++;
		}
	}

	if (m == 0) m = 1;

	std::gamma_distribution<double> dist(0.1+m-1, 1.0/(0.1-sum));
	alpha = dist(r);
}

void MCMC_state::sample_beta(size_t j, const mcmc_data &dat, ldmat_data &ldmat_dat) {
	size_t start_i = dat.boundary[j].first;
	size_t end_i = dat.boundary[j].second;

	vector<size_t> causal_list;
	for (size_t i=start_i; i<end_i; i++) {
		if (cls_assgn[i] != 0) {
			causal_list.push_back(i);
		}
	}

	arma::vec beta_j = beta.subvec(start_i, end_i-1);

	beta_j.zeros();

	if (causal_list.size() == 0) {
		ldmat_dat.num[j] = 0;
		ldmat_dat.denom[j] = 0;
		return;
	}
	else if (causal_list.size() == 1) {
		double var_k = cluster_var[cls_assgn[causal_list[0]]];
		double bj = b(causal_list[0]);
		double Bjj = ldmat_dat.B[j](causal_list[0]-start_i, causal_list[0]-start_i);
		double C = var_k / (N*var_k*square(eta)*Bjj + 1.0);
		std::normal_distribution<double> dist(0.0, sqrt(C));
		double rv = dist(r) + C*N*bj;
		beta_j(causal_list[0]-start_i) = rv;
		ldmat_dat.num[j] = bj*rv;
		ldmat_dat.denom[j] = square(rv)*Bjj;
		return;
	}

	arma::vec A_vec(causal_list.size());
	arma::vec A_vec2(causal_list.size());

	double C = square(eta)*N;

	arma::mat B(causal_list.size(), causal_list.size());

	for (size_t i=0; i<causal_list.size(); i++) {
		A_vec(i) = N*eta*ldmat_dat.calc_b_tmp[j](causal_list[i]-start_i);

		for (size_t k=0; k<causal_list.size(); k++) {
			if (i != k) {
				B(i, k) = C * ldmat_dat.B[j](causal_list[i]-start_i, causal_list[k]-start_i);
			}
			else {
				B(i, k) = C * ldmat_dat.B[j](causal_list[i]-start_i, causal_list[i]-start_i) + 1.0/cluster_var[cls_assgn[causal_list[i]]];
			}
		}
	}

	A_vec2 = A_vec;

	arma::vec beta_c(causal_list.size());
	std::normal_distribution<double> dist(0.0, 1.0);
	for (size_t i=0; i<causal_list.size(); i++) {
		beta_c(i) = dist(r);
	}

	// (N B_gamma + \Sigma_0^-1) = L L^T
	arma::mat L = arma::chol(B, "lower");

	// \mu = L^{-1} A_vec
	arma::vec mu = arma::solve(arma::trimatl(L), A_vec);

	// N(\mu, I)
	beta_c += mu;

	// X ~ N(\mu, I), L^{-T} X ~ N( L^{-T} \mu, (L L^T)^{-1} )
	beta_c = arma::solve(arma::trimatu(L.t()), beta_c);

	// compute eta related terms
	for (size_t i=0; i<causal_list.size(); i++) {
		B(i, i) = C*ldmat_dat.B[j](causal_list[i]-start_i, causal_list[i]-start_i);
	}

	ldmat_dat.num[j] = arma::dot(A_vec2, beta_c);
	arma::vec tmp = B * beta_c;
	ldmat_dat.denom[j] = arma::dot(beta_c, tmp);
	ldmat_dat.denom[j] /= square(eta);
	ldmat_dat.num[j] /= eta;

	for (size_t i=0; i<causal_list.size(); i++) {
		beta_j(causal_list[i]-start_i) = beta_c(i);
	}
}

void MCMC_state::compute_h2(const mcmc_data &dat) {
	double h2_tmp = 0;
	h2 = 0;
	for (size_t j=0; j<dat.ref_ld_mat.size(); j++) {
		size_t start_i = dat.boundary[j].first;
		size_t end_i = dat.boundary[j].second;
		arma::vec tmp(end_i-start_i);
		arma::vec beta_j = beta.subvec(start_i, end_i-1);
		tmp = dat.ref_ld_mat[j] * beta_j;
		h2_tmp = arma::dot(tmp, beta_j);
		h2 += h2_tmp;
	}
}

void MCMC_state::sample_eta(const ldmat_data &ldmat_dat) {
	double num_sum = std::accumulate(ldmat_dat.num.begin(), ldmat_dat.num.end(), 0.0);
	double denom_sum = std::accumulate(ldmat_dat.denom.begin(), ldmat_dat.denom.end(), 0.0);
	denom_sum += 1e-6;

	std::normal_distribution<double> dist(num_sum/denom_sum, sqrt(1.0/denom_sum));
	eta = dist(r);
}

void solve_ldmat(const mcmc_data &dat, ldmat_data &ldmat_dat, const double a, unsigned sz, int opt_llk) {
	for (size_t i=0; i<dat.ref_ld_mat.size(); i++) {
		size_t size = dat.boundary[i].second - dat.boundary[i].first;
		arma::mat A = dat.ref_ld_mat[i];
		arma::mat B = dat.ref_ld_mat[i];
		arma::mat L = dat.ref_ld_mat[i];

		if (opt_llk == 1) {
			// (R + aNI) / N A = R via cholesky decomp
			// Changed May 21 2021 to divide by N
			// replace aN with a
			B.diag() += a;
		}
		else {
			// R_ij N_s,ij / N_i N_j
			// Added May 24 2021
			for (size_t j=0; j<size; j++) {
				for (size_t k=0; k<size; k++) {
					double tmp = B(j, k);
					size_t idx1 = j + dat.boundary[i].first;
					size_t idx2 = k + dat.boundary[i].first;
					if ((dat.array[idx1] == 1 && dat.array[idx2] == 2) || (dat.array[idx1] == 2 && dat.array[idx2] == 1)) {
						tmp = 0;
					}
					else {
						tmp *= min(dat.sz[idx1], dat.sz[idx2]) / (1.1 * dat.sz[idx1] * dat.sz[idx2]);
					}
					B(j, k) = tmp;
				}
			}
			// force positive definite
			// B = Q \Lambda Q^T
			arma::vec eval;
			arma::mat evec;
			arma::eig_sym(eval, evec, B);
			double eval_min = eval.min();

			// restore lower half of B
			B = arma::symmatu(B);

			// if min eigen value < 0, add -1.1 * eval to diagonal
			for (size_t j=0; j<size; j++) {
				if (eval_min < 0) {
					B(j, j) = 1.0/dat.sz[j+dat.boundary[i].first] - 1.1*eval_min;
				}
				else {
					B(j, j) = 1.0/dat.sz[j+dat.boundary[i].first];
				}
			}
		}

		arma::mat L_B = arma::chol(B, "lower");
		A = arma::solve(arma::trimatl(L_B), A);
		A = arma::solve(arma::trimatu(L_B.t()), A);

		if (opt_llk == 1) {
			A *= sz;
		}

		B = L * A;
		L %= L;

		arma::vec beta_mrg(size);
		for (size_t j=0; j<size; j++) {
			beta_mrg(j) = dat.beta_mrg[j+dat.boundary[i].first];
		}
		arma::vec b_tmp = A.t() * beta_mrg;

		ldmat_dat.A.push_back(A);
		ldmat_dat.B.push_back(B);
		ldmat_dat.L.push_back(L);
		ldmat_dat.calc_b_tmp.push_back(b_tmp);
		ldmat_dat.beta_mrg.push_back(beta_mrg);
		ldmat_dat.denom.push_back(0);
		ldmat_dat.num.push_back(0);
	}
}

std::unordered_map<std::string, arma::vec> mcmc(
	mcmc_data& data,
	unsigned   sz,
	double     a = 0.1,
	double     c = 1.0,
	size_t     M = 1000,
	double     a0k = 0.5,
	double     b0k = 0.5,
	int        iter = 1000,
	int        burn = 200,
	int        thin = 5,
	unsigned   n_threads = 1,
	int        opt_llk = 1,
	bool       verbose = true
	) {

	int n_pst = (iter-burn) / thin;

	ldmat_data ldmat_dat;

	MCMC_state state(data.beta_mrg.size(), M, a0k, b0k, sz);

	for (size_t i=0; i<data.beta_mrg.size(); i++) {
		data.beta_mrg[i] /= c;
	}

	MCMC_samples samples(data.beta_mrg.size());

	solve_ldmat(data, ldmat_dat, a, sz, opt_llk);
	state.update_suffstats();

	Function_pool func_pool(n_threads);

	for (int j=1; j<iter+1; j++) {
		state.sample_sigma2();

		for (size_t i=0; i<data.ref_ld_mat.size(); i++) {
			state.calc_b(i, data, ldmat_dat);
		}

		for (size_t i=0; i<data.ref_ld_mat.size(); i++) {
			func_pool.push(std::bind(&MCMC_state::sample_assignment, &state, i, ref(data), ref(ldmat_dat)));
		}

		func_pool.waitFinished();

		state.update_suffstats();

		state.sample_V();
		state.update_p();
		state.sample_alpha();

		for (size_t i=0; i<data.ref_ld_mat.size(); i++) {
			state.sample_beta(i, data, ldmat_dat);
		}

		state.sample_eta(ldmat_dat);

		if ((j>burn) && (j%thin == 0)) {
			state.compute_h2(data);
			samples.h2 += state.h2*square(state.eta) / n_pst;
			samples.beta += state.eta/n_pst * state.beta;
		}

		if (verbose && j % 100 == 0) {
			state.compute_h2(data);
			cout << j << " iter. h2: " << state.h2*square(state.eta) << " max beta: " << arma::max(state.beta)*state.eta << endl;
		}
	}

	if (verbose) {
		cout << "h2: " << samples.h2 << " max: " << arma::max(samples.beta) << endl;
	}

	std::unordered_map<std::string, arma::vec> results;
	results["beta"] = samples.beta;
	results["h2"] = arma::vec(1, arma::fill::value(samples.h2));

	return results;
}