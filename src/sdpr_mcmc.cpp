// sdpr_mcmc.cpp -- SDPR MCMC sampler using Armadillo
//
// Clean Armadillo port of the original SDPR by Zhou et al.
// (https://github.com/eldronzhou/SDPR), which used GSL + x86 SSE intrinsics.
//
// Translation strategy:
//   GSL gsl_blas_dsymv  -> M * v, with M symmetrized once when it is built
//   GSL gsl_blas_dgemv  -> A.t() * v
//   GSL gsl_blas_dgemm  -> A * B
//   GSL gsl_blas_dtrsv  -> arma::solve(trimatl(L), v)
//   GSL gsl_blas_dtrsm  -> arma::solve(trimatl(L), A)
//   GSL gsl_blas_ddot   -> arma::dot(x, y)
//   GSL gsl_blas_daxpy  -> y += alpha * x
//   GSL gsl_linalg_cholesky_decomp1 -> arma::chol(A, "lower")
//   GSL gsl_ran_gamma   -> std::gamma_distribution
//   GSL gsl_ran_beta    -> beta_distribution (ratio of gammas)
//   GSL gsl_ran_ugaussian -> std::normal_distribution(0,1)
//   GSL gsl_rng_uniform -> std::uniform_real_distribution(0,1)
//
// The original SSE intrinsics (log_ps, exp_ps, _mm_max_ps, _mm_hadd_ps)
// in sample_assignment() are replaced with Armadillo vectorized operations
// (arma::log, arma::exp, arma::max, arma::accu) which delegate to
// platform-optimal SIMD (NEON on ARM, SSE/AVX on x86) through the
// underlying BLAS/compiler auto-vectorization.
//
// The sampler reads each LD block through a block type (sdpr_mcmc.h): the
// dense matrices SDPR itself precomputes (sdpr_dense_blocks()), or their
// low-rank form from a factor of the LD (sdpr_factor_block()).

#include <cpp11.hpp>
#include <cpp11armadillo.hpp>
#include <algorithm>
#include <cmath>
#include <thread>
#include <chrono>
#include <fstream>
#include <numeric>
#include <random>
#include "function_pool.h"
#include "sdpr_mcmc.h"

using std::cout; using std::endl;
using std::thread; using std::ref;
using std::vector; using std::ofstream;
using std::string; using std::min;

// A function rather than the original macro, which would also expand inside
// arma::square().
static inline double square(double x) { return x * x; }

// ---------------------------------------------------------------------------
// sample_sigma2: Sample cluster variances from inverse-gamma posterior
//   var_k ~ InvGamma(suff_stats[k]/2 + a0k, sumsq[k]/2 + b0k)
// Original: mcmc.cpp lines 25-40
// ---------------------------------------------------------------------------
void MCMC_state::sample_sigma2() {
    std::gamma_distribution<double> dist;
    for (size_t i = 1; i < M; i++) {
        double a = suff_stats[i] / 2.0 + a0k;
        double b = 1.0 / (sumsq[i] / 2.0 + b0k);
        dist = std::gamma_distribution<double>(a, b);
        cluster_var[i] = 1.0 / dist(r);
        if (std::isinf(cluster_var[i])) {
            cluster_var[i] = 1e5;
            REprintf("Cluster variance is infinite.\n");
        }
        else if (cluster_var[i] == 0) {
            cluster_var[i] = 1e-10;
            REprintf("Cluster variance is zero.\n");
        }
    }
}

// ---------------------------------------------------------------------------
// calc_b: Compute the b vector for block j
//   b_j = eta^2 * (diag(B)*beta - B*beta) + eta * A^T * beta_mrg
// Original: mcmc.cpp lines 42-62
//   gsl_blas_dsymv(CblasUpper, -eta*eta, B, beta, eta*eta, b)
//   gsl_blas_daxpy(eta, calc_b_tmp, b)
// ---------------------------------------------------------------------------
template <class Block>
void MCMC_state::calc_b(size_t j, const mcmc_data &dat,
                        const std::vector<Block> &blocks) {
    size_t start_i = dat.boundary[j].first;
    size_t end_i   = dat.boundary[j].second;

    arma::vec beta_j = beta.subvec(start_i, end_i - 1);

    // Original GSL: b = eta^2 * diag(B)*beta;
    //               b = eta^2*b - eta^2 * B*beta  (via dsymv with alpha=-eta^2, beta=eta^2)
    //               b += eta * calc_b_tmp
    arma::vec b_j = eta * eta * (blocks[j].diag_B % beta_j
                    - blocks[j].B_times(beta_j))
                    + eta * blocks[j].calc_b_tmp;

    b.subvec(start_i, end_i - 1) = b_j;
}

// ---------------------------------------------------------------------------
// sample_assignment: Sample cluster assignments for each SNP in block j
//   For each SNP i, compute log P(z_i = k) for k = 0..M-1, then sample.
//
// Original: mcmc.cpp lines 64-165
//   Used x86 SSE intrinsics (log_ps, exp_ps, _mm_max_ps, _mm_hadd_ps)
//   to vectorize log/exp over M clusters. We use Armadillo vectorized
//   ops instead, which auto-vectorize via NEON (ARM) or SSE/AVX (x86).
//
// Math for cluster k >= 1:
//   C_k = eta^2 * N * B[i,i] * var_k + 1
//   log P(z_i=k) = -0.5*log(C_k) + log(p_k) + (N*b_i)^2 * var_k / (2*C_k)
// ---------------------------------------------------------------------------
template <class Block>
void MCMC_state::sample_assignment(size_t j, const mcmc_data &dat,
                                   const std::vector<Block> &blocks) {
    size_t start_i   = dat.boundary[j].first;
    size_t end_i     = dat.boundary[j].second;
    size_t n_snp_blk = end_i - start_i;

    std::uniform_real_distribution<float> unif(0.0f, 1.0f);

    // N = 1.0 after May 21 2021 (absorbed into A, B in solve_ldmat)
    float C = static_cast<float>(eta * eta * N);

    // Pre-convert cluster variances and log-probabilities to float vectors
    // for Armadillo vectorized ops (matching original's float precision)
    arma::fvec cvar(M);
    arma::fvec log_p_fvec(M);
    for (size_t k = 0; k < M; k++) {
        cvar(k)       = static_cast<float>(cluster_var[k]);
        log_p_fvec(k) = static_cast<float>(log_p[k]);
    }

    for (size_t i = 0; i < n_snp_blk; i++) {
        float Bjj_i = static_cast<float>(blocks[j].diag_B(i));
        float bj_i  = static_cast<float>(b(start_i + i));
        float rnd_i = unif(r);

        // prob[0] = log_p[0] (null cluster)
        // prob[k] for k>=1: see math above
        arma::fvec prob(M);
        prob(0) = log_p_fvec(0);

        // Vectorized: Ck = C * Bjj * cvar[1..M-1] + 1
        arma::fvec Ck = C * Bjj_i * cvar.subvec(1, M - 1) + 1.0f;

        // prob[k] = -0.5*log(Ck) + log_p[k] + (N*bj)^2 * var_k / (2*Ck)
        // Original: tmp[k] = log(prob[k]); prob[k] = -0.5*tmp[k] + log_p[k] + ...
        prob.subvec(1, M - 1) = -0.5f * arma::log(Ck) + log_p_fvec.subvec(1, M - 1)
            + square(N * bj_i) * cvar.subvec(1, M - 1) / (2.0f * Ck);

        // Log-sum-exp for numerical stability (replaces SSE _mm_max_ps + exp_ps + _mm_hadd_ps)
        float max_elem    = prob.max();
        float log_exp_sum = max_elem
            + std::logf(arma::accu(arma::exp(prob - max_elem)));

        // Categorical sampling via inverse CDF
        // Original: mcmc.cpp lines 155-163
        cls_assgn[i + start_i] = M - 1;
        for (size_t k = 0; k < M - 1; k++) {
            rnd_i -= std::expf(prob(k) - log_exp_sum);
            if (rnd_i < 0) {
                cls_assgn[i + start_i] = k;
                break;
            }
        }
    }
}

// ---------------------------------------------------------------------------
// update_suffstats: Count SNPs per cluster and sum of squared effects
// Original: mcmc.cpp lines 167-175
// ---------------------------------------------------------------------------
void MCMC_state::update_suffstats() {
    std::fill(suff_stats.begin(), suff_stats.end(), 0);
    std::fill(sumsq.begin(), sumsq.end(), 0.0);
    for (size_t i = 0; i < n_snp; i++) {
        suff_stats[cls_assgn[i]]++;
        double tmp = beta(i);
        sumsq[cls_assgn[i]] += square(tmp);
    }
}

// ---------------------------------------------------------------------------
// sample_V: Sample stick-breaking fractions from Beta posterior
//   V[k] ~ Beta(1 + suff_stats[k], alpha + sum_{k'>k} suff_stats[k'])
// Original: mcmc.cpp lines 177-189
// ---------------------------------------------------------------------------
void MCMC_state::sample_V() {
    vector<double> a(M - 1);
    a[M - 2] = suff_stats[M - 1];
    for (int i = M - 3; i >= 0; i--) {
        a[i] = suff_stats[i + 1] + a[i + 1];
    }

    for (size_t i = 0; i < M - 1; i++) {
        beta_distribution dist(1 + suff_stats[i], alpha + a[i]);
        V[i] = dist(r);
    }
    V[M - 1] = 1;
}

// ---------------------------------------------------------------------------
// update_p: Convert stick-breaking fractions V to cluster probabilities p
//   p[0] = V[0]; p[k] = V[k] * prod_{j<k}(1-V[j])
// Original: mcmc.cpp lines 191-222
// ---------------------------------------------------------------------------
void MCMC_state::update_p() {
    vector<double> cumprod(M - 1);
    cumprod[0] = 1 - V[0];
    for (size_t i = 1; i < M - 1; i++) {
        cumprod[i] = cumprod[i - 1] * (1 - V[i]);
        if (V[i] == 1) {
            std::fill(cumprod.begin() + i + 1, cumprod.end(), 0.0);
            break;
        }
    }

    p[0] = V[0];
    for (size_t i = 1; i < M - 1; i++) {
        p[i] = cumprod[i - 1] * V[i];
    }

    double sum = std::accumulate(p.begin(), p.end() - 1, 0.0);
    p[M - 1] = (1 - sum > 0) ? (1 - sum) : 0;

    for (size_t i = 0; i < M; i++) {
        log_p[i] = std::logf(static_cast<float>(p[i]) + 1e-40f);
    }
}

// ---------------------------------------------------------------------------
// sample_alpha: Sample Dirichlet process concentration parameter
//   alpha ~ Gamma(0.1 + m - 1, 1/(0.1 - sum(log(1-V))))
// Original: mcmc.cpp lines 224-237
// ---------------------------------------------------------------------------
void MCMC_state::sample_alpha() {
    double sum = 0, m = 0;
    for (size_t i = 0; i < M; i++) {
        if (V[i] != 1) {
            sum += std::log(1 - V[i]);
            m++;
        }
    }
    if (m == 0) m = 1;

    std::gamma_distribution<double> dist(0.1 + m - 1, 1.0 / (0.1 - sum));
    alpha = dist(r);
}

// ---------------------------------------------------------------------------
// sample_beta: Sample effect sizes for causal SNPs in block j
//   Precision = eta^2*N*B_gamma + Sigma_0^{-1}
//   Mean via Cholesky: L*L^T = Precision; mu = L^{-T}*L^{-1}*A_vec
//   Sample: z ~ N(0,I); beta = L^{-T}*(L^{-1}*A_vec + z)
//
// Original: mcmc.cpp lines 239-329
//   Uses gsl_linalg_cholesky_decomp1, gsl_blas_dtrsv for forward/back-solve
// ---------------------------------------------------------------------------
template <class Block>
void MCMC_state::sample_beta(size_t j, const mcmc_data &dat,
                             const std::vector<Block> &blocks) {
    size_t start_i = dat.boundary[j].first;
    size_t end_i   = dat.boundary[j].second;

    // Build causal list (SNPs not in null cluster)
    vector<size_t> causal_list;
    for (size_t i = start_i; i < end_i; i++) {
        if (cls_assgn[i] != 0) {
            causal_list.push_back(i);
        }
    }

    // Zero out this block's betas
    beta.subvec(start_i, end_i - 1).zeros();

    if (causal_list.empty()) {
        num[j]   = 0;
        denom[j] = 0;
        return;
    }

    // Single causal SNP: closed-form sampling
    // Original: mcmc.cpp lines 259-270
    if (causal_list.size() == 1) {
        double var_k = cluster_var[cls_assgn[causal_list[0]]];
        double bj    = b(causal_list[0]);
        double Bjj   = blocks[j].diag_B(causal_list[0] - start_i);
        // C = var_k / (N * var_k * eta^2 * Bjj + 1)
        double C_val = var_k / (N * var_k * square(eta) * Bjj + 1.0);
        std::normal_distribution<double> dist(0.0, std::sqrt(C_val));
        double rv = dist(r) + C_val * N * bj;
        beta(causal_list[0]) = rv;
        num[j]   = bj * rv;
        denom[j] = square(rv) * Bjj;
        return;
    }

    // Multiple causal SNPs: multivariate normal sampling via Cholesky
    size_t nc = causal_list.size();

    // A_vec = N * eta * A^T * beta_mrg (restricted to causal indices)
    // Original: mcmc.cpp line 279: N*eta*gsl_vector_get(calc_b_tmp, idx)
    arma::vec A_vec(nc);
    double C_coeff = square(eta) * N;  // Original: C = square(eta)*N (line 273)

    // Build precision matrix B_gamma and A_vec
    // Original: mcmc.cpp lines 278-289
    arma::uvec causal_idx(nc);
    for (size_t i = 0; i < nc; i++) {
        causal_idx(i) = causal_list[i] - start_i;
        A_vec(i) = N * eta * blocks[j].calc_b_tmp(causal_idx(i));
    }
    // Off-diagonal: eta^2 * N * B[i,k]; diagonal: eta^2 * N * B[i,i] + 1/var_k
    arma::mat B_gamma = C_coeff * blocks[j].B_sub(causal_idx);
    for (size_t i = 0; i < nc; i++) {
        B_gamma(i, i) += 1.0 / cluster_var[cls_assgn[causal_list[i]]];
    }

    arma::vec A_vec2 = A_vec;  // Save for eta computation

    // Sample z ~ N(0, I)
    arma::vec beta_c(nc);
    std::normal_distribution<double> dist(0.0, 1.0);
    for (size_t i = 0; i < nc; i++) {
        beta_c(i) = dist(r);
    }

    // Cholesky: B_gamma = L * L^T
    // Original: gsl_linalg_cholesky_decomp1(&B.matrix)
    arma::mat L = arma::chol(B_gamma, "lower");

    // mu = L^{-1} * A_vec  (forward-solve)
    // Original: gsl_blas_dtrsv(CblasLower, CblasNoTrans, CblasNonUnit, &B, A_vec)
    A_vec = arma::solve(arma::trimatl(L), A_vec);

    // beta_c = mu + z ~ N(mu, I)
    beta_c += A_vec;

    // beta_c = L^{-T} * beta_c ~ N(L^{-T}*mu, (L*L^T)^{-1})
    // Original: gsl_blas_dtrsv(CblasLower, CblasTrans, CblasNonUnit, &B, beta_c)
    beta_c = arma::solve(arma::trimatu(L.t()), beta_c);

    // Compute eta-related terms for sample_eta()
    // Original: mcmc.cpp lines 312-321
    // Restore diagonal of B_gamma to eta^2 * N * B[i,i] (without 1/var_k)
    for (size_t i = 0; i < nc; i++) {
        B_gamma(i, i) = C_coeff * blocks[j].diag_B(causal_idx(i));
    }

    // num = A_vec2^T * beta_c / eta
    num[j] = arma::dot(A_vec2, beta_c) / eta;

    // denom = beta_c^T * B_gamma * beta_c / eta^2
    // Original: gsl_blas_dsymv(CblasUpper, 1.0, &B, beta_c, 0, A_vec)
    arma::vec tmp = B_gamma * beta_c;
    denom[j] = arma::dot(beta_c, tmp) / square(eta);

    // Write sampled betas back
    for (size_t i = 0; i < nc; i++) {
        beta(causal_list[i]) = beta_c(i);
    }
}

// ---------------------------------------------------------------------------
// compute_h2: Compute heritability h2 = sum_j beta_j^T * R_j * beta_j
// Original: mcmc.cpp lines 331-343
// ---------------------------------------------------------------------------
template <class Block>
void MCMC_state::compute_h2(const mcmc_data &dat,
                            const std::vector<Block> &blocks) {
    h2 = 0;
    for (size_t j = 0; j < blocks.size(); j++) {
        size_t start_i = dat.boundary[j].first;
        size_t end_i   = dat.boundary[j].second;
        arma::vec beta_j = beta.subvec(start_i, end_i - 1);
        // Original: gsl_blas_dsymv(CblasUpper, 1.0, ref_ld_mat, beta, 0, tmp)
        h2 += blocks[j].R_quad(beta_j);
    }
}

// ---------------------------------------------------------------------------
// sample_eta: Sample global scaling parameter
//   eta ~ N(num_sum/denom_sum, 1/denom_sum)
// Original: mcmc.cpp lines 345-352
// ---------------------------------------------------------------------------
void MCMC_state::sample_eta() {
    double num_sum   = std::accumulate(num.begin(), num.end(), 0.0);
    double denom_sum = std::accumulate(denom.begin(), denom.end(), 0.0);
    denom_sum += 1e-6;

    std::normal_distribution<double> dist(num_sum / denom_sum,
                                          std::sqrt(1.0 / denom_sum));
    eta = dist(r);
}

// ---------------------------------------------------------------------------
// sdpr_llk2_matrix: The optLlk = 2 matrix B2 from one block's LD
//   Each pair's LD is scaled by min(N_i, N_j) / (1.1 * N_i * N_j) and zeroed
//   between arrays 1 and 2; the diagonal is then set to 1 / N_i, raised by
//   1.1 times the scaled matrix's most negative eigenvalue, if any, to keep
//   B2 positive definite.
//
// Original: mcmc.cpp lines 368-402
// ---------------------------------------------------------------------------
static arma::mat sdpr_llk2_matrix(arma::mat B, const mcmc_data &dat,
                                  size_t start) {
    size_t blk_size = B.n_rows;
    for (size_t j = 0; j < blk_size; j++) {
        for (size_t k = 0; k < blk_size; k++) {
            size_t idx1 = j + start;
            size_t idx2 = k + start;
            if ((dat.array[idx1] == 1 && dat.array[idx2] == 2) ||
                (dat.array[idx1] == 2 && dat.array[idx2] == 1)) {
                B(j, k) = 0;
            } else {
                B(j, k) *= min(dat.sz[idx1], dat.sz[idx2])
                           / (1.1 * dat.sz[idx1] * dat.sz[idx2]);
            }
        }
    }

    // Force positive definite via eigenvalue correction. Only the smallest
    // eigenvalue is read, so no eigenvectors are computed.
    double eval_min = arma::eig_sym(B).min();

    // Restore exact symmetry before the diagonal is reset
    B = arma::symmatu(B);

    // Add to diagonal to ensure PD
    for (size_t j = 0; j < blk_size; j++) {
        double diag_val = 1.0 / dat.sz[j + start];
        if (eval_min < 0) diag_val -= 1.1 * eval_min;
        B(j, j) = diag_val;
    }
    return B;
}

// ---------------------------------------------------------------------------
// sdpr_dense_blocks: Precompute LD-derived matrices B and A' beta_mrg for
// each block
//   opt_llk=1: B0 = R + aI; A = B0^{-1} R * sz; B = R * A
//   opt_llk=2: B0 = sdpr_llk2_matrix(R); A = B0^{-1} R; B = R * A
//
// Original: mcmc.cpp lines 354-431 (solve_ldmat)
//   Uses gsl_linalg_cholesky_decomp1, gsl_blas_dtrsm for matrix solves. The
//   original also kept A and R .* R per block, which nothing reads after
//   this, so neither is kept here.
// ---------------------------------------------------------------------------
std::vector<SdprDenseBlock> sdpr_dense_blocks(const mcmc_data &dat,
                                              std::vector<arma::mat> ref_ld_mat,
                                              double a, unsigned sz,
                                              int opt_llk) {
    std::vector<SdprDenseBlock> blocks;
    for (size_t i = 0; i < ref_ld_mat.size(); i++) {
        size_t start    = dat.boundary[i].first;
        size_t blk_size = dat.boundary[i].second - start;
        const arma::mat &R = ref_ld_mat[i];

        arma::mat B0;
        if (opt_llk == 1) {
            // B0 = R + a*I  (diagonal shrinkage)
            // Original: gsl_vector_add_constant(&diag.vector, a)
            B0 = R;
            B0.diag() += a;
        } else {
            B0 = sdpr_llk2_matrix(R, dat, start);
        }

        // Cholesky: B0 = L_chol * L_chol^T
        // Original: gsl_linalg_cholesky_decomp1(B)
        arma::mat L_chol = arma::chol(B0, "lower");

        // A = B0^{-1} * R  via forward then back-solve
        // Original: gsl_blas_dtrsm(Left, Lower, NoTrans, NonUnit, 1.0, B, A)
        //           gsl_blas_dtrsm(Left, Lower, Trans,   NonUnit, 1.0, B, A)
        arma::mat A = arma::solve(arma::trimatl(L_chol), R);
        A = arma::solve(arma::trimatu(L_chol.t()), A);

        // Scale by sample size for opt_llk=1
        if (opt_llk == 1) {
            A *= sz;
        }

        // Precompute A^T * beta_mrg for calc_b()
        // Original: gsl_blas_dgemv(CblasTrans, 1.0, A, beta_mrg, 0, b_tmp)
        arma::vec beta_mrg(blk_size);
        for (size_t j = 0; j < blk_size; j++) {
            beta_mrg(j) = dat.beta_mrg[j + start];
        }

        SdprDenseBlock blk;
        // B = R * A. The sampler reads it as symmetric (the original through
        // dsymv), so it is symmetrized once here rather than on every product.
        // Original: gsl_blas_dgemm(NoTrans, NoTrans, 1.0, L, A, 0, B)
        blk.B          = arma::symmatu(R * A);
        blk.diag_B     = blk.B.diag();
        blk.calc_b_tmp = A.t() * beta_mrg;
        blk.R          = arma::symmatu(R);
        blocks.push_back(std::move(blk));
    }
    return blocks;
}

// ---------------------------------------------------------------------------
// sdpr_factor_root: R = Z Z' for the factor F (m x p) of R = F'F / (m - 1)
//   Z (p x r) has orthogonal columns with Z'Z = diag(lambda), lambda the
//   nonzero eigenvalues of R, from the eigendecomposition of the smaller of
//   G G' and G'G (G = F / sqrt(m - 1)). Eigenvalues at the rounding level of
//   that decomposition are dropped.
// ---------------------------------------------------------------------------
static arma::mat sdpr_factor_root(const arma::mat &F) {
    arma::mat G = F / std::sqrt(static_cast<double>(F.n_rows) - 1.0);
    bool wide = G.n_rows <= G.n_cols;
    arma::mat gram = wide ? arma::mat(G * G.t()) : arma::mat(G.t() * G);
    arma::vec lam;
    arma::mat Q;
    arma::eig_sym(lam, Q, gram);
    double tol = lam.is_empty() ? 0.0
        : lam.max() * gram.n_rows * arma::datum::eps;
    arma::uvec keep = arma::find(lam > tol);
    if (wide) {
        return G.t() * Q.cols(keep);
    }
    return Q.cols(keep) * arma::diagmat(arma::sqrt(lam(keep)));
}

// ---------------------------------------------------------------------------
// sdpr_llk2_solve: B2^{-1} rhs for the block R = Z Z', with B2 as
// sdpr_llk2_matrix() builds it, without forming B2 while it has a narrow
// low-rank part
//   B2 = T - diag(T) + D0 with T = R .* S .* Mask, D0 = diag(1 / N) + fix.
//   S_ij = 1 / (1.1 max(N_i, N_j)) is a sum over the K distinct sample sizes
//   n_1 < ... < n_K,
//     S = sum_l Delta_l 1_l 1_l' / 1.1,   1_l = [N_i <= n_l],
//     Delta_l = 1 / n_l - 1 / n_(l+1)  (1 / n_K for l = K),
//   so R .* S = U_S U_S' with U_S = [sqrt(Delta_l / 1.1) diag(1_l) Z]_l.
//   The array mask keeps the pairs within arrays {0, 1} and within {0, 2},
//   which counts the pairs within array 0 twice, so T = U Sigma U' with
//   U = [P01 U_S, P02 U_S, P0 U_S] and Sigma = diag(I, I, -I).
//   While U is narrower than p, T is singular, its most negative eigenvalue
//   is that of the small (U'U)^(1/2) Sigma (U'U)^(1/2), and the Woodbury
//   identity solves B2 = U Sigma U' + D through an inner matrix as wide as U.
//   Otherwise B2 is formed exactly as the dense route forms it.
// ---------------------------------------------------------------------------
static arma::mat sdpr_llk2_solve(const mcmc_data &dat, const arma::mat &Z,
                                 const arma::mat &rhs) {
    size_t p = Z.n_rows;
    size_t r = Z.n_cols;
    arma::vec sz(dat.sz);
    arma::ivec arr = arma::conv_to<arma::ivec>::from(dat.array);
    arma::vec levels = arma::unique(sz);
    size_t K = levels.n_elem;
    bool masked = arma::any(arr == 1) && arma::any(arr == 2);
    size_t width = (masked ? 3 : 1) * K * r;

    if (width >= p) {
        arma::mat L_chol = arma::chol(sdpr_llk2_matrix(Z * Z.t(), dat, 0),
                                      "lower");
        arma::mat Y = arma::solve(arma::trimatl(L_chol), rhs);
        return arma::solve(arma::trimatu(L_chol.t()), Y);
    }

    arma::mat U_S(p, K * r);
    for (size_t l = 0; l < K; l++) {
        double next = (l + 1 < K) ? 1.0 / levels(l + 1) : 0.0;
        double delta = 1.0 / levels(l) - next;
        arma::vec in_l = arma::conv_to<arma::vec>::from(sz <= levels(l));
        U_S.cols(l * r, (l + 1) * r - 1) =
            std::sqrt(delta / 1.1) * (Z.each_col() % in_l);
    }

    arma::mat U = U_S;
    arma::vec sigma(width, arma::fill::ones);
    // Unmasked, T = U_S U_S' is positive semidefinite (and singular), so its
    // smallest eigenvalue is 0 and B2 needs no fix.
    double eval_min = 0.0;
    if (masked) {
        arma::vec p01 = arma::conv_to<arma::vec>::from(arr != 2);
        arma::vec p02 = arma::conv_to<arma::vec>::from(arr != 1);
        arma::vec p0  = arma::conv_to<arma::vec>::from(arr == 0);
        U = arma::join_rows(U_S.each_col() % p01, U_S.each_col() % p02,
                            U_S.each_col() % p0);
        sigma.tail(K * r).fill(-1.0);
        arma::vec omega;
        arma::mat P;
        arma::eig_sym(omega, P, U.t() * U);
        arma::mat root = P
            * arma::diagmat(arma::sqrt(arma::clamp(omega, 0.0,
                                                   arma::datum::inf)))
            * P.t();
        arma::vec small = arma::eig_sym(root * arma::diagmat(sigma) * root);
        eval_min = std::min(0.0, small.min());
    }

    // D: 1 / N_i in place of T's own diagonal R_ii / (1.1 N_i), plus the fix.
    arma::vec d = 1.0 / sz - arma::sum(arma::square(Z), 1) / (1.1 * sz);
    if (eval_min < 0) d -= 1.1 * eval_min;

    arma::mat UD    = U.each_col() / d;
    arma::mat D_rhs = rhs.each_col() / d;
    arma::mat inner = arma::diagmat(sigma) + U.t() * UD;
    return D_rhs - UD * arma::solve(inner, U.t() * D_rhs);
}

// ---------------------------------------------------------------------------
// sdpr_factor_block: The one block's B and A' beta_mrg from a factor of R
//   With R = Z Z' = V diag(lambda) V':
//   opt_llk=1: A = n (R + aI)^{-1} R = n V diag(lambda / (lambda + a)) V' and
//              B = R A = W W', W = sqrt(n) V diag(lambda / sqrt(lambda + a)).
//   opt_llk=2: A = B2^{-1} R and B = R A = Z H Z', H = Z' B2^{-1} Z (r x r),
//              so W = Z H^(1/2) and B2^{-1} is needed only on Z and beta_mrg.
// ---------------------------------------------------------------------------
SdprFactorBlock sdpr_factor_block(const mcmc_data &dat, const arma::mat &F,
                                  double a, unsigned sz, int opt_llk) {
    SdprFactorBlock blk;
    blk.Z = sdpr_factor_root(F);
    const arma::mat &Z = blk.Z;
    size_t p = Z.n_rows;
    size_t r = Z.n_cols;
    arma::vec beta_mrg(dat.beta_mrg);

    if (r == 0) {
        // No LD at all (every variant monomorphic): B and A vanish.
        blk.W.zeros(p, 0);
        blk.calc_b_tmp.zeros(p);
    } else if (opt_llk == 1) {
        double n = static_cast<double>(sz);
        arma::vec lam = arma::sum(arma::square(Z), 0).t();
        arma::mat V = Z.each_row() / arma::sqrt(lam).t();
        blk.W = std::sqrt(n)
            * (V.each_row() % (lam / arma::sqrt(lam + a)).t());
        blk.calc_b_tmp = n * (V * ((lam / (lam + a)) % (V.t() * beta_mrg)));
    } else {
        arma::mat Y = sdpr_llk2_solve(dat, Z, arma::join_rows(Z, beta_mrg));
        arma::mat H = Z.t() * Y.head_cols(r);
        H = 0.5 * (H + H.t());
        arma::vec h;
        arma::mat Q;
        arma::eig_sym(h, Q, H);
        blk.W = Z * Q
            * arma::diagmat(arma::sqrt(arma::clamp(h, 0.0, arma::datum::inf)));
        blk.calc_b_tmp = Z * (Z.t() * Y.col(r));
    }
    blk.diag_B = arma::sum(arma::square(blk.W), 1);
    return blk;
}

// ---------------------------------------------------------------------------
// mcmc: Main MCMC loop
// Original: mcmc.cpp lines 434-513
// ---------------------------------------------------------------------------
template <class Block>
std::unordered_map<std::string, arma::vec> mcmc(
    const mcmc_data&          data,
    const std::vector<Block>& blocks,
    size_t                    M,
    double                    a0k,
    double                    b0k,
    int                       iter,
    int                       burn,
    int                       thin,
    unsigned                  n_threads,
    bool                      verbose,
    unsigned int              seed
    ) {

    int n_pst = (iter - burn) / thin;

    MCMC_state state(data.beta_mrg.size(), M, blocks.size(), a0k, b0k, seed);

    MCMC_samples samples(data.beta_mrg.size());

    state.update_suffstats();

    Function_pool func_pool(n_threads);

    for (int j = 1; j < iter + 1; j++) {
        state.sample_sigma2();

        for (size_t i = 0; i < blocks.size(); i++) {
            state.calc_b(i, data, blocks);
        }

        // sample_assignment dispatched to thread pool across LD blocks.
        // NOTE: std::mt19937 is not thread-safe. With n_threads > 1 there
        // is a data race on MCMC_state::r. This matches the original SDPR
        // (eldronzhou/SDPR) which has the same issue with gsl_rng. The
        // default n_threads=1 avoids the race. For safe parallelism, each
        // block would need its own RNG seeded from the shared one.
        for (size_t i = 0; i < blocks.size(); i++) {
            func_pool.push(std::bind(&MCMC_state::sample_assignment<Block>,
                                     &state, i, std::cref(data),
                                     std::cref(blocks)));
        }
        func_pool.waitFinished();

        state.update_suffstats();
        state.sample_V();
        state.update_p();
        state.sample_alpha();

        for (size_t i = 0; i < blocks.size(); i++) {
            state.sample_beta(i, data, blocks);
        }

        state.sample_eta();

        // Collect posterior samples
        if ((j > burn) && (j % thin == 0)) {
            state.compute_h2(data, blocks);
            samples.h2   += state.h2 * square(state.eta) / n_pst;
            samples.beta += state.eta / n_pst * state.beta;
        }

        if (verbose && j % 100 == 0) {
            state.compute_h2(data, blocks);
            Rprintf("%d iter. h2: %g max beta: %g\n",
                    j, state.h2 * square(state.eta),
                    arma::max(state.beta) * state.eta);
        }
    }

    if (verbose) {
        Rprintf("h2: %g max: %g\n", samples.h2, arma::max(samples.beta));
    }

    std::unordered_map<std::string, arma::vec> results;
    results["beta"] = samples.beta;
    results["h2"]   = arma::vec(1, arma::fill::value(samples.h2));

    return results;
}

template std::unordered_map<std::string, arma::vec> mcmc<SdprDenseBlock>(
    const mcmc_data&, const std::vector<SdprDenseBlock>&, size_t, double,
    double, int, int, int, unsigned, bool, unsigned int);
template std::unordered_map<std::string, arma::vec> mcmc<SdprFactorBlock>(
    const mcmc_data&, const std::vector<SdprFactorBlock>&, size_t, double,
    double, int, int, int, unsigned, bool, unsigned int);
