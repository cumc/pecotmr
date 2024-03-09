/**
 * @file mcmc.hpp
 * @brief Markov Chain Monte Carlo (MCMC) sampler for polygenic prediction with continuous shrinkage (CS) priors.
 */

#ifndef MCMC_HPP
#define MCMC_HPP

#include <armadillo>
#include <vector>
#include <string>
#include <random>
#include <cmath>
#include <fstream>

/**
 * @brief Generate random variates from the generalized inverse Gaussian distribution.
 *
 * @param p Shape parameter.
 * @param a Scale parameter.
 * @param b Scale parameter.
 * @return Random variate from the generalized inverse Gaussian distribution.
 */
double gigrnd(double p, double a, double b) {
    double lam = p;
    double omega = std::sqrt(a * b);

    bool swap = false;
    if (lam < 0) {
        lam = -lam;
        swap = true;
    }

    double alpha = std::sqrt(std::pow(omega, 2) + std::pow(lam, 2)) - lam;

    double x = -psi(1.0, alpha, lam);
    double t = 1.0;
    if (x > 2.0) {
        if (alpha == 0 && lam == 0) {
            t = 1.0;
        } else {
            t = std::sqrt(2.0 / (alpha + lam));
        }
    } else if (x < 0.5) {
        if (alpha == 0 && lam == 0) {
            t = 1.0;
        } else {
            t = std::log(4.0 / (alpha + 2.0 * lam));
        }
    }

    x = -psi(-1.0, alpha, lam);
    double s = 1.0;
    if (x > 2.0) {
        if (alpha == 0 && lam == 0) {
            s = 1.0;
        } else {
            s = std::sqrt(4.0 / (alpha * std::cosh(1) + lam));
        }
    } else if (x < 0.5) {
        if (alpha == 0 && lam == 0) {
            s = 1.0;
        } else if (alpha == 0) {
            s = 1.0 / lam;
        } else if (lam == 0) {
            s = std::log(1.0 + 1.0 / alpha + std::sqrt(1.0 / std::pow(alpha, 2) + 2.0 / alpha));
        } else {
            s = std::min(1.0 / lam, std::log(1.0 + 1.0 / alpha + std::sqrt(1.0 / std::pow(alpha, 2) + 2.0 / alpha)));
        }
    }

    double eta = -psi(t, alpha, lam);
    double zeta = -dpsi(t, alpha, lam);
    double theta = -psi(-s, alpha, lam);
    double xi = dpsi(-s, alpha, lam);

    double p_r = 1.0 / xi;
    double r = 1.0 / zeta;

    double td = t - r * eta;
    double sd = s - p_r * theta;
    double q = td + sd;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);

    while (true) {
        double U = dis(gen);
        double V = dis(gen);
        double W = dis(gen);

        double rnd;
        if (U < q / (p_r + q + r)) {
            rnd = -sd + q * V;
        } else if (U < (q + r) / (p_r + q + r)) {
            rnd = td - r * std::log(V);
        } else {
            rnd = -sd + p_r * std::log(V);
        }

        double f1 = std::exp(-eta - zeta * (rnd - t));
        double f2 = std::exp(-theta + xi * (rnd + s));
        if (W * g(rnd, sd, td, f1, f2) <= std::exp(psi(rnd, alpha, lam))) {
            break;
        }
    }

    rnd = std::exp(rnd) * (lam / omega + std::sqrt(1.0 + std::pow(lam, 2) / std::pow(omega, 2)));
    if (swap) {
        rnd = 1.0 / rnd;
    }

    rnd = rnd / std::sqrt(a / b);
    return rnd;
}

/**
 * @brief Markov Chain Monte Carlo (MCMC) sampler for polygenic prediction with continuous shrinkage (CS) priors.
 *
 * @param a Shape parameter for the prior distribution of psi.
 * @param b Scale parameter for the prior distribution of psi.
 * @param phi Global shrinkage parameter. If nullptr, it will be estimated automatically.
 * @param sst_dict Dictionary containing summary statistics.
 * @param n Sample size.
 * @param ld_blk List of LD blocks.
 * @param blk_size List of block sizes.
 * @param n_iter Number of MCMC iterations.
 * @param n_burnin Number of burn-in iterations.
 * @param thin Thinning interval.
 * @param chrom Chromosome number.
 * @param out_dir Output directory.
 * @param beta_std Whether to standardize the effect sizes.
 * @param write_psi Whether to write the posterior estimates of psi to a file.
 * @param seed Random seed. If nullptr, no seed is set.
 */
void prs_cs(double a, double b, double* phi, const std::vector<std::vector<double>>& sst_dict,
            int n, const std::vector<arma::mat>& ld_blk, const std::vector<int>& blk_size,
            int n_iter, int n_burnin, int thin, int chrom, const std::string& out_dir,
            bool beta_std, bool write_psi, int* seed) {
    std::cout << "... MCMC ..." << std::endl;

    // Seed the random number generator
    if (seed != nullptr) {
        std::srand(*seed);
    }

    // Derived statistics
    arma::vec beta_mrg(sst_dict[0].size());
    arma::vec maf(sst_dict[0].size());
    for (size_t i = 0; i < sst_dict[0].size(); ++i) {
        beta_mrg(i) = sst_dict[1][i];
        maf(i) = sst_dict[2][i];
    }
    int n_pst = (n_iter - n_burnin) / thin;
    int p = sst_dict[0].size();
    int n_blk = ld_blk.size();

    // Initialization
    arma::vec beta(p, arma::fill::zeros);
    arma::vec psi(p, arma::fill::ones);
    double sigma = 1.0;
    bool phi_updt = (phi == nullptr);
    if (!phi_updt) {
        *phi = 1.0;
    }

    arma::vec beta_est(p, arma::fill::zeros);
    arma::vec psi_est(p, arma::fill::zeros);
    double sigma_est = 0.0;
    double phi_est = 0.0;

    // MCMC
    for (int itr = 1; itr <= n_iter; ++itr) {
        if (itr % 100 == 0) {
            std::cout << "--- iter-" << itr << " ---" << std::endl;
        }

        int mm = 0;
        double quad = 0.0;
        for (int kk = 0; kk < n_blk; ++kk) {
            if (blk_size[kk] == 0) {
                continue;
            }

            arma::uvec idx_blk = arma::regspace<arma::uvec>(mm, mm + blk_size[kk] - 1);
            arma::mat dinvt = ld_blk[kk] + arma::diagmat(1.0 / psi(idx_blk));
            arma::mat dinvt_chol = arma::chol(dinvt);
            arma::vec beta_tmp = arma::solve(arma::trimatl(dinvt_chol.t()), beta_mrg(idx_blk), arma::solve_opts::fast) +
                                 arma::randn<arma::vec>(blk_size[kk]) * std::sqrt(sigma / n);
            beta(idx_blk) = arma::solve(arma::trimatu(dinvt_chol), beta_tmp, arma::solve_opts::fast);
            quad += arma::as_scalar(beta(idx_blk).t() * dinvt * beta(idx_blk));
            mm += blk_size[kk];
        }

        double err = std::max(n / 2.0 * (1.0 - 2.0 * arma::sum(beta % beta_mrg) + quad),
                              n / 2.0 * arma::sum(arma::pow(beta, 2) / psi));
        sigma = 1.0 / std::gamma_distribution<double>((n + p) / 2.0, 1.0 / err)(std::mt19937(std::random_device{}()));

        arma::vec delta = arma::randg<arma::vec>(p, arma::distr_param(a + b, 1.0 / (psi + *phi)));

        for (int jj = 0; jj < p; ++jj) {
            psi(jj) = gigrnd(a - 0.5, 2.0 * delta(jj), n * std::pow(beta(jj), 2) / sigma);
        }
        psi.elem(arma::find(psi > 1)).fill(1.0);

        if (phi_updt) {
            double w = std::gamma_distribution<double>(1.0, 1.0 / (*phi + 1.0))(std::mt19937(std::random_device{}()));
            *phi = std::gamma_distribution<double>(p * b + 0.5, 1.0 / (arma::sum(delta) + w))(std::mt19937(std::random_device{}()));
        }

        // Posterior
        if (itr > n_burnin && (itr % thin == 0)) {
            beta_est += beta / n_pst;
            psi_est += psi / n_pst;
            sigma_est += sigma / n_pst;
            phi_est += *phi / n_pst;
        }
    }

    // Convert standardized beta to per-allele beta
    if (!beta_std) {
        beta_est /= arma::sqrt(2.0 * maf % (1.0 - maf));
    }

    // Write posterior effect sizes
    std::string eff_file;
    if (phi_updt) {
        eff_file = out_dir + "_pst_eff_a" + std::to_string(a) + "_b" + std::to_string(b) + "_phiauto_chr" + std::to_string(chrom) + ".txt";
    } else {
        eff_file = out_dir + "_pst_eff_a" + std::to_string(a) + "_b" + std::to_string(b) + "_phi" + std::to_string(*phi) + "_chr" + std::to_string(chrom) + ".txt";
    }

    std::ofstream eff_out(eff_file);
    for (size_t i = 0; i < sst_dict[0].size(); ++i) {
        eff_out << chrom << "\t" << sst_dict[0][i] << "\t" << sst_dict[3][i] << "\t"
                << sst_dict[4][i] << "\t" << sst_dict[5][i] << "\t" << beta_est(i) << "\n";
    }
    eff_out.close();

    // Write posterior estimates of psi
    if (write_psi) {
        std::string psi_file;
        if (phi_updt) {
            psi_file = out_dir + "_pst_psi_a" + std::to_string(a) + "_b" + std::to_string(b) + "_phiauto_chr" + std::to_string(chrom) + ".txt";
        } else {
            psi_file = out_dir + "_pst_psi_a" + std::to_string(a) + "_b" + std::to_string(b) + "_phi" + std::to_string(*phi) + "_chr" + std::to_string(chrom) + ".txt";
        }

        std::ofstream psi_out(psi_file);
        for (size_t i = 0; i < sst_dict[0].size(); ++i) {
            psi_out << sst_dict[0][i] << "\t" << psi_est(i) << "\n";
        }
        psi_out.close();
    }

    // Print estimated phi
    if (phi_updt) {
        std::cout << "... Estimated global shrinkage parameter: " << phi_est << " ..." << std::endl;
    }

    std::cout << "... Done ..." << std::endl;
}

#endif // MCMC_HPP