#include <iostream>
#include <vector>
#include <armadillo>
#include <cmath>

typedef struct {
    std::vector<std::string> id;
    std::vector<std::string> A1;
    std::vector<std::string> A2;
    std::vector<double> beta_mrg;
    std::vector<std::pair<size_t, size_t>> boundary;
    std::vector<arma::mat> ref_ld_mat;
    std::vector<double> sz;
    std::vector<int> array;
} mcmc_data;

typedef struct {
    std::vector<arma::mat> A;
    std::vector<arma::mat> B;
    std::vector<arma::mat> L;
    std::vector<arma::vec> beta_mrg;
    std::vector<arma::vec> calc_b_tmp;
    std::vector<double> num;
    std::vector<double> denom;
} ldmat_data;

class MCMC_state {
    public:
        double alpha;
        double eta;
        double h2;
        double N;
        arma::vec beta;
        arma::vec b;
        std::vector<int> cls_assgn;
        std::vector<double> V;
        std::vector<double> p;
        std::vector<double> log_p;
        std::vector<double> cluster_var;
        std::vector<unsigned> suff_stats;
        std::vector<double> sumsq;
        MCMC_state(size_t num_snp, size_t max_cluster, \
            double a0, double b0, double sz) {
            a0k = a0; b0k = b0; N = sz;
            // Changed May 20 2021
            // Now N (sz) is absorbed into A, B; so set to 1.
            N = 1.0;

            n_snp = num_snp;
            M = max_cluster;
            alpha = 1;
            eta = 1;
            beta = arma::vec(num_snp, arma::fill::zeros);
            b = arma::vec(num_snp, arma::fill::zeros);
            p.assign(max_cluster, 1.0/max_cluster);
            log_p.assign(max_cluster, 0);
            for (size_t i=0; i<max_cluster; i++) {
                log_p[i] = logf(p[i] + 1e-40);
            }

            cluster_var.assign(max_cluster, 0.0);
            suff_stats.assign(max_cluster, 0);
            sumsq.assign(max_cluster, 0.0);
            V.assign(max_cluster, 0.0);
            cls_assgn.assign(num_snp, 0);
            r.seed(arma::randi<arma::uvec>(1));
            for (size_t i=0; i<num_snp; i++) {
                cls_assgn[i] = arma::randi<int>(arma::distr_param(0, M-1));
            }
        }

        void sample_sigma2();
        void calc_b(size_t j, const mcmc_data &dat, const ldmat_data &ldmat_dat);
        void sample_assignment(size_t j, const mcmc_data &dat, \
                    const ldmat_data &ldmat_dat);
        void update_suffstats();
        void sample_V();
        void update_p();
        void sample_alpha();
        void sample_beta(size_t j, const mcmc_data &dat, \
                   ldmat_data &ldmat_dat);
        void compute_h2(const mcmc_data &dat);
        void sample_eta(const ldmat_data &ldmat_dat);

    private:
        double a0k;
        double b0k;
        size_t M, n_snp;
        arma::arma_rng::set_seed_random r;
};

class MCMC_samples {
    public:
        arma::vec beta;
        double h2;

        MCMC_samples(size_t num_snps) {
            beta = arma::vec(num_snps, arma::fill::zeros);
            h2 = 0;
        }
};

void mcmc(const std::string &ref_path, const std::string &ss_path, \
    const std::string &valid_path, const std::string &ldmat_path, \
    const std::string &out_path, unsigned sz, double a, double c, \
    size_t M, double a0k, double b0k, \
    int iter, int burn, int thin, unsigned n_threads, int opt_llk);

