#include <RcppArmadillo.h>
#include <omp.h> // Required for parallel processing
#include <algorithm>
#include <random>
#include <vector>
#include <numeric> // For std::iota

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;
using namespace arma;


// Assuming sort_indexes is defined as provided
std::vector<size_t> sort_indexes(const std::vector<int>& v, unsigned int theSize) {
    std::vector<size_t> idx(theSize);
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});
    return idx;
}

// Improved generateSetOfNumbers function using C++11 random
std::vector<size_t> generateSetOfNumbers(int SIZE, int seed) {
    std::vector<int> numbers(SIZE, 0);
    std::mt19937 rng(seed); // Mersenne Twister: Good quality random number generator
    std::uniform_int_distribution<int> dist(0, INT_MAX);

    // Generate the first random number
    numbers[0] = dist(rng);
    for (int index = 1; index < SIZE; index++) {
        int tempNum;
        do {
            tempNum = dist(rng); // Generate a new random number
            // Check for uniqueness in the current list of generated numbers
            bool isUnique = true;
            for (int index2 = 0; index2 < index; index2++) {
                if (tempNum == numbers[index2]) {
                    isUnique = false;
                    break;
                }
            }
            // If the number is not unique, force the loop to try again
            if (!isUnique) tempNum = -1;
        } while (tempNum == -1);
        // Assign the unique number to the list
        numbers[index] = tempNum;
    }

    // Sort the indices of 'numbers' based on their values
    return sort_indexes(numbers, SIZE);
}

// Get a quantile value
double getQuantile(const std::vector<double>& dat, double whichQuantile) {
    std::vector<double> sortedData = dat;
    std::sort(sortedData.begin(), sortedData.end());
    size_t pos = ceil(sortedData.size() * whichQuantile) - 1;
    return sortedData.at(pos);
}

// Get a quantile value based on grouping
double getQuantile2(const std::vector<double>& dat, const std::vector<uint>& grouping, double whichQuantile) {
    std::vector<double> filteredData;
    for (size_t i = 0; i < dat.size(); ++i) if (grouping[i] == 1) filteredData.push_back(dat[i]);
    if (filteredData.size() < 50) return 0;
    return getQuantile(filteredData, whichQuantile);
}

// Calculate minus log p-value of chi-squared statistic
double minusLogPvalueChisq2(double stat) {
    double p = 1.0 - arma::chi2cdf(stat, 1.0);
    return -log10(p);
}

// Perform one iteration of the algorithm, assuming LDmat is an arma::mat
void oneIteration(const arma::mat& LDmat, const std::vector<uint>& idx, const std::vector<uint>& idx2,
                  arma::vec& zScore, arma::vec& imputedZ, arma::vec& rsqList, arma::vec& zScore_e,
                  uint nSample, float probSVD, int ncpus) {
    omp_set_num_threads(ncpus);

    uint K = std::min(static_cast<uint>(idx.size()), nSample) * probSVD;

    arma::mat LD_it(idx2.size(), idx.size());
    arma::vec zScore_eigen(idx.size());
    arma::mat VV(idx.size(), idx.size());

    // Fill LD_it and VV matrices using direct indexing
    #pragma omp parallel for collapse(2)
    for (size_t i = 0; i < idx2.size(); i++) {
        for (size_t k = 0; k < idx.size(); k++) {
            LD_it(i, k) = LDmat(idx2[i], idx[k]);
        }
    }

    #pragma omp parallel for
    for (size_t i = 0; i < idx.size(); i++) {
        zScore_eigen(i) = zScore[idx[i]];
        for (size_t j = 0; j < idx.size(); j++) {
            VV(i, j) = LDmat(idx[i], idx[j]);
        }
    }

    // Eigen decomposition
    arma::vec eigval;
    arma::mat eigvec;
    arma::eig_sym(eigval, eigvec, VV);

    int nRank = eigvec.n_rows;
    int nZeros = arma::sum(eigval < 0.0001);
    nRank -= nZeros;
    K = std::min(K, static_cast<uint>(nRank));
    if (K <= 1) {
        Rcpp::stop("Rank of eigen matrix <= 1");
    }

    arma::mat ui = arma::mat(eigvec.n_rows, K, arma::fill::zeros);
    arma::mat wi = arma::mat(K, K, arma::fill::zeros);
    for (uint m = 0; m < K; ++m) {
        int j = eigvec.n_rows - m - 1;
        ui.col(m) = eigvec.col(j);
        wi(m, m) = 1.0 / eigval(j);
    }

    // Calculate imputed Z scores and R squared values
    arma::mat beta = LD_it * ui * wi;
    arma::vec zScore_eigen_imp = beta * (ui.t() * zScore_eigen);
    arma::vec rsq_eigen = (beta * (ui.t() * LD_it.t())).diag();

    #pragma omp parallel for
    for (size_t i = 0; i < idx2.size(); ++i) {
        imputedZ[idx2[i]] = zScore_eigen_imp(i);
        rsqList[idx2[i]] = rsq_eigen(i);
        if (rsq_eigen(i) >= 1) {
            Rcpp::stop("Dividing zero: Rsq = " + std::to_string(rsq_eigen(i)));
        }
        uint j = idx2[i];
        zScore_e[j] = (zScore[j] - imputedZ[j]) / std::sqrt(LDmat(j, j) - rsqList[j]);
    }
}



// Adapted DENTIST function to RcppArmadillo, including all necessary steps and logic
// [[Rcpp::export]]
List DENTIST(arma::mat& LDmat, unsigned int markerSize, unsigned int nSample, arma::vec& zScore,
             arma::vec imputedZ, arma::vec rsq, arma::vec zScore_e, arma::ivec iterID,
             double pValueThreshold, arma::ivec interested, float propSVD, bool gcControl, int nIter,
             double groupingPvalue_thresh, int ncpus)

    omp_set_num_threads(ncpus);

    std::vector<size_t> randOrder = generateSetOfNumbers(markerSize, 10); // Example seed
    std::vector<uint> idx, idx2, fullIdx = randOrder;
    for (uint i = 0; i < markerSize; i++) {
        if (randOrder[i] > markerSize / 2) idx.push_back(i);
        else idx2.push_back(i);
    }

    std::vector<uint> groupingGWAS(markerSize, 0);
    for (uint i = 0; i < markerSize; i++) {
        if (minusLogPvalueChisq2(zScore[i] * zScore[i]) > -log10(groupingPvalue_thresh)) {
            groupingGWAS[i] = 1;
        }
    }

    for (int t = 0; t < nIter; t++) {
        std::vector<uint> idx2_QCed;
        std::vector<size_t> fullIdx_tmp;
        std::vector<double> diff, chisq;
        std::vector<uint> grouping_tmp;

        oneIteration(LDmat, idx, idx2, zScore, imputedZ, rsq, zScore_e, nSample, propSVD, ncpus);

        diff.resize(idx2.size());
        grouping_tmp.resize(idx2.size());
        for(size_t i = 0; i < idx2.size(); i++) {
            diff[i] = std::abs(zScore_e[idx2[i]]);
            grouping_tmp[i] = groupingGWAS[idx2[i]];
        }

        double threshold = getQuantile(diff, 0.995);
        double threshold1 = getQuantile2(diff, grouping_tmp, 0.995);
        double threshold0 = getQuantile2(diff, grouping_tmp, 0.995, false); // Assuming a modified version of getQuantile2 that accepts a flag for inverse grouping

        for(size_t i = 0; i < diff.size(); i++) {
            if ((grouping_tmp[i] == 1 && diff[i] <= threshold1) || (grouping_tmp[i] == 0 && diff[i] <= threshold0)) {
                idx2_QCed.push_back(idx2[i]);
            }
        }

        oneIteration(LDmat, idx, idx2_QCed, zScore, imputedZ, rsq, zScore_e, nSample, propSVD, ncpus);

        // Recalculate diff and grouping_tmp for the fullIdx
        diff.resize(fullIdx.size());
        grouping_tmp.resize(fullIdx.size());
        for(size_t i = 0; i < fullIdx.size(); i++) {
            diff[i] = std::abs(zScore_e[fullIdx[i]]);
            grouping_tmp[i] = groupingGWAS[fullIdx[i]];
        }

        threshold = getQuantile(diff, 0.995);
        threshold1 = getQuantile2(diff, grouping_tmp, 0.995);
        threshold0 = getQuantile2(diff, grouping_tmp, 0.995, false);

        // Adjust for inflation factor if gcControl is true
        if(gcControl) {
            double inflationFactor = calculateInflationFactor(chisq, pValueThreshold); // Assuming a function for inflation factor calculation

            for(size_t i = 0; i < diff.size(); i++) {
                if (!(diff[i] > threshold && minusLogPvalueChisq2(diff[i] * diff[i] / inflationFactor) > -log10(pValueThreshold))) {
                    fullIdx_tmp.push_back(fullIdx[i]);
                }
            }
        } else {
            for(size_t i = 0; i < diff.size(); i++) {
                if (minusLogPvalueChisq2(diff[i] * diff[i]) < -log10(pValueThreshold)) {
                    fullIdx_tmp.push_back(fullIdx[i]);
                }
            }
        }

        fullIdx = fullIdx_tmp;
        randOrder = generateSetOfNumbers(fullIdx.size(), 20000 + t * 20000); // Update seed for randomness
        idx.clear();
        idx2.clear();
        for(size_t i = 0; i < fullIdx.size(); i++) {
            if (randOrder[i] > fullIdx.size() / 2) idx.push_back(fullIdx[i]);
            else idx2.push_back(fullIdx[i]);
        }
    }
        return List::create(Named("imputedZ") = imputedZ,
                        Named("rsq") = rsq,
                        Named("zScore_e") = zScore_e,
                        Named("iterID") = iterID,
                        Named("groupingGWAS") = wrap(groupingGWAS)); // Convert std::vector<uint> to Rcpp::NumericVector
}