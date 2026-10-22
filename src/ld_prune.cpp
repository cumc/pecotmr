// Connected components of an LD threshold graph, for LD pruning.
//
// Single-linkage clustering cut at height 1 - threshold groups exactly the
// variants joined by a path of pairs with |r| >= threshold: the connected
// components of that graph. They are found here with a union-find over R's
// rows, computed a stripe at a time from the standardized dosages, so the
// p x p correlation matrix is never held.
#include <cpp11.hpp>
#include <cpp11armadillo.hpp>
#include <algorithm>
#include <numeric>
#include <vector>

using namespace cpp11;
using namespace arma;

namespace {

// The root of i's component, halving the path on the way up.
uword findRoot(std::vector<uword>& parent, uword i) {
	while (parent[i] != i) {
		parent[i] = parent[parent[i]];
		i = parent[i];
	}
	return i;
}

// Join the components of i and k under the smaller root, so every
// component's root is its smallest index.
void joinComponents(std::vector<uword>& parent, uword i, uword k) {
	uword ri = findRoot(parent, i);
	uword rk = findRoot(parent, k);
	if (ri < rk) {
		parent[rk] = ri;
	} else if (rk < ri) {
		parent[ri] = rk;
	}
}

}  // namespace

// Each variant's component in the graph joining variants whose |r| reaches
// `threshold`, labelled by the component's smallest (1-based) index. The
// correlation of columns j and k of `standardized` (samples x variants, each
// column centred and scaled to unit variance) is F_j'F_k / (m - 1). A pair is
// joined when 1 - |r| <= 1 - threshold: the test single-linkage clustering
// makes when it is cut at that height.
[[cpp11::register]]
cpp11::writable::integers ldThresholdComponents(const doubles_matrix<>& standardized, double threshold,
                                                int stripe) {
	mat F = as_Mat(standardized);
	uword p = F.n_cols;
	double scale = 1.0 / (static_cast<double>(F.n_rows) - 1.0);
	double cut = 1.0 - threshold;
	std::vector<uword> parent(p);
	std::iota(parent.begin(), parent.end(), 0);
	for (uword start = 0; start < p; start += stripe) {
		uword end = std::min(start + static_cast<uword>(stripe), p) - 1;
		// Rows start..end of R against columns start..p-1: every pair (j, k)
		// with k > j is met once, in the stripe holding j.
		mat rows = F.cols(start, end).t() * F.cols(start, p - 1) * scale;
		for (uword a = 0; a < rows.n_rows; ++a) {
			for (uword b = a + 1; b < rows.n_cols; ++b) {
				if (1.0 - std::abs(rows(a, b)) <= cut) {
					joinComponents(parent, start + a, start + b);
				}
			}
		}
	}
	cpp11::writable::integers labels(static_cast<R_xlen_t>(p));
	for (uword i = 0; i < p; ++i) {
		labels[static_cast<R_xlen_t>(i)] = static_cast<int>(findRoot(parent, i)) + 1;
	}
	return labels;
}
