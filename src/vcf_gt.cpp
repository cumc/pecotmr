/**
 * @file vcf_gt.cpp
 * @brief Allele counts from VCF GT strings.
 *
 * The R side (R/GenotypeDosageSeed.R) reads a block of records' GT field
 * through VariantAnnotation, which hands it back as strings ("0/1", "1|2",
 * "./.", "1"), and asks here how many copies of one ALT allele each call
 * carries. A row of the block may ask about any of its record's ALT alleles,
 * which is how a record with several ALT alleles is read as one variant per
 * ALT allele.
 */

#include <cpp11.hpp>

namespace {

// The copies of allele `allele` in one GT string: its alleles, of any ploidy,
// separated by '/' or '|'. A call with a missing allele ("./.", "0|.", ".")
// is missing as a whole, as PLINK reads a half-call, and so is anything that
// is not a GT call.
double gt_copies(SEXP call, int allele) {
	if (call == NA_STRING) return NA_REAL;
	const char* p = CHAR(call);
	if (*p == '\0') return NA_REAL;
	int copies = 0;
	while (true) {
		if (*p < '0' || *p > '9') return NA_REAL;
		long value = 0;
		while (*p >= '0' && *p <= '9') {
			value = value * 10 + (*p - '0');
			++p;
		}
		if (value == allele) ++copies;
		if (*p == '\0') break;
		if (*p != '/' && *p != '|') return NA_REAL;
		++p;
	}
	return copies;
}

}  // namespace

/**
 * @brief The dosage of each row's ALT allele in a block of GT calls.
 *
 * @param gt A records x samples character matrix of GT calls, as
 *   VariantAnnotation::geno() returns it.
 * @param alleles For each row, the 1-based index of the ALT allele it
 *   counts (1 for the first ALT allele).
 * @return A numeric matrix the shape of `gt`: each call's copies of its
 *   row's ALT allele, NA where the call is missing.
 */
[[cpp11::register]]
cpp11::writable::doubles_matrix<> vcfGtDosageCpp(cpp11::strings gt,
                                                  cpp11::integers alleles) {
	const R_xlen_t nRow = alleles.size();
	const R_xlen_t nCell = gt.size();
	if (nRow == 0) return cpp11::writable::doubles_matrix<>(0, 0);
	if (nCell % nRow != 0)
		cpp11::stop("The GT block does not have one row per allele.");
	const R_xlen_t nCol = nCell / nRow;
	cpp11::writable::doubles_matrix<> out(nRow, nCol);
	for (R_xlen_t c = 0; c < nCol; ++c) {
		for (R_xlen_t r = 0; r < nRow; ++r) {
			out(r, c) = gt_copies(STRING_ELT(gt, r + c * nRow), alleles[r]);
		}
	}
	return out;
}
