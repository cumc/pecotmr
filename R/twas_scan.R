#' TWAS Analysis
#'
#' Performs TWAS analysis using the provided weights matrix, GWAS summary statistics database,
#' and LD matrix. It extracts the necessary GWAS summary statistics and LD matrix based on the
#' specified variants and computes the z-score and p-value for each gene.
#'
#' @param weights_matrix A matrix containing weights for all methods.
#' @param gwas_sumstats_db A data frame containing the GWAS summary statistics.
#' @param LD_matrix A matrix representing linkage disequilibrium between variants.
#' @param extract_variants_objs A vector of variant identifiers to extract from the GWAS and LD matrix.
#'
#' @return A list with TWAS z-scores and p-values across four methods for each gene.
#' @export
twas_analysis <- function(weights_matrix, gwas_sumstats_db, LD_matrix, extract_variants_objs) {
    #
    # Extract gwas_sumstats
    gwas_sumstats_subset <- gwas_sumstats_db[match(extract_variants_objs,gwas_sumstats_db$variant_id),]
    # Validate that the GWAS subset is not empty
    if (nrow(gwas_sumstats_subset) == 0) {
    stop("No GWAS summary statistics found for the specified variants.")
    }
    # Extract LD_matrix
    LD_matrix_subset <- LD_matrix[extract_variants_objs,extract_variants_objs]
    # Validate the LD matrix subset
    if (is.null(LD_matrix_subset)) {
     stop("LD matrix subset extraction failed.")
    }
    # Caculate the z score and pvalue of each gene
    twas_z_pval <- apply(as.matrix(weights_matrix), 2, 
                     function(x) twas_z(x, gwas_sumstats_subset$z, R = LD_matrix_subset))
    return(twas_z_pval)
}