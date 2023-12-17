#' Calculate the pvalue and zscore of each gene using precomputed weights
#' @description there are five steps:1) use allele_qc function to QC the GWAS summary statistics with xQTL weights;
#' 2) load the corresponding LD block within the region (start and end) of gene;
#' 3) use twas_z function to calculate TWAS results (pvalue and zscore) from multiple methods of calculating xQTL weights.
#' @param weights_path a data frame with columns "ID" and "path", "ID" is the gene name, and "path" is the susie result path of the corresponding gene.
#' @param region a data frame with columns "chrom", "start", "end", and "ID", "chrom" is the chromosome of gene, "start" and "end" are the position, "ID" is the gene name.
#' @param GWAS_data a data frame of GWAS summary statistics with columns "chr","pos","A1","A2","beta","se" and "z".
#' @param LD_meta_file a data frame of LD block matrix path with columns "chrom","start","end" and "path"
#'
#' @return A list 
#' \describe{
#' \item{twas_z_format}{the input format of twas_z function with the columns "variants_name","gene_name","chr","beta"."se","z","susie_weights","lasso_weights","enet_weights" and "mrash_weights".}
#' \item{gene_weights_pq}{twas_z function results with pvalue, qvalue and zscore}
#' \item{outcome_QC}{the GWAS summary statistics of the outcome after QC (mainly accounting for allele flip)}
#' }
#' @importFrom Matrix bdiag
#' @importFrom stringr str_split
#' @importFrom stats setNames
#' @importFrom utils read.table
#' @export
twas_scan <- function(weights_path, region, GWAS_data, LD_meta_file) {
  # Load weights
  gene_name <- weights_path$ID
  qtl_weights <- readRDS(weights_path$path)
  
  if (!is.null(names(qtl_weights[[1]]$sets))) {
    qtl_reference <- str_split(qtl_weights[[1]]$variant_names, ":", simplify = TRUE) %>% 
      data.frame() %>% 
      setNames(c("chr", "pos", "A1", "A2"))
    
    # Load GWAS summary statistics
    outcome_QC <- allele_qc(GWAS_data, qtl_reference, match.min.prop = 0.2, 
                            remove_dups = TRUE, flip = TRUE, remove = TRUE) %>% 
      mutate(variant_allele_flip = paste(chr, pos, A1.sumstats, A2.sumstats, sep = ":"))
    
    LD.matrix <- load_LD_matrix(LD_meta_file, region)
    
    # Generate the twas_z format input
    twas_z_format <- data.frame(LD.matrix$variants_id_all) %>% 
      mutate(gene_name = gene_name) %>% 
      mutate(chr = region$chr) %>% 
      mutate(outcome_QC[match(LD.matrix$variants_id_all, outcome_QC$variant_allele_flip), ] %>% 
               select(beta, se, z)) %>% 
      mutate(susie_weights = qtl_weights[[1]]$susie_weights[match(LD.matrix$variants_id_all, qtl_weights[[1]]$variant_names)]) %>% 
      mutate(enet_weights = qtl_weights[[1]]$enet_weights[match(LD.matrix$variants_id_all, qtl_weights[[1]]$variant_names)]) %>% 
      mutate(lasso_weights = qtl_weights[[1]]$lasso_weights[match(LD.matrix$variants_id_all, qtl_weights[[1]]$variant_names)]) %>% 
      mutate(mrash_weights = qtl_weights[[1]]$mrash_weights[match(LD.matrix$variants_id_all, qtl_weights[[1]]$variant_names)]) %>% 
      rename("variants_name" = "LD.files.variants")
    
    weights <- apply(twas_z_format[, c("susie_weights", "enet_weights", "lasso_weights", "mrash_weights")], 2, 
                     function(x) twas_z(x, twas_z_format$z, R = LD.matrix$LD))
    
    twas_weights <- data.frame(gene_name = gene_name, chr = region$chr, 
                               weights$susie_weights$pval, weights$susie_weights$z, 
                               weights$lasso_weights$pval, weights$lasso_weights$z, 
                               weights$enet_weights$pval, weights$enet_weights$z, 
                               weights$mrash_weights$pval, weights$mrash_weights$z)
    
    names(twas_weights) <- c("gene_name", "chr", "susie_pval", "susie_z", "lasso_pval", 
                             "lasso_z", "enet_pval", "enet_z", "mrash_pval", "mrash_z")
    
    p_values <- twas_weights[, c("susie_pval", "lasso_pval", "enet_pval", "mrash_pval")]
    p_values[is.na(p_values)] <- 1
    q_values <- qvalue(p_values, lambda = 0)$qvalues
    
    gene_weights_pq <- data.frame(twas_weights, qvalue = q_values)
    names(gene_weights_pq)[11:14] <- c("susie_qval", "lasso_qval", "enet_qval", "mrash_qval")
    
    # Calculate the pvalue and zscore using twas_z function
    return(list(twas_z_format = twas_z_format,
                gene_weights_pq = gene_weights_pq,
                outcome_QC = outcome_QC))
  } else {
    cat("The 'qtl_weights' is NULL, so no output is generated.\n")
  }
}
