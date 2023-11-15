# heterogeneity:  calculate I2 statistics based on the Cochran's Q statistic
calc_I2 = function(Q, Est) {
    Q = Q[[1]]
    Est = length(unique(Est))
    I2 = if (Q > 1e-3) (Q - Est + 1)/Q else 0
    return(if (I2 < 0) 0 else I2)
}

#' Generate input format for MR
#'
#' @param target_genes A data frame of significant candidate genes after PTWAS step with columns "chr", "gene_id" and "gene_name"
#' @param susie_path A data frame of the susie results path with columns "ID" and "path", "ID" is the gene names of target_genes, and "path" is their susie results path.
#' @param weights_file_path A data frame of the output files path of twas_scan() with the columns "ID" and "path".
#' @return A dataframe  of the mr function input format with the columns "snp", "bhat_x",	"sbhat_x", "pip", "cs", "X_ID", "gene_id", "bhat_y" and "sbhat_y", 
#' "bhat_x" and "sbhat_x" are effect size and standard error of exposure, "bhat_y" and "sbhat_y" are effect size and standard error of outcome, "pip" is the posterior inclusion probability (PIP) and "cs" is the credible set index
#' of susie result. "X_ID" and "gene_id" are ensemble ID and gene name respectively.
#' 
#' @importFrom qvalue qvalue
#' @import dplyr
#' @export
twas_mr_format_input <- function(target_genes, susie_path, weights_file_path) {
  mr_format_input <- NULL
  for (k in 1:dim(target_genes)[1]) {
    gene_name <- target_genes$gene_name[k]
    gene_id <- target_genes$gene_id[k]
    qtl_susie_res <- readRDS(susie_path$path[susie_path$ID == gene_name])
    
    if (!is.null(names(qtl_susie_res[[1]]$sets$cs))) {
      qtl_finemap <- qtl_susie_res[[1]]$top_loci %>% 
        filter(cs_index_primary >= 1) %>% 
        select(-cs_index_secondary) %>% 
        mutate(X_ID = gene_name) %>% 
        rename("snp" = "variant_id", "bhat_x" = "bhat", "sbhat_x" = "sbhat", "cs" = "cs_index_primary")
      
      outcome_QC <- readRDS(weights_file_path$path[weights_file_path$ID == gene_name])$outcome_QC
      merge_snp <- intersect(qtl_finemap$snp, outcome_QC$variant_allele_flip)
      
      format_input <- qtl_finemap[match(merge_snp, qtl_finemap$snp), ] %>% 
        mutate(gene_id = gene_id) %>% 
        cbind(., outcome_QC[match(merge_snp, outcome_QC$variant_allele_flip), ] %>% 
                select(beta, se) %>% 
                rename("bhat_y" = "beta", "sbhat_y" = "se"))
      
      mr_format_input <- rbind(mr_format_input, format_input)
    }
  }
  return(mr_format_input)
}

#' Mendelian Randomization (MR)
#' 
#' @param formatted_input the output of twas_mr_format_input function
#' @param cpip_cutoff the threshold of cumulative posterior inclusion probability, default is 0.5
#' @return A single data frame of output with columns "X_ID", "gene_id", "chr", "num_CS", "num_IV",
#' "meta_eff", "se_meta_eff", "meta_pval", "Q", "Q_pval" and "I2". "X_ID" and "gene_id" are ensemble ID and gene name, respectively. "num_CS" is the number of credible sets
#' contained in each gene, "num_IV" is the number of variants contained in each gene. "meta_eff", "se_meta_eff" and "meta_pval" are the MR estimate, standard error and pvalue.
#' "Q" is Cochran’s Q statistic, "I2" quantifies the heterogeneity, range from 0 to 1.
#' @import dplyr
#' @importFrom qvalue qvalue
#' @importFrom stats pnorm
#' @export
fine_mr <- function(formatted_input, cpip_cutoff=0.5) {
output = formatted_input %>%
    mutate(
        bhat_x = bhat_x/sbhat_x,
        sbhat_x = 1) %>%
    group_by(X_ID, cs) %>%
    mutate(cpip = sum(pip)) %>%
    filter(cpip >= cpip_cutoff) %>% # Cumulative pip greater than a user defined cumulative pip threshold
    group_by(X_ID, cs) %>%
    mutate(
        beta_yx = bhat_y/bhat_x,
        se_yx = sqrt(
            (sbhat_y^2/bhat_x^2) + ((bhat_y^2*sbhat_x^2)/bhat_x^4)),
        composite_bhat = sum((beta_yx*pip)/cpip),
        composite_sbhat = sum((beta_yx^2 + se_yx^2)*pip/cpip)) %>%
    mutate(
        composite_sbhat = sqrt(composite_sbhat - composite_bhat^2),
        wv = composite_sbhat^-2) %>%
    ungroup() %>%
    group_by(X_ID) %>%
    mutate(
        meta_eff = sum(unique(wv) * unique(composite_bhat)),
        sum_w = sum(unique(wv)),
        se_meta_eff = sqrt(sum_w^-1),
        num_CS = length(unique(cs))) %>%
    mutate(
        num_IV = length(snp),
        meta_eff = meta_eff/sum_w,
        meta_pval = 2 * pnorm(abs(meta_eff) / se_meta_eff, lower.tail=FALSE),
        ######sum(unique(wv)*(unique(composite_bhat) - unique(meta_eff))^2)
        Q = sum(unique(wv)*(unique(composite_bhat) - unique(meta_eff))^2),
        I2 = calc_I2(Q, composite_bhat),
        Q_pval = pchisq(Q, df = length(unique(composite_bhat))-1, lower = F)) %>%
        ungroup() %>%
    distinct(X_ID, .keep_all = TRUE) %>%
    mutate(
        meta_qval = qvalue(meta_pval,lambda=0)$qvalues
        )%>%
    mutate(
        cpip = round(cpip, 3),
        meta_qval = round(meta_qval,3),
        meta_pval= round(meta_pval,3),
        #composite_bhat = round(composite_bhat, 3),
        meta_eff = round(meta_eff, 3),
        se_meta_eff = round(se_meta_eff, 3),
        Q = round(Q, 3),
        Q_pval = round(Q_pval,3),
        I2 = round(I2, 3)) %>%arrange(meta_pval)%>%
    select(X_ID, num_CS, num_IV, cpip, gene_id, meta_eff, se_meta_eff, meta_pval,meta_qval,Q, Q_pval,I2)
    return(output)
}