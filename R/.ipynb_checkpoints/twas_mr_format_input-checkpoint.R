#' Generate input format of mr.R
#'
#' @param cand_genes A data frame of significant candidate genes after PTWAS step with columns "chr", "gene_id" and "gene_name"
#' @param susie_path A data frame of the susie results path with columns "ID" and "path", "ID" is the gene names of cand_genes, and "path" is their susie results path.
#' @param weights_file_path A data frame of the output files path of ptwas.R with the columns "ID" and "path".
#'
#' @return A dataframe  of the mr function input format with the columns "snp", "bhat_x",	"sbhat_x", "pip", "cs", "X_ID", "gene_id", "bhat_y" and "sbhat_y", 
#' "bhat_x" and "sbhat_x" are effect size and standard error of exposure, "bhat_y" and "sbhat_y" are effect size and standard error of outcome, "pip" is the posterior inclusion probability (PIP) and "cs" is the credible set index
#' of susie result. "X_ID" and "gene_id" are ensemble ID and gene name respectively.
#'
twas_mr_format_input = function(cand_genes,susie_path,weights_file_path){
mr_format_input = NULL
for (k in 1:dim(cand_genes)[1]){
gene_name = cand_genes$gene_name[k]
gene_id = cand_genes$gene_id[k]
qtl_susie_res = readRDS(susie_path$path[susie_path$ID==gene_name])
if(!is.null(names(qtl_susie_res[[1]]$sets$cs))){
qtl_finemap = qtl_susie_res[[1]]$top_loci%>%filter(cs_index_primary>=1)%>%select(-cs_index_secondary)%>%mutate(X_ID=gene_name)%>%
           rename("snp"="variant_id","bhat_x" = "bhat","sbhat_x"="sbhat","cs" = "cs_index_primary")
outcome_QC = readRDS(weights_file_path$path[weights_file_path$ID==gene_name])$outcome_QC
merge_snp = intersect(qtl_finemap$snp,outcome_QC$variant_allele_flip)
format_input = qtl_finemap[match(merge_snp,qtl_finemap$snp),]%>%
    mutate(gene_id = gene_id)%>%
    cbind(.,outcome_QC[match(merge_snp,outcome_QC$variant_allele_flip),]%>%
    select(beta,se)%>%
    rename("bhat_y"="beta","sbhat_y"="se"))
mr_format_input = rbind(mr_format_input, format_input)
    }
}
return(mr_format_input)
}