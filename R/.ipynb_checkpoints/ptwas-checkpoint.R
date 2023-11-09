#' Calculate the pvalue and zscore of each gene using susie based
#' @description there are five steps:1) use allele_qc function to QC the GWAS summary statistics with xQTL susie results;
#' 2) load the corresponding LD block within the region (start and end) of gene;
#' 3) use twas_z function to calculate TWAS results (pvalue and zscore) from multiple methods of calculating xQTL weights.
#'
#' @param susie_path a data frame with columns "ID" and "path", "ID" is the gene name, and "path" is the susie result path of the corresponding gene.
#' @param region a data frame with columns "chr", "start", "end", and "ID", "chr" is the chromosome of gene, "start" and "end" are the position, "ID" is the gene name.
#' @param GWAS_data a data frame of GWAS summary statistics with columns "chr","pos","A1","A2","beta","se" and "z".
#' @param LD_path a data frame of LD block matrix path with columns "chr","start","end" and "path"
#'
#' @return A list 
#' \describe{
#' \item{twas_z_format}{the input format of twas_z function with the columns "variants_name","gene_name","chr","beta"."se","z","susie_weights","lasso_weights","enet_weights" and "mr_ash_weights".}
#' \item{gene_weights_pq}{twas_z function results with pvalue, qvalue and zscore}
#' \item{outcome_QC}{the GWAS summary statistics of the outcome after QC (mainly accounting for allele flip)}
#' }

ptwas = function(susie_path,region,GWAS_data,LD_block_path){
#####load susie_res_path
#susie = str_split(susie_path,"\\.",simplify=T)[,2]%>%cbind(.,susie_path)%>%data.frame()%>%setNames(c("ID","path"))
gene_name = susie_path$ID
susie_res = readRDS(susie_path$path)
if(!is.null(names(susie_res[[1]]$sets))){
qtl_reference = str_split(susie_res[[1]]$variant_names,":",simplify = T)%>%data.frame()%>%setNames(c("chr","pos","A1","A2"))
######load GWAS summary statistics
outcome_QC = allele_qc(GWAS_data,qtl_reference,match.min.prop=0.2,remove_dups=TRUE,flip=TRUE,remove=TRUE)%>%
              mutate(variant_allele_flip = paste(chr,pos,A1.sumstats,A2.sumstats,sep=":"))
####load LD matrix
LD.files <- bedtoolsr::bt.intersect(a = LD_block_path, b = region)
LD.files.name = unique(LD.files$V5)
LD.list = list()
LD.matrix.names=NULL
for (k in 1:length(LD.files.name)){
 npz = np$load(LD.files.name[k])
 LD.matrix = npz$f[["arr_0"]]
 LD.snps = str_split(LD.files.name[k],"[.]",simplify = T)%>%.[,-c(length(.),(length(.)-1))]%>%paste(.,collapse=".")%>%paste0(.,".bim",sep="")%>%read.table(.)
 #head(LD.snps)
 LD_names = colnames(LD.matrix) = rownames(LD.matrix) = gsub("_",":",LD.snps$V2)
 snp_merge = intersect(LD_names,outcome_QC$variant_allele_flip)
 LD.select = as.matrix(LD.matrix[snp_merge,snp_merge])
 LD.list[[k]] = LD.select
 LD.matrix.names = append(LD.matrix.names,snp_merge)
}
   LD.block = as.matrix(bdiag(LD.list))
   upperTriangle(LD.block,byrow=TRUE) = lowerTriangle(LD.block)
   colnames(LD.block) = rownames(LD.block) = LD.matrix.names
####generate the twas_z format input
twas_z_format = data.frame(LD.matrix.names)%>%mutate(gene_name = gene_name)%>%
    mutate(chr = region$chr)%>%
    mutate(outcome_QC[match(LD.matrix.names,outcome_QC$variant_allele_flip),]%>%select(beta,se,z))%>%
    mutate(susie_weights = susie_res[[1]]$susie_weights[match(LD.matrix.names,susie_res[[1]]$variant_names)])%>%
    mutate(enet_weights =  susie_res[[1]]$enet_weights[match(LD.matrix.names,susie_res[[1]]$variant_names)])%>%
    mutate(lasso_weights = susie_res[[1]]$lasso_weights[match(LD.matrix.names,susie_res[[1]]$variant_names)])%>%
    mutate(mr_ash_weights = susie_res[[1]]$mr_ash_weights[match(LD.matrix.names,susie_res[[1]]$variant_names)])%>%
    rename("variants_name"="LD.matrix.names")
weights = apply(twas_z_format[,c("susie_weights","enet_weights","lasso_weights","mr_ash_weights")],2,function(x) twas_z(x,twas_z_format$z,R = LD.block))
twas_weights = data.frame(gene_name=gene_name,chr = region$chr,weights$susie_weights$pval,weights$susie_weights$z,
                          weights$lasso_weights$pval,weights$lasso_weights$z,weights$enet_weights$pval,weights$enet_weights$z,
                          weights$mr_ash_weights$pval,weights$mr_ash_weights$z)
  names(twas_weights) = c("gene_name","chr","susie_pval","susie_z","lasso_pval","lasso_z","enet_pval","enet_z","mr_ash_pval","mr_ash_z")
  p_values = twas_weights[,c("susie_pval","lasso_pval","enet_pval","mr_ash_pval")]
  p_values[is.na(p_values)] = 1
  q_values = qvalue(p_values,lambda=0)$qvalues
  gene_weights_pq = data.frame(twas_weights,qvalue = q_values)
  names(gene_weights_pq)[11:14]=c("susie_qval","lasso_qval","enet_qval","mr_ash_qval")
####calculate the pvalue and zscore using twas_z function    
return(list(twas_z_format = twas_z_format,
            gene_weights_pq = gene_weights_pq,
            outcome_QC = outcome_QC))
 }
else {
    cat("The 'susie_result' is NULL, so no output is generated.\n")
}
}