#######generate the input format of mr.R
#cand_genes: the output of the twas_z
#susie_res_path: a list of files of susie results
#AD_GWAS_path: GWAS summary statistics path
twas_mr_format_input = function(cand_genes,susie_res_path,AD_GWAS_path){
mr_format_input = NULL
for (k in 1:dim(cand_genes)[1]){
gene_name = cand_genes$gene_name[k]
chr  = cand_genes$chr[k]
susie = fread(susie_res_path)
qtl_susie_res = readRDS(susie$path[susie$ID==gene_name])
if(!is.null(names(qtl_susie_res[[1]]$sets$cs))){
qtl_finemap = qtl_susie_res[[1]]$top_loci%>%filter(cs_index_primary>=1)%>%select(-cs_index_secondary)%>%mutate(X_ID=gene_name)%>%
           rename("snp"="variant_id","bhat_x" = "bhat","sbhat_x"="sbhat","cs" = "cs_index_primary")
qtl_reference = str_split(qtl_finemap$snp,":",simplify = T)%>%data.frame()%>%setNames(c("chr","pos","A1","A2"))
AD_dataset = fread(paste0(AD_GWAS_path,"/ADGWAS_Bellenguez_2022.",chr,"/ADGWAS2022.chr",chr,".sumstat.tsv",sep=""))
##transform the AD_dataset to the form of allele_qc format
AD_data = AD_dataset%>%mutate(chr = paste0("chr",chromosome))%>%mutate(pos = as.character(position))%>%select(-chromosome)%>%rename("A1"="ref","A2"="alt")%>%mutate(z=beta/se)%>%select(-position)
AD_allele_flip = allele_qc(AD_data,qtl_reference,match.min.prop=0.2,remove_dups=TRUE,flip=TRUE,remove=TRUE)%>%
              mutate(variant_allele_flip = paste(chr,pos,A1.sumstats,A2.sumstats,sep=":"))
format_input = qtl_finemap%>%cbind(.,AD_allele_flip[match(qtl_finemap$snp,AD_allele_flip$variant_allele_flip),]%>%select(beta,se)%>%rename("bhat_y"="beta","sbhat_y"="se"))
mr_format_input = rbind(mr_format_input, format_input)

      }
    }
return(mr_format_input)
}