######calculate the pvalue and zscore of each gene
#susie_res_path: a list of path of susie results
#qtl_TAD_region: a data frame of the TAD region of each gene
#AD_GWAS_path: AD GWAS summary statistics path
#LD path: LD matrix path
ptwas = function(susie_res_path,qtl_TAD_region,AD_GWAS_path,LD_path){
#####load LD path list
LD_list = read.table(LD_path,header=F,sep="\t")
LD_list_pos.bed = str_split(LD_list$V1,"_",simplify=T)%>%cbind(.,LD_list)

gene_name = susie_res_path$ID
######this one should be replaced by the rds file list, load qtl_select, load the qtl_variants
susie_res = readRDS(susie_res_path$path)
qtl_select = qtl_TAD_region%>%filter(ID==gene_name)
chr = str_sub(qtl_select$`#chr`,4)
#####after runing the susie fine mapping again, we will remove this line of code (this one just consider the nchar(ref)>1|nchar(alt)>1)
qtl_reference = str_split(susie_res[[1]]$variant_names,":",simplify = T)%>%data.frame()%>%setNames(c("chr","pos","A1","A2"))%>%
           mutate(remove = ifelse(nchar(A1)>1|nchar(A2)>1,1, 0))%>%
           filter(remove==0)
######load GWAS summary statistics
AD_dataset = fread(paste0(AD_GWAS_path,"/ADGWAS_Bellenguez_2022.",chr,"/ADGWAS2022.chr",chr,".sumstat.tsv",sep=""))
##transform the AD_dataset to the form of allele_qc format
AD_data = AD_dataset%>%mutate(chr = paste0("chr",chromosome))%>%mutate(pos = as.character(position))%>%select(-chromosome)%>%rename("A1"="ref","A2"="alt")%>%mutate(z=beta/se)%>%select(-position)
AD_allele_flip = allele_qc(AD_data,qtl_reference,match.min.prop=0.2,remove_dups=TRUE,flip=TRUE,remove=TRUE)%>%
              mutate(variant_allele_flip = paste(chr,pos,A1.sumstats,A2.sumstats,sep=":"))
####load LD matrix
LD.files <- bedtoolsr::bt.intersect(a = LD_list_pos.bed, b = qtl_select)
LD.files.name = unique(LD.files$V5)
LD.files.name
LD.list = list()
LD.matrix.names=NULL
for (k in 1:length(LD.files.name)){
 npz = np$load(LD.files.name[k])
 LD.matrix = npz$f[["arr_0"]]
 LD.snps = str_split(LD.files.name[k],"[.]",simplify = T)%>%.[,-c(length(.),(length(.)-1))]%>%paste(.,collapse=".")%>%paste0(.,".bim",sep="")%>%read.table(.)
 #head(LD.snps)
 LD_names = colnames(LD.matrix) = rownames(LD.matrix) = gsub("_",":",LD.snps$V2)
 snp_merge = intersect(LD_names,AD_allele_flip$variant_allele_flip)
 LD.select = as.matrix(LD.matrix[snp_merge,snp_merge])
 LD.list[[k]] = LD.select
 LD.matrix.names = append(LD.matrix.names,snp_merge)
}
   LD.block = as.matrix(bdiag(LD.list))
   upperTriangle(LD.block,byrow=TRUE) = lowerTriangle(LD.block)
   colnames(LD.block) = rownames(LD.block) = LD.matrix.names
####generate the twas_z format input
twas_z_format = data.frame(LD.matrix.names)%>%mutate(gene_name = gene_name)%>%
mutate(AD_allele_flip[match(LD.matrix.names,AD_allele_flip$variant_allele_flip),]%>%select(beta,se,z))%>%
mutate(susie_weights = susie_res[[1]]$susie_weights[match(LD.matrix.names,susie_res[[1]]$variant_names)])%>%
mutate(enet_weights =  susie_res[[1]]$enet_weights[match(LD.matrix.names,susie_res[[1]]$variant_names)])%>%
mutate(lasso_weights = susie_res[[1]]$lasso_weights[match(LD.matrix.names,susie_res[[1]]$variant_names)])%>%
mutate(mr_ash_weights = susie_res[[1]]$mr_ash_weights[match(LD.matrix.names,susie_res[[1]]$variant_names)])
####calculate the pvalue and zscore using twas_z function    
weights = apply(twas_z_format[,c("susie_weights","enet_weights","lasso_weights","mr_ash_weights")],2,function(x) twas_z(x,twas_z_format$z,R = LD.block))
twas_weights = data.frame(gene_name,weights$susie_weights$pval,weights$susie_weights$z,weights$lasso_weights$pval,weights$lasso_weights$z,
                       weights$enet_weights$pval,weights$enet_weights$z,weights$mr_ash_weights$pval,weights$mr_ash_weights$z)
names(twas_weights) = c("gene_name","susie_pval","susie_z","lasso_pval","lasso_z","enet_pval","enet_z","mr_ash_pval","mr_ash_z")
return(twas_weights)
}