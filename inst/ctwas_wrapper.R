# requires ctwas package.
# remotes::install_github("xinhe-lab/ctwas",ref = "main")
#
#' Causal inference for TWAS with Summary Statistics
#'
#' @title cTWAS (causal-TWAS) wrapper
#' 
#' @description 
#' 
#' @param weight_matrices a list of weight vectors, one weight vector corresponds
#' to one gene, name of each list is gene name.
#' 
#' @param gwas_sumstat a data frame of gwas summary statistics with the columns of:
#' "id", "A1", "A2", "z". 
#' 
#' @param ld_dir a string, pointing to a directory containing all LD matrix files 
#' (*.cor.xz) and variant information (*.bim). The .bim file has 5 required columns: 
#' "chrom", "id", "pos", "alt", "ref" . 
#'
#' @param wd is a path as working directory that will later create folders for output
#' files and formatted LDs. 
#'
#' @param outputdir a string, the directory to store output
#' 
#' @param harmonize_z TRUE/FALSE. If TRUE, GWAS and eQTL genotype alleles are 
#' harmonized
#'
#' @param harmonize_wgt TRUE/FALSE. If TRUE, GWAS and eQTL genotype alleles are 
#' harmonized
#'
#' @param recover_strand_ambig_wgt TRUE/FALSE. If TRUE, a procedure is used to 
#' recover strand ambiguous variants. If FALSE, these variants are dropped from 
#' the prediction model. 

#' @param strand_ambig_action_z the action to take to harmonize strand ambiguous variants 
#' (A/T, G/C) between the z scores and LD reference. "drop" removes the ambiguous variant 
#' from the z scores. "none" treats the variant as unambiguous, flipping the z score to 
#' match the LD reference and then taking no additional action. "recover" imputes the sign 
#' of ambiguous z scores using unambiguous z scores and the LD reference and flips the z 
#' scores if there is a mismatch between the imputed sign and the observed sign of the z 
#' score. This option is computationally intensive.
#'
#' @param outname a string, the output name. 
#'
#' @return 
#'
#' 
#'
#'
# @importFrom ctwas impute_expr_z ctwas_rss
# @importFrom biomaRt useMart getBM
# @importFrom RSQLite dbConnect dbDisconnect dbWriteTable
#' @importFrom data.table fread read_delim
#' @importFrom readr read_delim
#' @importFrom reshape2 melt 
#' @importFrom foreach foreach

# work in progress 
# library("biomaRt")
# library("DBI")
# library("RSQLite")
# library("data.table")
# library("ctwas")
# library("reshape2")
# library("readr")
# library("foreach")
#
## raw example 
# z_snp <- data.table::fread("/mnt/vast/hpc/csg/cl4215/mrmash/workflow/ctwas_gwas.tsv", 
#     data.table=FALSE, header=TRUE)
#
# wgt_all <- readRDS("/mnt/vast/hpc/csg/cl4215/mrmash/workflow/weight_all.rds")
#
# region_list <- read.table("/mnt/vast/hpc/csg/snuc_pseudo_bulk/eight_celltypes_analysis/output/data_preprocessing/ALL/phenotype_data/ALL.log2cpm.region_list", 
#                          sep="\t", header=FALSE)
# colnames(region_list) <- c("chrom", "start", "stop", "id", "gene_name")
#
#
# k <- ctwas_wrapper(weight_matrices=wgt_all, 
#               gwas_sumstat=z_snp,
#               ld_dir="/mnt/vast/hpc/csg/cl4215/mrmash/workflow/LD_WGS/",
#               wd="/mnt/vast/hpc/csg/cl4215/mrmash/workflow/",
#               region_list=region_list,
#               plot=FALSE,
#               harmonize_z = T,
#               harmonize_wgt = T,
#               outname="test",
#               strand_ambig_action_z = "none",
#               recover_strand_ambig_wgt = T)
# 

ctwas_wrapper <- function(weight_matrices, 
                          gwas_sumstat,
                          ld_dir,
                          wd,
                          region_list,
                          plot=FALSE,
                          harmonize_z = T,
                          harmonize_wgt = T,
                          outname="test",
                          strand_ambig_action_z = "none",
                          recover_strand_ambig_wgt = T,
                          estimate_group_prior = T,
                          estimate_group_prior_var = T,
                          ncore = 16
                          ){
    
    
    # FORMAT WEIGHTS - > make db file
    all_snps <- c()
    for ( i in names(wgt_all)){
        all_snps <- c(all_snps, names(wgt_all[[i]]))
        }
    length(all_snps)

    wgt_table <- do.call(rbind, lapply(names(wgt_all), get_wgt_table, wgt_list=wgt_all)) 
    extra <- get_extra(wgt_list=wgt_all, region_list)
    
    extra <-  extra[extra$gene_type=="protein_coding",,drop=F]
    wgt_table <- wgt_table[wgt_table$gene %in% extra$gene,]
    
    mydb <- dbConnect(RSQLite::SQLite(), paste0(wd, "/weights.db"))
    dbWriteTable(mydb, "weights", wgt_table, overwrite=TRUE)
    dbWriteTable(mydb,"extra", extra, overwrite=TRUE)                            
    dbDisconnect(mydb)
    
    # FORMAT LD
    ifelse(!dir.exists(file.path(paste0(wd), "LD_format")), dir.create(file.path(paste0(wd), "LD_format")), FALSE)
    ld_R_dir <- paste0(wd, "/LD_format/") # output of formatted ld
    
    # get original LD
    lds <- list.files(path=ld_dir, pattern="*.cor.xz.bim", all.files=TRUE, full.names=FALSE)
    lds <- strsplit(lds, "*.cor.xz.bim")

    foreach (i = lds) %do% {
        rvar_file <- data.table::fread(paste0(ld_dir, "/" ,i,".cor.xz.bim"), header=FALSE)
        colnames(rvar_file) <- c("chrom", "id","posg", "pos", "alt", "ref") 
        rvar_file <- rvar_file[,c("chrom", "id", "pos", "alt", "ref")]
        indx <- which(rvar_file$id %in% all_snps)
    
        if (length(indx) >=2){
            rvar <- rvar_file[indx, ]
            rvar$variance <- 1 
            rvar$chrom <- as.integer(rvar$chrom)
            
            ld_fname <- paste0(ld_dir, "/", i,".cor")
            if (!file.exists(ld_fname)){
                system(paste("xz -dk", paste0(ld_fname, ".xz")))
                }
            ld <- data.table::fread(paste0(ld_dir, "/", i,".cor"), sep=" ")
            ld <- as.matrix(ld)
            ld <- ld[indx, indx]
            ld[upper.tri(ld)] <- t(ld)[upper.tri(ld)] 
            start <- min(rvar$pos)
            stop <-  max(rvar$pos)
            CHR <- as.integer(unique(rvar$chrom))
            if(length(CHR)>=2) {stop("has more than one chromosome provided in the LD file")}

            saveRDS(ld, paste0(ld_R_dir, "/ld_chr",CHR, ".R_snp.",start,"_",stop,".RDS"))
            data.table::fwrite(rvar, file = paste0(ld_R_dir, "/ld_chr",CHR, ".R_snp.",start,"_",stop, ".Rvar"), 
                               sep="\t", quote = F)
            }
    }
    
    # get pair-wise covariance /LD for harmonization
    ls <- list.files(path=ld_R_dir, pattern="*.RDS", all.files=TRUE, full.names=FALSE)
    ls <- strsplit(ls, ".RDS")
    
    covariance <- do.call(rbind, lapply(ls ,function(x){
            myld <- readRDS(paste0(ld_R_dir, "/", x, ".RDS"))
            Rvar <- data.table::fread(paste0(ld_R_dir, x, ".Rvar"), header = T)
    
            rownames(myld) <- paste0("chr", Rvar$chrom, "_", sub("^.*?\\:", "", Rvar$id), "_b38")
            colnames(myld) <- rownames(myld)
    
            cov <- reshape2::melt(myld)
            cov <- data.frame(t(apply(cov, 1, sort)))
            cov <- cov[duplicated(cov[, 1 : 2], MARGIN = 1), ]
            cov <- cov[, 3 : 1]
            colnames(cov) <- c("RSID1", "RSID2", "VALUE")
            return(cov)}
       )
    )
    
    data.table::fwrite(covariance, paste0(wd, "/weights.txt.gz"), sep="\t", 
                   quote = FALSE, row.names=FALSE)
    data.table::fwrite(region_list, paste0(wd, "/region_list.txt"), sep="\t", 
                   quote = FALSE, row.names=FALSE)
    
    ifelse(!dir.exists(file.path(paste0(wd), "output")), dir.create(file.path(paste0(wd), "output")), FALSE)
    
    res <- impute_expr_z(z_snp = z_snp,
                     weight = paste0(wd, "/weights.db"), 
                     ld_R_dir = ld_R_dir, 
                     outputdir = paste0(wd, "/output/"), 
                     outname = outname,
                     harmonize_z = harmonize_z,
                     harmonize_wgt = harmonize_wgt,
                     strand_ambig_action_z = strand_ambig_action_z,
                     recover_strand_ambig_wgt = T)

     
    # pars <- ctwas_rss(z_gene=res$z_gene,
    #              z_snp = res$z_snp, 
    #              ld_R_dir = ld_R_dir, 
    #              ld_regions_custom = paste0(wd, "/region_list.txt"), 
    #              outputdir = paste0(wd, "/output/"), 
    #              outname = outname,
    #              ncore = ncore)

    return(res)
    
}



get_wgt_table <- function(gene, wgt_list){
    wgt_table <- data.frame(gene=gene, 
                     rsid=names(wgt_list[[gene]]), 
                     varID=paste0(gsub( ":.*$", "", names(wgt_list[[gene]])) , "_" ,sub("^.*?\\:", "", names(wgt_list[[gene]])),"_b38"),
                     ref_allele=gsub( "_.*$", "", sub("^.*?\\_", "", names(wgt_list[[gene]]))), 
                     eff_allele=sub("^.*?\\_", "", sub("^.*?\\_", "", names(wgt_list[[gene]]))), 
                     weight=wgt_list[[gene]])
    return(wgt_table)
}


get_extra <- function(wgt_list, region_list){
    mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
    gtype_table <- getBM(values = names(wgt_all), attributes = c("ensembl_gene_id", "gene_biotype"), mart=mart)
    
    extra <- data.frame(
        gene=names(wgt_list),
        genename=unlist(lapply(names(wgt_list), function(x) region_list$gene_name[region_list$id == x])),
        n.snps.in.model=unlist(lapply(names(wgt_list), function(x) length(wgt_list[[x]]))),
        gene_type=gtype_table$gene_biotype[gtype_table$ensembl_gene_id %in% names(wgt_list)], 
        pred.perf.R2=NA,
        pred.perf.pval=NA,
        pred.perf.qval=NA
        )
    return(extra)  
}
