
#FIXME: need to be updated with roxygen documentation
source("~/githubrepo/twas_rosmap/workflow/mrmash_utility_functions.R")


#input data processing
mrmash_input <- function(gene_rds, sample_id_list, data_dir, wd, k_fold, n_con){
    sample_folds <- sample_partition(sample_id_list, wd, k_fold)
    #cv_processing # package: MBSP - run without container
    sumstat_ls <- do.call(c, lapply(gene_rds, cv_processing, 
                                    partition=sample_folds, 
                                    data_dir=data_dir, wd=wd)) #list of summary stats file name 
    
    #compute prior_grids for all folds
    prior_grid_files <- do.call(c, lapply(1:k_fold, compute_grid, 
                                          sumstat_ls=sumstat_ls, 
                                          wd=wd, n_con=n_con))
    return(list(sample_folds=sample_folds,prior_grid_files=prior_grid_files, 
                sumstat_ls=sumstat_ls))  
}  



#fit mr.mash, compute accuracy, and filter  
mrmash_run <- function(gene, sample_folds, wd, data_dir, k_fold, prior_grid_files, mixture_pr_file){
    #run for all folds
    lapply(1:k_fold, mrmash_prediction, gene=gene, prior_grid_files=prior_grid_files, 
               wd=wd, data_dir=data_dir, sample_folds=sample_folds, 
               mixture_pr_file=mixture_pr_file)
                #save mr.mash weights
    
    prediction_accuracy(gene=gene, wd=wd, data_dir=data_dir, k_fold=k_fold)
                #save prediction scores
}


#select genes based on prediction accuracy results
select_genes <- function(gene_rds, wd){
    pred_scors <-  do.call(rbind, lapply(gene_rds,summary_pred_score, wd=wd))
    pred_scors <- pred_scors[pred_scors$r2 >= 0.01 & pred_scors$f.pval <=0.05, ]
    gene_selected <- basename(unique(gsub("_score.rds", "", pred_scors$fname, perl = TRUE)))
    return(gene_selected)
    }

