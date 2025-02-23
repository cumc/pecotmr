

#' Load colocboost regional data: this function will go to file_utils.R
#' 
#' This function loads a mixture data sets for a specific region, including individual-level data (genotype, phenotype, covariate data) 
#' or summary statistics (sumstats, LD). Run \code{load_regional_univariate_data} and \code{load_rss_data} multiple times for different datasets
#' 
#' @section Loading individual level data from multiple corhorts
#' @param region A string of chr:start-end for the phenotype region.
#' @param genotype_list a vector of PLINK bed file containing genotype data.
#' @param phenotype_list A vector of phenotype file names.
#' @param covariate_list A vector of covariate file names corresponding to the phenotype file vector.
#' @param conditions_list_individual A vector of strings representing different conditions or groups.
#' @param match_geno_pheno A vector of index of phentoypes matched to genotype if mulitple genotype PLINK files
#' @param maf_cutoff Minimum minor allele frequency (MAF) cutoff. Default is 0.
#' @param mac_cutoff Minimum minor allele count (MAC) cutoff. Default is 0.
#' @param xvar_cutoff Minimum variance cutoff. Default is 0.
#' @param imiss_cutoff Maximum individual missingness cutoff. Default is 0.
#' @param association_window A string of chr:start-end for the association analysis window (cis or trans). If not provided, all genotype data will be loaded.
#' @param extract_region_name A list of vectors of strings (e.g., gene ID ENSG00000269699) to subset the information when there are multiple regions available. Default is NULL.
#' @param region_name_col Column name containing the region name. Default is NULL.
#' @param keep_indel Logical indicating whether to keep insertions/deletions (INDELs). Default is TRUE.
#' @param keep_samples A vector of sample names to keep. Default is NULL.
#' @param phenotype_header Number of rows to skip at the beginning of the transposed phenotype file (default is 4 for chr, start, end, and ID).
#' @param scale_residuals Logical indicating whether to scale residuals. Default is FALSE.
#' @param tabix_header Logical indicating whether the tabix file has a header. Default is TRUE.
#' 
#' @section Loading summary statistics from multiple corhorts or data set
#' @param sumstat_path_list A vector of file path to the summary statistics.
#' @param column_file_path_list A vector of file path to the column file for mapping.
#' @param LD_meta_file_path_list A vector of path of LD_metadata, LD_metadata is a data frame specifying LD blocks with columns "chrom", "start", "end", and "path". "start" and "end" denote the positions of LD blocks. "path" is the path of each LD block, optionally including bim file paths.
#' @param match_LD_sumstat A vector of index of sumstat matched to LD if mulitple sumstat files
#' @param conditions_list_sumstat A vector of strings representing different sumstats.
#' @param conditions_list_LD A vector of strings representing different LD reference for different sumstat
#' @param n_sample User-specified sample size. If unknown, set as 0 to retrieve from the sumstat file.
#' @param n_case User-specified number of cases.
#' @param n_control User-specified number of controls.
#' @param region The region where tabix use to subset the input dataset.
#' @param target User-specified gene/phenotype name used to further subset the phenotype data.
#' @param target_column_index Filter this specific column for the target.
#' @param comment_string comment sign in the column_mapping file, default is #
#' @param extract_coordinates Optional data frame with columns "chrom" and "pos" for specific coordinates extraction.
#'
#'
#'
#' @return A list containing the individual_data and sumstat_data:
#' individual_data contains the following components if exist
#' \itemize{
#'   \item residual_Y: A list of residualized phenotype values (either a vector or a matrix).
#'   \item residual_X: A list of residualized genotype matrices for each condition.
#'   \item residual_Y_scalar: Scaling factor for residualized phenotype values.
#'   \item residual_X_scalar: Scaling factor for residualized genotype values.
#'   \item dropped_sample: A list of dropped samples for X, Y, and covariates.
#'   \item covar: Covariate data.
#'   \item Y: Original phenotype data.
#'   \item X_data: Original genotype data.
#'   \item X: Filtered genotype matrix.
#'   \item maf: Minor allele frequency (MAF) for each variant.
#'   \item chrom: Chromosome of the region.
#'   \item grange: Genomic range of the region (start and end positions).
#'   \item Y_coordinates: Phenotype coordinates if a region is specified.
#' }
#' sumstat_data contains the following components if exist
#' \itemize{
#'   \item sumstats: A list of summary statistics for the matched LD_info, each sublist contains sumstats, n, var_y from \code{load_rss_data}.
#'   \item LD_info: A list of LD information, each sublist contains combined_LD_variants, combined_LD_matrix, ref_panel  \code{load_LD_matrix}.
#' }
#'
#' @export
load_colocboost_regional_data <-  function(region, # a string of chr:start-end for phenotype region
                                           genotype_list = NULL, # PLINK file
                                           phenotype_list = NULL, # a vector of phenotype file names
                                           covariate_list = NULL, # a vector of covariate file names corresponding to the phenotype file vector
                                           conditions_list_individual = NULL, # a vector of strings
                                           match_geno_pheno = NULL, # a vector of index of phentoypes matched to genotype if mulitple genotype files
                                           maf_cutoff = 0,
                                           mac_cutoff = 0,
                                           xvar_cutoff = 0,
                                           imiss_cutoff = 0,
                                           association_window = NULL,
                                           extract_region_name = NULL,
                                           region_name_col = NULL,
                                           keep_indel = TRUE,
                                           keep_samples = NULL,
                                           keep_variants = NULL,
                                           phenotype_header = 4, # skip first 4 rows of transposed phenotype for chr, start, end and ID
                                           scale_residuals = FALSE,
                                           tabix_header = TRUE,

                                           # sumstat if need
                                           sumstat_path_list = NULL, 
                                           column_file_path_list = NULL, 
                                           LD_meta_file_path_list = NULL,
                                           match_LD_sumstat = NULL, # a vector of index of sumstat matched to LD if mulitple sumstat files
                                           conditions_list_sumstat = NULL,
                                           conditions_list_LD = NULL,
                                           subset = TRUE, 
                                           n_sample = 0, 
                                           n_case = 0, 
                                           n_control = 0, 
                                           target = "",
                                           target_column_index = "", 
                                           comment_string = "#",
                                           extract_coordinates = NULL
                                           ) {
    
    
    if (is.null(genotype_list)&is.null(sumstat_path_list)){
        stop("Data load error. Please make sure at least one data set (sumstat_path_list or genotype_list) exists.")
    }
    
    # - if exist individual level data
    individual_data <- NULL
    if (!is.null(genotype_list)){
        
        if (length(genotype_list)!=1 & is.null(match_geno_pheno)){
            stop("Data load error. Please make sure 'match_geno_pheno' exists if you load data from multiple individual-level data.")
        } else if (length(genotype_list)==1 & is.null(match_geno_pheno)){
            match_geno_pheno <- rep(1, length(phenotype_list))
        }

        # - load individual data from multiple datasets
        n_dataset <- unique(match_geno_pheno)
        for (i_data in 1:n_dataset){
            # extract genotype file name
            genotype <- genotype_list[i_data]
            # extract phenotype and covariate file names
            pos <- which(match_geno_pheno == i_data)
            phenotype <- phenotype_list[pos]
            covariate <- covariate_list[pos]
            conditions <- conditions_list_individual[pos]
            dat <- load_regional_univariate_data(genotype = genotype, phenotype = phenotype,
                                                 covariate = covariate, region = region, 
                                                 association_window = association_window,
                                                 conditions = conditions, xvar_cutoff = xvar_cutoff,
                                                 maf_cutoff = maf_cutoff, mac_cutoff = mac_cutoff,
                                                 imiss_cutoff = imiss_cutoff, keep_indel = keep_indel,
                                                 keep_samples = keep_samples, keep_variants = keep_variants,
                                                 extract_region_name = extract_region_name,
                                                 phenotype_header = phenotype_header,
                                                 region_name_col = region_name_col,
                                                 scale_residuals = scale_residuals)
            if (is.null(individual_data)) {
                individual_data <- dat
            } else {
                individual_data <- lapply(names(dat), function(k) {
                    c(individual_data[[k]], dat[[k]])
                })
                individual_data$chrom <- dat$chrom
                individual_data$grange <- dat$grange
            }
        }
        
    }
    
    # - if exist summstat data
    sumstat_data <- NULL
    if (!is.null(sumstat_path_list)){
        
        if (length(LD_meta_file_path_list)!=1 & is.null(match_LD_sumstat)){
            stop("Data load error. Please make sure 'match_LD_sumstat' exists if you load data from multiple sumstats.")
        } else if (length(LD_meta_file_path_list)==1 & is.null(match_LD_sumstat)){
            match_LD_sumstat <- rep(1, length(sumstat_path_list))
        }
        if (length(sumstat_path_list)==1&is.null(conditions_list_sumstat)) {conditions_list_sumstat = "LD_sumstat"}
        if (length(LD_meta_file_path_list)==1&is.null(conditions_list_LD)) {conditions_list_LD = "LD"}
        
        
        # - load sumstat data from multiple datasets
        n_LD <- unique(match_LD_sumstat)
        for (i_ld in 1:n_LD){
            
            # extract LD meta file path name
            LD_meta_file_path <- LD_meta_file_path_list[i_ld]
            LD_info <- load_LD_matrix(LD_meta_file_path, 
                                      region = association_window, 
                                      extract_coordinates = extract_coordinates)
            # extract sumstat information
            pos <- which(match_LD_sumstat == i_ld)
            conditions <- conditions_list_sumstat[pos]
            sumstats <- lapply(pos, function(ii){
                sumstat_path <- sumstat_path_list[ii]
                column_file_path <- column_file_path_list[ii]
                tmp <- load_rss_data(
                    sumstat_path = sumstat_path, column_file_path = column_file_path,
                    n_sample = n_sample, n_case = n_case, n_control = n_control, target = target, region = association_window,
                    target_column_index = target_column_index, comment_string = comment_string)
                if (!("variant_id" %in% colnames(tmp$sumstats))){
                    tmp$sumstats <- tmp$sumstats %>% rowwise() %>% mutate(variant_id = paste0(c(chrom,pos,A1,A2), collapse = ":"))
                }
                return(tmp)
            })
            names(sumstats) <- conditions
            sumstat_data$sumstats <- c(sumstat_data$sumstats, list(sumstats))
            sumstat_data$LD_info <- c(sumstat_data$LD_info, list(LD_info))
        }
        names(sumstat_data$sumstats) <- names(sumstat_data$LD_info) <- conditions_list_LD
    }
    
    return(list(individual_data = individual_data, 
                sumstat_data = sumstat_data))
    
}



#' The following functions will go to colocboost_pipline.R
#' 
#' 
#' Initial QC for the region data loading from \code{load_regional_data}
#' 
#' This function do the initial QC including: check PIP; check maf for individual_data; check QC and impute for sumstat_data
#' 
#' @section Loading individual level data from multiple corhorts
#' @param region_data A region data loaded from \code{load_regional_data}.
#' @param maf_cutoff A scalar to remove variants with maf < maf_cutoff, dafault is 0.005.
#' @param pip_cutoff_to_skip_ind A vector of cutoff values for skipping analysis based on PIP values for each context. Default is 0.
#' @param skip_region A character vector specifying regions to be skipped in the analysis (optional).
#'                    Each region should be in the format "chrom:start-end" (e.g., "1:1000000-2000000").
#' @param pip_cutoff_to_skip_sumstat A vector of cutoff values for skipping analysis based on PIP values for each sumstat Default is 0.
#' @param qc_method Quality control method to use. Options are "rss_qc", "dentist", or "slalom" (default: "rss_qc").
#' @param impute Logical; if TRUE, performs imputation for outliers identified in the analysis (default: TRUE).
#' @param impute_opts A list of imputation options including rcond, R2_threshold, and minimum_ld (default: list(rcond = 0.01, R2_threshold = 0.6, minimum_ld = 5)).
#'
#'
#' @return A list containing the individual_data and sumstat_data after QC:
#' individual_data contains the following components if exist
#' \itemize{
#'   \item residual_Y: A list of residualized phenotype values (either a vector or a matrix).
#'   \item residual_X: A list of residualized genotype matrices for each condition.
#'   \item residual_Y_scalar: Scaling factor for residualized phenotype values.
#'   \item residual_X_scalar: Scaling factor for residualized genotype values.
#'   \item dropped_sample: A list of dropped samples for X, Y, and covariates.
#'   \item covar: Covariate data.
#'   \item Y: Original phenotype data.
#'   \item X_data: Original genotype data.
#'   \item X: Filtered genotype matrix.
#'   \item maf: Minor allele frequency (MAF) for each variant.
#'   \item chrom: Chromosome of the region.
#'   \item grange: Genomic range of the region (start and end positions).
#'   \item Y_coordinates: Phenotype coordinates if a region is specified.
#' }
#' sumstat_data contains the following components if exist
#' \itemize{
#'   \item sumstats: A list of summary statistics for the matched LD_info, each sublist contains sumstats, n, var_y from \code{load_rss_data}.
#'   \item LD_info: A list of LD information, each sublist contains combined_LD_variants, combined_LD_matrix, ref_panel  \code{load_LD_matrix}.
#' }
#'
#' @export
qc_regional_data <- function(region_data,
                             # - individual
                             maf_cutoff = 0.0005, 
                             pip_cutoff_to_skip_ind = 0,
                             # - sumstat
                             skip_region = NULL,
                             remove_indels = FALSE,
                             pip_cutoff_to_skip_sumstat = 0,
                             qc_method = c("rss_qc", "dentist", "slalom"),
                             impute = TRUE, 
                             impute_opts = list(rcond = 0.01, R2_threshold = 0.6, minimum_ld = 5, lamb = 0.01)){
    
    # - individual level data QC
    individual_data <- region_data$individual_data
    if (!is.null(individual_data)){
        X <- individual_data$residual_X
        Y <- individual_data$residual_Y
        MAF <- individual_data$maf
        # 1. remove maf < maf_cutoff
        # 2. initial check PIP
        individual_data <- data_initial_screen_individual(X = X, Y = Y, MAF = MAF, 
                                                          maf_cutoff = maf_cutoff, 
                                                          pip_cutoff_to_skip = pip_cutoff_to_skip_ind )
    }
    
    # - summary statistics QC
    sumstat_data = region_data$sumstat_data
    if (!is.null(sumstat_data)){
        # - 1. initial check PIP
        sumstat_data <- data_initial_screen_sumstat(sumstat_data, skip_region = skip_region,
                                                    remove_indels = remove_indels, 
                                                    pip_cutoff_to_skip = pip_cutoff_to_skip_sumstat)
        
        # - qc or impute
        sumstat_data <- sumstat_qc(sumstat_data, qc_method = qc_method,
                                   impute = impute, impute_opts = impute_opts)
        
    }
    return(list(individual_data = individual_data, 
                sumstat_data = sumstat_data))

}





############# Individual - level data
#' Add context names to colname of Y if missing
#' @return A list of residual Y
#' @noRd
add_context_to_Y <- function(res_Y){
    res <- lapply(1:length(res_Y), function(iy){
        y <- res_Y[[iy]]
        if (is.null(colnames(y))){ 
            colnames(y) <- names(res_Y)[iy] 
        } else {
            colnames(y) <- paste0(names(res_Y)[iy], "_", colnames(y))
        }
        return(y)
    })
    names(res) <- names(res_Y)
    return(res)
}

#' Remove variants with maf < maf_cutoff
#' @return A list of residual X and maf
#' @export
filter_resX_maf <- function(res_X, maf, maf_cutoff = 0.0005){
    pos <- which(maf < maf_cutoff)
    if (length(pos)!=0){
        res_X <- res_X[, -pos, drop = FALSE] %>% .[, order(colnames(.)), drop = FALSE]
        maf <- maf[-pos]
    }
    return(list(res_X = res_X, maf = maf))
}

#' Initial PIP check for individual level data
#' @return A list of residual Y
#' @export
filter_resY_pip <- function(res_X, res_Y, pip_cutoff_to_skip = 0, context = NULL){
    
    # Initial PIP check
    if (pip_cutoff_to_skip != 0){
        if (pip_cutoff_to_skip < 0){
            # automatically determine the cutoff to use
            pip_cutoff_to_skip <- 3 * 1 / ncol(res_X)
        }
        top_model_pip <- lapply(1:ncol(res_Y), function(i) susieR::susie(res_X, res_Y[,i], L = 1)$pip)
        check_model_pip <- sapply( top_model_pip, function(pip) any(pip > pip_cutoff_to_skip) )
        include_idx <- which(check_model_pip)
        if (length(include_idx)==0){
            message(paste("Skipping follow-up analysis for individual-context", context,  
                          ". No signals above PIP threshold", pip_cutoff_to_skip, "in initial model screening."))
            return(list())
        } else if (length(include_idx) == ncol(res_Y)) {
            message(paste("Keep all individual-phenotypes in context", context,  "."))
        } else {
            exclude_idx <- setdiff(1:ncol(res_Y), include_idx)
            exclude_pheno <- paste(colnames(res_Y)[exclude_idx], collapse = ";")
            message(paste("Skipping follow-up analysis for individual-phenotypes", exclude_pheno, "in context", context,  
                          ". No signals above PIP threshold", pip_cutoff_to_skip, "in initial model screening."))
            res_Y <- res_Y[, include_idx, drop = FALSE] %>% .[, order(colnames(.)), drop = FALSE]
        }
    }
    return(res_Y)
    
}

#' Initial check for all contexts with individual-level data
#' @return A list of individual_data after initial checking and QC:
#' \itemize{
#'   \item X: A list of residual X and ready to do analysis.
#'   \item Y: A list of residual Y and ready to do analysis.
#'   \item MAF: A list of MAF for the remaining genotype.
#' }
#' @export
data_initial_screen_individual <- function(
    X,
    Y,
    MAF,
    maf_cutoff = 0.0005,
    pip_cutoff_to_skip = 0 ){
    
    # - add context to colname of Y
    Y <- add_context_to_Y(Y)
    n_context <- length(X)
    residual_X <- residual_Y <- list()
    keep_contexts <- c()
    for (i_context in 1:n_context){
        resX <- X[[i_context]]
        resY <- Y[[i_context]]
        maf <- MAF[[i_context]]
        context <- names(Y)[i_context]
        # - remove variants with maf < maf_cutoff
        tmp <- filter_resX_maf(resX, maf, maf_cutoff = maf_cutoff)
        resX <- tmp$res_X
        MAF[[i_context]] <- tmp$maf
        # Initial PIP check
        resY <- filter_resY_pip(resX, resY, pip_cutoff_to_skip = pip_cutoff_to_skip[i_context], context = context)
        if (!is.null(resY)){
            residual_X <- c(residual_X, list(resX))
            residual_Y <- c(residual_Y, list(resY))
            keep_contexts <- c(keep_contexts, context)
        }
    }
    if (length(keep_contexts)==0){
        message(paste("Skipping follow-up analysis for all contexts."))
        return(list())
    } else {
        message(paste("Region includes the following contexts after inital screening:", paste(keep_contexts, collapse = ";"),  "."))
        names(residual_X) <- names(residual_Y) <- keep_contexts
        return(list(
            X = residual_X,
            Y = residual_Y,
            MAF = MAF
        ))
    }
    
}
                                  
################# summary statistics
#' Initial PIP check for summary statistics
#' @return A list of sumstat including sumstat and LD passing PIP check.
#' @export            
data_initial_screen_sumstat <- function(
    sumstat_data,
    skip_region = NULL,
    remove_indels = FALSE,
    pip_cutoff_to_skip = 0){
    
    n_LD <- length(sumstat_data$LD_info)
    conditions_sumstat <- names(pip_cutoff_to_skip_sumstat)
    for (i in 1:n_LD){
        LD_data <- sumstat_data$LD_info[[i]]
        sumstats <- sumstat_data$sumstats[[i]]
        # Initial PIP check
        if (all(pip_cutoff_to_skip != 0)){
            ld_condition <- names(sumstat_data$LD_info)[i]
            pos <- grepl(ld_condition, conditions_sumstat)
            pip_cutoff_to_skip_ld <- pip_cutoff_to_skip[pos]

            check_model_pip <- sapply(1:length(sumstats), function(ii){
                sumstat <- sumstats[[ii]]
                n <- sumstat$n
                var_y = sumstat$var_y
                preprocess_results <- rss_basic_qc(sumstat$sumstats, LD_data, skip_region = skip_region, remove_indels = remove_indels)
                sumstat <- preprocess_results$sumstats
                LD_mat <- preprocess_results$LD_mat   
                pip <- susie_rss_wrapper(z = sumstat$z, R = LD_mat, L = 1,
                                         n = n, var_y = var_y)$pip
                if (pip_cutoff_to_skip_ld[ii] < 0){
                    # automatically determine the cutoff to use
                    pip_cutoff_to_skip_ld[ii] <- 3 * 1 / length(LD_data$combined_LD_variants)
                }
                any(pip > pip_cutoff_to_skip_ld[ii])
            })
            include_idx <- which(check_model_pip)
            if (length(include_idx)==0){
                message(paste("Skipping follow-up analysis for all summary statistic based on LD -", ld_condition,  
                              ". No signals above PIP threshold", paste0(pip_cutoff_to_skip_ld,";"), "in initial model screening."))
                sumstat_data$sumstats[i] <- list(NULL)
            } else if (length(include_idx) == length(sumstats)) {
                message(paste("Keep all summary statistics based on LD -", ld_condition,  "."))
            } else {
                exclude_idx <- setdiff(1:length(sumstats), include_idx)
                exclude_pheno <- paste(names(sumstats)[exclude_idx], collapse = ";")
                message(paste("Skipping follow-up analysis for summary statistics", exclude_pheno, "based on LD -", ld_condition,  
                              ". No signals above PIP threshold", paste0(pip_cutoff_to_skip_ld[exclude_idx],";"), "in initial model screening."))
                sumstat_data$sumstats[[i]] <- sumstats[include_idx]
            }
        }
    }
    return(sumstat_data)
    
}
                                      
                           
#' sumstat_data QC, imputation using raiss
#' @return A list of sumstat_data after initial checking and QC:
#' \itemize{
#'   \item sumstats: A list of summary statistics and ready to do analysis.
#'   \item LD_mat: A list of LD matrix and ready to do analysis.
#'   \item LD_match: A vector of strings to indicating sumstats and LD matching (save space since multiple sumstats may link to the same LD matrix).
#' }
#' @export     
sumstat_qc <- function(sumstat_data,
                       qc_method = c("rss_qc", "dentist", "slalom"),
                       impute = TRUE, impute_opts = list(rcond = 0.01, R2_threshold = 0.6, minimum_ld = 5, lamb = 0.01)){

    n_LD <- length(sumstat_data$LD_info)
    conditions_sumstat <- names(pip_cutoff_to_skip_sumstat)
    final_sumstats <- final_LD <- list()
    LD_match <- c()
    for (i in 1:n_LD){
        LD_data <- sumstat_data$LD_info[[i]]
        sumstats <- sumstat_data$sumstats[[i]]
        ld_condition <- names(sumstat_data$LD_info)[i]
        pos <- grepl(ld_condition, conditions_sumstat)
        conditions_sumstat_ld <- conditions_sumstat[pos]
        for (ii in 1:length(sumstats)){
            sumstat <- sumstats[[ii]]
            n <- sumstat$n
            var_y = sumstat$var_y
            
            # Perform quality control
            if (!is.null(qc_method)) {
                qc_results <- summary_stats_qc(sumstat$sumstats, LD_data, n = n, var_y = var_y, method = qc_method)
                sumstat$sumstats <- qc_results$sumstats
                LD_mat <- qc_results$LD_mat
            }
            # Perform imputation
            if (impute) {
                impute_results <- raiss(LD_data$ref_panel, 
                                        sumstat$sumstats, 
                                        LD_data$combined_LD_matrix, 
                                        rcond = impute_opts$rcond, 
                                        R2_threshold = impute_opts$R2_threshold, 
                                        minimum_ld = impute_opts$minimum_ld, 
                                        lamb = impute_opts$lamb)
                sumstat$sumstats <- impute_results$result_filter
                LD_mat <- impute_results$LD_mat
            }
            
            # - check if LD exist
            variants <- colnames(LD_mat)
            if (length(final_LD)==0){
                final_LD <- c(final_LD, list(LD_mat)%>%setNames(conditions_sumstat_ld[ii]))
                final_sumstats <- c(final_sumstats, list(sumstat)%>%setNames(conditions_sumstat_ld[ii]))
                LD_match <- c(LD_match, conditions_sumstat_ld[ii])
            } else {
                final_sumstats <- c(final_sumstats, list(sumstat)%>%setNames(conditions_sumstat_ld[ii]))
                exist_variants <- lapply(final_LD, colnames)
                if_exist <- sapply(exist_variants, function(v) all(variants==v) )
                pos <- which(if_exist)
                if (length(pos)==0){
                   final_LD <- c(final_LD, list(LD_mat)%>%setNames(conditions_sumstat_ld[ii]))
                   LD_match <- c(LD_match, conditions_sumstat_ld[ii])
                } else {
                   LD_match <- c(LD_match, LD_match[pos[1]])
                }
            }
        }
    }
    return(list(sumstats = final_sumstats, LD_mat = final_LD, LD_match = LD_match))
    
}

                                  

#' QUESTION: INPUT `region_data` or `X, Y, MAF, sumstats, LD`?           
colocboost_analysis_pipline <- function(region_data, 
                                        target_trait = NULL,
                                        keep_PR = TRUE, 
                                        rm_IN = TRUE,
                                        # - individual QC
                                        maf_cutoff = 0.0005, 
                                        pip_cutoff_to_skip_ind = 0,
                                        # - sumstat QC
                                        skip_region = NULL,
                                        remove_indels = FALSE,
                                        pip_cutoff_to_skip_sumstat = 0,
                                        qc_method = c("rss_qc", "dentist", "slalom"),
                                        impute = TRUE, 
                                        impute_opts = list(rcond = 0.01, R2_threshold = 0.6, minimum_ld = 5, lamb = 0.01),
                                        ...){
    # - QC for the region_data
    region_data <- qc_regional_data(region_data, maf_cutoff = maf_cutoff, 
                                    pip_cutoff_to_skip_ind = pip_cutoff_to_skip_ind,
                                    skip_region = skip_region,
                                    remove_indels = remove_indels,
                                    pip_cutoff_to_skip_sumstat = pip_cutoff_to_skip_sumstat,
                                    qc_method =qc_method,
                                    impute = impute, 
                                    impute_opts = impute_opts)
    
    # - individual level 
    individual_data <- region_data$individual_data
    if (!is.null(individual_data)){
        X <- individual_data$X
        Y <- individual_data$Y
        null_Y <- which(sapply(Y, is.null))
        if (length(null_Y)!=0&length(null_Y)!=length(Y)){
            message(paste("Skipping follow-up analysis for individual traits", paste(names(Y)[null_Y], collapse = ";"), "."))
            X <- X[-null_Y]
            Y <- Y[-null_Y]
        } else if (length(null_Y)==length(Y)){
            message(paste("Skipping follow-up analysis for all individual traits", paste(names(Y)[null_Y], collapse = ";"), "."))
            X <- NULL
            Y <- NULL
        }
        if (!is.null(Y)){
            Y <- lapply(1:length(Y), function(i){
                y <- Y[[i]]
                events <- colnames(y)
                condition <- names(Y)[i]
                if (all(grepl("sQTL", events))){
                    events_keep <- filter_sQTL_events(events, keep_PR = keep_PR, rm_IN = rm_IN)
                    if (length(events_keep)==length(events)){
                        message(paste("All sQTL events in", condition, "inclued in following analysis based on keep_PR and rm_IN filtering."))
                    } else if (length(events_keep)==0){
                        message(paste("No sQTL events in", condition, "pass the filtering based on keep_PR =", keep_PR,
                                      "and rm_IN", rm_IN, "."))
                        return(NULL)
                    } else {
                        exclude_events <- paste0(setdiff(events, events_keep), collapse = ";")
                        message(paste("Some sQTL events,", exclude_events, "in", condition, "are removed based on keep_PR =", keep_PR,
                                      "and rm_IN =", rm_IN, "."))
                        y <- y[,events_keep,drop=FALSE]
                    }
                }
                lapply(seq_len(ncol(y)), function(j) y[, j, drop=FALSE])%>%setNames(condition)
            })
            dict_YX <- cbind(seq_along(Reduce("c", Y)), rep(seq_along(Y), sapply(Y, length)))
            Y <- Reduce("c", Y)
        }
        
    } else {X <- Y <- dict_YX <- NULL}
    
    # - summary statistics
    sumstat_data = region_data$sumstat_data
    if (!is.null(sumstat_data)){
        sumstats <- lapply(sumstat_data$sumstats, function(ss){
            z <- ss$sumstats$z
            variant <- paste0("chr", ss$sumstats$variant_id)
            n <- ss$n
            data.frame("z" = z, "n" = n, "variant" = variant)
        })
        names(sumstats) <- names(sumstat_data$sumstats)
        LD_mat <- lapply(sumstat_data$LD_mat, function(ld){
            colnames(ld)<-rownames(ld)<-paste0("chr",colnames(ld))
            return(ld)
        })
        LD_match <- sumstat_data$LD_match
        dict_sumstatLD <- cbind(seq_along(sumstats), match(LD_match, names(sumstat_data$LD_mat)))
    } else {sumstats <- LD_mat <- dict_sumstatLD <- NULL}
                                                     
    if (is.null(X)&is.null(sumstats)){
        stop("No data in this region!")
    }
    
    traits <- c(names(Y), names(sumstats))
    if (!is.null(target_trait)){
        if (!(target_trait%in%traits)){
            stop(paste("Trying to run target version for target trait", target_trait, 
                       ". But this target trait is not in analysis traits", paste0(traits, collapse = ";"), "."))
        } else {
            target_idx <- which(traits == target_trait)
        }
    } else { target_idx = NULL }
                                                     
    res_cb <- colocboost(X = X, Y = Y, sumstat = sumstats, LD = LD_mat,
                         dict_YX = dict_YX, dict_sumstatLD = dict_sumstatLD, 
                         traits_names = traits, target_idx = target_idx, ...)
                                                     
    return(res_cb)                        

}
                                                     
                                                     
# sQTL filtering: for sQTL, we keep two strategy before analyzing: i) keep the clusters including PR. ii) remove the clusters including IN.
filter_sQTL_events <- function(events, keep_PR = TRUE, rm_IN = TRUE){
    
    # function to extract cluster number and type
    extract_clu <- function(x) { strsplit(x, ":")[[1]][4] }
    extract_type <- function(x) { strsplit(x, ":")[[1]][5] }
    
    df <- data.frame(events, clu_types = sapply(events, extract_clu), types = sapply(events, extract_type))
    df_summary <- df %>% group_by(clu_types) %>% 
            summarize(cluster_types = paste(sort(unique(types)), collapse = ","),
                      include_IN = ("IN"%in%types),
                      include_PR = ("PR"%in%types))
    df <- df %>% inner_join(df_summary, by = "clu_types")
    df_filter <- df %>% rowwise() %>%
            filter(ifelse(rm_IN, !include_IN, TRUE)&ifelse(keep_PR, include_PR, TRUE))
    return(df_filter$events)

}