## genotype matrix preprocesing
compute_maf <- function(geno){
  f <- mean(geno,na.rm = TRUE)/2
  return(min(f, 1-f))
}

compute_missing <- function(geno){
  miss <- sum(is.na(geno))/length(geno)
  return(miss)
}

mean_impute <- function(geno){
  f <- apply(geno, 2, function(x) mean(x,na.rm = TRUE))
  for (i in 1:length(f)) geno[,i][which(is.na(geno[,i]))] <- f[i]
  return(geno)
}

is_zero_variance <- function(x) {
  if (length(unique(x))==1) return(T)
  else return(F)
}

filter_X <- function(X, missing_rate_thresh, maf_thresh) {
    rm_col <- which(apply(X, 2, compute_missing) > missing_rate_thresh)
    if (length(rm_col)) X <- X[, -rm_col]
    rm_col <- which(apply(X, 2, compute_maf) <= maf_thresh)
    if (length(rm_col)) X <- X[, -rm_col]
    rm_col <- which(apply(X, 2, is_zero_variance))
    if (length(rm_col)) X <- X[, -rm_col]
    return(mean_impute(X))
}

#' @importFrom plink2R read_plink
#' @importFrom readr read_delim
load_regional_association_data <- function(genotype, # PLINK file
                                           phenotype, # a vector of phenotype file names 
                                           covariate, # a vector of covariate file names corresponding to the phenotype file vector
                                           region, # a string of chr:start-end
                                           conditions, # a vector of strings
                                           maf_cutoff = 0,
                                           mac_cutoff = 0,
                                           imiss_cutoff = 0,
                                           y_as_matrix = FALSE,
                                           keep_indel = TRUE) {
    library(plink2R)
    library(dplyr)
    library(readr)
    library(stringr)
    library(purrr)

    ## Load genotype
    geno = read_plink(genotype)
    rownames(geno$bed) = read.table(text = rownames(geno$bed), sep= ":")$V2
    ## if indel is false, remove the indel in the genotype
    if (!keep_indel){
    geno_bim = geno$bim%>%rename("chrom" = "V1","variant_id" = "V2","alt" = "V5","ref"="V6")%>%mutate(indel = ifelse(grepl("[^ATCG]",alt)=="TRUE"|grepl("[^ATCG]",ref)=="TRUE"|nchar(alt)>1|nchar(ref)>1,1, 0))
    geno_bed = geno$bed[,geno_bim$indel==0]}
    else {
    geno_bed = geno$bed
    }
    ## Load phenotype and covariates and perform some pre-processing
    ### including Y ( cov ) and specific X and covar match, filter X variants based on the overlapped samples.
    data_list = tibble(covariate_path = covariate, phenotype_path =phenotype) %>%
        mutate(covar = map(covariate_path, ~read_delim(.x,"\t")%>%select(-1)%>%na.omit%>%t()),
        Y = map2(phenotype_path,covar, ~{
          y_data <- tabix_region(.x, region)%>%select(-4)%>%select(rownames(.y))%>%t()%>%as.matrix
          return(y_data)
          }),
        Y = map(Y, ~.x%>%na.omit),    # remove na where Y raw data has na which block regression
        dropped_sample = map2(covar, Y , ~rownames(.x)[!rownames(.x) %in% rownames(.y)]),
        covar = map2(covar, Y , ~.x[intersect(.x%>%rownames,rownames(.y)),]), # remove the dropped samples from Y
        X_data = map(covar,~ filter_X( geno_bed[intersect(rownames(.x),rownames(geno_bed)),], imiss_cutoff, max(maf_cutoff, mac_cutoff/(2*length(intersect(rownames(.x),rownames(geno_bed))) ) ))   ))
              
    ## Get residue Y for each of condition and its mean and sd
    data_list = data_list%>%mutate(Y_resid_mean = map2(Y,covar,~.lm.fit(x = cbind(1,.y), y = .x)$residuals%>%mean),
                               Y_resid_sd = map2(Y,covar,~.lm.fit(x = cbind(1,.y), y = .x)$residuals%>%sd),
                               Y_resid = map2(Y,covar,~.lm.fit(x = cbind(1,.y), y = .x)$residuals%>%scale%>%t%>%as_tibble)) ## T so that it can be unnest
    if(y_as_matrix) {
        Y_resid = data_list%>%select(Y_resid)%>%tidyr::unnest(Y_resid)%>%t%>%as.matrix
        colnames(Y_resid) = conditions
        print(paste("Dimension of Y matrix:", nrow(Y_resid), ncol(Y_resid)))
    } else {
        Y_resid = map(data_list$Y_resid,~.x%>%t) # Transpose back 
        names(Y_resid) = conditions
    }
    # Get X matrix for union of samples
    all_samples = map(data_list$covar, ~rownames(.x))%>%unlist%>%unique()
    maf_cutoff = max(maf_cutoff,mac_cutoff/(2*length(all_samples)))
    X = filter_X(geno_bed[all_samples,], imiss_cutoff, maf_cutoff) ## Filter X for mvSuSiE
    #
    maf_list = lapply(data_list$X_data, function(x) apply(x, 2, compute_maf))
    ## Get residue X for each of condition and its mean and sd
    print(paste0("Dimension of input genotype data is row:", nrow(X), " column: ", ncol(X) ))
    X_list = data_list%>%mutate(X_resid_mean= map2(X_data,covar,~.lm.fit(x = cbind(1,.y), y = .x)$residuals%>%data.frame()%>%apply(.,2,mean)),
                               X_resid_sd= map2(X_data,covar,~.lm.fit(x = cbind(1,.y), y = .x)$residuals%>%data.frame()%>%apply(.,2,sd)),
                               X_resid = map2(X_data,covar,~.lm.fit(x = cbind(1,.y), y = .x)$residuals%>%scale))%>%select(X_resid_mean,X_resid_sd,X_resid)
    ## residual_Y_scaled: if y_as_matrix is true, then return a matrix of R conditions, with column names being the names of the conditions (phenotypes) and row names being sample names. Even for one condition it has to be a matrix with just one column. if y_as_matrix is false, then return a list of y either vector or matrix (CpG for example), and they need to match with residual_X_scaled in terms of which samples are missing.
    ## residual_X_scaled: is a list of R conditions each is a matrix, with list names being the names of conditions, column names being SNP names and row names being sample names.
    ## X: is the somewhat original genotype matrix output from `filter_X`, with column names being SNP names and row names being sample names. Sample names of X should match example sample names of residual_Y_scaled matrix form (not list); but the matrices inside residual_X_scaled would be subsets of sample name of residual_Y_scaled matrix form (not list).
    return (list(
            residual_Y_scaled = Y_resid,
            residual_X_scaled = X_list$X_resid,
            residual_Y_sd = data_list$Y_resid_sd,
            residual_X_sd = X_list$X_resid_sd,
            dropped_sample = data_list$dropped_sample,
            covar = data_list$covar,
            Y = data_list$Y,
            X = X,
            X_data = data_list$X_data,
            maf = maf_list
            ))
}

load_regional_finemapping_data <- function(...) {
  dat <- load_regional_association_data(...)
  return (list(
          residual_Y_scaled = dat$residual_Y_scaled,
          residual_X_scaled = dat$residual_X_scaled,
          residual_Y_sd = dat$residual_Y_sd,
          residual_X_sd = dat$residual_X_sd,
          X = dat$X,
          dropped_sample = dat$dropped_sample,
          maf = dat$maf
          ))
}

load_regional_regression_data <- function(...) {
  dat <- load_regional_association_data(...)
  return (list(
          Y = dat$Y,
          X_data = dat$X_data,
          covar = dat$covar,
          dropped_sample = dat$dropped_sample,
          maf = dat$maf
          ))
}