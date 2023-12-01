library("mr.mash.alpha")
library("stringr")
library("flashier")
library("mashr")
library("foreach")
library("MBSP")
#library("udr")


###Set some parameter variables (These should be set in the SoS script)
standardize <- TRUE
standardize_response <- FALSE
nthreads <- 2
missing_rate_cutoff <- 0.05
maf_cutoff <- 0.05
var_cutoff <- 0.05
options(stringsAsFactors = FALSE)
                     
#parameter for fitting mr.mash                     
###Set some parameter variables (These should be set in the SoS script)
n_nonmiss_Y <- 100
canonical_mats <- FALSE
standardize <- TRUE
update_w0 <- TRUE
w0_threshold <- 1e-08
update_V <- TRUE
update_V_method <- "full"
B_init_method <- "enet"
max_iter <- 5000
tol <- 0.01
verbose <- FALSE
save_model <- FALSE
glmnet_pred <- TRUE
analysis_stage <- "first_pass"



###Functions to compute MAF, missing genotype rate, impute missing, and filter X accordingly 
compute_maf <- function(geno){
  f <- mean(geno,na.rm = TRUE)/2
  return(min(f, 1-f))
}

compute_missing <- function(geno){
  miss <- sum(is.na(geno))/length(geno)
  return(miss)
}

compute_non_missing_y <- function(y){
  nonmiss <- sum(!is.na(y))
  return(nonmiss)
}

mean_impute <- function(geno){
  f <- apply(geno, 2, function(x) mean(x,na.rm = TRUE))
  for (i in 1:length(f)) geno[,i][which(is.na(geno[,i]))] <- f[i]
  return(geno)
}

filter_X <- function(X, missing_rate_thresh, maf_thresh, var_thresh) {
  rm_col <- which(apply(X, 2, compute_missing) > missing_rate_thresh)
  if (length(rm_col)) X <- X[, -rm_col]
  rm_col <- which(apply(X, 2, compute_maf) < maf_thresh)
  if (length(rm_col)) X <- X[, -rm_col]
  X <- mean_impute(X)
  rm_col <- which(matrixStats::colVars(X) < var_thresh)
  if (length(rm_col)) X <- X[, -rm_col]
  return(X)
}

                     
filter_Y <- function(Y, n_nonmiss){
  rm_col <- which(apply(Y, 2, compute_non_missing_y) < n_nonmiss)
  if (length(rm_col)) Y <- Y[, -rm_col]
  if(is.matrix(Y)){
    rm_rows <- which(apply(Y, 1, compute_all_missing_y))
    if (length(rm_rows)) Y <- Y[-rm_rows, ]
  } else {
    Y <- Y[which(!is.na(Y))]
  }
  return(Y)
}
                
             
get_center <- function(k,n) {
  ## For given number k, get the range k surrounding n/2
  ## but have to make sure it does not go over the bounds
  if (is.null(k)) {
      return(1:n)
  }
  start = floor(n/2 - k/2)
  end = floor(n/2 + k/2)
  if (start<1) start = 1
  if (end>n) end = n
  return(start:end)
}
             
             
sample_partition <- function(sampleid, wd, k_fold){    
    set.seed(100)
    sampleid_kfold <- sampleid
    sampleid_kfold$fold <-sample(rep(1:k_fold, length.out=nrow(sampleid)))

    write.table(sampleid_kfold, file=paste0(wd, "/sample_partition.txt"), 
           sep = "\t", row.names=F, quote=FALSE)
    return(sampleid_kfold)
}
             
             
             

matxMax <- function(mtx) {
  return(arrayInd(which.max(mtx), dim(mtx)))
}
             
remove_rownames = function(x) {
    for (name in names(x)) rownames(x[[name]]) = NULL
    return(x)
}
             
handle_nan_etc = function(x) {
  x$bhat[which(is.nan(x$bhat))] = 0
  x$sbhat[which(is.nan(x$sbhat) | is.infinite(x$sbhat))] = 1E3
  return(x)
}
             
extract_one_data = function(dat, n_random, n_null, infile) {
    if (is.null(dat)) return(NULL)
    z = abs(dat$Bhat/dat$Shat)
    max_idx = matxMax(z)
    # strong effect samples
    strong = list(bhat = dat$Bhat[max_idx[1],,drop=F], sbhat = dat$Shat[max_idx[1],,drop=F])
    # random samples excluding the top one
    if (max_idx[1] == 1) {
        sample_idx = 2:nrow(z)
    } else if (max_idx[1] == nrow(z)) {
        sample_idx = 1:(max_idx[1]-1)
    } else {
        sample_idx = c(1:(max_idx[1]-1), (max_idx[1]+1):nrow(z))
    }
    random_idx = sample(sample_idx, min(n_random, length(sample_idx)), replace = F)
    random = list(bhat = dat$Bhat[random_idx,,drop=F], sbhat = dat$Shat[random_idx,,drop=F])
    # null samples defined as |z| < 2
    null.id = which(apply(abs(z), 1, max) < 2)
    if (length(null.id) == 0) {
      warning(paste("Null data is empty for input file", infile))
      null = list()
    } else {
      null_idx = sample(null.id, min(n_null, length(null.id)), replace = F)
      null = list(bhat = dat$Bhat[null_idx,,drop=F], sbhat = dat$Shat[null_idx,,drop=F])
    }
    dat = (list(random = remove_rownames(random), null = remove_rownames(null), strong = remove_rownames(strong)))
    dat$random = handle_nan_etc(dat$random)
    dat$null = handle_nan_etc(dat$null)
    dat$strong = handle_nan_etc(dat$strong)
    return(dat)
}
             
reformat_data = function(dat, z_only = TRUE) {
    # make output consistent in format with 
    # https://github.com/stephenslab/gtexresults/blob/master/workflows/mashr_flashr_workflow.ipynb      
    res = list(random.z = dat$random$bhat/dat$random$sbhat, 
              strong.z = dat$strong$bhat/dat$strong$sbhat,  
              null.z = dat$null$bhat/dat$null$sbhat)
    if (!z_only) {
      res = c(res, list(random.b = dat$random$bhat,
       strong.b = dat$strong$bhat,
       null.b = dat$null$bhat,
       null.s = dat$null$sbhat,
       random.s = dat$random$sbhat,
       strong.s = dat$strong$sbhat))
  }
  return(res)
}
             
merge_data = function(res, one_data) {
  if (length(res) == 0) {
      return(one_data)
  } else if (is.null(one_data)) {
      return(res)
  } else {
      for (d in names(one_data)) {
        if (is.null(one_data[[d]])) {
          next
        } else {
            res[[d]] = rbind(res[[d]], one_data[[d]])
        }
      }
      return(res)
  }
}
             

merge_data_cache = function(res, one_data) {
      if (length(res) == 0) {
          return(one_data)
      } else {
          for (d in names(one_data)) {
            res[[d]] = rbind(res[[d]], one_data[[d]])
          }
          return(res)
      }
    }

             
cv_processing <- function(gene, partition, data_dir, wd){
    if(!dir.exists(file.path(wd, "sumstats"))){dir.create(file.path(wd, "sumstats"))}
    
    ###Get fold names
    folds <- sort(unique(partition$fold))
    
    ###List to store results
    res_all_folds <- vector("list", length(folds))
    names(res_all_folds) <- paste0("fold_", folds)
    
    ###Extract Y and (filter) X
    dat <- readRDS(paste0(data_dir, "/", gene))
    Y <- dat$y_res
    X <- filter_X(dat$X, missing_rate_cutoff, maf_cutoff, var_cutoff)
    X <- as.matrix(X[,get_center(NULL, ncol(X))])

    for(i in folds){
      test_ids <- partition[which(partition$fold == i), 1]
      Xtrain <- X[!(rownames(X) %in% test_ids), ]
      Ytrain <- Y[!(rownames(Y) %in% test_ids), ]

      univ_sumstats <- mr.mash.alpha::compute_univariate_sumstats(X=Xtrain, Y=Ytrain, standardize=standardize,
                                                                  standardize.response=standardize_response,
                                                                  mc.cores=nthreads)
      res_all_folds[[i]] <- univ_sumstats
    }
    
    #save summary statistics
    saveRDS(res_all_folds, paste0(wd, "/sumstats/", paste0(strsplit(basename(gene), ".rds")[[1]], "_sumstats_cv.rds")))
    return(paste0(strsplit(basename(gene), ".rds")[[1]], "_sumstats_cv.rds"))
}
             
             
extr_eff <- function(fold, wd, sumstat_ls, n_con){
       set.seed(999)
       res = list() #store temporary result for one fold
       #for one fold for all genes
       for (f in paste0(wd, "/sumstats/", sumstat_ls)) {
          # If cannot read the input for some reason then we just skip it, assuming we have other enough data-sets to use.
          dat = tryCatch(readRDS(f), error = function(e) return(NULL))[[paste0("fold_",fold)]]
          if (is.null(dat)) {
              message(paste("Skip loading file", f, "due to load failure."))
              next
          }
          if (n_con > 0 && (ncol(dat$Bhat) != n_con || ncol(dat$Shat) != n_con)) {
              message(paste("Skip loading file", f, "because it has", ncol(dat$Bhat), "columns different from required", n_con))
               next
              }
              res = merge_data(res, reformat_data(extract_one_data(dat, 4, 4, f), TRUE)) #res has all gene's info now      
            }#end of for - genes  
                         
        #merge and update for one fold for all genes
        dat = list()
        dat = merge_data_cache(dat, res)

        # compute empirical covariance XtX                   
        X = dat$strong.z
        X[is.na(X)] = 0
        dat$XtX = t(as.matrix(X)) %*% as.matrix(X) / nrow(X)
        
        saveRDS(dat, paste0(wd,"/data_driven_matrices/fold_", fold, ".rds"))                
        return(dat) #for one fold
        
    }# extr_eff function end

        

#ed methopd                         
mixture_model <- function(fold, effects, wd){
    #flash
    bhat = effects[[paste0("fold_", fold)]][["strong.z"]]
    sbhat = bhat
    sbhat[!is.na(sbhat)] = 1
    dat = mashr::mash_set_data(bhat,sbhat)
    flash = mashr::cov_flash(dat, factors="default", remove_singleton=FALSE, 
                       output_model=paste0(wd, "/data_driven_matrices/fold_", fold, ".flash.model.rds"))
    saveRDS(flash, paste0(wd, "/data_driven_matrices/fold_", fold,".flash.rds"))
    
    #flash_nonneg
    flash_nonneg = mashr::cov_flash(dat, factors="nonneg", remove_singleton=FALSE, 
                           output_model=paste0(wd, "data_driven_matrices/fold_", fold, ".flash_nonneg.model.rds"))
    saveRDS(flash_nonneg, paste0(wd, "/data_driven_matrices/fold_", fold,".flash_nonneg.rds")) 
    
    #pca
    pca = mashr::cov_pca(dat, 3)
    saveRDS(pca, paste0(wd, "/data_driven_matrices/fold_", fold,".pca.rds"))
    
    
    rds_files <- paste0(wd, "/data_driven_matrices/fold_", fold, c(".rds", ".flash.rds", ".flash_nonneg.rds", ".pca.rds"))
    
    ##ED method
    dat = readRDS(rds_files[1])
    U = list(XtX = dat$XtX)
    for (f in rds_files[2:length(rds_files)]) U = c(U, readRDS(f))
    if (FALSE) {
      V = readRDS('.')
    } else {
      V = cor(dat$null.z)
    }
    # Fit mixture model using ED code by J. Bovy
    mash_data = mashr::mash_set_data(dat$strong.z, V=V)
    message(paste("Running ED via J. Bovy's code for", length(U), "mixture components"))
    res = mashr:::bovy_wrapper(mash_data, U, logfile=paste0(wd, "/data_driven_matrices/fold_", fold, ".ed_bovy"), tol = 1e-06)
    saveRDS(list(U=res$Ulist, w=res$pi, loglik=scan(paste0(wd, "/data_driven_matrices/fold_", fold, ".ed_bovy_loglike.log"))), 
            paste0(wd, "/data_driven_matrices/fold_", fold, ".ed_bovy.rds"))
    return(list(U=res$Ulist, w=res$pi))
}
                         
                         

                         
#ud method                        
mixture_model_ud <- function(fold, effects, wd){
    set.seed(999)
    #flash
    bhat = effects[[paste0("fold_", fold)]][["strong.z"]]
    sbhat = bhat
    sbhat[!is.na(sbhat)] = 1
    dat = mashr::mash_set_data(bhat,sbhat)
    flash = mashr::cov_flash(dat, factors="default", remove_singleton=FALSE, 
                       output_model=paste0(wd, "data_driven_matrices/fold_", fold, ".flash.model.rds"))
    saveRDS(flash, paste0(wd, "/data_driven_matrices/fold_", fold,".flash.rds"))
    
    #flash_nonneg
    flash_nonneg = mashr::cov_flash(dat, factors="nonneg", remove_singleton=FALSE, 
                           output_model=paste0(wd, "data_driven_matrices/fold_", fold, ".flash_nonneg.model.rds"))
    saveRDS(flash_nonneg, paste0(wd, "/data_driven_matrices/fold_", fold,".flash_nonneg.rds")) 
    
    #pca
    pca = mashr::cov_pca(dat, 3)
    saveRDS(pca, paste0(wd, "/data_driven_matrices/fold_", fold,".pca.rds"))
    
    
    rds_files <- paste0(wd, "/data_driven_matrices/fold_", fold, c(".rds", ".flash.rds", ".flash_nonneg.rds", ".pca.rds"))

    dat = readRDS(rds_files[1])
    U = list(XtX = dat$XtX)
    U_scaled = list()
    mixture_components =  c('flash','flash_nonneg','pca')
    scale_only =  c()
    scale_idx = which(mixture_components %in% scale_only )
    for (f in 2:length(rds_files) ) {
        if ((f - 1) %in% scale_idx ) {
          U_scaled = c(U_scaled, readRDS(rds_files[f]))
        } else {
          U = c(U, readRDS(rds_files[f]))
        }
    }
    #
    if (FALSE) {
      V = readRDS('.')
    } else {
      V = cor(dat$null.z)
    }
    # Fit mixture model using udr package
    message(paste("Running TED via udr package for", length(U), "mixture components"))
    f0 = ud_init(X = as.matrix(dat$strong.z), V = V, U_scaled = U_scaled, U_unconstrained = U, n_rank1=0)
    res = ud_fit(f0, X = na.omit(f0$X), control = list(unconstrained.update = "ted", resid.update = 'none', scaled.update = "fa", maxiter=5000, tol.lik = 0.001), verbose=TRUE)
    saveRDS(list(U=res$U, w=res$w, loglik=res$loglik), paste0(wd, "/data_driven_matrices/fold_", fold, ".ted_unconstrained.rds"))
     
}       
                         
                         
###Function to compute the grid
autoselect_mixsd <- function(gmin, gmax, mult=2){
  if (mult == 0) {
    return(c(0, gmax/2))
  }
  else {
    npoint = ceiling(log2(gmax/gmin)/log2(mult))
    return(mult^((-npoint):0) * gmax)
  }
}

###Compute endpoints
compute_grid_endpoints = function(data){
  include = !(data$Shat==0 | !is.finite(data$Shat) | is.na(data$Bhat))
  gmax = grid_max(data$Bhat[include], data$Shat[include])
  gmin = grid_min(data$Bhat[include], data$Shat[include])

  return(list(gmin=gmin, gmax=gmax))
}


###Compute the minimum value for the grid
grid_min = function(Bhat,Shat){
  min(Shat)
}

###Compute the maximum value for the grid
grid_max = function(Bhat,Shat){
  if (all(Bhat^2 <= Shat^2)) {
    8 * grid_min(Bhat,Shat) # the unusual case where we don't need much grid
  } else {
    2 * sqrt(max(Bhat^2 - Shat^2))
  }
}
   
                         
                         
compute_grid <- function(fold, sumstat_ls, wd, n_con){
    set.seed(999)
    grid_mins = c()
    grid_maxs = c()

    for (f in paste0(wd, "/sumstats/", sumstat_ls)) {
      # If cannot read the input for some reason then we just skip it, assuming we have other enough data-sets to use.
      dat = tryCatch(readRDS(f), error = function(e) return(NULL))[[paste0("fold_", fold)]]
      if (is.null(dat)) {
          message(paste("Skip loading file", f, "due to load failure."))
          next
      }
      if (n_con > 0 && (ncol(dat$Bhat) != n_con || ncol(dat$Shat) != n_con)) {
          message(paste("Skip loading file", f, "because it has", ncol(dat$Bhat), "columns different from required", n_con))
          next
      }
      endpoints = compute_grid_endpoints(dat)
      grid_mins = c(grid_mins, endpoints$gmin)
      grid_maxs = c(grid_maxs, endpoints$gmax)

    }

    gmin_tot = min(grid_mins)
    gmax_tot = max(grid_maxs)
    grid = autoselect_mixsd(gmin_tot, gmax_tot, mult=sqrt(2))^2  

    saveRDS(grid, paste0(wd, "/grid/fold_", fold, "_grid.rds"))
    return(paste0(wd, "/grid/fold_", fold, "_grid.rds"))        
}
                                          
                     
                  
                     
###Function to compute initial estimates of the coefficients from group-lasso
compute_coefficients_glasso <- function(X, Y, standardize, nthreads, Xnew=NULL, version=c("Rcpp", "R")){

  version <- match.arg(version)

  n <- nrow(X)
  p <- ncol(X)
  r <- ncol(Y)
  Y_has_missing <- any(is.na(Y))
  tissue_names <- colnames(Y)

  if(Y_has_missing){
    ###Extract per-individual Y missingness patterns
    Y_miss_patterns <- mr.mash.alpha:::extract_missing_Y_pattern(Y)

    ###Compute V and its inverse
    V <- mr.mash.alpha:::compute_V_init(X, Y, matrix(0, p, r), method="flash")
    Vinv <- chol2inv(chol(V))

    ###Initialize missing Ys
    muy <- colMeans(Y, na.rm=TRUE)
    for(l in 1:r){
      Y[is.na(Y[, l]), l] <- muy[l]
    }

    ###Compute expected Y (assuming B=0)
    mu <- matrix(rep(muy, each=n), n, r)

    ###Impute missing Ys
    Y <- mr.mash.alpha:::impute_missing_Y(Y=Y, mu=mu, Vinv=Vinv, miss=Y_miss_patterns$miss, non_miss=Y_miss_patterns$non_miss,
                                          version=version)$Y
  }
  ##Fit group-lasso
  if(nthreads>1){
    doMC::registerDoMC(nthreads)
    paral <- TRUE
  } else {
    paral <- FALSE
  }

  cvfit_glmnet <- glmnet::cv.glmnet(x=X, y=Y, family="mgaussian", alpha=1, standardize=standardize, parallel=paral)
  coeff_glmnet <- coef(cvfit_glmnet, s="lambda.min")

  ##Build matrix of initial estimates for mr.mash
  B <- matrix(as.numeric(NA), nrow=p, ncol=r)

  for(i in 1:length(coeff_glmnet)){
    B[, i] <- as.vector(coeff_glmnet[[i]])[-1]
  }

  ##Make predictions if requested
  if(!is.null(Xnew)){
    Yhat_glmnet <- drop(predict(cvfit_glmnet, newx=Xnew, s="lambda.min"))
    colnames(Yhat_glmnet) <- tissue_names
    res <- list(Bhat=B, Ytrain=Y, Yhat_new=Yhat_glmnet)
  } else {
    res <- list(Bhat=B, Ytrain=Y)
  }
  return(res)
    
} #end of compute_coefficients_glasso

                     


compute_coefficients_univ_glmnet <- function(X, Y, alpha, standardize, nthreads, Xnew=NULL){

  r <- ncol(Y)

  linreg <- function(i, X, Y, alpha, standardize, nthreads, Xnew){
    if(nthreads>1){
      doMC::registerDoMC(nthreads)
      paral <- TRUE
    } else {
      paral <- FALSE
    }

    samples_kept <- which(!is.na(Y[, i]))
    Ynomiss <- Y[samples_kept, i]
    Xnomiss <- X[samples_kept, ]

    cvfit <- glmnet::cv.glmnet(x=Xnomiss, y=Ynomiss, family="gaussian", alpha=alpha, standardize=standardize, parallel=paral)
    coeffic <- as.vector(coef(cvfit, s="lambda.min"))
    lambda_seq <- cvfit$lambda

    ##Make predictions if requested
    if(!is.null(Xnew)){
      yhat_glmnet <- drop(predict(cvfit, newx=Xnew, s="lambda.min"))
      res <- list(bhat=coeffic, lambda_seq=lambda_seq, yhat_new=yhat_glmnet)
    } else {
      res <- list(bhat=coeffic, lambda_seq=lambda_seq)
    }

    return(res)
  }
  out <- lapply(1:r, linreg, X, Y, alpha, standardize, nthreads, Xnew)

  Bhat <- sapply(out,"[[","bhat")

  if(!is.null(Xnew)){
    Yhat_new <- sapply(out,"[[","yhat_new")
    colnames(Yhat_new) <- colnames(Y)
    results <- list(Bhat=Bhat[-1, ], intercept=Bhat[1, ], Yhat_new=Yhat_new)
  } else {
    results <- list(Bhat=Bhat[-1, ], intercept=Bhat[1, ])
  }
  return(results)
    
}  #end of compute_coefficients_univ_glmnet 

    
compute_all_missing_y <- function(y){
  allmiss <- all(is.na(y))
  return(allmiss)
}
                     
load_prior_grid = function(prior_grid, sumstats=NULL) {
  res <- tryCatch(readRDS(prior_grid), 
                   error = function(e) {
                     return(NULL)
                   },
                   warning = function(w) {
                     return(NULL)
                 }
    )
   ###Compute prior covariance
  if(is.null(res) && !is.null(sumstats)){
    res <- autoselect_mixsd(sumstats, mult=sqrt(2))^2
  }
   return(res)
} 

###Filter S0 and w0Drop mixture components with weight equal to 0
filter_S0_w0 <- function(S0, w0, thresh=.Machine$double.eps){
  comps_to_keep <- which(w0 > thresh)
  S0 <- S0[comps_to_keep]
  w0 <- w0[comps_to_keep]

  return(list(S0=S0, w0=w0))
}

###Filter data-driven matrices and summary stats based on tissues used
filter_datadriven_mats_and_sumstats <- function(Y, datadriven_mats, sumstats){
  tissues_to_keep <- colnames(Y)
  #Handle different data structure between udr and Bovy's ed
  if(!is.list(datadriven_mats$U[[1]])){
    datadriven_mats_filt <- lapply(datadriven_mats$U, function(x, to_keep){x[to_keep, to_keep]}, tissues_to_keep)
  } else {
    datadriven_mats_filt <- lapply(datadriven_mats$U, function(x, to_keep){x$mat[to_keep, to_keep]}, tissues_to_keep)
  }
  if(!is.null(sumstats)){
    sumstats_filt <- lapply(sumstats[[1]], function(x, to_keep){x[, to_keep]}, tissues_to_keep)
  } else {
    sumstats_filt <- sumstats
  }

  return(list(datadriven_mats_filt=datadriven_mats_filt, sumstats_filt=sumstats_filt))
}

###Split the data in training and test
split_data <- function(X, Y, gtex_ids_folds, fold){
  test_ids <- gtex_ids_folds[which(gtex_ids_folds$fold == fold), "id"]
  Xtrain <- X[!(rownames(X) %in% test_ids), ]
  Ytrain <- Y[!(rownames(Y) %in% test_ids), ]
  Xtest <- X[rownames(X) %in% test_ids, ]
  Ytest <- Y[rownames(Y) %in% test_ids, ]

  return(list(Xtrain=Xtrain, Ytrain=Ytrain, Xtest=Xtest, Ytest=Ytest))
}

###Compute prior weights from coefficients estimates
compute_w0 <- function(Bhat, ncomps){
  prop_nonzero <- sum(rowSums(abs(Bhat))>0)/nrow(Bhat)
  w0 <- c((1-prop_nonzero), rep(prop_nonzero/(ncomps-1), (ncomps-1)))

  if(sum(w0 != 0)<2){
    w0 <- rep(1/ncomps, ncomps)
  }

  return(w0)
}

###Compute column sums of the posterior assignment probabilities
compute_posterior_weight_colsum <- function(w, S0_labels){
  w1_colsums <- colSums(w)
  if(!is.null(S0_labels)){
    tmp <- rep(0, length(S0_labels))
    names(tmp) <- S0_labels
    tmp[which(S0_labels %in% names(w1_colsums))] <- w1_colsums
    w1_colsums <- tmp
  }
  return(w1_colsums)
}


#load precomputed grid                  
load_prior_grid = function(prior_grid, sumstats=NULL) {
  res <- tryCatch(readRDS(prior_grid), 
                   error = function(e) {
                     return(NULL)
                   },
                   warning = function(w) {
                     return(NULL)
                 }
    )
   ###Compute prior covariance
  if(is.null(res) && !is.null(sumstats)){
    res <- autoselect_mixsd(sumstats, mult=sqrt(2))^2
  }
   return(res)
} 


                     
                     
                     
#function for mrmash prediction 
mrmash_prediction <- function(fold, gene, prior_grid_files, wd, data_dir, sample_folds, mixture_pr_file){
    if(!dir.exists(file.path(wd, "prediction"))){dir.create(file.path(wd, "prediction"))}
    myfold <- paste0("fold_", fold)
    resdir <- paste0(wd, "/prediction")
    if(!dir.exists(file.path(resdir, myfold))){dir.create(file.path(resdir, myfold))}  
    
    ###Read in the data
    dat <- readRDS(paste0(data_dir, "/", gene))
    Y <- filter_Y(dat$y_res, n_nonmiss_Y)
    
    sumstats <- tryCatch(readRDS(paste0("./", strsplit(basename(gene), ".rds")[[1]],"_sumstats_cv.rds")), 
                     error = function(e) {
                       return(NULL)
                     },
                     warning = function(w) {
                       return(NULL)
                     }
    )
    sumstats <- sumstats[myfold]
    
   #read Mixture prior
    tryCatch({
        datadriven_mats <- readRDS(mixture_pr_file)
      }, error = function(e) {
        # FIXME: we can implement it and provide a warning instead
        stop("Default prior is not yet implemented. Please provide a prior to use.")
    })
    w0_init <- tryCatch(readRDS('.'), 
                error = function(e) {
                  message("Prior weights not provided. Computing them from initial estimates of the coefficients.")
                  return(NULL)
                },
                warning=function(w){
                  message("Prior weights not provided. Computing them from initial estimates of the coefficients.")
                  return(NULL)
                }
    ) 


    
    if(is.matrix(Y)){
        X <- filter_X(dat$X, missing_rate_cutoff, maf_cutoff, var_cutoff)
        X <- X[rownames(Y), ]
        w0_init <- NULL # Prior weights not provided. Computing them from initial estimates of the coefficients

        
        ###Drop tissues with < n_nonmiss_Y in data-driven matrices and sumstats
        datadriven_mats_sumstats_filt <- filter_datadriven_mats_and_sumstats(Y=Y, datadriven_mats=datadriven_mats, sumstats=NULL)
        
        S0_data <- datadriven_mats_sumstats_filt$datadriven_mats_filt
        sumstats <- datadriven_mats_sumstats_filt$sumstats_filt
        rm(datadriven_mats_sumstats_filt)
        
        #Read grid
        prior_grid_file <- prior_grid_files[fold]
        prior_grid = load_prior_grid(prior_grid_file, sumstats)
          if(is.null(sumstats) && is.null(prior_grid)){
            # FIXME: we can implement it and provide a warning instead
            stop("Computing summary stats and grid on the fly is not yet implemented. Please provide either proper summary stats path or prior grid file.")
          }
        
        ###Split the data in training and test sets
        dat_split <- split_data(X=X, Y=Y, gtex_ids_folds=sample_folds, fold=fold)
        Xtrain <- dat_split$Xtrain
        Ytrain <- dat_split$Ytrain
        Xtest <- dat_split$Xtest
        Ytest <- dat_split$Ytest
        rm(dat_split)
        
        ###Compute canonical matrices, if requested
        if(canonical_mats){
        S0_can <- mr.mash.alpha::compute_canonical_covs(ncol(Ytrain), singletons=TRUE, hetgrid=c(0, 0.25, 0.5, 0.75, 1))
        S0_raw <- c(S0_can, S0_data)
          } else {
            S0_raw <- S0_data
          } 
        ###Compute prior covariance
        S0 <- mr.mash.alpha::expand_covs(S0_raw, prior_grid, zeromat=TRUE)
        time1 <- proc.time()
        
        ###Compute initial estimates of regression coefficients and prior weights (if not provided)
        if(glmnet_pred){
          Xnew <- Xtest
        } else {
          Xnew <- NULL
        }
        
        if(B_init_method == "enet"){
                  out <- compute_coefficients_univ_glmnet(Xtrain, Ytrain, alpha=0.5, standardize=standardize, nthreads=nthreads, Xnew=Xnew)
                } else if(B_init_method == "glasso"){
                  out <- compute_coefficients_glasso(Xtrain, Ytrain, standardize=standardize, nthreads=nthreads, Xnew=Xnew)
                }
        B_init <- out$Bhat
          
        
        if(is.null(w0_init)){
            external_w0 <- FALSE
            w0 <- compute_w0(B_init, length(S0))
          } else {
            external_w0 <- TRUE
            w0 <- w0_init
          } 
        
        ###Filter prior components based on weights
        comps_filtered <- filter_S0_w0(S0=S0, w0=w0)
        S0 <- comps_filtered$S0
        w0 <- comps_filtered$w0
        rm(comps_filtered)
        
        
        ###Fit mr.mash
        fit_mrmash <- tryCatch({mr.mash.alpha::mr.mash(X=Xtrain, Y=Ytrain, S0=S0, w0=w0, update_w0=update_w0, tol=tol,
                                                 max_iter=max_iter, convergence_criterion="ELBO", compute_ELBO=TRUE,
                                                 standardize=standardize, verbose=verbose, update_V=update_V, update_V_method=update_V_method,
                                                 w0_threshold=w0_threshold, nthreads=nthreads, mu1_init=B_init)
                          },
                         error=function(e) {
                              message("Original mr.mash error message:")
                              message(e)
                              return(NULL)
                          })



        if(!is.null(fit_mrmash)){
        time2 <- proc.time()
        elapsed_time <- time2["elapsed"] - time1["elapsed"]
        ###Make predictions
        Yhat_test <- predict(fit_mrmash, Xtest)
            
        ###Save results
        if(external_w0) {
          resu <- list(Ytest=Ytest, Yhat_test=Yhat_test, elapsed_time=elapsed_time)
        } else {
          if(w0_threshold > 0){
            S0_labels = names(S0)
          } else {
            S0_labels = NULL
          } 
          w1_colsums <- compute_posterior_weight_colsum(fit_mrmash$w1, S0_labels)
          resu <- list(w1_colsums=w1_colsums, Bhat=fit_mrmash$mu1, Ytest=Ytest, Yhat_test=Yhat_test, Ybar_train=colMeans(Ytrain, na.rm=TRUE), elapsed_time=elapsed_time)
        }

        if(save_model){
          resu$model <-  fit_mrmash
        }

        if(glmnet_pred){
          resu$Yhat_test_glmnet <- out$Yhat_new
            }
        } #end if(!is.null(fit_mrmash))
        saveRDS(resu, paste0(resdir, "/", myfold, "/", strsplit(basename(gene), ".rds")[[1]], "_", myfold ,"_mrmash.first_pass.rds")) 
      
    }else{ #else if(is.matrix(Y))
        message("Filtered Y has fewer than 2 tissues. Gene will not be analyzed.") 
        }
        
    return(NULL)
    #}, error=function(e){cat("ERROR :",conditionMessage(e), paste0(gene), "\n")}) #end of tryCatch which start with ###Read in the data 
} 
                     
                     
                     
                     
###Function to load the data
load_data <- function(path, nfolds, gene, model_suffix){
  dat_list <- list()
  for(i in 1:nfolds){
    dat <- readRDS(paste0(path, "fold_", i, "/", gene, "_fold_", i, "_", model_suffix, ".rds"))
    if(is.null(dat)){
      dat_list[[i]] <- NA
    } else {
      dat_list[[i]] <- dat
    }
  }
  return(dat_list)
}                     


###Function to compute accuracy
compute_accuracy <- function(Y, Yhat) {
  bias <- rep(as.numeric(NA), ncol(Y))
  names(bias) <- colnames(Y)
  r2 <- rep(as.numeric(NA), ncol(Y))
  names(r2) <- colnames(Y)
  pval <- rep(as.numeric(NA), ncol(Y))
  names(pval) <- colnames(Y)
  mse <- rep(as.numeric(NA), ncol(Y))
  names(mse) <- colnames(Y)
  rmse <- rep(as.numeric(NA), ncol(Y))
  names(rmse) <- colnames(Y)

  for(i in 1:ncol(Y)){
    dat <- na.omit(data.frame(Y[, i], Yhat[, i]))
    colnames(dat) <- c("Y", "Yhat")
    fit  <- lm(Y ~ Yhat, data=dat)
    bias[i] <- coef(fit)[2]
    r2[i] <- summary(fit)$r.squared
    pval[i] <-  ifelse(sd(dat$Yhat)!=0, summary(fit)$coefficients[2,4], runif(1))
    mse[i] <- mean((dat$Y - dat$Yhat)^2)
    rmse[i] <- sqrt(mse[i])
  }

  return(list(bias=bias, r2=r2, pval=pval, mse=mse, rmse=rmse))
}



###Function to compute accuracy of glmnet for all the folds
compute_accuracy_glmnet <- function(dat, sample_size, thresh){
  r2 <- vector("list", length(dat))
  pval <- vector("list", length(dat))
  scaled_mse <- vector("list", length(dat))
  scaled_rmse <- vector("list", length(dat))

  for(i in 1:length(dat)){
    if(length(dat[[i]]) == 1 && is.na(dat[[i]])){
      r2[[i]] <- rep(as.numeric(NA), sum(sample_size[, 2] > thresh))
      pval[[i]] <- rep(as.numeric(NA), sum(sample_size[, 2] > thresh))
      scaled_mse[[i]] <- rep(as.numeric(NA), sum(sample_size[, 2] > thresh))
      scaled_rmse[[i]] <- rep(as.numeric(NA), sum(sample_size[, 2] > thresh))
    } else {
      acc_Yhat_test <- compute_accuracy(dat[[i]]$Ytest, dat[[i]]$Yhat_test_glmnet)
      acc_Ybar_train <- compute_accuracy(dat[[i]]$Ytest, matrix(dat[[i]]$Ybar_train, nrow=nrow(dat[[i]]$Ytest), ncol=ncol(dat[[i]]$Ytest), byrow=TRUE))
      r2[[i]] <- acc_Yhat_test$r2
      pval[[i]] <- acc_Yhat_test$pval
      scaled_mse[[i]] <- acc_Yhat_test$mse/acc_Ybar_train$mse
      scaled_rmse[[i]] <- acc_Yhat_test$rmse/acc_Ybar_train$rmse
    }
  }
  return(list(r2=r2, pval=pval, scaled_mse=scaled_mse, scaled_rmse=scaled_rmse))
}

###Function to compute accuracy of mr.mash/mtlasso for all the folds
compute_accuracy_general <- function(dat, sample_size, thresh){
  r2 <- vector("list", length(dat))
  pval <- vector("list", length(dat))
  scaled_mse <- vector("list", length(dat))
  scaled_rmse <- vector("list", length(dat))

  for(i in 1:length(dat)){
    if(length(dat[[i]]) == 1 && is.na(dat[[i]])){
      r2[[i]] <- rep(as.numeric(NA), sum(sample_size[, 2] > thresh))
      pval[[i]] <- rep(as.numeric(NA), sum(sample_size[, 2] > thresh))
      scaled_mse[[i]] <- rep(as.numeric(NA), sum(sample_size[, 2] > thresh))
      scaled_rmse[[i]] <- rep(as.numeric(NA), sum(sample_size[, 2] > thresh))
    } else {
      acc_Yhat_test <- compute_accuracy(dat[[i]]$Ytest, dat[[i]]$Yhat_test)
      acc_Ybar_train <- compute_accuracy(dat[[i]]$Ytest, matrix(dat[[i]]$Ybar_train, nrow=nrow(dat[[i]]$Ytest), ncol=ncol(dat[[i]]$Ytest), byrow=TRUE))
      r2[[i]] <- acc_Yhat_test$r2
      pval[[i]] <- acc_Yhat_test$pval
      scaled_mse[[i]] <- acc_Yhat_test$mse/acc_Ybar_train$mse
      scaled_rmse[[i]] <- acc_Yhat_test$rmse/acc_Ybar_train$rmse
    }
  }
  return(list(r2=r2, pval=pval, scaled_mse=scaled_mse, scaled_rmse=scaled_rmse))
}

                     


#Function of Fisher's method to combine pvalues
fishr_pval <- function(pvals){
    p <- pchisq(-2*sum(log(pvals)), 2*length(pvals), lower.tail = FALSE)
    return(p)
}
                     



                     
#Functions to compute mr.mash prediction accuracy
prediction_accuracy <- function(gene, wd, data_dir, k_fold){
    thresh <- 100
    if(!dir.exists(file.path(wd, "pred_score"))){dir.create(file.path(wd, "pred_score"))}
    
    ###Load the data
    dat_first <- load_data(paste0(wd, "/prediction/"), k_fold, strsplit(basename(gene), ".rds")[[1]], "mrmash.first_pass")
    
    if('NULL' != "NULL"){
        dat_mtlasso <- load_data(paste0(wd, "/prediction/"), k_fold, strsplit(basename(gene), ".rds")[[1]], "mtlasso")
    }
    
    dat_input <- readRDS(paste0(data_dir,"/", gene))
    
    ###Extract sample size
    sample_size <- apply(dat_input$y_res, 2, function(x){sum(!is.na(x))})
    sample_size <- data.frame(tissue=names(sample_size), sample_size)

    ###enet accuracy
    accuracy_enet <- compute_accuracy_glmnet(dat_first, sample_size, thresh)
    r2_enet <- do.call(cbind, accuracy_enet$r2)
    pval_enet <- do.call(cbind, accuracy_enet$pval)
    scaled_mse_enet <- do.call(cbind, accuracy_enet$scaled_mse)
    scaled_rmse_enet <- do.call(cbind, accuracy_enet$scaled_rmse)

    ###mr.mash accuracy
    ##First pass
    accuracy_mrmash_first <- compute_accuracy_general(dat_first, sample_size, thresh)
    r2_mrmash_first <- do.call(cbind, accuracy_mrmash_first$r2)
    pval_mrmash_first <- do.call(cbind, accuracy_mrmash_first$pval)
    scaled_mse_mrmash_first <- do.call(cbind, accuracy_mrmash_first$scaled_mse)
    scaled_rmse_mrmash_first <- do.call(cbind, accuracy_mrmash_first$scaled_rmse)

    ###mtlasso accuracy
    if('NULL' != "NULL"){
      accuracy_mtlasso <- compute_accuracy_general(dat_mtlasso, sample_size, thresh)
      r2_mtlasso <- do.call(cbind, accuracy_mtlasso$r2)
      pval_mtlasso <- do.call(cbind, accuracy_mtlasso$pval)
      scaled_mse_mtlasso <- do.call(cbind, accuracy_mtlasso$scaled_mse)
      scaled_rmse_mtlasso <- do.call(cbind, accuracy_mtlasso$scaled_rmse)
    }
    
    ###Combined accuracy
    mean_r2 <- cbind(enet=rowMeans(r2_enet, na.rm=TRUE), mrmash_first=rowMeans(r2_mrmash_first, na.rm=TRUE))
    f.pval <-  cbind(enet = apply(pval_enet, 1, fishr_pval), mrmash_first = apply(pval_mrmash_first, 1, fishr_pval))


    se_r2 <- cbind(enet=matrixStats::rowSds(r2_enet, na.rm=TRUE)/apply(r2_enet, 1, function(x){sum(is.finite(x))}),
                   mrmash_first=matrixStats::rowSds(r2_mrmash_first, na.rm=TRUE)/apply(r2_mrmash_first, 1, function(x){sum(is.finite(x))}))

    sd_r2 <- cbind(enet=matrixStats::rowSds(r2_enet, na.rm=TRUE),
                   mrmash_first=matrixStats::rowSds(r2_mrmash_first, na.rm=TRUE))

    mean_scaled_rmse <- cbind(enet=rowMeans(scaled_rmse_enet, na.rm=TRUE), mrmash_first=rowMeans(scaled_rmse_mrmash_first, na.rm=TRUE))
    se_scaled_rmse <- cbind(enet=matrixStats::rowSds(scaled_rmse_enet, na.rm=TRUE)/apply(scaled_rmse_enet, 1, function(x){sum(is.finite(x))}),
                            mrmash_first=matrixStats::rowSds(scaled_rmse_mrmash_first, na.rm=TRUE)/apply(scaled_rmse_mrmash_first, 1, function(x){sum(is.finite(x))}))

    mean_scaled_mse <- cbind(enet=rowMeans(scaled_mse_enet, na.rm=TRUE), mrmash_first=rowMeans(scaled_mse_mrmash_first, na.rm=TRUE))
    se_scaled_mse <- cbind(enet=matrixStats::rowSds(scaled_mse_enet, na.rm=TRUE)/apply(scaled_mse_enet, 1, function(x){sum(is.finite(x))}),
                            mrmash_first=matrixStats::rowSds(scaled_mse_mrmash_first, na.rm=TRUE)/apply(scaled_mse_mrmash_first, 1, function(x){sum(is.finite(x))}))
    
    
    if('NULL' != "NULL"){
      mean_r2 <- cbind(mean_r2, mtlasso=rowMeans(r2_mtlasso, na.rm=TRUE))
      f.pval <- cbind(f.pval, mrmash_second=apply(pval_mrmash_sec, 1, fishr_pval))

      se_r2 <- cbind(se_r2, mtlasso=matrixStats::rowSds(r2_mtlasso, na.rm=TRUE)/apply(r2_mtlasso, 1, function(x){sum(is.finite(x))}))
      sd_r2 <- cbind(sd_r2, mtlasso=matrixStats::rowSds(r2_mtlasso, na.rm=TRUE))

      mean_scaled_rmse <- cbind(mean_scaled_rmse, mtlasso=rowMeans(scaled_rmse_mtlasso, na.rm=TRUE))
      se_scaled_rmse <- cbind(se_scaled_rmse, mtlasso=matrixStats::rowSds(scaled_rmse_mtlasso, na.rm=TRUE)/apply(scaled_rmse_mtlasso, 1, function(x){sum(is.finite(x))}))

      mean_scaled_mse <- cbind(mean_scaled_mse, mtlasso=rowMeans(scaled_mse_mtlasso, na.rm=TRUE))
      se_scaled_mse <- cbind(se_scaled_mse, mtlasso=matrixStats::rowSds(scaled_mse_mtlasso, na.rm=TRUE)/apply(scaled_mse_mtlasso, 1, function(x){sum(is.finite(x))}))
    }

    if(!all(is.nan(mean_r2))){
       mean_r2 <- data.frame(tissue=rownames(mean_r2), mean_r2)
       mean_r2_sample_size <- merge(mean_r2, sample_size, by="tissue", sort=FALSE, all.x=TRUE)

       f.pval <- data.frame(tissue=rownames(f.pval), f.pval)
       f.pval_sample_size <- merge(f.pval, sample_size, by="tissue", sort=FALSE, all.x=TRUE)

       se_r2 <- data.frame(tissue=rownames(se_r2), se_r2)
       se_r2_sample_size <- merge(se_r2, sample_size, by="tissue", sort=FALSE, all.x=TRUE)

       sd_r2 <- data.frame(tissue=rownames(se_r2), sd_r2)
       sd_r2_sample_size <- merge(sd_r2, sample_size, by="tissue", sort=FALSE, all.x=TRUE)

       ci_r2_enet <- cbind(lower=mean_r2_sample_size$enet-2*se_r2_sample_size$enet,
                           upper=mean_r2_sample_size$enet+2*se_r2_sample_size$enet)

       ci_r2_mrmash_first <- cbind(lower=mean_r2_sample_size$mrmash_first-2*se_r2_sample_size$mrmash_first,
                                   upper=mean_r2_sample_size$mrmash_first+2*se_r2_sample_size$mrmash_first)

       res <- list(mean_r2=mean_r2_sample_size, f.pval=f.pval_sample_size, se_r2=se_r2_sample_size, sd_r2=sd_r2_sample_size, ci_mean_r2_enet=ci_r2_enet, ci_mean_r2_mrmash_first=ci_r2_mrmash_first)


       if('NULL' != "NULL"){
         ci_r2_mtlasso <- cbind(lower=mean_r2_sample_size$mtlasso-2*se_r2_sample_size$mtlasso,
                                upper=mean_r2_sample_size$mtlasso+2*se_r2_sample_size$mtlasso)

         res$ci_mean_r2_mtlasso <- ci_r2_mtlasso
       }
    } else { #else (!all(is.nan(mean_r2)))
      res <- NULL
    }
    
    if(!all(is.nan(mean_scaled_rmse))){
      mean_scaled_rmse <- data.frame(tissue=rownames(mean_scaled_rmse), mean_scaled_rmse)
      mean_scaled_rmse_sample_size <- merge(mean_scaled_rmse, sample_size, by="tissue", sort=FALSE, all.x=TRUE)

      se_scaled_rmse <- data.frame(tissue=rownames(se_scaled_rmse), se_scaled_rmse)
      se_scaled_rmse_sample_size <- merge(se_scaled_rmse, sample_size, by="tissue", sort=FALSE, all.x=TRUE)

      ci_scaled_rmse_enet <- cbind(lower=mean_scaled_rmse_sample_size$enet-2*se_scaled_rmse_sample_size$enet,
                                   upper=mean_scaled_rmse_sample_size$enet+2*se_scaled_rmse_sample_size$enet)

      ci_scaled_rmse_mrmash_first <- cbind(lower=mean_scaled_rmse_sample_size$mrmash_first-2*se_scaled_rmse_sample_size$mrmash_first,
                                           upper=mean_scaled_rmse_sample_size$mrmash_first+2*se_scaled_rmse_sample_size$mrmash_first)

      res$mean_scaled_rmse <- mean_scaled_rmse_sample_size
      res$se_scaled_rmse <- se_scaled_rmse_sample_size
      res$ci_mean_scaled_rmse_enet <- ci_scaled_rmse_enet
      res$ci_mean_scaled_rmse_mrmash_first <- ci_scaled_rmse_mrmash_first


      if('NULL' != "NULL"){
        ci_scaled_rmse_mtlasso <- cbind(lower=mean_scaled_rmse_sample_size$mtlasso-2*se_scaled_rmse_sample_size$mtlasso,
                                        upper=mean_scaled_rmse_sample_size$mtlasso+2*se_scaled_rmse_sample_size$mtlasso)

        res$ci_mean_scaled_rmse_mtlasso <- ci_scaled_rmse_mtlasso
      }
    }
    

    if(!all(is.nan(mean_scaled_mse))){
      mean_scaled_mse <- data.frame(tissue=rownames(mean_scaled_mse), mean_scaled_mse)
      mean_scaled_mse_sample_size <- merge(mean_scaled_mse, sample_size, by="tissue", sort=FALSE, all.x=TRUE)

      se_scaled_mse <- data.frame(tissue=rownames(se_scaled_mse), se_scaled_mse)
      se_scaled_mse_sample_size <- merge(se_scaled_mse, sample_size, by="tissue", sort=FALSE, all.x=TRUE)

      ci_scaled_mse_enet <- cbind(lower=mean_scaled_mse_sample_size$enet-2*se_scaled_mse_sample_size$enet,
                                   upper=mean_scaled_mse_sample_size$enet+2*se_scaled_mse_sample_size$enet)

      ci_scaled_mse_mrmash_first <- cbind(lower=mean_scaled_mse_sample_size$mrmash_first-2*se_scaled_mse_sample_size$mrmash_first,
                                           upper=mean_scaled_mse_sample_size$mrmash_first+2*se_scaled_mse_sample_size$mrmash_first)

      res$mean_scaled_mse <- mean_scaled_mse_sample_size
      res$se_scaled_mse <- se_scaled_mse_sample_size
      res$ci_mean_scaled_mse_enet <- ci_scaled_mse_enet
      res$ci_mean_scaled_mse_mrmash_first <- ci_scaled_mse_mrmash_first

      if("first_pass" == "second_pass"){
        ci_scaled_mse_mrmash_sec <- cbind(lower=mean_scaled_mse_sample_size$mrmash_sec-2*se_scaled_mse_sample_size$mrmash_sec,
                                           upper=mean_scaled_mse_sample_size$mrmash_sec+2*se_scaled_mse_sample_size$mrmash_sec)

        res$ci_mean_scaled_mse_mrmash_sec <- ci_scaled_mse_mrmash_sec
      }

      if('NULL' != "NULL"){
        ci_scaled_mse_mtlasso <- cbind(lower=mean_scaled_mse_sample_size$mtlasso-2*se_scaled_mse_sample_size$mtlasso,
                                        upper=mean_scaled_mse_sample_size$mtlasso+2*se_scaled_mse_sample_size$mtlasso)

        res$ci_mean_scaled_mse_mtlasso <- ci_scaled_mse_mtlasso
      }
    }
    saveRDS(res, paste0(wd,"/pred_score/", strsplit(basename(gene), ".rds")[[1]],"_score.rds"))
    
}

                     

#function to load prediction score
summary_pred_score <- function(gene, wd) {
    file <- paste0(wd, "/pred_score/", strsplit(basename(gene), ".rds")[[1]],"_score.rds")
    dat <- readRDS(file)
    pval <- dat$f.pval
    
    rs <- dat$mean_r2
    rs$ratio <- (rs$mrmash_first-rs$enet)/rs$enet
    rs$fname <- file
    rs$f.pval <- pval$mrmash_first
    colnames(rs)[which(colnames(rs)=="mrmash_first")] <- "r2"
    return(rs)
}
                     
