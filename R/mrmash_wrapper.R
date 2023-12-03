
# requires mr.mash.alpha package. 
# remotes::install_github("stephenslab/mr.mash.alpha")
#
# mr.mash pipeline 
 


# Load regional association data
# for each gene from analysis unit 
fdat <- load_regional_multivariate_data <- function(gene,
                                            maf_cutoff = 0.05,
                                            mac_cutoff = 0,
                                            var_cutoff = 0.05,
                                            imiss_cutoff = 0.05,  
                                            matrix_y_min_complete = 100, 
                                            keep_indel = TRUE,
                                            data_dir){
            X <- readRDS(paste0(data_dir, "/", gene, ".rds"))$X
            y_res <- readRDS(paste0(data_dir, "/", gene, ".rds"))$y_res
            Y <- filter_Y(y_res, matrix_y_min_complete)      
            X <- filter_X(X, imiss_cutoff, maf_cutoff, var_cutoff)
            X <- X[rownames(Y), ]
            dropped_samples <- setdiff(colnames(Y), colnames(y_res)) 
            return(list(X=X, residual_Y_scaled=Y, dropped_samples=dropped_samples))
    }

#fit mr.mash
#function for mrmash prediction 
#' @importFrom mr.mash.alpha compute_V_init extract_missing_Y_pattern impute_missing_Y compute_canonical_cov mr.mash expand_covs
#' @importFrom matrixStats colVars rowSds
#' @importFrom glmnet cv.glmnet
#' @importFrom doMC registerDoMC
#' @export
mr.mash.pipeline  <- function(X, 
                              Y, 
                              sumstats_file, 
                              prior_matrices, 
                              prior_weights,
                              prior_grid,
                              sample_partition=NULL, 
                              analysis_stage ="first_pass", 
                              fold, 
                              nthreads=2,
                              canonical_mats = FALSE,
                              standardize=TRUE, 
                              update_w0 = TRUE,
                              w0_threshold = 1e-8,
                              update_V=TRUE, 
                              update_V_method = "full", 
                              B_init_method="enet",
                              max_iter=5000, 
                              tol = 0.01, 
                              verbose = FALSE, 
                              save_model=FALSE, 
                              glmnet_pred = TRUE){
    set.seed(999)                      
    if(!dir.exists(file.path(wd, "prediction"))){dir.create(file.path(wd, "prediction"))}
    resdir <- paste0(wd, "/prediction")
    
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
    
    # cross validation parameters
    if(!is.null(sample_partition)) {                
        myfold <- paste0("fold_", fold)
        }else{
        myfold <- "all"
        message("sample_partition not provided, fitting mr.mash without cross validation. ")
        }

    if(!dir.exists(file.path(resdir, myfold))){dir.create(file.path(resdir, myfold))}  
        

    sumstats <- NULL
    datadriven_mats <- prior_matrices
    

    ###Drop tissues with < n_nonmiss_Y in data-driven matrices and sumstats
    datadriven_mats_sumstats_filt <- filter_datadriven_mats_and_sumstats(Y=Y, datadriven_mats=datadriven_mats, sumstats=NULL)
    S0_data <- datadriven_mats_sumstats_filt$datadriven_mats_filt
    sumstats <- datadriven_mats_sumstats_filt$sumstats_filt
    rm(datadriven_mats_sumstats_filt)


    if(is.null(prior_grid) && !is.null(sumstats)){
         res <- autoselect_mixsd(sumstats, mult=sqrt(2))^2
      }
    
    if(is.null(sumstats) && is.null(prior_grid)){
    # FIXME: we can implement it and provide a warning instead
    stop("Computing summary stats and grid on the fly is not yet implemented. Please provide either proper summary stats path or prior grid file.")
      }
        
     if(!is.null(sample_partition)) {     
     ###Split the data in training and test sets
        dat_split <- split_data(X=X, Y=Y, gtex_ids_folds=sample_folds, fold=fold)
        Xtrain <- dat_split$Xtrain
        Ytrain <- dat_split$Ytrain
        Xtest <- dat_split$Xtest
        Ytest <- dat_split$Ytest
        rm(dat_split)
     }
        
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
      if(!is.null(sample_partition)) {
           Xnew <- Xtest
        } else{
        Xnew <- X
        } #end glmnet_pred =TRUE
    } else {
      Xnew <- NULL
    }
        
    if(B_init_method == "enet"){
          if(!is.null(sample_partition)) {
          out <- compute_coefficients_univ_glmnet(Xtrain, Ytrain, alpha=0.5, standardize=standardize, nthreads=nthreads, Xnew=Xnew)
           } else {
              out <- compute_coefficients_univ_glmnet(X, Y, alpha=0.5, standardize=standardize, nthreads=nthreads, Xnew=Xnew)
              }
        } else if(B_init_method == "glasso"){
            if(!is.null(sample_partition)) {
              out <- compute_coefficients_glasso(Xtrain, Ytrain, standardize=standardize, nthreads=nthreads, Xnew=Xnew)
                } else {
                 out <- compute_coefficients_glasso(X, Y, standardize=standardize, nthreads=nthreads, Xnew=Xnew)
                }
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
        if(is.null(sample_partition)){
            Xtrain=X 
            Ytrain=Y
            } #For no CV
    
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
            if(!is.null(sample_partition)){
            Yhat_test <- predict(fit_mrmash, Xtest)
            }
            ###Save results
            if(external_w0) {
              if(!is.null(sample_partition)){
                resu <- list(Ytest=Ytest, Yhat_test=Yhat_test, analysis_time=elapsed_time)
                }
            } else {
              if(w0_threshold > 0){
                S0_labels = names(S0)
              } else {
                S0_labels = NULL
              } 
              w1_colsums <- compute_posterior_weight_colsum(fit_mrmash$w1, S0_labels)
              if(!is.null(sample_partition)){
                  resu <- list(w1_colsums=w1_colsums, Bhat=fit_mrmash$mu1, Ytest=Ytest, Yhat_test=Yhat_test, 
                           Ybar_train=colMeans(Ytrain, na.rm=TRUE), analysis_time=elapsed_time)
                  } else {
                  resu <- fit_mrmash
                  resu$w1_colsums <- w1_colsums
                  resu$analysis_time<- elapsed_time
                  }
            }  
            
            if(save_model){
              resu$model <-  fit_mrmash
            }
            
            if(!is.null(sample_partition)){
                if(glmnet_pred){
                  resu$Yhat_test_glmnet <- out$Yhat_new
                    }
            }
        return(resu)
        } else {return(NULL)} # if(is.null(fit_mrmash)) save nothing
} 
 


# get prediction accuracy score results and summaries          
dat_first <- load_data(pre_wd, k_fold, strsplit(basename(gene), ".rds")[[1]], "mrmash.first_pass")
pred_scores <- prediction_accuracy(dat_first, fdat)


