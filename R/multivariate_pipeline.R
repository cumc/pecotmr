
# FIXME: use a single prior for input into mvsusie and mr.mash 
# FIXME: variant selection: add in common variants into max_cv_variants selection
#' 
#' Multivariate Pipeline 
#' 
#' This function performs weights computation for Transcriptome-Wide Association Study (TWAS) with fitting 
#' models using mvSuSiE and mr.mash with the option of using limited number of variants selected from 
#' mvSuSiE fine-mapping for computing TWAS weights with cross validation.  
#'
#' @param X A matrix of genotype data where rows represent samples and columns represent genetic variants.
#' @param y A matrix of phenotype measurements, represent samples and columns represent condidtions.
#' @param maf A list of vectors for minor allele frequencies for each variant in X
#' @param mvsusie_fit An object returned by the mvSuSiE function, containing the mvSuSiE model fit.
#' @param ld_reference_meta_file An optional path to a file containing linkage disequilibrium reference data. If provided, variants in X are filtered based on this reference.
#' @param X_scalar Scalar for the genotype data, used in residual scaling at the step of susie_post_processor. Defaults to 1 (no scaling).
#' @param y_scalar Scalar for the phenotype data, used in residual scaling at the step of susie_post_processor. Defaults to 1 (no scaling).
#' @param dropped_sample A list of list with X, Y, and covar, each list contains list of dropped samples as a list of vectors for each condition. Can be obtained from the output from load_regional_multivariate_data. 
#' @param mrmash_weights_prior_matrices A list of data-driven covariance matrices.
#' @param prior_canonical_matrices If set to TRUE, will compute canonical covariance matrices add into prior covariance matrice list in mrmash_wrapper. 
#' @param cv_folds The number of folds to use for cross-validation. Set to 0 to skip cross-validation.
#' @param max_cv_variants The maximum number of variants to be included in cross-validation. Defaults to 5000.
#' @param signal_cutoff Cutoff value for signal identification in PIP values for susie_post_processor. 
#' @param secondary_coverage A vector of secondary coverage probabilities for credible set refinement. Defaults to c(0.7, 0.5).
#' @param cv_seed The seed for random number generation in cross-validation. Defaults to 999.
#' @param cv_threads The number of threads to use for parallel computation in cross-validation. Defaults to 1.
#' 
#' @importFrom mvsusieR mvsusie create_mixture_prior 
#' @export
run_multivariate_pipeline <- function(X, y, maf, 
                                      dropped_sample,
                                      max_L=30,
                                      ld_reference_meta_file = NULL,  
                                      X_scalar = rep(1, ncol(y)), 
                                      y_scalar = rep(1, ncol(y)), 
                                      pip_cutoff_to_skip = rep(0, ncol(y)), 
                                      signal_cutoff = 0.025, #based on the susie_twas pipeline
                                      secondary_coverage = c(0.5, 0.7), 
                                      mrmash_weights_prior_matrices=NULL,
                                      mrmash_weights_prior_matrices_cv=NULL,
                                      prior_canonical_matrices=TRUE,
                                      sample_partition = NULL, 
                                      mrmash_max_iter=5000,
                                      mvsusie_max_iter=200, 
                                      max_cv_variants=5000, 
                                      cv_folds=5, 
                                      cv_threads = 1,
                                      cv_seed = 999,...){
    for(r in 1:ncol(y)){
      if (pip_cutoff_to_skip[r]>0) {
          top_model_pip <- susie(X, y[, r],L=1)$pip
          if (!any(top_model_pip>pip_cutoff_to_skip[r])){
             message(paste0("skipping condition ", colnames(y)[r], ", because all top_model_pip < pip_cutoff_to_skip=",
                  pip_cutoff_to_skip[r],  ", top loci model does not show any potentially significant variants. "))     
             y[, r] <- NA
           }
        }
     }
    # if all columns of y is NA, return empty set
    if (all(unlist(lapply(1:ncol(y), function(x){is.null(y[, x])})))) return(list())
    drop_indx <- which(colSums(is.na(y)) == nrow(y))

    if(length(drop_indx)>=1){
        y <- y[, -c(drop_indx)]
        maf <- maf[-drop_indx]
        X_scalar <- X_scalar[-drop_indx]
        y_scalar <- y_scalar[-drop_indx]
        dropped_sample <- lapply(dropped_sample, function(x) x[-drop_indx])  
        # filter prior_matrices                        
        if(!is.null(mrmash_weights_prior_matrices)) {
          mrmash_weights_prior_matrices = lapply(prior_matrices, function(x) x[colnames(y), colnames(y)])
        }
        if(!is.null(mrmash_weights_prior_matrices_cv)){ 
          mrmash_weights_prior_matrices_cv=lapply(prior_matrices_cv, function(j){lapply(j, function(r){r[colnames(y), colnames(y)]})})
        }
                                                 
    }
                                 
    # convert into per condition dropped_sample list with $X, $Y, $covar
    dropped_samples <-  lapply(1:ncol(y), function(r){list(X=dropped_sample$X[[r]],
                                    y=dropped_sample$Y[[r]], covar=dropped_sample$covar[[r]])})
             
    
    prior = create_mixture_prior(R=ncol(y))
    mrmash_output <- mrmash_wrapper(X=X,Y=y,prior_data_driven_matrices=mrmash_weights_prior_matrices, prior_grid = NULL,
                                             prior_canonical_matrices=prior_canonical_matrices, 
                                             max_iter=mrmash_max_iter) 
    resid_Y = mrmash_output$V  #resid_Y = mr.mash.alpha:::compute_cov_flash(y)

    mvsusie_fitted <- mvsusie(X, 
                        Y=y, 
                        L=max_L,  
                        prior_variance=prior, 
                        residual_variance=resid_Y, 
                        precompute_covariances=F, 
                        compute_objective=T, 
                        estimate_residual_variance=F, 
                        estimate_prior_variance=T, 
                        estimate_prior_method='EM',
                        max_iter = mvsusie_max_iter, 
                        n_thread=1, 
                        approximate=F)

    mvsusie_prefit <- lapply(1:ncol(y), function(x){susie_post_processor(mvsusie_fitted, 
                                                    X, 
                                                    y[,x,drop=FALSE], 
                                                    X_scalar=X_scalar,  
                                                    y_scalar=y_scalar,  
                                                    maf[[x]],
                                                    secondary_coverage = secondary_coverage, 
                                                    signal_cutoff = signal_cutoff, 
                                                    other_quantities = list(dropped_samples = dropped_samples[[x]]))})
    names(mvsusie_prefit) <- colnames(y)                         
    mvsusie_prefit <- lapply(mvsusie_prefit, function(x) x$susie_result_trimmed) 

    result <-  twas_multivariate_weights_pipeline(X, y, maf, mvsusie_prefit, mvsusie_fitted, prior=prior, resid_Y=resid_Y, 
                                                  dropped_samples=dropped_samples, cv_folds=cv_folds, 
                                                  ld_reference_meta_file=ld_reference_meta_file, max_cv_variants=max_cv_variants, 
                                                  mvsusie_max_iter=mvsusie_max_iter, mrmash_max_iter=mrmash_max_iter,signal_cutoff=signal_cutoff,
                                                  pip_cutoff_to_skip =pip_cutoff_to_skip, secondary_coverage=secondary_coverage,
                                                  prior_canonical_matrices = prior_canonical_matrices, 
                                                  mrmash_weights_prior_matrices=mrmash_weights_prior_matrices,
                                                  mrmash_weights_prior_matrices_cv=mrmash_weights_prior_matrices_cv,
                                                  cv_seed=cv_seed)  
    return(result)         
                             
}
                             
#' TWAS Weights Multivariate Pipeline                            
#' @importFrom mvsusieR mvsusie 
#' @export                           
twas_multivariate_weights_pipeline <- function(X, y, maf, mvsusie_prefit, mvsusie_fitted, dropped_samples, prior, resid_Y,
                                               X_scalar=rep(1, ncol(y)),y_scalar=rep(1, ncol(y)), 
                                               max_cv_variants=5000, mvsusie_max_iter=200, mrmash_max_iter=5000, 
                                               pip_cutoff_to_skip = 0, signal_cutoff = 0.025,
                                               secondary_coverage = c(0.5, 0.7), ld_reference_meta_file=NULL,
                                               cv_folds=5, sample_partition = NULL, 
                                               mrmash_weights_prior_matrices=NULL,
                                               mrmash_weights_prior_matrices_cv=NULL,
                                               prior_canonical_matrices = FALSE,cv_seed = 999, cv_threads = 1,...){ 

    max_L <- c()
    for (r in 1:ncol(y)){ 
        if (!is.null(mvsusie_prefit[[r]])) {
            L <- length(which(mvsusie_prefit[[r]]$V > 1E-9))
            max_L <- c(max_L, L + 3)
        } else {
            # susie_fit did not detect anything significant
            max_L <- c(max_L, 2)
        }
    }
    max_L <- unique(max_L)     
    if(length(max_L)>=2)  max_L <- max(max_L) 

    res <- list()     

    #filter variants
    if (!is.null(ld_reference_meta_file)) {
        variants_kept <- filter_variants_by_ld_reference(colnames(X), ld_reference_meta_file)
        X <- X[, variants_kept$data, drop = FALSE]
        maf <- lapply(maf, function(x, idx) {x[idx]}, idx = variants_kept$idx)

        # mvsusie with updated max_L, filtered X, filtered Y
        mvsusie_fitted <- mvsusie(X=X, 
                        Y=y, 
                        L=max_L, # updated 
                        prior_variance=prior, 
                        residual_variance=resid_Y, 
                        precompute_covariances=F, 
                        compute_objective=T, 
                        estimate_residual_variance=F, 
                        estimate_prior_variance=T, 
                        estimate_prior_method='EM',
                        max_iter = mvsusie_max_iter, 
                        n_thread=1, 
                        approximate=F)
    } 

    # susie_post_processor with filtered variants 
    res <- lapply(1:ncol(y), function(x){list(preset_variants_result=list(susie_post_processor(mvsusie_fitted, 
                                            X, 
                                            y[,x,drop=FALSE], 
                                            X_scalar=X_scalar,  
                                            y_scalar=y_scalar, 
                                            maf[[x]],
                                            secondary_coverage = secondary_coverage, 
                                            signal_cutoff = signal_cutoff,
                                            other_quantities = list(dropped_samples = dropped_samples[[x]]))))})
    names(res) <- colnames(y)  

    # twas  
    weight_methods <- list(mrmash_weights = list(prior_data_driven_matrices=mrmash_weights_prior_matrices,
                                                 prior_canonical_matrices=prior_canonical_matrices, 
                                                 max_iter=mrmash_max_iter),
                           mvsusie_weights= list(mvsusie_fit=mvsusie_fitted, prior_variance=prior, residual_variance=resid_Y, 
                                                 L=max_L, max_iter=mvsusie_max_iter))
    twas_weight <- twas_weights(X=X,Y=y,weight_methods=weight_methods)
    twas_predictions <- twas_predict(X=X, twas_weight)

    # split twas results by condition 
    for (i in names(res)){
     res[[i]]$twas_weights <- lapply(twas_weight, function(wgts){wgts[, i]})
     res[[i]]$twas_predictions <- lapply(twas_predictions, function(pred){pred[, i]})
    } 

    # Twas with Cross Validation
    # Select variants for cv by pip values - top `max_cv_variants` number of high pip variants get selected 
    if (max_cv_variants>length(mvsusie_fitted$pip))  max_cv_variants=length(mvsusie_fitted$pip)
    top_sig_idx <- which(mvsusie_fitted$pip %in% mvsusie_fitted$pip[order(-mvsusie_fitted$pip)[1:max_cv_variants]])
    mvsusie_fitted$coef <- mvsusie_fitted$coef[c(1, top_sig_idx+1), ] #first row is intercept - which will later handled in mvsusie.coef()

    if (cv_folds > 0) {
        weight_methods <- list(mrmash_weights = list(), mvsusie_weights= list())
        message(paste0("Performing cross validation with ", length(top_sig_idx), " variants. "))
        
        twas_cv_result <- twas_weights_cv(X=X[, top_sig_idx],Y=y, 
                fold=cv_folds, 
                weight_methods=weight_methods, 
                sample_partition=sample_partition, 
                mrmash_weights_prior_matrices=mrmash_weights_prior_matrices_cv, 
                mrmash_max_iter=mrmash_max_iter, 
                prior_canonical_matrices=prior_canonical_matrices, 
                mvsusie_fit=mvsusie_fitted, residual_variance=resid_Y, L=max_L, 
                mvsusie_max_iter=mvsusie_max_iter, 
                num_threads = cv_threads, seed = cv_seed) 
                # if no mvsusie_fit were provided, mvsusie_weights will need to use all other parameters in the list 
                # to fit mvSuSiE to get mvsusie weight matrix. 

        for (i in names(res)){
             res[[i]]$twas_cv_result$sample_partition <- twas_cv_result$sample_partition
             res[[i]]$twas_cv_result$prediction <- lapply(twas_cv_result$prediction, function(predicted){predicted[, i]})
             res[[i]]$twas_cv_result$performance <- lapply(twas_cv_result$performance, function(perform){perform[i, ]})
             res[[i]]$twas_cv_result$time_elapsed <- twas_cv_result$time_elapsed
        } 
    } 
  
    return(res) 

}