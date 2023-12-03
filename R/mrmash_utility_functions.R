#' mr.mash pipeline utility functions
#'
#'
#' @importFrom mr.mash.alpha compute_V_init extract_missing_Y_pattern impute_missing_Y
#' @importFrom matrixStats colVars rowSds
#' @importFrom glmnet cv.glmnet
#' @importFrom doMC registerDoMC
#' @export
#'
#'
#'
#'



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

             
compute_all_missing_y <- function(y){
  allmiss <- all(is.na(y))
  return(allmiss)
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
    
} 

                     

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
prediction_accuracy <- function(dat_first, dat_input){
    thresh <- 100
    
    ###Extract sample size
    sample_size <- apply(dat_input$y_res, 2, function(x){sum(!is.na(x))})
    sample_size <- data.frame(tissue=names(sample_size), sample_size)

    ###enet accuracy
    accuracy_enet <- compute_accuracy_glmnet(dat_first, sample_size, thresh)
    r2_enet <- do.call(cbind, accuracy_enet$r2)
    pval_enet <- do.call(cbind, accuracy_enet$pval)
    scaled_mse_enet <- do.call(cbind, accuracy_enet$scaled_mse)
    scaled_rmse_enet <- do.call(cbind, accuracy_enet$scaled_rmse)
    
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
    
    pval <- res$f.pval
    summary <- res$mean_r2
    summary$ratio <- (summary$mrmash_first-summary$enet)/summary$enet
    summary$cv.pval <- pval$mrmash_first
    colnames(summary)[which(colnames(summary)=="mrmash_first")] <- "r2"
    
    return(list(summary=summary, res=res))
    # saveRDS(res, paste0(wd,"/pred_score/", strsplit(basename(gene), ".rds")[[1]],"_score.rds"))
    
}

                        

#select genes based on prediction accuracy results
select_genes <- function(gene_rds, wd){
    pred_scors <-  do.call(rbind, lapply(gene_rds,summary_pred_score, wd=wd))
    pred_scors <- pred_scors[pred_scors$r2 >= 0.01 & pred_scors$f.pval <=0.05, ]
    gene_selected <- basename(unique(gsub("_score.rds", "", pred_scors$fname, perl = TRUE)))
    return(gene_selected)
    }