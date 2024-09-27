#' @title Quantile TWAS Weight Calculation and QTL Analysis
#'
#' @description
#' This file contains functions for performing Quantile Transcriptome-Wide 
#' Association Studies (TWAS) weight calculations and Quantile QTL analysis.
#' It provides tools for screening quantile regression results, performing
#' LD clumping and pruning, and calculating TWAS weights.
#'
#' @details
#' The main function in this file is `quantile_twas_weight_pipeline`, which
#' orchestrates the entire analysis process. Other functions are helper
#' functions used within the main pipeline.
#'

#' Qrank Score Test Screen
#' @param X Matrix of predictors
#' @param Y Matrix or vector of response variables
#' @param Z Matrix of covariates (optional)
#' @param tau.list Vector of quantiles to be analyzed
#' @param threshold Significance threshold for adjusted p-values
#' @param method Method for p-value adjustment ('fdr' or 'qvalue')
#' @param top_count Number of top SNPs to select
#' @param top_percent Percentage of top SNPs to select
#' @return A list containing various results from the QR screen
#' @importFrom quantreg rq rq.fit.br
#' @importFrom tidyr separate
#' @importFrom dplyr %>% mutate select
#' @export
qr_screen <- function(X, Y, Z = NULL, tau.list, threshold = 0.05, method = 'qvalue', top_count = 10, top_percent = 15) {
    p = ncol(X)
    pvec = rep(NA, p)
    ltau = length(tau.list)
    quantile.pvalue = matrix(NA, nrow = p, ncol = ltau, dimnames = list(colnames(X), paste("p_qr", tau.list, sep = "_")))
    quantile.zscore = matrix(NA, nrow = p, ncol = ltau, dimnames = list(colnames(X), paste("zscore_qr", tau.list, sep = "_")))
    y = as.matrix(Y)
    
    if (!is.null(Z)) {
        zz = cbind(rep(1, nrow(y)), Z)
    } else {
        zz = matrix(1, nrow = nrow(y), ncol = 1)
    }

    ranks_list = lapply(tau.list, function(tau) {
        suppressWarnings(rq.fit.br(zz, y, tau = tau)$dual - (1 - tau))
    })
  
    for (ip in 1:p) {
        x = as.matrix(X[, ip])
        VN = matrix(0, nrow = ltau, ncol = ltau)
        for (i in 1:ltau) {
            for (j in 1:ltau) {
                VN[i, j] = min(tau.list[i], tau.list[j]) - tau.list[i] * tau.list[j]
            }
        }

        if (!is.null(Z)) {
            #xstar = lm(x ~ zz - 1)$residual
            xstar = .lm.fit(zz, x)$residual
        } else {
            xstar = x
        }

        SN = NULL
        for (itau in 1:ltau) {
            Sn = as.matrix(t(xstar) %*% ranks_list[[itau]])
            SN = c(SN, Sn)
        }
        VN2 = matrix(outer(VN, t(xstar) %*% xstar, "*"), nrow = ltau)
        z_score = SN / sqrt(diag(VN2))
        pvalue1 = pchisq(SN^2 / diag(VN2), 1, lower.tail = F)
        names(pvalue1) = tau.list
        quantile.pvalue[ip, ] <- pvalue1
        quantile.zscore[ip, ] <- z_score
        e = solve(chol(VN2))
        SN2 = t(e) %*% SN
        pvalue = pchisq(sum(SN2^2), ltau, lower.tail = F)
        pvec[ip] = pvalue
    }

    pvec <- apply(quantile.pvalue, 1, pval_cauchy)
    
    if (method == 'fdr') {
        adjusted_pvalues = p.adjust(pvec)
        method_col_name = "fdr_p_qr"
        method_quantile_names = paste0("fdr_p_qr_", tau.list)
        quantile_adjusted_pvalues = apply(quantile.pvalue, 2, p.adjust)
    } else if (method == 'qvalue') {
        adjusted_pvalues = compute_qvalues(pvec)
        method_col_name = "qvalue_qr"
        method_quantile_names = paste0("qvalue_qr_", tau.list)
        quantile_adjusted_pvalues = apply(quantile.pvalue, 2, compute_qvalues)
    } else {
        stop("Invalid method. Choose 'fdr' or 'qvalue'.")
    }

    sig_SNP_threshold = which(adjusted_pvalues < threshold)
    sig_SNP_top_count = order(adjusted_pvalues)[1:top_count]
    sig_SNP_top_percent = order(adjusted_pvalues)[1:max(1, round(length(adjusted_pvalues) * top_percent / 100))]

    sig.SNPs_names = colnames(X)[sig_SNP_threshold]
    sig.SNPs_names_top_count = colnames(X)[sig_SNP_top_count]
    sig.SNPs_names_top_percent = colnames(X)[sig_SNP_top_percent]
    phenotype_id = colnames(y)[1]

    df_result = data.frame(
        phenotype_id = phenotype_id,
        variant_id = colnames(X),
        p_qr = pvec
    )

    # Add quantile-specific p-values
    for (tau in tau.list) {
        df_result[[paste0("p_qr_", tau)]] = quantile.pvalue[, paste0("p_qr_", tau)]
    }

    # Add overall q-value
    df_result[[method_col_name]] = adjusted_pvalues

    # Add quantile-specific q-values
    for (tau in tau.list) {
        df_result[[paste0(method_col_name, "_", tau)]] = quantile_adjusted_pvalues[, paste0("p_qr_", tau)]
    }

    # Add quantile-specific z-scores
    for (tau in tau.list) {
        df_result[[paste0("zscore_qr_", tau)]] = quantile.zscore[, paste0("zscore_qr_", tau)]
    }

    # Split variant_id and reorder columns
    df_result <- df_result %>%
        separate(variant_id, into = c("chr", "pos", "ref", "alt"), sep = "[:_]", remove = FALSE) %>%
        mutate(chr = as.numeric(gsub("chr", "", chr)), pos = as.numeric(pos))

    # Define the column order
    col_order <- c("chr", "pos", "ref", "alt", "phenotype_id", "variant_id", "p_qr")
    col_order <- c(col_order, paste0("p_qr_", tau.list))
    col_order <- c(col_order, method_col_name)
    col_order <- c(col_order, paste0(method_col_name, "_", tau.list))
    col_order <- c(col_order, paste0("zscore_qr_", tau.list))

    # Reorder the columns
    df_result <- df_result %>% select(all_of(col_order))

    return(list(df_result = df_result,
                tau_list = tau.list, 
                quantile_pvalue = quantile.pvalue, 
                quantile_zscore = quantile.zscore,
                pvec = pvec, 
                adjusted_pvalues = adjusted_pvalues, 
                sig_SNP_threshold = sig_SNP_threshold, 
                sig.SNPs_names = sig.SNPs_names,
                sig_SNP_top_count = sig_SNP_top_count, 
                sig_SNP_top_percent = sig_SNP_top_percent))
}

#' Perform Clumping and Pruning
#' @param X Matrix of genotypes
#' @param qr_results Results from QR_screen
#' @param maf_list List of minor allele frequencies (optional)
#' @param ld_clump_r2 R-squared threshold for initial LD clumping based on pvalue
#' @param final_clump_r2 R-squared threshold for final LD clumping based on MAF
#' @return A list containing final SNPs and clumped SNPs
#' @importFrom bigstatsr FBM.code256 
#' @importFrom bigsnpr snp_clumping 
#' @export
multicontext_ld_clumping <- function(X, qr_results, maf_list = NULL, ld_clump_r2 = 0.2, final_clump_r2 = 0.8) {
    # Extract significant SNP names
    sig_SNPs_names <- qr_results$sig.SNPs_names
    
    # If no significant SNPs, return empty result
    if (length(sig_SNPs_names) == 0) {
        return(list(final_SNPs = NULL, message = "No significant SNPs found"))
    }
    
    # Check if X only contains one SNP (one column)
    if (ncol(X) == 1) {
        print("Only one SNP in X. Skipping LD clumping.")
        final_SNPs <- colnames(X)
        return(list(final_SNPs = final_SNPs, message = "Only one SNP, no LD clumping performed"))
    }
    
    # Extract genotype matrix for significant SNPs
    G_all <- X  # Only retain columns for significant SNPs
    
    # Convert the genotype matrix into FBM format
    code_vec <- c(0, 1, 2, rep(NA, 256 - 3))
    G_all <- FBM.code256(
        nrow = nrow(G_all),
        ncol = ncol(G_all),
        init = G_all,
        code = code_vec
    )
    
    # Parse SNP names to extract chromosome and position information
    parsed_snp_info <- do.call(rbind, strsplit(sig_SNPs_names, ":"))
    chr <- as.numeric(gsub("^chr", "", parsed_snp_info[, 1])) 
    pos <- as.numeric(parsed_snp_info[, 2])  # Extract position
    
    # Step 1: Perform LD clumping for each tau based on p-values
    clumped_snp_list <- list()
    
    for (itau in seq_along(qr_results$tau_list)) {
        tau_name <- paste0("p_qr_", qr_results$tau_list[itau])
        
        # Extract p-values for the given quantile
        p_values_quantile <- qr_results$quantile_pvalue[, tau_name][qr_results$sig_SNP_threshold]
        log_p_values <- -log10(p_values_quantile)  # Calculate log10 p-values      
        
        # Perform LD clumping using p-values as S
        ind_clumped <- snp_clumping(G = G_all, 
                                    infos.chr = chr, 
                                    infos.pos = pos, 
                                    S = log_p_values,  # Use log10 p-values for each quantile
                                    thr.r2 = ld_clump_r2, 
                                    size = 100 / ld_clump_r2)  # Window size of 500kb
        
        # Store clumping results for each quantile
        clumped_snp_list[[tau_name]] <- ind_clumped
    }
    
    # Step 2: Take the union of clumping results across all quantiles
    clumped_snp_union <- unique(unlist(clumped_snp_list))  # This is the SNP index, not the name
    print(paste("Number of SNPs after union of clumping:", length(clumped_snp_union)))
    
    # Step 3: Sort results from union
    sorted_indices <- order(chr[clumped_snp_union], pos[clumped_snp_union])
    chr_sorted <- chr[clumped_snp_union][sorted_indices]
    pos_sorted <- pos[clumped_snp_union][sorted_indices]
    
    # Step 4: Initialize maf_values to NULL, and update only if maf_list is provided
    maf_values <- NULL
    if (!is.null(maf_list)) {
        maf_values <- maf_list[qr_results$sig_SNP_threshold][clumped_snp_union][sorted_indices]
    }
    G_union = X[, clumped_snp_union, drop = FALSE][, sorted_indices]
    if (!inherits(G_union, "FBM")) {
        G_union <- FBM.code256(
        nrow = nrow(G_union),
        ncol = ncol(G_union),
        init = G_union,
        code = code_vec
    )}
    # Step 5: Perform final clumping using maf_values (if available), otherwise proceed with NULL S
    final_clumped <- snp_clumping(G = G_union, 
                                infos.chr = chr_sorted, 
                                infos.pos = pos_sorted, 
                                S = maf_values,  # Use MAF values if provided, otherwise NULL
                                thr.r2 = final_clump_r2, 
                                size = 100 / final_clump_r2)  # Final clumping
    
    # Restrict the genotype matrix to only the final clumped SNPs
    G_final_clumped <- G_union[, final_clumped, drop = FALSE]  # Limit G to the final clumped SNPs
    # Get the final SNP names
    final_SNPs <- sig_SNPs_names[clumped_snp_union][sorted_indices][final_clumped]
    print(paste("Number of final SNPs after MAF-based clumping (if applied):", length(final_SNPs)))
    
    return(list(final_SNPs = final_SNPs, clumped_SNPs = clumped_snp_union))
}

#' Perform Quantile Regression Analysis to get beta
#' @param X Matrix of predictors
#' @param Y Matrix or vector of response variables
#' @param Z Matrix of covariates (optional)
#' @param tau_values Vector of quantiles to be analyzed
#' @return A data frame with QR coefficients for each quantile
#' @importFrom quantreg rq rq.fit.br
#' @importFrom tidyr pivot_wider separate
#' @importFrom dplyr %>% mutate select
#' @export
perform_qr_analysis <- function(X, Y, Z = NULL, tau_values = seq(0.05, 0.95, by = 0.05)) {
# Convert Y and X to matrices if they aren't already
pheno.mat <- as.matrix(Y)
geno.mat <- as.matrix(X)

# Initialize an empty result table to store results
result_table <- data.frame(
    phenotype_id = character(),
    variant_id = character(),
    tau = numeric(),
    predictor_coef = numeric(),
    stringsAsFactors = FALSE
)

# Loop over each tau value to perform quantile regression
for (tau in tau_values) {
    # Loop over each SNP/variant in geno.mat (X)
    for (n in 1:ncol(geno.mat)) {
    response <- pheno.mat  # Y
    predictor <- geno.mat[, n]  # X
    phenotype_id <- colnames(pheno.mat)
    variant_id <- colnames(geno.mat)[n]
    
    # Construct the design matrix based on whether Z is provided
    if (is.null(Z)) {
        # If no covariates, include intercept and predictor
        X_design <- cbind(1, predictor)
    } else {
        # If covariates are provided, include them in the design matrix
        X_design <- cbind(1, predictor, as.matrix(Z))
    }
    
    # Fit the quantile regression model using rq.fit.br
    mod <- rq.fit.br(X_design, response, tau = tau)
    
    # Extract the coefficient for the predictor (second coefficient)
    predictor_coef <- mod$coefficients[2]  # Coefficient for predictor
    
    # Create a row with the results and append to the result table
    row <- data.frame(
        phenotype_id = phenotype_id, 
        variant_id = variant_id, 
        tau = tau, 
        predictor_coef = predictor_coef, 
        stringsAsFactors = FALSE
    )
    result_table <- rbind(result_table, row)
    }
}

# Reshape result_table to a wide format, so each tau's results are in separate columns
result_table_wide <- result_table %>%
    pivot_wider(
    id_cols = c(phenotype_id, variant_id),
    names_from = tau,
    values_from = predictor_coef,
    names_prefix = "coef_qr_"
    ) %>%
    separate(variant_id, into = c("chr", "pos", "ref", "alt"), sep = "[:_]", remove = FALSE, fill = "right") %>%
    mutate(
    chr = as.numeric(gsub("chr", "", chr)),
    pos = as.numeric(pos)
    ) %>%
    select(chr, pos, ref, alt, everything())

# Return the wide format result table
return(result_table_wide)
}

#' Filter Highly Correlated SNPs
#' @param X Matrix of genotypes
#' @param cor_thres Correlation threshold for filtering
#' @return A list containing filtered X matrix and filter IDs
#' @importFrom Rfast cora
corr_filter <- function(X, cor_thres = 0.8) {
    p <- ncol(X)
    
    # Use Rfast for faster correlation calculation if available
    if (requireNamespace("Rfast", quietly = TRUE)) {
        cor.X <- Rfast::cora(X, large = TRUE)
    } else {
        cor.X <- cor(X)
    }        
    Sigma.distance = as.dist(1 - abs(cor.X))
    fit = hclust(Sigma.distance, method="single")
    clusters = cutree(fit, h=1-cor_thres)
    groups = unique(clusters)
    ind.delete = NULL
    X.new = X
    filter.id = c(1:p)
    for (ig in 1:length(groups)) {
        temp.group = which(clusters == groups[ig])
        if (length(temp.group) > 1) {
            ind.delete = c(ind.delete, temp.group[-1])
        }
    }
    ind.delete = unique(ind.delete)
    if ((length(ind.delete) ) > 0) {
        X.new = as.matrix(X[,-ind.delete])
        filter.id = filter.id[-ind.delete]
    }
    return(list(X.new = X.new, filter.id = filter.id))
}


#' Calculate QR Coefficients and Pseudo R-squared
#' @param ExprData List containing X, Y, C, and X.filter
#' @param tau.list Vector of quantiles to be analyzed
#' @return A list containing beta matrix as twas weight and pseudo R-squared values
#' @importFrom quantreg rq rq.fit.br
calculate_qr_and_pseudo_R2 <- function(ExprData, tau.list) {
    # Fit models for all taus
    fit_full <- rq(Y ~ X.filter + C, tau = tau.list, data = ExprData)
    fit_intercept <- rq(Y ~ 1, tau = tau.list, data = ExprData)
    
    # Define rho function
    rho <- function(u, tau) {
    u * (tau - (u < 0))
    }
    
    # Prepare to store results
    pseudo_R2 <- numeric(length(tau.list))
    names(pseudo_R2) <- tau.list
    
    # Calculate pseudo R^2 for each tau
    for (i in seq_along(tau.list)) {
    tau <- tau.list[i]
    
    # Get residuals
    residuals0 <- residuals(fit_intercept, subset = i)
    residuals1 <- residuals(fit_full, subset = i)
    
    # Calculate and store pseudo R^2
    rho0 <- sum(rho(residuals0, tau))
    rho1 <- sum(rho(residuals1, tau))
    pseudo_R2[i] <- 1 - rho1 / rho0
    }
    
    # Extract coefficients of snps, removing intercept if included
    num_filter_vars <- ncol(ExprData$X.filter)
    beta_mat <- coef(fit_full)[2:(1 + num_filter_vars), , drop = FALSE]
    rownames_beta <- rownames(beta_mat)
    rownames(beta_mat) <- gsub("^X.filter", "", rownames_beta)
    list(beta_mat = beta_mat,  pseudo_R2 = pseudo_R2)
}

#' Quantile TWAS Weight Pipeline
#'
#' @param X Matrix of genotypes
#' @param Y Matrix or vector of phenotypes
#' @param Z Matrix of covariates (optional)
#' @param maf Vector of minor allele frequencies (optional)
#' @param extract_region_name Name of the region being analyzed
#' @param quantile_qtl_tau_list Vector of quantiles for QTL analysis
#' @param quantile_twas_tau_list Vector of quantiles for TWAS analysis
#'
#'
#' @return A list containing various results from the TWAS weight pipeline:
#' \itemize{
#'   \item qr_screen_pvalue_df: Data frame with QR screening results: pavlue, qvalue and zscore.
#'   \item message: Any informational or warning messages.
#'   \item twas_variant_names: Names of variants used in TWAS weight calculation.
#'   \item rq_coef_df: Data frame with quantile regression coefficients.
#'   \item twas_weight: Matrix of TWAS weights.
#'   \item pseudo_R2: Vector of pseudo R-squared values.
#'   \item quantile_twas_prediction: Matrix of TWAS predictions.
#' }
#'
#' @details
#' The function performs the following steps:
#' 1. QR screening to identify significant SNPs.
#' 2. Filtering of highly correlated SNPs.
#' 3. LD clumping and pruning(use filtered SNPs from step 1).
#' 4. Calculation of QR coefficients for selected SNPs(use filtered SNPs from step 3).
#' 5. Calculation of TWAS weights and pseudo R-squared values(use filtered SNPs from step 2).
#'
#' @examples
#' # Example usage:
#' # X <- matrix of genotypes
#' # Y <- vector of phenotypes
#' # Z <- matrix of covariates
#' # results <- quantile_twas_weight_pipeline(X, Y, Z, extract_region_name = "GeneA")
#'
#' @export
quantile_twas_weight_pipeline <- function(X, Y, Z = NULL, maf = NULL, extract_region_name, 
                                          quantile_qtl_tau_list = seq(0.05, 0.95, by = 0.05),
                                          quantile_twas_tau_list = seq(0.01, 0.99, by = 0.01)) {
    
    # Step 1: QR screen
    p.screen <- qr_screen(X = X, Y = Y, Z = Z, tau.list = quantile_qtl_tau_list, threshold = 0.05, method = "qvalue", top_count = 10,top_percent = 15)
    # Initialize results list
    results <- list(qr_screen_pvalue_df = p.screen$df_result)

    if (length(p.screen$sig_SNP_threshold) == 0) {
        results$message <- paste0("No significant SNPs detected in gene ", extract_region_name, names(fdat$residual_Y)[r])
        return(results)
    }
    # Step 2: Filter highly correlated SNPs
    if (length(p.screen$sig_SNP_threshold) > 1) {
        filtered <- corr_filter(X[, p.screen$sig_SNP_threshold, drop = FALSE], 0.8)
        X.filter <- filtered$X.new
    } else {
        X.filter <- X[, p.screen$sig_SNP_threshold, drop = FALSE]
        results$message <- paste0("Only one significant SNP in gene ", extract_region_name, names(fdat$residual_Y)[r], ", skipping correlation filter.")
    }
    # Step 3: LD clumping and pruning from results of QR_screen
    LD_SNPs <- multicontext_ld_clumping(X = X[, p.screen$sig_SNP_threshold, drop = FALSE], qr_results = p.screen, maf_list = maf)
    x_clumped <- X[, p.screen$sig_SNP_threshold, drop = FALSE][, LD_SNPs$final_SNPs, drop = FALSE]
    
    # Step 4: Only fit marginal QR to get beta with SNPs after LD pruning for quantile_qtl_tau_list values
    rq_coef_result <- perform_qr_analysis(X = x_clumped, Y = Y, Z = Z, tau_values = quantile_qtl_tau_list)
    
    # Step 5: Fit QR and get twas weight and R2 for all taus
    ExprData <- list(X = X, Y = Y, C = Z, X.filter = X.filter)
    qr_beta_R2_results <- calculate_qr_and_pseudo_R2(ExprData, quantile_twas_tau_list)
    
    # Add additional results
    results$twas_variant_names <- colnames(X.filter)
    results$rq_coef_df <- rq_coef_result
    results$twas_weight <- qr_beta_R2_results$beta_mat
    results$pseudo_R2 <- qr_beta_R2_results$pseudo_R2
    results$quantile_twas_prediction <- X.filter %*% results$twas_weight
    
    return(results)
}
