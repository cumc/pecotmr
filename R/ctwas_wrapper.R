# requires ctwas package.
# remotes::install_github("xinhe-lab/ctwas",ref = "multigroup_test")
#
#' Causal inference for TWAS with Summary Statistics with Multigroup Data
#'
#' @title cTWAS (causal-TWAS) wrapper
#
#

#' Determine Imputability and Variant Selection
#'
#' This function load TWAS weights and assesses the imputability of genes across different contexts, 
#' and select top variants for imputable gene-context pair, and output the extracted variant weights 
#' from best performing method model. Imputability of a gene-context pair are determined by having 
#' least one of the twas methods' model being imputable based on cross livadtion metrics. The imputable 
#' model is if the model surpasses the r-square and p-value threashold from cross validation metrics. 
#' If non of the contexts of a gene has at least one method being imputable, then this gene will be 
#' considered as unimputable, and do not return any weight results. For imputable gene-context pair, 
#' we select a subset of variants with high Posterior Inclusion Probabilities (PIPs) from SuSiE 
#' fine-mapping results for the imputable gene-context pair. After selecting variants, we extracts 
#' variant weights from the best performing model among imputable models with highest r-square or 
#' smallest p-value from cross validation metrics. This function is essential for preparing data 
#' for complex trait weighted analysis with sparse weight.
#' 
#' @param weight_db_files File paths to `.rds` files containing SuSiE-TWAS weights.
#' @param variable_name_obj Name of the element in SuSiE-TWAS output that lists variant names aligned with TWAS weights.
#' @param max_var_selection Maximum number of SNPs to be selected from each gene-context pair.
#' @param min_rsq_threshold Minimum \(R^2\) threshold for determining imputability of a gene-context pair.
#' @param p_val_cutoff Maximum allowable p-value for a gene to be considered imputable in a given context.
#' @return A list of elements containing fine-mapping and TWAS related data:
#' \itemize{
#'   \item{refined_twas_weights}{Context imputability, model selection, and selected variants, and extracted 
#' weights for the selected variants across all contexts for a imputable gene. }
#'   \item{susie_results}{SuSiE fine-mapping results, including variant selection from the SuSiE-TWAS pipeline.}
#'   \item{weights}{TWAS weights from the SuSiE-TWAS pipeline, to be harmonized for further analysis.}
#'   \item{gene}{Gene name of the twas db weight.}
#'   \item{cv_performance}{Cross validation performance metrics for all contexts.}
#' }
#' @export
#' @examples
#' results <- generate_twas_db(weight_db_file = "path/to/weights.rds", conditions=c("Mic", "Oli", "Exc"),
#'                                   max_var_selection = 10, min_rsq_threshold = 0.01, p_val_cutoff = 0.05)
# select variants for ctwas weights input
generate_twas_db <- function(weight_db_file, contexts = NULL, 
                                 variable_name_obj = c("preset_variants_result", "variant_names"),
                                 twas_weights_table = "twas_weights", max_var_selection,
                                 min_rsq_threshold = 0.01, p_val_cutoff = 0.05){
  
  # determine if the region is imputable and select the best model
  # Function to pick the best model based on adj_rsq and p-value
  pick_best_model <- function(twas_data_combined, contexts, min_rsq_threshold, p_val_cutoff) {
    best_adj_rsq <- min_rsq_threshold
    # Extract performance table
    performance_tables <- lapply(contexts, function(context) {
      get_nested_element(twas_data_combined, c("twas_cv_performance", context))
    })
    names(performance_tables) <- contexts
    # Determine if a gene/region is imputable and select the best model
    model_selection <- lapply(contexts, function(context) {
      selected_model <- NULL
      available_models <- do.call(c, lapply(names(performance_tables[[context]]), function(model) {
        if (!is.na(performance_tables[[context]][[model]][, "adj_rsq"])) {
          return(model)
        }
      }))
      if (length(available_models) <= 0) {
        print(paste0("No model provided TWAS cross validation performance metrics information at context ", context, ". "))
        return(NULL)
      }
      for (model in available_models) {
        model_data <- performance_tables[[context]][[model]]
        if (model_data[, "adj_rsq"] >= best_adj_rsq) {
          best_adj_rsq <- model_data[, "adj_rsq"]
          selected_model <- model
        }
      }
      region_names <- do.call(c, lapply(contexts, function(context) {
        twas_data_combined$susie_results[[context]]$region_info$region_name
      }))
      if (is.null(selected_model)) {
        print(paste0(
          "No model has p-value < ", p_val_cutoff, " and r2 >= ", min_rsq_threshold, ", skipping context ", context,
          " at region ", unique(region_names), ". "
        ))
        return(list(selected_model=c(NULL), imputable=FALSE)) # No significant model found
      } else {
        selected_model <- unlist(strsplit(selected_model, "_performance"))
        print(paste0("The selected best performing model for context ", context, " at region ", unique(region_names), " is ", selected_model, ". "))
        return(list(selected_model=selected_model, imputable=TRUE)) 
      }
    })
    names(model_selection) <- contexts
    return(model_selection)
  }
  
  # Based on selected best model and imputable contexts, select variants based on susie output
  twas_select <- function(twas_data_combined, contexts, imputable_status, max_var_selection, variable_name_obj) {
    # loop through context
    names(imputable_status) <- contexts
    var_selection_list <- lapply(contexts, function(context) {
      if (imputable_status[context] == "non_imputable") {
        variant_selected <- rep(NA, max_var_selection)
        return(variant_selected)
      } 
      
      # if the context is imputable 
      variant_names = get_nested_element(twas_data_combined, c("susie_results", context,"variant_names"))
      # susie_result_trimmed is from "preset_variants_result"-"susie_result_trimmed" from load_twas_weights()
      susie_result_trimmed = get_nested_element(twas_data_combined, c("susie_results", context, "susie_result"))

      # Go through cs to select top variants from each set until we select all variants
      select_cs_var <- function(set_cs_list, max_var_selection) {
        selected_indices <- c()
        while (length(selected_indices) < max_var_selection) {
          for (cs in set_cs_list$sets$cs) {
            if (length(selected_indices) < max_var_selection) {
              # Filter out any indices that have already been selected
              available_indices <- setdiff(cs, selected_indices)
              if (length(available_indices) > 0) {
                # Get the index with the highest PIP score
                top_index <- available_indices[which.max(set_cs_list$pip[available_indices])]
                selected_indices <- c(selected_indices, top_index)
              }
            }
          }
        }
        return(selected_indices)
      }
      # check if the context has top_loci table
      if ("top_loci" %in% names(twas_data_combined$susie_results[[context]])) {
        top_loci <- twas_data_combined$susie_results[[context]][["top_loci"]]
        if (length(top_loci[, "variant_id"]) <= max_var_selection) {
          variant_selected <- top_loci$variant_id
        } else {
          # if top_loci variants more than max_var_selection, we loop through CS set to select top pip variants from each CS set

          if (!is.null(susie_result_trimmed$sets$cs)) {
            # check total variant number in the CS sets
            totl_cs_indice <- do.call(sum, lapply(names(susie_result_trimmed$sets$cs), function(L) {
              length(susie_result_trimmed$sets$cs[[L]])
            }))
            if (max_var_selection <= totl_cs_indice) {
              selec_idx <- select_cs_var(susie_result_trimmed, max_var_selection)
              variant_selected <- variant_names[selec_idx]
            } else {
              # when top loci has larger number of variants than max_var_selection, but CS sets has less number than max_var_selection
              selec_idx <- unlist(susie_result_trimmed$sets$cs)
              variant_selected <- variant_names[selec_idx]
            }
          } else {
            top_idx <- order(susie_result_trimmed$pip, decreasing = TRUE)[1:max_var_selection]
            variant_selected <- variant_names[top_idx]
          }
        }
      } else {
        # if the context did not come with the top_loci table, we select top pip variants from [[context]]$preset_variants_result$susie_result_trimmed$pip
        print(paste0(context, " do not have top_loci table, select top pip variants. "))
        top_idx <- order(susie_result_trimmed$pip, decreasing = TRUE)[1:max_var_selection]
        variant_selected <- variant_names[top_idx]
      } 
      return(variant_selected)
    })
    names(var_selection_list) <- contexts
    return(var_selection_list)
  }
  
  
  ## load twas weights
  twas_data_combined <- load_twas_weights(weight_db_file, conditions=contexts, variable_name_obj=variable_name_obj,
                                                        twas_weights_table = twas_weights_table)
  
  weights <- twas_data_combined$weights
  if (is.null(contexts)) {contexts <- names(weights)}
  
  #get gene name for the weight file 
  gene <- unique(twas_data_combined$susie_results[[1]]$region_info$region_name)
  
  ## we first select best model, determine imputable contexts, then select variants based on susie obj output
  model_selection <- pick_best_model(twas_data_combined, contexts, min_rsq_threshold, p_val_cutoff) 
  # Determine Imputable Contexts  
  imputable_contexts <- c()
  for (contx in contexts){
    if (isTRUE(model_selection[[contx]][["imputable"]])) imputable_contexts <- c(imputable_contexts, contx)
  }
  
  # select variants 
  names(contexts) <- rep("non_imputable", length(contexts))
  names(contexts)[which(contexts %in% imputable_contexts)] <- "imputable"
  twas_select_result <- twas_select(twas_data_combined, contexts, names(contexts), max_var_selection, variable_name_obj)
  
  # output refined_twas_weights for imputable genes across all contexts 
  refined_twas_weights <- list()
  if (length(imputable_contexts)>0){
    for (context in contexts){
      #susie_result_name update
      names(twas_data_combined$susie_results[[context]])[which(names(twas_data_combined$susie_results[[context]]) %in% "susie_result")] <- "susie_result_trimmed"
      # construct refined_twas_weights
      refined_twas_weights[[context]] <- list(selected_model=model_selection[[context]][["selected_model"]])
      refined_twas_weights[[context]][["is_imputable"]] <- model_selection[[context]][["imputable"]]

      if (isTRUE(refined_twas_weights[[context]][["is_imputable"]])){
          all_weight_variants <- rownames(weights[[context]])
          refined_twas_weights[[context]][["selected_top_variants"]] <- twas_select_result[[context]]
          refined_twas_weights[[context]][["selected_model_weights"]] <- weights[[context]][match(twas_select_result[[context]],                         
                            all_weight_variants), paste0(refined_twas_weights[[context]][["selected_model"]], "_weights")]

          if (model_selection[[context]][["selected_model"]]=="susie"){
              refined_twas_weights[[context]][["susie_parameters"]] <- get_nested_element(twas_data_combined, c("susie_results", context, 
                                                        "susie_result_trimmed"))[c("X_column_scale_factors","lbf_variable", "mu")]
          }
      }
    }
    # return results
    return(list(
      refined_twas_weights=refined_twas_weights, 
      weights = weights, gene=gene, cv_performance=twas_data_combined$twas_cv_performance, 
      susie_results = twas_data_combined$susie_results))
  } else {
      print(paste0("Weight input ", weight_db_file, " is non-imputable for all contexts. "))
      return(NULL)
  }
}



#' Function to perform allele flip QC and harmonization on the weights and GWAS against LD for a region. 
#' FIXME: GWAS loading function from Haochen for both tabix & column-mapping yml application 
#' Main Steps: 
#' 1. allele QC for TWAS weights against the LD meta
#' 2. allele QC for GWA Ssummary stats against the LD meta
#' 3. adjust susie/mvsusie weights based on the overlap variants
#' @param refined_twas_weights_data List of list of twas weights output from twas pipeline with following list items. 
#' @param gwas_file A file path for a dataframe of GWAS summary statistics with column name "chrom" (or #chrom" if tabix-indexed),
#' "pos", "A2", "A1".
#' @param ld_meta_file_path A tab-delimited data frame with colname "#chrom", "start", "end", "path", where "path" column
#' contains file paths for both LD matrix and bim file and is separated by ",". Bim file input would expect no headers, while the
#' columns are aligned in the order of "chrom", "variants", "GD", "pos", "A1", "A2", "variance", "allele_freq", "n_nomiss". 
#' @param scale_weights TRUE/FALSE statement. If turn off, the post-qc/harmonized weights will not be scaled by SNP variance from LD. 
#' @return A list of list for harmonized weights and dataframe of gwas summary statistics that is add to the original input of 
#' refined_twas_weights_data under each context.
#' @import IRanges IRanges findOverlaps queryHits subjectHits start end
#' @import data.table fread
#' @export
update_weights_gwas <- function(refined_twas_weights_data, gwas_file, ld_meta_file_path, scale_weights=FALSE){
  # Get region info
  region <- unique(names(refined_twas_weights_data))
  genes <- names(refined_twas_weights_data[[region]])
  chrom <- as.integer(readr::parse_number(gsub( "_.*$", "", region)))
  
  # Function to group contexts based on start and end positions
  group_contexts_by_region <- function(refined_twas_weights_data, gene, chrom, tolerance = 10000) {
      contexts <- do.call(c, lapply(names(refined_twas_weights_data[[1]][[gene]]), function(context) {
                          if (isTRUE(refined_twas_weights_data[[1]][[gene]][[context]][["is_imputable"]])) context}))
      region_info_df <- do.call(rbind, lapply(contexts, function(context) {
                                    wgt_range <- sapply(get_nested_element(refined_twas_weights_data[[1]], c(gene, context, 
                                              "selected_top_variants")),function(variant_id) {strsplit(variant_id, "\\:")[[1]][2]}) 
                                    min <- min(wgt_range)
                                    max <- max(wgt_range)
                                    data.frame(context = context, start = min, end = max)}))
      # Calculate distance matrix and perform hierarchical clustering
      clusters <- cutree(hclust(dist(region_info_df[, c("start", "end")])), h = tolerance)
      # Group contexts and determine query regions
      region_groups <- split(region_info_df, clusters) %>%
          lapply(function(group){list(contexts = group$context, query_region = paste0(chrom, ":", min(group$start), 
                                 "-", max(group$end)))})
      # Create IRanges objects and merge overlapping intervals
      intervals <- IRanges(start = unlist(lapply(region_groups, function(context_group) as.numeric(gsub('^.*:\\s*|\\s*-.*$', '',
                       context_group$query_region)))), end = unlist(lapply(region_groups, 
                       function(context_group) as.numeric(sub("^.*?\\-", "", context_group$query_region)))))
      reduced_intervals <- reduce(intervals)

      # Find which original groups are merged, and update region_groups lists
      overlaps <- findOverlaps(intervals, reduced_intervals)
      # Create merged groups based on overlap mapping
      merged_groups <- lapply(seq_along(reduced_intervals), function(i) {
          context_indices <- queryHits(overlaps)[subjectHits(overlaps) == i]
          merged_contexts <- unlist(lapply(context_indices, function(idx) region_groups[[idx]]$contexts))
          list(contexts = merged_contexts, query_region = paste0(chrom, ":", start(reduced_intervals[i]), "-", end(reduced_intervals[i])))
      })
      names(merged_groups) <- paste0("context_group_", seq_along(merged_groups))
      return(merged_groups)
  }
  # Function to extract LD variance for the query region                                 
  query_variance <- function(ld_meta_file_path, region, query_region) {
    # Extract chromosome, start, and end from the region string
    region_parts <- strsplit(region, "_")[[1]]
    region_chr <- region_parts[1]
    region_start <- as.numeric(region_parts[2]) #this start and end correspond to the LD block meta file to locate the ld block
    region_end <- as.numeric(region_parts[3])

    # Find the row in the LD meta file that matches the region
    ld_meta <- fread(ld_meta_file_path, header = TRUE, sep = "\t", data.table = FALSE)
    ld_row <- ld_meta[ld_meta[,1] == region_chr & ld_meta$start <= region_start & ld_meta$end >= region_end,]
    bim_file_path <- paste0(dirname(ld_meta_file_path), "/", strsplit(ld_row$path, ",")[[1]][2])

    # Filter the BIM data for the specified query region
    query_parts <- strsplit(query_region, "[:-]")[[1]]
    query_start <- as.numeric(query_parts[2])
    query_end <- as.numeric(query_parts[3])
    bim_data <- fread(bim_file_path, header = FALSE, data.table = FALSE)
    bim_filtered <- bim_data[bim_data$V1 == as.integer(query_parts[1]) & bim_data$V4 >= query_start & bim_data$V4 <= query_end, ]

    # Extract the variance column (7th column)
    variance_df <- bim_filtered[, c(1,4,5:7)] 
    colnames(variance_df) <- c("chrom", "pos", "A1", "A2", "variance")
    return(variance_df)
  }                                                         
                                  
  for (gene in genes) {
    context_clusters <- group_contexts_by_region(refined_twas_weights_data, gene, chrom, tolerance = 10000)
    for (context_group in names(context_clusters)){
      contexts <- context_clusters[[context_group]]$contexts
      query_region <- context_clusters[[context_group]]$query_region
      region_of_interest <- data.frame(chrom = chrom, start = as.numeric(gsub('^.*:\\s*|\\s*-.*$', '', query_region)), 
                                       end = as.numeric(sub("^.*?\\-", "", query_region)))
      # ??? do we enable multiple gwas study loading & harmonization per function call? # for (s in seq_along(gwas_studies)){}
      gwas_data <- list()
      gwas_sumstats <- tabix_region(gwas_file, query_region) #extension for yml file for column name mapping
      if (colnames(gwas_sumstats)[1]=="#chrom") colnames(gwas_sumstats)[1] <- "chrom"
      # Load LD list containing LD matrix and corresponding variants
      gwas_LD_list <- load_LD_matrix(ld_meta_file_path, region_of_interest, gwas_sumstats)
      # remove duplicate variants
      dup_idx <- which(duplicated(gwas_LD_list$combined_LD_variants))
      if (length(dup_idx) >= 1) {
        gwas_LD_list$combined_LD_variants <- gwas_LD_list$combined_LD_variants[-dup_idx] 
        gwas_LD_list$combined_LD_matrix <- gwas_LD_list$combined_LD_matrix[-dup_idx, -dup_idx] 
        gwas_LD_list$ref_panel <- gwas_LD_list$ref_panel[-dup_idx,]
      }
      # Allele flip
      gwas_allele_flip <- allele_qc(gwas_sumstats[, c("chrom", "pos", "A2", "A1")], gwas_LD_list$combined_LD_variants, 
                                gwas_sumstats, c("beta", "se", "z"))
      # Load LD matrix and sumstats
      gwas_data[["LD"]] <- gwas_LD_list$combined_LD_matrix
      gwas_data[["sumstats"]] <- gwas_allele_flip$target_data_qced
      for (context in contexts){
        # Intersect with gwas summary statistics and adjust susie weights
        if (refined_twas_weights_data[[region]][[gene]][[context]]$selected_model=="susie"){
            adjusted_susie_weights <- adjust_susie_weights(refined_twas_weights_data, condition=context, 
                                              keep_variants = gwas_data[["sumstats"]]$variant_id, allele_qc = TRUE, 
                                              variable_name_obj = c(region, gene, context,"selected_top_variants"), 
                                              susie_obj = c(region, gene, context, "susie_parameters"),
                                              twas_weights_table = c(region, gene, context, "selected_model_weights"))
            refined_twas_weights_data[[region]][[gene]][[context]]$selected_model_weights <- adjusted_susie_weights$adjusted_susie_weights
        }
        weights <- matrix(get_nested_element(refined_twas_weights_data, c(region, gene, context, "selected_model_weights")), ncol=1, 
                          dimnames = list(refined_twas_weights_data[[region]][[gene]][[context]][["selected_top_variants"]], 
                          paste0(refined_twas_weights_data[[region]][[gene]][[context]]$selected_model,"_weights")))
        # Allele flip for weight variants and harmonize weights    
        weights_qced <- allele_qc(get_nested_element(refined_twas_weights_data, c(region, gene, context,"selected_top_variants")), 
                                  gwas_LD_list$combined_LD_variants, weights, 1,  target_gwas = FALSE)
        weights_subset <- weights_qced$target_data_qced[, !colnames(weights_qced$target_data_qced) %in% c("chrom","pos", "A2", "A1", "variant_id")]
        names(weights_subset) <- get_nested_element(weights_qced, c("target_data_qced", "variant_id"))
        refined_twas_weights_data[[region]][[gene]][[context]]$qced_weights <- weights_subset
        # scale weights by variance 
        if (isTRUE(scale_weights)){
          variance_df <- query_variance(ld_meta_file_path, region, query_region) %>%
                        mutate(key1 = paste(chrom, pos, A1, A2, sep = ":"), key2 = paste(chrom, pos, A2, A1, sep = ":"))
          match_idx1 <- match(names(weights_subset), variance_df$key1)
          match_idx2 <- match(names(weights_subset), variance_df$key2)
          var_idx <- ifelse(!is.na(match_idx1), match_idx1, match_idx2)# Combine indices and choose non-NA matches by order
          refined_twas_weights_data[[region]][[gene]][[context]]$qced_weights <- weights_subset * sqrt(variance_df$variance[var_idx])
        }
      }
    }
  }
  # return results
  return(list(refined_twas_weights_qced=refined_twas_weights_data, gwas_sumstats_qced = gwas_data[["sumstats"]]))
}




#' Function to generate weights_list from harmonized refined_twas_weights_data for cTWAS multigroup analysis
#' - need to qc for allele-flip for weights against LD 
#' - need to update susie_weights with susie_parameters 
#' - need to get variance to scale weights 
#' @param refined_twas_weights_data output from twas pipeline as refined_twas_weights_data that is QC/harmonizaed with LD. 
#' If the selected_model is susie, the input of weights are updated  
#' @param     
#' @param  
#' @return A list of list for weight information for each gene-context pair.  
#' @export
format_ctwas_weights <- function(refined_twas_weights_data, xqtl_meta, ld_meta_file) {

  if(x==0){
  weights_list <- lapply(genes, function(gene) {
    contexts <- names(refined_twas_weights_data[[region]][[gene]])
    for (context in contexts){
      
    }
    variants <- paste0("chr", colnames(twas_rs[[study]]$gwas_qced$LD))[var_idx] # get final variants name
    # get variant's variance to scale weights
    variance_idx <- match(variants, paste0("chr", twas_rs[[study]]$gwas_qced$sumstats$variant_id))
    variance_snp <- twas_rs[[study]]$gwas_qced$variance[variance_idx]
    # scale weights
    wgt <- wgt * sqrt(variance_snp) # wgt <- wgt*sqrt(ld_snpinfo_wgt$variance[ld_snpinfo_wgt.idx])
    wgt <- as.matrix(wgt[variants])
    colnames(wgt) <- "weight"
    # get LD - correlation matrix for selected variance
    ld_wgt <- as.matrix(twas_rs[[study]]$gwas_qced$LD[var_idx, var_idx])
    colnames(ld_wgt) <- rownames(ld_wgt) <- variants
    # p0 and p1 represent the min max position of weights variants
    p0 <- min(sapply(variants, function(x) {
      as.integer(unlist(strsplit(x, ":"))[2])
    }))
    p1 <- max(sapply(variants, function(x) {
      as.integer(unlist(strsplit(x, ":"))[2])
    }))
    return(list(
      chrom = chrom,
      p0 = p0,
      p1 = p1,
      wgt = wgt,
      R_wgt = ld_wgt,
      gene_name = gene,
      weight_name = gsub("_.*$", "", condition),
      type = type,
      context = gsub("_.*$", "", condition),
      n_wgt = length(var_idx)
    ))
  })
    names(weights_list) <- paste0(gene, "|", gsub("_.*$", "", imputable_conditions))
    return(weights_list)
  } else {
    print(paste0("Gene ", gene, " is not imputable, skipping. "))
  }
}


#' Convert Lower Triangle LD Matrix to Symmetrical Full Matrix fot cTWAS analysis. 
#'
#' This function converts a lower triangle LD (Linkage Disequilibrium) matrix into a symmetrical full matrix suitable for cTWAS analysis. It then saves the result in both RDS and Rvar formats in the specified output directory. This is particularly useful for genetic analysis workflows that require full matrix formats for computational processes.
#'
#' @param ld_paths A vector of full pathnames to LD matrix files, expected to follow the naming convention `/LD_path/chrN/chrN_start_end.cor.xz` with a corresponding bim file in the same directory, e.g., `/LD_path/chrN/chrN_start_end.cor.xz.bim`.
#' @param outdir A character string specifying the output directory where the processed files will be saved.
#' @return A data frame containing region information for each formatted LD file, detailing the locations and specifics of the regions processed.
#' @details On execution, the function writes two types of files to the `outdir`:
#' \itemize{
#'   \item An `.RDS` file, which contains the fully processed symmetrical LD matrix.
#'   \item A `.Rvar` file, which includes variance estimates for the LD regions processed.
#' }
#' @importFrom data.table fwrite
#' @export
#' @examples
#' convert_ld_matrix(ld_paths = c("/path/to/ld/chr1/chr1_12345_17755.cor.xz"),
#'                   outdir = "path/to/output")

format_ctwas_ld <- function(ld_paths, outdir) {
  region_info <- data.frame()
  for (ld_path in ld_paths) {
    if (!file.exists(paste0(outdir, "/LD_", basename(ld_path), ".RDS"))){
      ld <- as.matrix(read.table(ld_path, sep = " "))
      ld[lower.tri(ld)] <- t(ld)[lower.tri(ld)]
      print(paste0("writing ld file ", outdir, "/LD_", basename(ld_path), ".RDS"))
      saveRDS(ld, paste0(outdir, "/LD_", basename(ld_path), ".RDS"), compress='xz')
    }
    if (!file.exists(paste0(outdir, "/LD_", basename(ld_path), ".Rvar"))){
      bim <- read.table(paste0(ld_path, ".bim"), header = FALSE)[, -c(3, 9)] # remove posg and number missing sample column
      colnames(bim) <- c("chrom", "id", "pos", "alt", "ref", "variance", "allele_freq")
      bim$id <- gsub("_", ":", bim$id)
      print(paste0("writing ld file ", outdir, "/LD_", basename(ld_path), ".Rvar"))
      data.table::fwrite(bim, file = paste0(outdir, "/LD_", basename(ld_path), ".Rvar"), sep = "\t", quote = F)
    }
    ld_info <- data.frame(
      chrom = as.integer(sub(".*chr", "", gsub("_.*$", "", basename(ld_path)))),
      start = as.integer(unlist(strsplit(basename(ld_path), "_"))[2]),
      stop = as.integer(gsub("\\..*", "", gsub(".*\\_", "", basename(ld_path)))),
      region_tag = paste0(gsub("_.*$", "", basename(ld_path)), ":", unlist(strsplit(basename(ld_path), "_"))[2], "-", gsub("\\..*", "", gsub(".*\\_", "", basename(ld_path)))),
      LD_matrix = paste0(outdir, "/LD_", basename(ld_path), ".RDS"),
      SNP_info = paste0(outdir, "/LD_", basename(ld_path), ".Rvar")
    )
    region_info <- rbind(region_info, ld_info)
  }
  return(region_info)
}
