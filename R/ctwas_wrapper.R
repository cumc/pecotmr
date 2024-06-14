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
                                 susie_obj = c("preset_variants_result", "susie_result_trimmed"),
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
  twas_select <- function(twas_data_combined, contexts, imputable_status, max_var_selection, variable_name_obj, susie_obj) {
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
  twas_select_result <- twas_select(twas_data_combined, contexts, names(contexts), max_var_selection, variable_name_obj, susie_obj)
  
  # output refined_twas_weights for imputable genes across all contexts 
  refined_twas_weights <- list()
  if (length(imputable_contexts)>0){
    for (context in contexts){
      #susie_result_name update
      names(twas_data_combined$susie_results[[context]])[which(names(twas_data_combined$susie_results[[context]]) %in% "susie_result")] <- "susie_result_trimmed"

      refined_twas_weights[[context]] <- list(selected_model=model_selection[[context]][["selected_model"]])
      refined_twas_weights[[context]][["is_imputable"]] <- model_selection[[context]][["imputable"]]

      if (isTRUE(refined_twas_weights[[context]][["is_imputable"]])){
          all_weight_variants <- rownames(weights[[context]])
          refined_twas_weights[[context]][["selected_top_variants"]] <- twas_select_result[[context]]
          refined_twas_weights[[context]][["selected_model_weights"]] <- weights[[context]][match(twas_select_result[[context]],                                
                            all_weight_variants), paste0(refined_twas_weights[[context]][["selected_model"]], "_weights")]
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




#' Format Multigroup cTWAS Input Lists
#' Aggregates gene data from the 'summary_report' table and organizes it into nested lists. These lists are hierarchically structured by `gwas_study_name`, chromosome (chr`chr_number`), and then grouped under "weights_list" with each gene identified by `gene_name` and its `context` (`gene_name`|`context`).  
#' @param summary_report A data frame containing metadata about each gene, including imputable contexts and the best-performing TWAS methods.
#' @param outdir The directory where the output will be saved. If not specified, the current working directory is used.
#' @param weights_input_file Path to an existing weight list file. This file is used to append additional weights based on new entries in the `summary_report`.
#' @param auto_save Boolean flag indicating whether to automatically save the output to the `outdir` or use the current working directory if `outdir` is not specified. Default is `FALSE`.
#' @importFrom foreach foreach %do%
#' @export
# weights_input_file is an existing weight file that we want to add additional weights on top of it.  
get_ctwas_input <- function(summary_report, outdir=NULL, outname=NULL, weights_input_file = NULL, auto_save = TRUE) {
  if (!is.null(weights_input_file)) {
    out <- readRDS(weights_input_file)
  } else {
    out <- NULL
  }
  # only process new genes and add to pre-existing ctwas input file
  summary_report <- summary_report[summary_report$IsImputable, ]
  genes <- summary_report$gene
  if (!is.null(weights_input_file)) {
    processed_genes <- unique(do.call(c, lapply(names(out[[1]]), function(chr) {
      gsub("\\|.*$", "", names(out[[1]][[chr]]$weights))
    })))
    genes <- genes[!genes %in% processed_genes]
  }
  file_list <- na.omit(summary_report$path[summary_report$gene %in% genes])
  if (length(file_list) >= 1) {
    # loop through genes in a chromosome
    foreach(file = file_list) %do% {
      twas_rs <- readRDS(file) # file is the ctwas variant selection output
      gene <- summary_report$gene[summary_report$path == file]
      chr <- as.integer(summary_report$chrom[summary_report$path == file])
      data_type <- summary_report$type[summary_report$path == file]
      studies <- na.omit(names(twas_rs))

      for (study in studies) {
        # update colnames of gwas sumstat z_snp
        gwas <- twas_rs[[study]]$gwas_qced
        id_idx <- match("variant_id", colnames(gwas$sumstats))
        colnames(gwas$sumstats)[id_idx] <- "id"
        # merge gwas sumstats and weights info
        out[[study]][[paste0("chr", chr)]][["z_snp"]] <- rbind(out[[study]][[paste0("chr", chr)]][["z_snp"]], gwas$sumstats[, c("id", "A1", "A2", "z")])
        out[[study]][[paste0("chr", chr)]][["z_snp"]] <- out[[study]][[paste0("chr", chr)]][["z_snp"]][!duplicated(out[[study]][[paste0("chr", chr)]][["z_snp"]][, "id"]), ]
        out[[study]][[paste0("chr", chr)]][["weights_list"]] <- c(out[[study]][[paste0("chr", chr)]][["weights_list"]], format_ctwas_weights(twas_rs, study, data_type))
      }
    }
    if (!is.null(weights_input_file)) {
      if (auto_save){
        if (is.null(outdir)) outdir <- dirname(weights_input_file)
        print(paste0("Updating formated weights list of ", outdir, "/", basename(weights_input_file)))
        saveRDS(out, paste0(outdir, "/", basename(weights_input_file)), compress='xz')
      }
      return(out)
    } else {
      if (auto_save){
        if (is.null(outdir)) outdir <- getwd()
        print(paste0("saving formated output as ", outdir, "/ctwas_weights_", outname, ".rds. "))
        saveRDS(out, paste0(outdir, "/ctwas_weights_", outname, ".rds"), compress='xz')
      }
      return(out)
    }
  } else {
    print(paste0("No ctwas input file found. "))
    return(NULL)
  }
}


#' Function to extract selected weights from a gene across all contexts for cTWAS multigroup analysis
#' @param twas_rs twas_sparse pipeline selection results.
#' @param study Study name of the GWAS data that was incorporated for TWAS analysis.   
#' @param type xQTL data type such as expression or protein
#' @return A list of list for weight information for each gene-context pair.  
#' @export
format_ctwas_weights <- function(twas_rs, study, type) {
  # Register parallel backend to use multiple cores
  region_info <- twas_rs[[study]][["region_info"]]
  chrom <- as.integer(region_info$grange$chrom)
  gene <- unique(region_info$region_name)
  # check if this gene is imputable at study
  if (twas_rs[[study]]$model_selection[[1]]$imputable) {
    conditions <- names(twas_rs[[study]]$weights)
    imputable_conditions <- unlist(lapply(conditions, function(condition) {
      ifelse(is.null(twas_rs[[study]][["model_selection"]][[condition]][["method"]]), return(NULL), return(condition))
    }))
    weights_list <- lapply(imputable_conditions, function(condition) {
      # remove NA weights
      wgt <- twas_rs[[study]][["weights"]][[condition]][["selected_weights"]]
      wgt_idx <- !is.na(wgt)
      wgt <- wgt[wgt_idx]
      if (length(wgt) == 0) {
        return(NULL)
      }
      variants <- twas_rs[[study]][["weights"]][[condition]][["variant_selected"]][wgt_idx]
      var_idx <- na.omit(match(variants, paste0("chr", colnames(twas_rs[[study]]$gwas_qced$LD))))
      # return NULL if no variants found in common of weights and LD
      if (length(var_idx) == 0) {
        return(NULL)
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
