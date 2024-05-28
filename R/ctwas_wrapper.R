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
#' Assesses the imputability of genes across different contexts and selects a subset of SNPs with high Posterior Inclusion Probabilities (PIPs) from SuSiE fine-mapping results. It retrieves variant weights and identifies the best-performing methods based on cross-validation metrics. This function is essential for preparing data for complex trait weighted analysis.
#' @param weight_db_files File paths to `.rds` files containing SuSiE-TWAS weights.
#' @param variable_name_obj Name of the element in SuSiE-TWAS output that lists variant names aligned with TWAS weights.
#' @param max_var_selection Maximum number of SNPs to be selected from each gene-context pair.
#' @param min_rsq_threshold Minimum \(R^2\) threshold for determining imputability of a gene-context pair.
#' @param p_val_cutoff Maximum allowable p-value for a gene to be considered imputable in a given context.
#' @return A list of elements containing fine-mapping and TWAS related data:
#' \itemize{
#'   \item{susie_results}{SuSiE fine-mapping results, including variant selection from the SuSiE-TWAS pipeline.}
#'   \item{model_selection}{Identification of the best-performing TWAS method for each context, based on cross-validation.}
#'   \item{weights}{TWAS weights from the SuSiE-TWAS pipeline, to be harmonized for further analysis.}
#'   \item{region_info}{Region information corresponding to each gene, critical for interpreting SuSiE-TWAS weight inputs.}
#' }
#' @export
#' @examples
#' results <- determine_imputability(weight_db_files = "path/to/weights.rds", conditions=c("Mic", "Oli", "Exc"),
#'                                   max_var_selection = 10, min_rsq_threshold = 0.01, p_val_cutoff = 0.05)
# select variants for ctwas weights input
select_ctwas_weights <- function(weight_db_files, conditions = NULL,
                                 variable_name_obj = c("preset_variants_result", "variant_names"),
                                 susie_obj = c("preset_variants_result", "susie_result_trimmed"),
                                 twas_weights_table = "twas_weights", max_var_selection,
                                 min_rsq_threshold = 0.01, p_val_cutoff = 0.05) {
  ## Internal function to load and validate data from RDS files
  load_and_validate_data <- function(weight_db_files, conditions, variable_name_obj) {
    all_data <- lapply(weight_db_files, readRDS)
    unique_regions <- unique(unlist(lapply(all_data, function(data) names(data))))

    # Check if region from all RDS files are the same
    if (length(unique_regions) != 1) {
      stop("The RDS files do not refer to the same region.")
    } else {
      # Assuming all data refer to the same region, now combine data by conditions
      combined_all_data <- do.call("c", lapply(all_data, function(data) data[[1]]))
    }
    # Set default for 'conditions' if they are not specified
    if (is.null(conditions)) {
      conditions <- names(combined_all_data) # get only condition name
    }
    ## Check if the specified condition and variable_name_obj are available in all files
    if (!all(conditions %in% names(combined_all_data))) {
      stop("The specified condition is not available in all RDS files.")
    }
    return(combined_all_data)
  }
  # determine if the region is imputable and select the best model
  # Function to pick the best model based on adj_rsq and p-value
  pick_best_model <- function(combined_all_data, conditions, min_rsq_threshold, p_val_cutoff) {
    best_adj_rsq <- min_rsq_threshold
    # Extract performance table
    performance_tables <- lapply(conditions, function(condition) {
      get_nested_element(combined_all_data, c(condition, "twas_cv_result", "performance"))
    })
    names(performance_tables) <- conditions
    # Determine if a gene/region is imputable and select the best model
    model_selection <- lapply(conditions, function(condition) {
      best_model <- NULL
      available_models <- do.call(c, lapply(names(performance_tables[[condition]]), function(model) {
        if (!is.na(performance_tables[[condition]][[model]][, "adj_rsq"])) {
          return(model)
        }
      }))
      if (length(available_models) <= 0) {
        print(paste0("No model provided TWAS cross validation performance metrics information at condition ", condition, ". "))
        return(NULL)
      }
      for (model in available_models) {
        model_data <- performance_tables[[condition]][[model]]
        if (model_data[, "adj_rsq"] >= best_adj_rsq) {
          best_adj_rsq <- model_data[, "adj_rsq"]
          best_model <- model
        }
      }
      region_names <- do.call(c, lapply(conditions, function(condition) {
        combined_all_data[[condition]]$region_info$region_name
      }))
      if (is.null(best_model)) {
        print(paste0(
          "No model has p-value < ", p_val_cutoff, " and r2 >= ", min_rsq_threshold, ", skipping condition ", condition,
          " at genes/regions ", unique(region_names), ". "
        ))
        return(NULL) # No significant model found
      } else {
        best_model <- unlist(strsplit(best_model, "_performance"))
        print(paste0("The best model for condition ", condition, " at region ", unique(region_names), " is ", best_model, ". "))
        return(best_model)
      }
    })
    names(model_selection) <- conditions
    return(model_selection)
  }
  # Based on selected best model and imputable conditions, select variants based on susie output
  ctwas_select <- function(combined_all_data, conditions, imputable_status, max_var_selection, variable_name_obj, susie_obj) {
    names(imputable_status) <- conditions
    var_selection_list <- lapply(conditions, function(condition) {
      if (imputable_status[condition] == "non_imputable") {
        out <- list(
          variant_selection = rep(NA, max_var_selection),
          susie_result_trimmed = list(
            sets = list(cs = c(NULL)),
            pip = rep(NA, max_var_selection)
          )
        )
        return(out)
      }
      # if the condition is imputable
      out <- list(
        variant_names = get_nested_element(combined_all_data, c(condition, variable_name_obj)),
        susie_result_trimmed = get_nested_element(combined_all_data, c(condition, susie_obj))
      )
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
      # check if the condition has top_loci table
      if ("top_loci" %in% names(combined_all_data[[condition]][["preset_variants_result"]])) {
        out$top_loci <- get_nested_element(combined_all_data, c(condition, "preset_variants_result", "top_loci"))
        if (length(out$top_loci[, "variant_id"]) <= max_var_selection) {
          out$variant_selection <- out$top_loci$variant_id
        } else {
          # if top_loci variants more than max_var_selection, we loop through CS set to select top pip variants from each CS set
          # out$out$susie_result_trimmed is from "preset_variants_result"-"susie_result_trimmed"
          if (!is.null(out$susie_result_trimmed$sets$cs)) {
            # check total variant number in the CS sets
            totl_cs_indice <- do.call(sum, lapply(names(out$susie_result_trimmed$sets$cs), function(L) {
              length(out$susie_result_trimmed$sets$cs[[L]])
            }))
            if (max_var_selection <= totl_cs_indice) {
              selec_idx <- select_cs_var(out$susie_result_trimmed, max_var_selection)
              out$variant_selection <- out$variant_names[selec_idx]
            } else {
              # when top loci has larger number of variants than max_var_selection, but CS sets has less number than max_var_selection
              selec_idx <- unlist(out$susie_result_trimmed$sets$cs)
              out$variant_selection <- out$variant_names[selec_idx]
            }
          } else {
            top_idx <- order(out$susie_result_trimmed$pip, decreasing = TRUE)[1:max_var_selection]
            out$variant_selection <- out$variant_names[top_idx]
          }
        }
      } else {
        # if the condition did not come with the top_loci table, we select top pip variants from [[condition]]$preset_variants_result$susie_result_trimmed$pip
        print(paste0(condition, " do not have top_loci table, select top pip variants. "))
        top_idx <- order(out$susie_result_trimmed$pip, decreasing = TRUE)[1:max_var_selection]
        out$variant_selection <- out$variant_names[top_idx]
      }
      return(out)
    })
    names(var_selection_list) <- conditions
    return(var_selection_list)
  }

  # Internal function to align and merge weight matrices
  align_and_merge <- function(weights_list, variable_objs) {
    # Get the complete list of variant names across all files
    all_variants <- unique(unlist(variable_objs))
    consolidated_list <- list()
    # Fill the matrix with weights, aligning by variant names
    for (i in seq_along(weights_list)) {
      # Initialize the temp matrix with zeros
      existing_colnames <- character(0)
      temp_matrix <- matrix(0, nrow = length(all_variants), ncol = ncol(weights_list[[i]]))
      rownames(temp_matrix) <- all_variants
      idx <- match(variable_objs[[i]], all_variants)
      temp_matrix[idx, ] <- weights_list[[i]]
      # Ensure no duplicate column names
      new_colnames <- colnames(weights_list[[i]])
      dups <- duplicated(c(existing_colnames, new_colnames))
      if (any(dups)) {
        duplicated_names <- paste(c(existing_colnames, new_colnames)[dups], collapse = ", ")
        stop("Duplicate column names detected during merging process: ", duplicated_names, ".")
      }
      existing_colnames <- c(existing_colnames, new_colnames)
      consolidated_list[[i]] <- temp_matrix
      colnames(consolidated_list[[i]]) <- existing_colnames
    }
    return(consolidated_list)
  }

  # Internal function to consolidate weights for a given context
  consolidate_weights_list <- function(combined_all_data, conditions, variable_name_obj, twas_weights_table) {
    # Set default for 'conditions' if they are not specified
    if (is.null(conditions)) {
      conditions <- names(combined_all_data)
    }
    combined_weights_by_condition <- lapply(conditions, function(condition) {
      temp_list <- get_nested_element(combined_all_data, c(condition, twas_weights_table))
      temp_list <- temp_list[!names(temp_list) %in% "variant_names"]
      sapply(temp_list, cbind)
    })
    names(combined_weights_by_condition) <- conditions
    if (is.null(variable_name_obj)) {
      # Standard processing: Check for identical row numbers and consolidate
      row_numbers <- sapply(combined_weights_by_condition, function(data) nrow(data))
      if (length(unique(row_numbers)) > 1) {
        stop("Not all files have the same number of rows for the specified condition.")
      }
      weights <- combined_weights_by_condition
    } else {
      # Processing with variable_name_obj: Align and merge data, fill missing with zeros
      variable_objs <- lapply(conditions, function(condition) {
        get_nested_element(combined_all_data, c(condition, variable_name_obj))
      })
      weights <- align_and_merge(combined_weights_by_condition, variable_objs)
    }
    names(weights) <- conditions
    return(weights)
  }

  ## Load, validate, and consolidate data
  combined_all_data <- load_and_validate_data(weight_db_files, conditions, variable_name_obj)
  ## we first select best model, then select variants based on susie obj output, then extract weight value for these variants in selected model
  model_selection <- pick_best_model(combined_all_data, conditions, min_rsq_threshold, p_val_cutoff)
  region_info <- combined_all_data[[1]]$region_info
  if (is.null(unlist(model_selection))) {
    print("No model meets the p_value threshold and R-squared minimum in all conditions. Region is not imputable. ")
    model_selection$imputable <- FALSE
    return(list(susie_results = NULL, model_selection = model_selection, weights = NULL, region_info = region_info))
  } else {
    model_selection$imputable <- TRUE
  }
  # check for non-imputable conditions for this gene
  names(conditions) <- rep("impputable", length(conditions))
  non_imp_conditions <- na.omit(unlist(lapply(conditions, function(condition) {
    if (is.null(model_selection[[condition]])) {
      return(condition)
    }
  })))
  if (length(non_imp_conditions) >= 1) {
    names(conditions)[which(conditions %in% non_imp_conditions)] <- "non_imputable"
  }
  # at imputable conditions, sometimes susie_obj is not available.
  # if susie_obj is not available
  if (model_selection$imputable == TRUE) {
    for (k in conditions) {
      if (is.null(combined_all_data[[k]][[susie_obj[1]]][[susie_obj[2]]])) {
        model_selection[k] <- list(NULL)
        names(conditions)[match(k, conditions)] <- "non_imputable"
      }
    }
  }
  # run ctwas_select on all conditions first to keep list structure the same, regardless of imputable status.
  ctwas_select_result <- ctwas_select(combined_all_data, conditions, names(conditions), max_var_selection, variable_name_obj, susie_obj)
  weights <- consolidate_weights_list(combined_all_data, conditions, variable_name_obj, twas_weights_table)
  # assign the weights_selection as NULL in ctwas_select_result #
  if (length(non_imp_conditions) >= 1) {
    for (unimp_conds in non_imp_conditions) {
      ctwas_select_result[[unimp_conds]][["variant_selection"]] <- rep(NA, length(ctwas_select_result[[unimp_conds]][["variant_selection"]]))
    }
  }
  # return results
  return(list(
    susie_results = ctwas_select_result, model_selection = model_selection,
    weights = weights, region_info = region_info
  ))
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
      variants <- twas_rs[[study]][["weights"]][[condition]][["variant_selection"]][wgt_idx]
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


#' Scan Directory for LD Files and Generate Region Information Table
#'
#' Scans a specified directory for formatted LD (Linkage Disequilibrium) matrix files and corresponding bim files, generating a table with detailed region information. This function is designed to handle files named in a structured format, such as `LD_chr2_12345_17755.cor.xz.RDS` for RDS files and `LD_chr2_12345_17755.cor.xz.Rvar` for variance files.
#' @param ld_dir Path to the directory containing the LD matrix RDS files.
#' @return A data frame containing region information for all formatted LD matrices located in the specified directory. This includes details such as chromosome number, start and end positions, and file types.
#' @export
#' @examples
#' region_info <- scan_ld_files("path/to/ld_directory")                          
get_dir_region_info <- function(ld_dir) {
  ld_file_names <- list.files(path = ld_dir, pattern = "*.RDS", all.files = FALSE, full.names = FALSE)

  # Initialize vectors to store the positions, a2, and a1 elements
  chrom <- start <- stop <- vector("numeric", length(ld_file_names))
  # Extract the elements
  for (i in seq_along(ld_file_names)) {
    parts <- strsplit(ld_file_names[i], "_")[[1]]
    chrom[i] <- parts[2]
    start[i] <- parts[3]
    stop[i] <- gsub("\\..*", "", parts[4])
  }
  chrom <- as.integer(readr::parse_number(chrom))
  start <- as.integer(start)
  stop <- as.integer(stop)
  return(data.frame(
    chrom = chrom, start = start, stop = stop,
    region_tag = paste0("chr", chrom, ":", start, "-", stop),
    LD_matrix = paste0(ld_dir, "/", ld_file_names),
    SNP_info = paste0(ld_dir, "/", sub(".[^.]+$", "", ld_file_names), ".Rvar")
  ))
}
