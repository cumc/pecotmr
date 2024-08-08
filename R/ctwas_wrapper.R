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
#' results <- generate_twas_db(
#'   weight_db_file = "path/to/weights.rds", conditions = c("Mic", "Oli", "Exc"),
#'   max_var_selection = 10, min_rsq_threshold = 0.01, p_val_cutoff = 0.05
#' )
# select variants for ctwas weights input
generate_twas_db <- function(weight_db_file, contexts = NULL, variable_name_obj = c("preset_variants_result", "variant_names"),
                             twas_weights_table = "twas_weights", max_var_selection, min_rsq_threshold = 0.01,
                             p_val_cutoff = 0.05, data_type = NULL) {
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
        return(list(selected_model = c(NULL), imputable = FALSE)) # No significant model found
      } else {
        selected_model <- unlist(strsplit(selected_model, "_performance"))
        print(paste0("The selected best performing model for context ", context, " at region ", unique(region_names), " is ", selected_model, ". "))
        return(list(selected_model = selected_model, imputable = TRUE))
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
      variant_names <- get_nested_element(twas_data_combined, c("susie_results", context, "variant_names"))
      # susie_result_trimmed is from "preset_variants_result"-"susie_result_trimmed" from load_twas_weights()
      susie_result_trimmed <- get_nested_element(twas_data_combined, c("susie_results", context, "susie_result"))
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
  twas_data_combined <- load_twas_weights(weight_db_file, conditions = contexts, variable_name_obj = variable_name_obj, twas_weights_table = twas_weights_table)
  weights <- twas_data_combined$weights
  if (is.null(contexts)) {
    contexts <- names(weights)
  }
  # get gene name for the weight file
  gene <- unique(twas_data_combined$susie_results[[1]]$region_info$region_name)
  ## we first select best model, determine imputable contexts, then select variants based on susie obj output
  model_selection <- pick_best_model(twas_data_combined, contexts, min_rsq_threshold, p_val_cutoff)

  # Determine Imputable Contexts
  imputable_contexts <- c()
  for (contx in contexts) {
    if (isTRUE(model_selection[[contx]][["imputable"]])) imputable_contexts <- c(imputable_contexts, contx)
  }
  # select variants
  names(contexts) <- rep("non_imputable", length(contexts))
  names(contexts)[which(contexts %in% imputable_contexts)] <- "imputable"
  twas_select_result <- twas_select(twas_data_combined, contexts, names(contexts), max_var_selection, variable_name_obj)

  # output refined_twas_weights for imputable genes across all contexts
  refined_twas_weights <- list()
  if (length(imputable_contexts) > 0) {
    for (context in contexts) {
      # construct refined_twas_weights
      refined_twas_weights[[context]] <- list(selected_model = model_selection[[context]][["selected_model"]])
      refined_twas_weights[[context]][["is_imputable"]] <- model_selection[[context]][["imputable"]]
      if (!is.null(data_type)) refined_twas_weights[[context]][["data_type"]] <- data_type
      if (isTRUE(refined_twas_weights[[context]][["is_imputable"]])) {
        all_weight_variants <- rownames(weights[[context]])
        refined_twas_weights[[context]][["selected_top_variants"]] <- twas_select_result[[context]]
        refined_twas_weights[[context]][["selected_model_weights"]] <- weights[[context]][match(
          twas_select_result[[context]],
          all_weight_variants
        ), paste0(refined_twas_weights[[context]][["selected_model"]], "_weights")]
        if (model_selection[[context]][["selected_model"]] == "susie") {
          refined_twas_weights[[context]][["susie_parameters"]] <- get_nested_element(twas_data_combined, c(
            "susie_results", context,
            "susie_result"
          ))[c("lbf_variable", "X_column_scale_factors", "mu")]
          idx <- match(refined_twas_weights[[context]]$selected_top_variants, twas_data_combined$susie_results[[context]]$variant_names)
          # intersect susie_obj variants with selected weights variants
          refined_twas_weights[[context]][["susie_parameters"]] <- lapply(refined_twas_weights[[context]][["susie_parameters"]], function(obj) {
            if (is.vector(obj)) {
              return(obj[idx])
            } else {
              return(obj[, idx, drop = FALSE])
            }
          })
        }
      }
    }
    # return results
    return(list(
      refined_twas_weights = refined_twas_weights,
      weights = weights, gene = gene, cv_performance = twas_data_combined$twas_cv_performance,
      susie_results = twas_data_combined$susie_results
    ))
  } else {
    print(paste0("Weight input ", weight_db_file, " is non-imputable for all contexts. "))
    return(NULL)
  }
}



#' Utility function to specify the path to access the target list item in a nested list, especially when some list layers
#' in between are dynamic or uncertain.
#' @export
find_data <- function(x, depth_obj, show_path = FALSE, rm_null = TRUE, rm_dup = FALSE, docall = c) {
  depth <- as.integer(depth_obj[1])
  list_name <- if (length(depth_obj) > 1) depth_obj[2:length(depth_obj)] else NULL
  if (depth == 1 | depth == 0) {
    if (!is.null(list_name)) {
      if (list_name[1] %in% names(x)) {
        return(get_nested_element(x, list_name))
      }
    } else {
      return(x)
    }
  } else if (is.list(x)) {
    result <- lapply(x, find_data,
      depth_obj = c(depth - 1, list_name), show_path = show_path,
      rm_null = rm_null, rm_dup = rm_dup
    )
    shared_list_names <- list()
    if (isTRUE(rm_null)) {
      result <- result[!sapply(result, is.null)]
      result <- result[!sapply(result, function(x) length(x) == 0)]
    }
    if (isTRUE(rm_dup)) {
      unique_result <- list()
      unique_counter <- 1
      for (i in seq_along(result)) {
        duplicate_found <- FALSE
        for (j in seq_along(unique_result)) {
          if (identical(result[[i]], unique_result[[j]])) {
            duplicate_found <- TRUE
            shared_list_names[[paste0("unique_list_", j)]] <- c(shared_list_names[[paste0("unique_list_", j)]], names(result)[i])
            break
          }
        }
        if (!duplicate_found) {
          unique_name <- paste0("unique_list_", unique_counter)
          unique_result[[names(result)[i]]] <- result[[i]]
          shared_list_names[[unique_name]] <- names(result)[i]
          unique_counter <- unique_counter + 1
        }
      }
      result <- unique_result
    }

    if (isTRUE(show_path)) {
      if (length(shared_list_names) > 0 & depth == 2) result$shared_list_names <- shared_list_names
      return(result) # Carry original list structure
    } else {
      flat_result <- do.call(docall, unname(result))
      if (length(shared_list_names) > 0 & depth == 2) {
        names(result) <- paste0("unique_list_", 1:length(result))
        result$shared_list_names <- shared_list_names
        return(result)
      } else {
        return(flat_result) # Only return values
      }
    }
  } else {
    stop("Please input correct depth number. ")
  }
}

#' Utility function to convert LD region_ids to `region of interest` dataframe
#' @param ld_region_id A string of region in the format of chrom_start_end.
#' @importFrom data.table fread
#' @export
region_to_df <- function(ld_region_id, colnames = c("chrom", "start", "end")) {
  region_of_interest <- as.data.frame(do.call(rbind, lapply(strsplit(ld_region_id, "[_:-]"), function(x) as.integer(sub("chr", "", x)))))
  colnames(region_of_interest) <- colnames
  return(region_of_interest)
}
#' Utility function to load LD in ctwas analyses
#' @param ld_matrix_file_path A string of file path to the LD matrix.
#' @export
ctwas_ld_loader <- function(ld_matrix_file_path) {
  ld_loaded <- process_LD_matrix(ld_matrix_file_path, paste0(ld_matrix_file_path, ".bim"))
  ld_loaded <- ld_loaded$LD_matrix
  return(ld_loaded)
}

#' Utility function to format meta data dataframe for cTWAS analyses
#' @importFrom data.table fread
#' @export
get_ctwas_meta_data <- function(ld_meta_data_file, regions_table) {
  LD_info <- fread(ld_meta_data_file, header = TRUE, data.table = FALSE)
  colnames(LD_info)[1] <- "chrom"
  LD_info$region_id <- gsub("chr", "", paste(LD_info$chrom, LD_info$start, LD_info$end, sep = "_"))
  LD_info$LD_matrix <- paste0(dirname(ld_meta_data_file), "/", gsub(",.*$", "", LD_info$path))
  LD_info <- LD_info[, c("region_id", "LD_matrix")]
  region_info <- read.table(regions_table, sep = "\t", header = TRUE) # to get exact LD bim file without over-including neighboring LD's info.
  colnames(region_info)[1] <- "chrom"
  region_info$chrom <- as.integer(gsub("chr", "", region_info$chrom))
  region_info$region_id <- paste(region_info$chrom, region_info$start, region_info$stop, sep = "_")
  return(list(LD_info = LD_info, region_info = region_info))
}


#' Data Loader Function to load twas weights and variant names from output of `load_twas_weights` function.
#' @param twas_weights_data twas weights data generated from `load_twas_weights`
#' @importFrom R6 R6Class
#' @export
# Data Loader Function for importing twas weights from `load_twas_weights` function output.
twas_weights_loader <- R6Class("twas_weights_loader",
  public = list(
    data = NULL,
    initialize = function(twas_weights_data) {
      self$data <- twas_weights_data # This is the output from load_twas_weights()
    },
    get_weights = function() {
      return(weights = self$data$weights)
    },
    get_susie_obj = function() {
      return(lapply(self$data$susie_results, function(context_data) {
        context_data$susie_result[c("lbf_variable", "X_column_scale_factors", "mu")]
      }))
    },
    get_variant_names = function() {
      return(lapply(self$data$susie_results, function(context_data) context_data$variant_names))
    },
    get_gene_name = function() {
      return(unique(self$data$susie_results[[1]]$region_info$region_name))
    },
    get_chrom_num = function() {
      return(as.integer(unique(self$data$susie_results[[1]]$region_info$grange$chrom)))
    },
    get_data_type = function() {
      return(find_data(self$data$refined_twas_weights, c(2, "data_type"), show_path = TRUE))
    },
    get_data = function() {
      weights <- self$get_weights()
      susie_obj <- self$get_susie_obj()
      gene_name <- self$get_gene_name()
      variant_names <- self$get_variant_names()
      chrom <- self$get_chrom_num()
      data <- list(list(chrom = chrom, variant_names = variant_names, weights = weights, susie_obj = susie_obj))
      if ("refined_twas_weights" %in% names(self$data)) data[[1]] <- c(list(data_type = self$get_data_type()), data[[1]])
      names(data) <- gene_name
      return(data)
    }
  )
)

#' Data Loader Function to load twas weights and variant names from output of `generate_twas_db` function.
#' @importFrom R6 R6Class
#' @export
refined_twas_weights_loader <- R6Class("refined_twas_weights_loader",
  public = list(
    data = NULL,
    #' @param twas_weights_data Data from `generate_twas_db`
    initialize = function(twas_weights_data) {
      self$data <- twas_weights_data # This is the output from load_twas_weights()
    },
    get_weights = function(gene_data) {
      weights <- lapply(names(gene_data), function(context) {
        if (isTRUE(gene_data[[context]]$is_imputable)) {
          weight <- gene_data[[context]]$selected_model_weights
          weight <- matrix(weight, ncol = 1, dimnames = list(
            names(gene_data[[context]]$selected_model_weights),
            paste0(gene_data[[context]]$selected_model, "_weights")
          ))
          return(weight)
        }
      })
      names(weights) <- names(gene_data)
      return(weights)
    },
    get_susie_obj = function(gene_data) {
      susie_list <- lapply(names(gene_data), function(context) {
        if (isTRUE(gene_data[[context]]$selected_model == "susie")) {
          return(gene_data[[context]]$susie_parameters)
        }
      })
      names(susie_list) <- names(gene_data)
      return(susie_list)
    },
    get_variant_names = function(gene_data) {
      return(lapply(gene_data, function(context_data) {
        context_data$selected_top_variants
      }))
    },
    get_gene_name = function() {
      return(names(self$data[[1]]))
    },
    get_chrom_num = function() {
      return(as.integer(gsub("chr", "", (gsub("_.*$", "", names(self$data))))))
    },
    get_imputable_data = function(gene, gene_data) {
      imputable_context <- do.call(c, lapply(
        names(gene_data),
        function(context) if (isTRUE(self$data[[1]][[gene]][[context]]$is_imputable)) context
      ))
      return(gene_data[imputable_context])
    },
    get_data_type = function(gene_data) {
      return(lapply(gene_data, function(context_data) {
        context_data$data_type
      }))
    },
    get_data = function() {
      genes <- self$get_gene_name()
      chrom <- self$get_chrom_num()
      data <- lapply(genes, function(gene) {
        list(
          chrom = chrom,
          variant_names = self$get_imputable_data(gene, self$get_variant_names(self$data[[1]][[gene]])),
          weights = self$get_imputable_data(gene, self$get_weights(self$data[[1]][[gene]])),
          susie_obj = self$get_imputable_data(gene, self$get_susie_obj(self$data[[1]][[gene]])),
          data_type = self$get_imputable_data(gene, self$get_data_type(self$data[[1]][[gene]]))
        )
      })
      names(data) <- genes
      return(data)
    }
  )
)


#' Function to perform allele flip QC and harmonization on the weights and GWAS against LD for a region.
#' FIXME: GWAS loading function from Haochen for both tabix & column-mapping yml application
#'
#' Function Conditions:
#' - processes data in the format of either the output from load_twas_weights/generate_weight_db or
#'   refined_twas_weights_data from twas pipeline.
#' - For the first format, we expect there is only one gene's information, that can be accessed through `region_info_obj`
#'   and refined_twas_weights_data contains per region multiple gene's refined weights data.
#'
#' Main Steps:
#' 1. allele QC for TWAS weights against the LD meta
#' 2. allele QC for GWA Ssummary stats against the LD meta
#' 3. adjust susie/mvsusie weights based on the overlap variants
#'
#' @param twas_weights_data List of list of twas weights output from twas pipeline with following list items, with was weights
#' and variant name object specified by "variant_name_obj" and "twas_weights_table".
#' @param gwas_meta_file A file path for a dataframe table with column of "study_id", "chrom" (integer), "file_path",
#' "column_mapping_file". Each file in "file_path" column is tab-delimited dataframe of GWAS summary statistics with column name
#' "chrom" (or #chrom" if tabix-indexed), "pos", "A2", "A1".
#' @param ld_meta_file_path A tab-delimited data frame with colname "#chrom", "start", "end", "path", where "path" column
#' contains file paths for both LD matrix and bim file and is separated by ",". Bim file input would expect no headers, while the
#' columns are aligned in the order of "chrom", "variants", "GD", "pos", "A1", "A2", "variance", "allele_freq", "n_nomiss".
#' @param scale_weights TRUE/FALSE statement. If turn off, the post-qc/harmonized weights will not be scaled by SNP variance from LD.
#' @return A list of list for harmonized weights and dataframe of gwas summary statistics that is add to the original input of
#' twas_weights_data under each context.
#' @importFrom data.table fread
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom IRanges IRanges findOverlaps start end reduce
#' @export
harmonize_twas <- function(twas_weights_data, ld_meta_file_path, gwas_meta_file, twas_data_loader, scale_weights = TRUE) {
  # Function to group contexts based on start and end positions
  group_contexts_by_region <- function(twas_weights_data, gene, chrom, tolerance = 5000) {
    region_info_df <- do.call(rbind, lapply(names(twas_weights_data[[gene]]$variant_names), function(context) {
      wgt_range <- as.integer(sapply(
        twas_weights_data[[gene]][["variant_names"]][[context]],
        function(variant_id) {
          strsplit(variant_id, "\\:")[[1]][2]
        }
      ))
      min <- min(wgt_range)
      max <- max(wgt_range)
      data.frame(context = context, start = min, end = max)
    }))
    if (nrow(region_info_df) == 1) {
      # Handle case with only one context
      single_context_group <- list(
        context_group_1 = list(
          contexts = region_info_df$context,
          query_region = paste0(chrom, ":", region_info_df$start, "-", region_info_df$end),
          all_variants = unique(twas_weights_data[[gene]][["variant_names"]][[region_info_df$context]])
        )
      )
      return(single_context_group)
    }
    # Calculate distance matrix and perform hierarchical clustering
    clusters <- cutree(hclust(dist(region_info_df[, c("start", "end")])), h = tolerance)
    # Group contexts and determine query regions
    region_groups <- split(region_info_df, clusters) %>%
      lapply(function(group) {
        list(contexts = group$context, query_region = paste0(
          chrom, ":", min(group$start),
          "-", max(group$end)
        ))
      })
    # Create IRanges objects and merge overlapping intervals
    intervals <- IRanges(start = unlist(lapply(region_groups, function(context_group) {
      as.numeric(gsub(
        "^.*:\\s*|\\s*-.*$", "",
        context_group$query_region
      ))
    })), end = unlist(lapply(
      region_groups,
      function(context_group) as.numeric(sub("^.*?\\-", "", context_group$query_region))
    )))
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
    # add varinat names for coordinate extraction
    for (group in names(merged_groups)) {
      contexts <- merged_groups[[group]]$contexts
      merged_groups[[group]]$all_variants <- unique(do.call(c, lapply(
        contexts,
        function(context) {
          twas_weights_data[[gene]][["variant_names"]][[context]]
        }
      )))
    }
    return(merged_groups)
  }

  # Function to extract LD variance for the query region
  query_variance <- function(ld_meta_file_path, chrom, query_region, extract_coordinates) {
    # Filter the BIM data for the specified query region
    query_region_df <- region_to_df(query_region)
    query_region_df <- transform(query_region_df, start = start + 1, end = end - 1)
    # Find bim file in LD meta file that overlaps with the query region
    bim_file_path <- get_regional_ld_meta(ld_meta_file_path, query_region_df)$intersections$bim_file_paths
    bim_data <- do.call(rbind, lapply(bim_file_path, function(bim_file) fread(bim_file, header = FALSE, data.table = FALSE)))
    bim_filtered <- bim_data[bim_data$V4 %in% extract_coordinates$pos, , drop = FALSE]
    # Extract the variance column (7th column)
    variance_df <- bim_filtered[, c(1, 4, 5:7)]
    colnames(variance_df) <- c("chrom", "pos", "A1", "A2", "variance")
    return(variance_df)
  }

  # Step1: load TWAS weights data
  loader <- twas_data_loader$new(twas_weights_data)
  twas_weights_data <- loader$get_data()
  genes <- names(twas_weights_data)
  chrom <- as.integer(twas_weights_data[[1]]$chrom)

  gwas_meta_df <- fread(gwas_meta_file, header = TRUE, sep = "\t", data.table = FALSE)
  gwas_files <- unique(gwas_meta_df$file_path[gwas_meta_df$chrom == chrom])
  names(gwas_files) <- unique(gwas_meta_df$study_id[gwas_meta_df$chrom == chrom])
  results <- list()

  # Step2: Load LD for all genes by clustered context region
  region_variants <- variant_id_to_df(unique(unlist(find_data(twas_weights_data, c(2, "variant_names")))))
  region_of_interest <- data.frame(chrom = chrom, start = min(region_variants$pos), end = max(region_variants$pos))
  LD_list <- load_LD_matrix(ld_meta_file_path, region_of_interest, region_variants)

  # remove duplicate variants
  dup_idx <- which(duplicated(LD_list$combined_LD_variants))
  if (length(dup_idx) >= 1) {
    LD_list$combined_LD_variants <- LD_list$combined_LD_variants[-dup_idx]
    LD_list$combined_LD_matrix <- LD_list$combined_LD_matrix[-dup_idx, -dup_idx]
    LD_list$ref_panel <- LD_list$ref_panel[-dup_idx, ]
  }

  #########
  for (gene in genes) {
    results[[gene]] <- twas_weights_data[[gene]][c("chrom", "data_type")]
    # group contexts based on the variant position
    context_clusters <- group_contexts_by_region(twas_weights_data, gene, chrom, tolerance = 5000)

    ############
    for (context_group in names(context_clusters)) {
      contexts <- context_clusters[[context_group]]$contexts
      query_region <- context_clusters[[context_group]]$query_region
      region_of_interest <- region_to_df(query_region)
      all_variants <- context_clusters[[context_group]]$all_variants # all variants for this context group
      all_variants <- variant_id_to_df(all_variants)

      ############## Step3: load GWAS data for clustered context groups
      for (study in names(gwas_files)) {
        gwas_file <- gwas_files[study]
        gwas_sumstats <- as.data.frame(tabix_region(gwas_file, query_region)) # extension for yml file for column name mapping
        if (colnames(gwas_sumstats)[1] == "#chrom") colnames(gwas_sumstats)[1] <- "chrom" # colname update for tabix
        gwas_sumstats$chrom <- as.integer(gwas_sumstats$chrom)
        gwas_allele_flip <- allele_qc(gwas_sumstats[, c("chrom", "pos", "A1", "A2")], LD_list$combined_LD_variants, gwas_sumstats, c("beta", "z"),
          match_min_prop = 0.001
        )
        gwas_data_sumstats <- gwas_allele_flip$target_data_qced

        ############### context in the context group
        for (context in contexts) {
          weights_matrix <- twas_weights_data[[gene]][["weights"]][[context]]
          if (is.null(rownames(weights_matrix))) rownames(weights_matrix) <- twas_weights_data[[gene]][["variant_names"]][[context]]
          # Step4: adjust susie weights for susie_weights
          if ("susie_weights" %in% colnames(twas_weights_data[[gene]][["weights"]][[context]])) {
            adjusted_susie_weights <- adjust_susie_weights(twas_weights_data[[gene]],
              keep_variants = gwas_data_sumstats$variant_id, allele_qc = TRUE,
              variable_name_obj = c("variant_names", context),
              susie_obj = c("susie_obj", context),
              twas_weights_table = c("weights", context), LD_list$combined_LD_variants, match_min_prop = 0.001
            )
            weights_matrix_subset <- cbind(
              susie_weights = adjusted_susie_weights$adjusted_susie_weights,
              weights_matrix[adjusted_susie_weights$remained_variants_ids, !colnames(weights_matrix) %in% "susie_weights"]
            )
          } else {
            weights_matrix_subset <- weights_matrix
          }

          # Step5: harmonize weights
          if (!all(c("chrom", "pos", "A2", "A1") %in% colnames(weights_matrix_subset))) {
            weights_matrix_subset <- cbind(variant_id_to_df(rownames(weights_matrix_subset)), weights_matrix_subset)
          }
          weights_matrix_qced <- allele_qc(rownames(weights_matrix_subset), LD_list$combined_LD_variants, weights_matrix_subset,
            colnames(weights_matrix_subset)[!colnames(weights_matrix_subset) %in% c("chrom", "pos", "A2", "A1")],
            match_min_prop = 0.001, target_gwas = FALSE
          )
          weights_matrix_subset <- as.matrix(weights_matrix_qced$target_data_qced[, !colnames(weights_matrix_qced$target_data_qced) %in% c(
            "chrom",
            "pos", "A2", "A1", "variant_id"
          ), drop = FALSE])
          rownames(weights_matrix_subset) <- get_nested_element(weights_matrix_qced, c("target_data_qced", "variant_id"))

          # intersect post-qc gwas and post-qc weight variants
          gwas_LD_variants <- intersect(gwas_data_sumstats$variant_id, LD_list$combined_LD_variants)
          weights_matrix_subset <- weights_matrix_subset[which(rownames(weights_matrix_subset) %in% gwas_LD_variants), , drop = FALSE]
          results[[gene]][["variant_names"]][[context]] <- rownames(weights_matrix_subset)

          # Step6: scale weights by variance
          if (isTRUE(scale_weights)) {
            variance_df <- query_variance(ld_meta_file_path, chrom, query_region, all_variants) %>%
              mutate(variants = paste(chrom, pos, A2, A1, sep = ":"))
            weight_variants <- rownames(weights_matrix_subset)
            variance <- variance_df[match(weight_variants, variance_df$variants), "variance"]
            results[[gene]][["weights_qced"]][[context]] <- weights_matrix_subset * sqrt(variance)
          } else {
            results[[gene]][["weights_qced"]][[context]] <- weights_matrix_subset
          }
        }
        # Combine gwas sumstat across different context for a single context group
        gwas_data_sumstats <- gwas_data_sumstats[gwas_data_sumstats$variant_id %in% unique(unlist(results[[gene]][["variant_names"]])), ]
        results[[gene]][["gwas_qced"]][[study]] <- rbind(results[[gene]][["gwas_qced"]][[study]], gwas_data_sumstats)
        results[[gene]][["gwas_qced"]][[study]] <- results[[gene]][["gwas_qced"]][[study]][!duplicated(results[[gene]][["gwas_qced"]][[study]][, c("variant_id", "z")]), ]
      }
    }
    # extract LD matrix for variants intersect with gwas and twas weights at gene level
    all_gene_variants <- unique(find_data(results[[gene]][["gwas_qced"]], c(2, "variant_id")))
    var_indx <- match(all_gene_variants, LD_list$combined_LD_variants)
    results[[gene]][["LD"]] <- LD_list$combined_LD_matrix[var_indx, var_indx]
  }
  # return results
  return(results)
}


#' Function to generate weights_list from harmonized twas_weights_data for cTWAS multigroup analysis
#'
#' @param post_qc_twas_data output from harmonize_twas function for twas results
#' @return A list of list for weight information for each gene-context pair.
#' @export
get_ctwas_weights <- function(post_qc_twas_data, LD_meta_file_path) {
  chrom <- unique(find_data(post_qc_twas_data, c(2, "chrom")))
  if (length(chrom) != 1) stop("Data provided contains more than one chromosome. ")
  genes <- names(post_qc_twas_data)
  # Compile weights_list
  weights_list <- list()
  for (gene in genes) {
    contexts <- names(post_qc_twas_data[[gene]][["weights_qced"]])
    for (context in contexts) {
      standardized_context <- clean_context_names(context, gene)
      data_type <- post_qc_twas_data[[gene]][["data_type"]][[context]]
      colnames(post_qc_twas_data[[gene]][["weights_qced"]][[context]]) <- "weight"
      context_variants <- post_qc_twas_data[[gene]][["variant_names"]][[context]]
      context_range <- sapply(context_variants, function(variant) as.integer(strsplit(variant, "\\:")[[1]][2]))
      weights_list[[paste0(gene, "|", data_type, "_", standardized_context)]] <- list(
        chrom = chrom,
        p0 = min(context_range),
        p1 = max(context_range),
        wgt = post_qc_twas_data[[gene]][["weights_qced"]][[context]],
        R_wgt = post_qc_twas_data[[gene]]$LD[context_variants, context_variants, drop = FALSE], # ld_list$combined_LD_matrix[context_variants, context_variants],
        gene_name = gene,
        weight_name = paste0(data_type, "_", standardized_context),
        type = data_type,
        context = standardized_context,
        n_wgt = length(context_variants)
      )
    }
  }
  return(weights_list)
}
