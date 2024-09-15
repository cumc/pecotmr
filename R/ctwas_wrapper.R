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
#'   min_rsq_threshold = 0.01, p_val_cutoff = 0.05
#' )
# select variants for ctwas weights input
generate_twas_db <- function(weight_db_file, contexts = NULL, variable_name_obj = c("variant_names"),
                             susie_obj="susie_weights_intermediate", twas_weights_table = "twas_weights", 
                             min_rsq_threshold = 0.01, p_val_cutoff = 0.05, 
                             data_type_table = NULL, region_name) {
  # determine if the region is imputable and select the best model
  # Function to pick the best model based on adj_rsq and p-value
  pick_best_model <- function(twas_data_combined, contexts, min_rsq_threshold, p_val_cutoff) {
    best_rsq <- min_rsq_threshold
    # Extract performance table
    performance_tables <- lapply(contexts, function(context) {
      get_nested_element(twas_data_combined, c("twas_cv_performance", context))
    })
    names(performance_tables) <- contexts
    # Determine if a gene/region is imputable and select the best model
    model_selection <- lapply(contexts, function(context) {
      selected_model <- NULL
      available_models <- do.call(c, lapply(names(performance_tables[[context]]), function(model) {
        if (!is.na(performance_tables[[context]][[model]][, "rsq"])) {
          return(model)
        }
      }))
      if (length(available_models) <= 0) {
        print(paste0("No model provided TWAS cross validation performance metrics information at context ", context, ". "))
        return(NULL)
      }
      for (model in available_models) {
        model_data <- performance_tables[[context]][[model]]
        if (model_data[, "rsq"] >= best_rsq) {
          best_rsq <- model_data[, "rsq"]
          selected_model <- model
        }
      }
      if (is.null(selected_model)) { 
        print(paste0(
          "No model has p-value < ", p_val_cutoff, " and r2 >= ", min_rsq_threshold, ", skipping context ", context,
          " at region ", unique(region_name), ". "
        ))
        return(list(selected_model = c(NULL), imputable = FALSE)) # No significant model found
      } else {
        selected_model <- unlist(strsplit(selected_model, "_performance"))
        print(paste0("The selected best performing model for context ", context, " at region ", region_name, " is ", selected_model, ". "))
        return(list(selected_model = selected_model, imputable = TRUE))
      }
    })
    names(model_selection) <- contexts
    return(model_selection)
  }

  ## load twas weights
  twas_data_combined <- load_twas_weights(weight_db_file, conditions = contexts, variable_name_obj = variable_name_obj, susie_obj=susie_obj,
                                  twas_weights_table = twas_weights_table)
  if(!"weights" %in% names(twas_data_combined)) stop("TWAS weights not loaded. ")
  weights <- setNames(lapply(names(twas_data_combined$weights), 
                      function(context) twas_data_combined$weights[[context]][names(twas_data_combined$susie_results[[context]]$X_column_scale_factors),,drop=FALSE]), 
                      names(twas_data_combined$weights))
  if (is.null(contexts)) contexts <- names(weights)
  
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

  # output export_twas_weights_db for imputable genes across all contexts
  export_twas_weights_db <- list()
  if (length(imputable_contexts) > 0) {
    for (context in contexts) {
      # construct export_twas_weights_db
      export_twas_weights_db[[context]][["is_imputable"]] <- model_selection[[context]][["imputable"]]
      export_twas_weights_db[[context]][["variant_names"]] <- rownames(weights[[context]])
      export_twas_weights_db[[context]][["selected_model"]] <- model_selection[[context]][["selected_model"]]
      export_twas_weights_db[[context]][["susie_weights_intermediate"]] <- twas_data_combined$susie_result[[context]]
      if (model_selection[[context]]$imputable){
        export_twas_weights_db[[context]][["model_weights"]] <- weights[[context]][, paste0(model_selection[[context]][["selected_model"]], "_weights"), drop=FALSE] 
      } else {
        export_twas_weights_db[[context]][["model_weights"]] <- rep(NA, length(rownames(weights[[context]])))
      }
      if (!is.null(data_type_table)) export_twas_weights_db[[context]][["data_type"]] <- data_type_table$type[sapply(data_type_table$context, function(x) grepl(x, context))]
    }
    # return results
    return(list(
      export_twas_weights_db = export_twas_weights_db,
      weights = weights, gene = region_name, cv_performance = twas_data_combined$twas_cv_performance,
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


#' Data Loader Function to load twas weights and variant names from output of `generate_twas_db` function.
#' @param twas_weights_data twas weights data generated from `generate_twas_db`
#' @importFrom R6 R6Class
#' @export
# Data Loader Function for importing twas weights from `generate_twas_db` function output.
twas_weights_loader <- R6Class("twas_weights_loader",
  public = list(
    data = NULL,
    twas_weights_table = NULL, 
    susie_obj = NULL,           
    variable_name_obj = NULL,  
    initialize = function(twas_weights_data, variable_name_obj="variant_names", susie_obj="susie_weights",
                                     twas_weights_table = "twas_weights") {
      self$data <- twas_weights_data # This is the output from generate_twas_db()
      self$twas_weights_table <- twas_weights_table # Store the twas_weights_table name
      self$susie_obj <- susie_obj # Store the susie_obj name
      self$variable_name_obj <- variable_name_obj # Store the variable_name_obj name
    },
    get_weights = function() {
      return(weights = get_nested_element(self$data, self$twas_weights_table))
    },
    get_susie_obj = function() {
      return(lapply(self$data$export_twas_weights_db, function(context_data) get_nested_element(context_data, self$susie_obj)))
    },
    get_variant_names = function() {
      return(lapply(self$data$export_twas_weights_db, function(context_data) get_nested_element(context_data, self$variable_name_obj)))
    },
    get_gene_name = function() {
      return(self$data$gene)
    },
    get_chrom_num = function() {
      return(as.integer(readr::parse_number(strsplit(self$data$export_twas_weights_db[[1]]$variant_names[1], "[_:-]")[[1]][1])))
    },
    get_data_type = function() {
      return(find_data(self$data$export_twas_weights_db, c(2, "data_type"), show_path = TRUE))
    },
    get_data = function() {
      weights <- self$get_weights()
      susie_weights_intermediate <- self$get_susie_obj()
      gene_name <- self$get_gene_name()
      variant_names <- self$get_variant_names()
      chrom <- self$get_chrom_num()
      data <- list(list(chrom = chrom, variant_names = variant_names, weights = weights, susie_weights_intermediate = susie_weights_intermediate))
      if ("export_twas_weights_db" %in% names(self$data)) data[[1]] <- c(list(data_type = self$get_data_type()), data[[1]])
      names(data) <- gene_name
      return(data)
    }
  )
)

#' Data Loader Function to load top model's twas weights and variant from export_twas_weights_db generated from `generate_twas_db` function.
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
#' - processes data in the format of either the output from load_twas_weights/generate_twas_db or
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
harmonize_twas <- function(twas_weights_data, ld_meta_file_path, gwas_meta_file, twas_data_loader, scale_weights = FALSE) {
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
  loader <- twas_data_loader$new(twas_weights_data, variable_name_obj="variant_names", susie_obj="susie_weights_intermediate",
                                    twas_weights_table = "weights")
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

  # loop through genes: 
  for (gene in genes) {
    results[[gene]] <- twas_weights_data[[gene]][c("chrom", "data_type")]
    # group contexts based on the variant position
    context_clusters <- group_contexts_by_region(twas_weights_data, gene, chrom, tolerance = 5000)

    # loop through contexts: grouping contexts can be useful during ctwas data harmonization to stratify variants for LD loading
    for (context_group in names(context_clusters)) {
      contexts <- context_clusters[[context_group]]$contexts
      query_region <- context_clusters[[context_group]]$query_region
      region_of_interest <- region_to_df(query_region)
      all_variants <- context_clusters[[context_group]]$all_variants # all variants for this context group
      all_variants <- variant_id_to_df(all_variants)

      # Step3: load GWAS data for clustered context groups
      for (study in names(gwas_files)) {
        gwas_file <- gwas_files[study]
        gwas_sumstats <- as.data.frame(tabix_region(gwas_file, query_region)) # extension for yml file for column name mapping
        if (colnames(gwas_sumstats)[1] == "#chrom") colnames(gwas_sumstats)[1] <- "chrom" # colname update for tabix
        gwas_sumstats$chrom <- as.integer(gwas_sumstats$chrom)
        gwas_allele_flip <- allele_qc(gwas_sumstats[, c("chrom", "pos", "A1", "A2")], LD_list$combined_LD_variants, gwas_sumstats, c("beta", "z"),
          match_min_prop = 0.001
        )
        gwas_data_sumstats <- gwas_allele_flip$target_data_qced # post-qc gwas data that is flipped and corrected - gwas study level

        # loop through context within the context group:
        for (context in contexts) {
          weights_matrix <- twas_weights_data[[gene]][["weights"]][[context]]
          if (is.null(rownames(weights_matrix))) rownames(weights_matrix) <- twas_weights_data[[gene]][["variant_names"]][[context]]

          # Step 4: harmonize weights, flip allele 
          weights_matrix <- cbind(variant_id_to_df(rownames(weights_matrix)), weights_matrix)
          weights_matrix_qced <- allele_qc(rownames(weights_matrix), LD_list$combined_LD_variants, weights_matrix,
            colnames(weights_matrix)[!colnames(weights_matrix) %in% c("chrom", "pos", "A2", "A1")],
            match_min_prop = 0.001, target_gwas = FALSE
          )
          weights_matrix_subset <- as.matrix(weights_matrix_qced$target_data_qced[, !colnames(weights_matrix_qced$target_data_qced) %in% c(
            "chrom",
            "pos", "A2", "A1", "variant_id"
          ), drop = FALSE])
          rownames(weights_matrix_subset) <- weights_matrix_qced$target_data_qced$variant_id # weight variant names are flipped/corrected

          # intersect post-qc gwas and post-qc weight variants
          gwas_LD_variants <- intersect(gwas_data_sumstats$variant_id, LD_list$combined_LD_variants)
          weights_matrix_subset <- weights_matrix_subset[rownames(weights_matrix_subset) %in% gwas_LD_variants, , drop = FALSE]
          rownames(weights_matrix_subset) <- rownames(weights_matrix_subset)
          postqc_weight_variants <- rownames(weights_matrix_subset)

          # Step 5: adjust susie weights
          if ("susie_weights" %in% colnames(twas_weights_data[[gene]][["weights"]][[context]])) {
            adjusted_susie_weights <- adjust_susie_weights(twas_weights_data[[gene]],
              keep_variants = postqc_weight_variants, allele_qc = TRUE,
              variable_name_obj = c("variant_names", context),
              susie_obj = c("susie_weights_intermediate", context),
              twas_weights_table = c("weights", context), postqc_weight_variants, match_min_prop = 0.001
            )
            weights_matrix_subset <- cbind(
              susie_weights = setNames(adjusted_susie_weights$adjusted_susie_weights, adjusted_susie_weights$remained_variants_ids),
              weights_matrix_subset[gsub("chr", "", adjusted_susie_weights$remained_variants_ids), !colnames(weights_matrix_subset) %in% "susie_weights"]
            )
          }
          results[[gene]][["variant_names"]][[context]] <- rownames(weights_matrix_subset)

          # Step6: scale weights by variance
          if (isTRUE(scale_weights)) {
            variance_df <- query_variance(ld_meta_file_path, chrom, query_region, all_variants) %>%
              mutate(variants = paste(chrom, pos, A2, A1, sep = ":"))
            variance <- variance_df[match(rownames(weights_matrix_subset), variance_df$variants), "variance"]
            results[[gene]][["weights_qced"]][[context]] <- weights_matrix_subset * sqrt(variance)
          } else {
            results[[gene]][["weights_qced"]][[context]] <- weights_matrix_subset
          }
        }
        # Combine gwas sumstat across different context for a single context group based on all variants included in this gene/region
        gwas_data_sumstats$variant_id <- paste0("chr", gwas_data_sumstats$variant_id)
        gwas_data_sumstats <- gwas_data_sumstats[gwas_data_sumstats$variant_id %in% unique(unlist(results[[gene]][["variant_names"]])), ,drop=FALSE]
        results[[gene]][["gwas_qced"]][[study]] <- rbind(results[[gene]][["gwas_qced"]][[study]], gwas_data_sumstats)
        results[[gene]][["gwas_qced"]][[study]] <- results[[gene]][["gwas_qced"]][[study]][!duplicated(results[[gene]][["gwas_qced"]][[study]][, c("variant_id", "z")]), ]
      }
    }
    # extract LD matrix for variants intersect with gwas and twas weights at gene level
    all_gene_variants <- unique(find_data(results[[gene]][["gwas_qced"]], c(2, "variant_id")))
    var_indx <- match(all_gene_variants, paste0("chr", LD_list$combined_LD_variants))
    results[[gene]][["LD"]] <- LD_list$combined_LD_matrix[var_indx, var_indx]
    rownames(results[[gene]][["LD"]]) <- colnames(results[[gene]][["LD"]]) <- paste0("chr", colnames(results[[gene]][["LD"]]))
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
