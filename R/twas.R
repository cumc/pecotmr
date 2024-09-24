#' Determine Imputability and Variant Selection
#'
#' This function load TWAS weights and assesses the imputability of genes across different contexts,
#' and select top variants for imputable gene-context pair, and output the extracted variant weights
#' from best performing method model. Imputability of a gene-context pair are determined by having
#' least one of the twas methods' model being imputable based on cross livadtion metrics. The imputable
#' model is if the model surpasses the r-square and p-value threashold from cross validation metrics.
#' If non of the contexts of a gene has at least one method being imputable, then this gene will be
#' considered as unimputable, and do not return any weight results. This function is essential for preparing data
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
generate_twas_db <- function(weight_db_file, contexts = NULL, variable_name_obj = "variant_names",
                             susie_obj = "susie_weights_intermediate", twas_weights_table = "twas_weights",
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
        message(paste0("No model provided TWAS cross validation performance metrics information at context ", context, ". "))
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
        message(paste0(
          "No model has p-value < ", p_val_cutoff, " and r2 >= ", min_rsq_threshold, ", skipping context ", context,
          " at region ", unique(region_name), ". "
        ))
        return(list(selected_model = c(NULL), imputable = FALSE)) # No significant model found
      } else {
        selected_model <- unlist(strsplit(selected_model, "_performance"))
        message(paste0("The selected best performing model for context ", context, " at region ", region_name, " is ", selected_model, ". "))
        return(list(selected_model = selected_model, imputable = TRUE))
      }
    })
    names(model_selection) <- contexts
    return(model_selection)
  }

  ## load twas weights
  twas_data_combined <- load_twas_weights(weight_db_file,
    conditions = contexts, variable_name_obj = variable_name_obj, susie_obj = susie_obj,
    twas_weights_table = twas_weights_table
  )
  if (!"weights" %in% names(twas_data_combined)) stop("TWAS weights not loaded. ")
  weights <- twas_data_combined$weights
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
      if (model_selection[[context]]$imputable) {
        export_twas_weights_db[[context]][["model_weights"]] <- weights[[context]][, paste0(model_selection[[context]][["selected_model"]], "_weights")]
      } else {
        export_twas_weights_db[[context]][["model_weights"]] <- rep(NA, length(rownames(weights[[context]])))
      }
      if (!is.null(data_type_table)) export_twas_weights_db[[context]][["data_type"]] <- data_type_table$type[sapply(data_type_table$context, function(x) grepl(x, context))]
    }
    # return results
    return(list(
      export_twas_weights_db = export_twas_weights_db,
      weights = weights, molecular_id = region_name, cv_performance = twas_data_combined$twas_cv_performance,
      susie_results = twas_data_combined$susie_results
    ))
  } else {
    message(paste0("Weight input ", weight_db_file, " is non-imputable for all contexts. \n "))
    return(NULL)
  }
}

#' Function to perform allele flip QC and harmonization on the weights and GWAS against LD for a region.
#' FIXME: GWAS loading function from Haochen for both tabix & column-mapping yml application
#'
#' Function Conditions:
#' - processes data in the format of either the output from load_twas_weights/generate_twas_db or
#'   refined_twas_weights_data from twas pipeline.
#' - For the first format, we expect there is only one gene/events's information, that can be accessed through `region_info_obj`
#'   and refined_twas_weights_data contains per region multiple gene/event's refined weights data.
#'
#' Main Steps:
#' 1. allele QC for TWAS weights against the LD meta
#' 2. allele QC for GWA summary stats against the LD meta
#' 3. adjust susie/mvsusie weights based on the overlap variants
#'
#' @param twas_weights_data List of list of twas weights output from from generate_twas_db function.
#' @param gwas_meta_file A file path for a dataframe table with column of "study_id", "chrom" (integer), "file_path",
#' "column_mapping_file". Each file in "file_path" column is tab-delimited dataframe of GWAS summary statistics with column name
#' "chrom" (or #chrom" if tabix-indexed), "pos", "A2", "A1".
#' @param ld_meta_file_path A tab-delimited data frame with colname "#chrom", "start", "end", "path", where "path" column
#' contains file paths for both LD matrix and bim file and is separated by ",". Bim file input would expect no headers, while the
#' columns are aligned in the order of "chrom", "variants", "GD", "pos", "A1", "A2", "variance", "allele_freq", "n_nomiss".
#' @return A list of list for harmonized weights and dataframe of gwas summary statistics that is add to the original input of
#' twas_weights_data under each context.
#' @importFrom data.table fread
#' @importFrom readr parse_number
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom IRanges IRanges findOverlaps start end reduce
#' @export
harmonize_twas <- function(twas_weights_data, ld_meta_file_path, gwas_meta_file) {
  # Function to group contexts based on start and end positions
  group_contexts_by_region <- function(twas_weights_data, molecular_id, chrom, tolerance = 5000) {
    region_info_df <- do.call(rbind, lapply(names(twas_weights_data$export_twas_weights_db), function(context) {
      wgt_range <- as.integer(sapply(
        twas_weights_data[["export_twas_weights_db"]][[context]][["variant_names"]],
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
          all_variants = unique(twas_weights_data[["export_twas_weights_db"]][[region_info_df$context]][["variant_names"]])
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
          twas_weights_data[["export_twas_weights_db"]][[context]][["variant_names"]]
        }
      )))
    }
    return(merged_groups)
  }
  # function to load bim file variants
  load_bim_file_info <- function(ld_meta_file_path, region_of_interest) {
    bim_file_path <- get_regional_ld_meta(ld_meta_file_path, region_of_interest)$intersections$bim_file_paths
    bim_data <- lapply(bim_file_path, function(bim_file) fread(bim_file, header = FALSE, data.table = FALSE))
    return(bim_data)
  }
  # function to extract LD variance for the query region
  query_variance <- function(ld_variant_info_data, extract_coordinates) {
    ld_info_data <- do.call(rbind, ld_variant_info_data)
    ld_info_data_filtered <- ld_info_data[ld_info_data$V4 %in% extract_coordinates$pos, , drop = FALSE]
    variance_df <- ld_info_data_filtered[, c(1, 4, 5:7)] # Extract the variance column (7th column)
    colnames(variance_df) <- c("chrom", "pos", "A1", "A2", "variance")
    return(variance_df)
  }

  # Step 1: load TWAS weights data
  molecular_ids <- twas_weights_data$molecular_id
  chrom <- as.integer(parse_number(gsub(":.*$", "", twas_weights_data[["export_twas_weights_db"]][[1]][["variant_names"]][1])))
  gwas_meta_df <- fread(gwas_meta_file, header = TRUE, sep = "\t", data.table = FALSE)
  gwas_files <- unique(gwas_meta_df$file_path[gwas_meta_df$chrom == chrom])
  names(gwas_files) <- unique(gwas_meta_df$study_id[gwas_meta_df$chrom == chrom])
  results <- list()

  # Step 2: Load LD for all events/genes by clustered context region
  region_variants <- variant_id_to_df(unique(unlist(find_data(twas_weights_data$export_twas_weights_db, c(2, "variant_names")))))
  region_of_interest <- data.frame(chrom = chrom, start = min(region_variants$pos), end = max(region_variants$pos))
  LD_list <- load_LD_matrix(ld_meta_file_path, region_of_interest, region_variants)
  # load snp info once
  ld_variant_info <- load_bim_file_info(ld_meta_file_path, region_of_interest)
  snp_info <- setNames(lapply(ld_variant_info, function(info_table) {
    info_table <- info_table[, c(1, 2, 4:8)]
    info_table$V2 <- gsub("chr", "", gsub("_", ":", info_table$V2))
    colnames(info_table) <- c("chrom", "id", "pos", "alt", "ref", "variance", "allele_freq") # A1:alt, A2: ref
    return(info_table)
  }), sapply(names(ld_variant_info), function(x) gsub("chr", "", paste(strsplit(basename(x), "[_:/.]")[[1]][1:3], collapse = "_"))))

  # remove duplicate variants
  dup_idx <- which(duplicated(LD_list$combined_LD_variants))
  if (length(dup_idx) >= 1) {
    LD_list$combined_LD_variants <- LD_list$combined_LD_variants[-dup_idx]
    LD_list$combined_LD_matrix <- LD_list$combined_LD_matrix[-dup_idx, -dup_idx]
    LD_list$ref_panel <- LD_list$ref_panel[-dup_idx, ]
  }

  # loop through genes/events:
  for (molecular_id in molecular_ids) {
    results[[molecular_id]][["chrom"]] <- chrom
    results[[molecular_id]][["data_type"]] <- if ("data_type" %in% names(twas_weights_data$export_twas_weights_db[[1]])) lapply(twas_weights_data$export_twas_weights_db, function(x) x$data_type)
    # group contexts based on the variant position
    context_clusters <- group_contexts_by_region(twas_weights_data, molecular_id, chrom, tolerance = 5000)

    # loop through contexts: grouping contexts can be useful during TWAS data harmonization to stratify variants for LD loading
    for (context_group in names(context_clusters)) {
      contexts <- context_clusters[[context_group]]$contexts
      query_region <- context_clusters[[context_group]]$query_region
      region_of_interest <- region_to_df(query_region)
      all_variants <- context_clusters[[context_group]]$all_variants # all variants for this context group
      all_variants <- variant_id_to_df(all_variants)

      # Step 3: load GWAS data for clustered context groups
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
          results[[molecular_id]][["model_selection"]][[context]] <- twas_weights_data$export_twas_weights_db[[context]][["selected_model"]]
          weights_matrix <- twas_weights_data[["weights"]][[context]]
          if (is.null(rownames(weights_matrix))) rownames(weights_matrix) <- twas_weights_data$export_twas_weights_db[[context]][["variant_names"]]

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

          # Step 5: adjust SuSiE weights based on available variants
          if ("susie_weights" %in% colnames(twas_weights_data[["weights"]][[context]])) {
            adjusted_susie_weights <- adjust_susie_weights(twas_weights_data,
              keep_variants = postqc_weight_variants, allele_qc = TRUE,
              variable_name_obj = c("export_twas_weights_db", context, "variant_names"),
              susie_obj = c("susie_results", context),
              twas_weights_table = c("weights", context), postqc_weight_variants, match_min_prop = 0.001
            )
            weights_matrix_subset <- cbind(
              susie_weights = setNames(adjusted_susie_weights$adjusted_susie_weights, adjusted_susie_weights$remained_variants_ids),
              weights_matrix_subset[gsub("chr", "", adjusted_susie_weights$remained_variants_ids), !colnames(weights_matrix_subset) %in% "susie_weights"]
            )
          }
          rownames(weights_matrix_subset) <- if (!grepl("^chr", rownames(weights_matrix_subset)[1])) paste0("chr", rownames(weights_matrix_subset)) else rownames(weights_matrix_subset)
          results[[molecular_id]][["variant_names"]][[context]] <- rownames(weights_matrix_subset)

          # Step 6: scale weights by variance
          variance_df <- query_variance(ld_variant_info, all_variants) %>%
            mutate(variants = paste(chrom, pos, A2, A1, sep = ":"))
          variance <- variance_df[match(rownames(weights_matrix_subset), paste0("chr", variance_df$variants)), "variance"]
          results[[molecular_id]][["weights_qced"]][[context]] <- list(scaled_weights = weights_matrix_subset * sqrt(variance), weights = weights_matrix_subset)
        }

        # Combine gwas sumstat across different context for a single context group based on all variants included in this molecular_id/gene/region
        gwas_data_sumstats$variant_id <- paste0("chr", gwas_data_sumstats$variant_id)
        gwas_data_sumstats <- gwas_data_sumstats[gwas_data_sumstats$variant_id %in% unique(unlist(results[[molecular_id]][["variant_names"]])), , drop = FALSE]
        results[[molecular_id]][["gwas_qced"]][[study]] <- rbind(results[[molecular_id]][["gwas_qced"]][[study]], gwas_data_sumstats)
        results[[molecular_id]][["gwas_qced"]][[study]] <- results[[molecular_id]][["gwas_qced"]][[study]][!duplicated(results[[molecular_id]][["gwas_qced"]][[study]][, c("variant_id", "z")]), ]
      }
    }
    # extract LD matrix for variants intersect with gwas and twas weights at molecular_id level
    all_molecular_variants <- unique(find_data(results[[molecular_id]][["gwas_qced"]], c(2, "variant_id")))
    var_indx <- match(all_molecular_variants, paste0("chr", LD_list$combined_LD_variants))
    results[[molecular_id]][["LD"]] <- LD_list$combined_LD_matrix[var_indx, var_indx]
    rownames(results[[molecular_id]][["LD"]]) <- colnames(results[[molecular_id]][["LD"]]) <- paste0("chr", colnames(results[[molecular_id]][["LD"]]))
  }
  # return results
  return(list(twas_data_qced = results, snp_info = snp_info))
}

#' Function to perform TWAS analysis for across multiple contexts.
#' This function peforms TWAS analysis for multiple contexts for imputable genes within an LD region and summarize the twas results.
#' @param twas_weights_data List of list of twas weights output from generate_twas_db function.
#' @param region_block A string with informaiton of chromosome number, startind position, and ending position of LD block conneced with "_".
#' @return A list of list containing twas result table and formatted TWAS data compatible with ctwas_sumstats() function.
#' \itemize{
#'   \item{twas_table}{ A dataframe of twas results summary is generated for each gene-contexts-method pair of all methods for imputable genes.}
#'   \item{twas_data}{ A list of list containing formatted TWAS data.}
#' }
#' @export
twas_pipeline <- function(twas_weights_data,
                          ld_meta_file_path,
                          gwas_meta_file,
                          region_block,
                          mr_pval_cutoff = 0.05) {
  # internal function to format TWAS output
  format_twas_data <- function(post_qc_twas_data, twas_table) {
    weights_list <- do.call(c, lapply(names(post_qc_twas_data), function(molecular_id) {
      contexts <- names(post_qc_twas_data[[molecular_id]][["weights_qced"]])
      do.call(c, lapply(contexts, function(context) {
        weight <- list()
        data_type <- post_qc_twas_data[[molecular_id]][["data_type"]][[context]]
        model_selected <- post_qc_twas_data[[molecular_id]][["model_selection"]][[context]]
        postqc_scaled_weight <- post_qc_twas_data[[molecular_id]][["weights_qced"]][[context]][["scaled_weights"]][, paste0(model_selected, "_weights"), drop = FALSE]
        colnames(postqc_scaled_weight) <- "weight"
        context_variants <- rownames(post_qc_twas_data[[molecular_id]][["weights_qced"]][[context]][["scaled_weights"]])
        context_range <- as.integer(sapply(context_variants, function(variant) strsplit(variant, "\\:")[[1]][2]))
        weight[[paste0(molecular_id, "|", data_type, "_", context)]] <- list(
          chrom = chrom, p0 = min(context_range), p1 = max(context_range),
          wgt = postqc_scaled_weight, molecular_id = molecular_id, weight_name = paste0(data_type, "_", context), type = data_type,
          context = context, n_wgt = length(context_variants)
        )
        return(weight)
      }))
    }))
    weights <- weights_list[!sapply(weights_list, is.null)]
    # gene_z table
    twas_table <- twas_table[na.omit(twas_table$is_selected_method), , drop = FALSE]
    twas_table$id <- paste0(twas_table$molecular_id, "|", twas_table$type, "_", twas_table$context)
    twas_table$z <- twas_table$twas_z
    twas_table$group <- paste0(twas_table$context, "|", twas_table$type)
    twas_table <- twas_table[, c("id", "z", "type", "context", "group", "gwas_study"), drop = FALSE]
    studies <- unique(twas_table$gwas_study)
    z_gene_list <- list()
    z_snp <- list()
    for (study in studies) {
      z_gene_list[[study]] <- twas_table[twas_table$gwas_study == study, , drop = FALSE]
      z_snp[[study]] <- do.call(rbind, lapply(post_qc_twas_data, function(x) find_data(x, c(1, "gwas_qced", study), docall = rbind)))
      colnames(z_snp[[study]])[which(colnames(z_snp[[study]]) == "variant_id")] <- "id"
      z_snp[[study]] <- z_snp[[study]][, c("id", "A1", "A2", "z")]
      z_snp[[study]] <- z_snp[[study]][!duplicated(z_snp[[study]]$id), , drop = FALSE]
    }
    return(list(weights = weights, z_gene = z_gene_list, z_snp = z_snp))
  }

  # Step 1: TWAS and MR analysis for all methods for imputable gene
  twas_results_db <- lapply(names(twas_weights_data), function(weight_db) {
    # harmonize twas weights and gwas sumstats against LD
    twas_data_qced_result <- harmonize_twas(twas_weights_data[[weight_db]], ld_meta_file_path, gwas_meta_file)
    twas_data_qced <- twas_data_qced_result$twas_data_qced
    molecular_id <- names(twas_data_qced)
    if (length(molecular_id) < 1) stop(paste0("No gene's data processed at harmonization process for ", weight_db, ". "))
    contexts <- names(twas_data_qced[[molecular_id]][["weights_qced"]])
    gwas_studies <- names(twas_data_qced[[molecular_id]][["gwas_qced"]])

    # Combined loop for TWAS and MR analysis
    mr_cols <- c("gene_name", "num_CS", "num_IV", "cpip", "meta_eff", "se_meta_eff", "meta_pval", "Q", "Q_pval", "I2")

    # Nested lapply for contexts and gwas studies
    twas_gene_results <- lapply(contexts, function(context) {
      study_results <- lapply(gwas_studies, function(study) {
        # twas analysis
        twas_rs <- twas_analysis(
          twas_data_qced[[molecular_id]][["weights_qced"]][[context]][["weights"]], twas_data_qced[[molecular_id]][["gwas_qced"]][[study]],
          twas_data_qced[[molecular_id]][["LD"]], rownames(twas_data_qced[[molecular_id]][["weights_qced"]][[context]][["weights"]])
        )
        twas_rs_df <- data.frame(
          gwas_study = study, method = sub("_[^_]+$", "", names(twas_rs)), twas_z = find_data(twas_rs, c(2, "z")),
          twas_pval = find_data(twas_rs, c(2, "pval")), context = context, molecular_id = molecular_id
        )
        # MR analysis
        if (any(na.omit(twas_rs_df$twas_pval) < mr_pval_cutoff) & "top_loci" %in% names(twas_weights_data[[weight_db]]$susie_results[[context]])) {
          mr_formatted_input <- mr_format(twas_weights_data[[weight_db]], context, twas_data_qced[[molecular_id]][["gwas_qced"]][[study]],
            coverage = "cs_coverage_0.95", allele_qc = TRUE, molecular_name_obj = c("molecular_id")
          )
          if (all(is.na(mr_formatted_input$bhat_y))) {
            # FIXME: after updating gwas beta and se NA problem, mr analysis will be restored
            mr_rs_df <- as.data.frame(matrix(rep(NA, length(mr_cols)), nrow = 1))
            colnames(mr_rs_df) <- mr_cols
          } else {
            mr_rs_df <- as.data.frame(mr_analysis(mr_formatted_input, cpip_cutoff = 0.5))
          }
        } else {
          mr_rs_df <- as.data.frame(matrix(rep(NA, length(mr_cols)), nrow = 1))
          colnames(mr_rs_df) <- mr_cols
        }
        mr_rs_df$context <- context
        mr_rs_df$gwas_study <- study
        return(list(twas_rs_df = twas_rs_df, mr_rs_df = mr_rs_df))
      })
      twas_context_table <- do.call(rbind, lapply(study_results, function(x) x$twas_rs_df))
      mr_context_table <- do.call(rbind, lapply(study_results, function(x) x$mr_rs_df))
      return(list(twas_context_table = twas_context_table, mr_context_table = mr_context_table))
    })
    twas_gene_table <- do.call(rbind, lapply(twas_gene_results, function(x) x$twas_context_table))
    mr_gene_table <- do.call(rbind, lapply(twas_gene_results, function(x) x$mr_context_table))
    return(list(twas_table = twas_gene_table, twas_data_qced = twas_data_qced, mr_result = mr_gene_table, snp_info = twas_data_qced_result$snp_info))
  })
  twas_results_table <- do.call(rbind, lapply(twas_results_db, function(x) x$twas_table))
  mr_results <- do.call(rbind, lapply(twas_results_db, function(x) x$mr_result))
  twas_data <- do.call(c, lapply(twas_results_db, function(x) x$twas_data_qced))
  snp_info <- do.call(c, lapply(twas_results_db, function(x) x$snp_info))

  # Step 2: Summarize and merge twas cv results and region information for all methods for all contexts for imputable genes.
  twas_table <- do.call(rbind, lapply(names(twas_weights_data), function(molecular_id) {
    contexts <- names(twas_weights_data[[molecular_id]][["export_twas_weights_db"]])
    # merge twas_cv information for same gene across all weight db files
    cv_data <- do.call(c, lapply(names(twas_weights_data), function(file) {
      if (twas_weights_data[[molecular_id]]$molecular_id == molecular_id) twas_weights_data[[molecular_id]]$cv_performance
    }))
    # loop through each context for all methods
    gene_table <- do.call(rbind, lapply(contexts, function(context) {
      methods <- sub("_[^_]+$", "", names(cv_data[[context]]))
      is_imputable <- twas_weights_data[[molecular_id]][["export_twas_weights_db"]][[context]]$is_imputable
      selected_method <- twas_weights_data[[molecular_id]][["export_twas_weights_db"]][[context]]$selected_model
      if (is.null(selected_method)) selected_method <- NA
      is_selected_method <- ifelse(methods == selected_method, TRUE, FALSE)
      cv_rsqs <- sapply(cv_data[[context]], function(x) x[, "rsq"])
      cv_pvals <- sapply(cv_data[[context]], function(x) x[, "pval"])
      context_table <- data.frame(
        context = context, method = methods, is_imputable = is_imputable, is_selected_method = is_selected_method,
        rsq_cv = cv_rsqs, pval_cv = cv_pvals, type = twas_weights_data[[molecular_id]][["export_twas_weights_db"]][[context]]$data_type
      )
      return(context_table)
    }))
    gene_table$molecular_id <- molecular_id
    gene_table$start <- xqtl_meta_df$start[xqtl_meta_df$region_id == molecular_id]
    gene_table$end <- xqtl_meta_df$end[xqtl_meta_df$region_id == molecular_id]
    gene_table$TSS <- xqtl_meta_df$TSS[xqtl_meta_df$region_id == molecular_id]
    return(gene_table)
  }))
  twas_table$chr <- chrom
  twas_table$block <- region_block

  # Step 3. merge twas result table
  colname_ordered <- c(
    "chr", "start", "end", "molecular_id", "TSS", "context", "gwas_study", "method", "is_imputable", "is_selected_method",
    "rsq_cv", "pval_cv", "twas_z", "twas_pval", "type", "block"
  )
  twas_table <- merge(twas_table, twas_results_table, by = c("molecular_id", "context", "method"))
  twas_data_subset <- format_twas_data(twas_data, twas_table)
  twas_data_subset$snp_info <- snp_info
  return(list(twas_result = twas_table[, colname_ordered], twas_data = twas_data_subset, mr_result = mr_results))
}

#' Calculate TWAS z-score and p-value
#'
#' This function calculates the TWAS z-score and p-value given the weights, z-scores,
#' and optionally the correlation matrix (R) or the genotype matrix (X).
#'
#' @param weights A numeric vector of weights.
#' @param z A numeric vector of z-scores.
#' @param R An optional correlation matrix. If not provided, it will be calculated from the genotype matrix X.
#' @param X An optional genotype matrix. If R is not provided, X must be supplied to calculate the correlation matrix.
#'
#' @return A list containing the following elements:
#' \itemize{
#'   \item z: The TWAS z-score.
#'   \item pval: The corresponding p-value.
#' }
#'
#' @importFrom Rfast cora
#' @importFrom stats cor pchisq
#'
#' @export
twas_z <- function(weights, z, R = NULL, X = NULL) {
  # Check that weights and z-scores have the same length
  if (length(weights) != length(z)) {
    stop("Weights and z-scores must have the same length.")
  }

  if (is.null(R)) R <- compute_LD(X)

  stat <- t(weights) %*% z
  denom <- t(weights) %*% R %*% weights
  zscore <- stat / sqrt(denom)
  pval <- pchisq(zscore * zscore, 1, lower.tail = FALSE)

  return(list(z = zscore, pval = pval))
}

#' Multi-condition TWAS joint test
#'
#' This function performs a multi-condition TWAS joint test using the GBJ method.
#' It assumes that the input genotype matrix (X) is standardized.
#'
#' @param R An optional correlation matrix. If not provided, it will be calculated from the genotype matrix X.
#' @param X An optional genotype matrix. If R is not provided, X must be supplied to calculate the correlation matrix.
#' @param weights A matrix of weights, where each column corresponds to a different condition.
#' @param z A vector of GWAS z-scores.
#'
#' @return A list containing the following elements:
#' \itemize{
#'   \item Z: A matrix of TWAS z-scores and p-values for each condition.
#'   \item GBJ: The result of the GBJ test.
#' }
#'
#' @importFrom GBJ GBJ
#' @importFrom stats cor pnorm
#'
#' @export
twas_joint_z <- function(weights, z, R = NULL, X = NULL) {
  # Check that weights and z-scores have the same number of rows
  if (nrow(weights) != length(z)) {
    stop("Number of rows in weights must match the length of z-scores.")
  }
  if (is.null(R)) R <- compute_LD(X)
  idx <- which(rownames(R) %in% rownames(weights))
  D <- R[idx, idx]
  cov_y <- crossprod(weights, D) %*% weights
  y_sd <- sqrt(diag(cov_y))
  x_sd <- rep(1, nrow(weights)) # Assuming X is standardized

  # Get gamma matrix MxM (snp x snp)
  g <- lapply(colnames(weights), function(x) {
    gm <- diag(x_sd / y_sd[x], length(x_sd), length(x_sd))
    return(gm)
  })
  names(g) <- colnames(weights)

  ######### Get TWAS - Z statistics & P-value, GBJ test ########
  z_matrix <- do.call(rbind, lapply(colnames(weights), function(x) {
    Zi <- crossprod(weights[, x], g[[x]]) %*% as.numeric(z)
    pval <- 2 * pnorm(abs(Zi), lower.tail = FALSE)
    Zp <- c(Zi, pval)
    names(Zp) <- c("Z", "pval")
    return(Zp)
  }))
  rownames(z_matrix) <- colnames(weights)

  # GBJ test
  lam <- matrix(rep(NA, ncol(weights) * nrow(weights)), nrow = ncol(weights))
  rownames(lam) <- colnames(weights)
  for (p in colnames(weights)) {
    la <- as.matrix(weights[, p] %*% g[[p]])
    lam[p, ] <- la
  }
  sig <- tcrossprod((lam %*% D), lam)
  gbj <- GBJ(test_stats = z_matrix[, 1], cor_mat = sig)

  rs <- list("Z" = z_matrix, "GBJ" = gbj)
  return(rs)
}

#' TWAS Analysis
#'
#' Performs TWAS analysis using the provided weights matrix, GWAS summary statistics database,
#' and LD matrix. It extracts the necessary GWAS summary statistics and LD matrix based on the
#' specified variants and computes the z-score and p-value for each gene.
#'
#' @param weights_matrix A matrix containing weights for all methods.
#' @param gwas_sumstats_db A data frame containing the GWAS summary statistics.
#' @param LD_matrix A matrix representing linkage disequilibrium between variants.
#' @param extract_variants_objs A vector of variant identifiers to extract from the GWAS and LD matrix.
#'
#' @return A list with TWAS z-scores and p-values across four methods for each gene.
#' @export
twas_analysis <- function(weights_matrix, gwas_sumstats_db, LD_matrix, extract_variants_objs) {
  # Extract gwas_sumstats
  gwas_sumstats_subset <- gwas_sumstats_db[match(extract_variants_objs, gwas_sumstats_db$variant_id), ]
  # Validate that the GWAS subset is not empty
  if (nrow(gwas_sumstats_subset) == 0 | all(is.na(gwas_sumstats_subset))) {
    stop("No GWAS summary statistics found for the specified variants.")
  }
  # Check if extract_variants_objs are in the rownames of LD_matrix
  valid_indices <- extract_variants_objs %in% rownames(LD_matrix)
  if (!any(valid_indices)) {
    stop("None of the specified variants are present in the LD matrix.")
  }
  # Extract only the valid indices from extract_variants_objs
  valid_variants_objs <- extract_variants_objs[valid_indices]
  # Extract LD_matrix subset using valid indices
  LD_matrix_subset <- LD_matrix[valid_variants_objs, valid_variants_objs]
  # Caculate the z score and pvalue of each gene
  twas_z_pval <- apply(
    as.matrix(weights_matrix), 2,
    function(x) twas_z(x, gwas_sumstats_subset$z, R = LD_matrix_subset)
  )
  return(twas_z_pval)
}
