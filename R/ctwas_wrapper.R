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
  weights <- setNames(
    lapply(
      names(twas_data_combined$weights),
      function(context) twas_data_combined$weights[[context]][names(twas_data_combined$susie_results[[context]]$X_column_scale_factors), , drop = FALSE]
    ),
    names(twas_data_combined$weights)
  )
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
        export_twas_weights_db[[context]][["model_weights"]] <- weights[[context]][, paste0(model_selection[[context]][["selected_model"]], "_weights"), drop = FALSE]
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
    message(paste0("Weight input ", weight_db_file, " is non-imputable for all contexts. "))
    return(NULL)
  }
}

#' Utility function to load LD in ctwas analyses, to interface with cTWAS package
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
    initialize = function(twas_weights_data, variable_name_obj = "variant_names", susie_obj = "susie_weights",
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
