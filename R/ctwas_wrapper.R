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
  LD_info$LD_file <- paste0(dirname(ld_meta_data_file), "/", gsub(",.*$", "", LD_info$path))
  LD_info$SNP_file <- paste0(dirname(ld_meta_data_file), "/", gsub(",.*$", "", LD_info$path), ".bim")
  LD_info <- LD_info[, c("region_id", "LD_file", "SNP_file")]
  region_info <- read.table(regions_table, sep = "\t", header = TRUE) # to get exact LD bim file without over-including neighboring LD's info.
  colnames(region_info)[1] <- "chrom"
  region_info$chrom <- as.integer(gsub("chr", "", region_info$chrom))
  region_info$region_id <- paste(region_info$chrom, region_info$start, region_info$stop, sep = "_")
  return(list(LD_info = LD_info, region_info = region_info))
}


#' Function to extract information from harmonize_twas output as cTWAS input. For each imputable context within a gene, 
#' we select strong variants based on cs set and pip value.  
#' required ctwas package: remotes::install_github("xinhe-lab/ctwas",ref = "multigroup")
#' @importFrom ctwas compute_weight_LD_from_ref 
get_ctwas_input <- function(post_qc_data, twas_weights_data, LD_meta_file_path, cs_min_cor=0.8, min_pip_cutoff=0.1, 
                      min_var_selection=10, region_info=NULL, LD_info = NULL){
  # select variants based on cs and size of scaled weights
  select_variants <- function(post_qc_data, max_var_selection, twas_weights_data, cs_min_cor, min_pip_cutoff, min_var_selection){
    variant_names_list <- find_data(post_qc_data, c(2, "variant_names"))
    contexts <- names(variant_names_list)
    selected_variants_by_context <- list()
    for (context in contexts){
      if (find_data(twas_weights_data, c(3, context, "is_imputable"))){
        LD_ref_variants <- colnames(post_qc_data[[1]][["LD"]])
        scaled_weights <- post_qc_data[[1]][["weights_qced"]][[context]][["scaled_weights"]][,1]
        if (is.null(names(scaled_weights))) names(scaled_weights) <- post_qc_data[[1]][["variant_names"]][[context]]
        if ("cs_variants" %in% names(find_data(twas_weights_data, c(3, context,"susie_weights_intermediate")))){
          cs_min_abs_cor <- find_data(twas_weights_data, c(3, context,"susie_weights_intermediate", "cs_purity", "min.abs.corr"))
          for (L in 1:length(find_data(twas_weights_data, c(3, context,"susie_weights_intermediate","cs_variants")))){
            # we includ all variants in $cs_variant if min_abs_corr > cs_min_cor for the set 
            if (cs_min_abs_cor[L]>=cs_min_cor){
              cs_variants <-  find_data(twas_weights_data, c(3, context,"susie_weights_intermediate", "cs_variants"))[[L]]
              cs_variant_qced <- allele_qc(cs_variants, LD_ref_variants, cs_variants, match_min_prop = 0)
              selected_variants_by_context[[context]] <- paste0("chr", cs_variant_qced$target_data_qced$variant_id)
            }
          }
        }
        context_pip <- find_data(twas_weights_data, c(3, context, "susie_weights_intermediate", "pip"))
        susie_variants_qced <- allele_qc(names(context_pip), LD_ref_variants, names(context_pip), match_min_prop = 0)
        names(context_pip) <- paste0("chr",susie_variants_qced$target_data_qced$variant_id)
        available_variants <- names(context_pip[context_pip> min_pip_cutoff])[names(context_pip[context_pip> min_pip_cutoff]) %in% names(scaled_weights)]
        selected_variants_by_context[[context]] <- c(selected_variants_by_context[[context]], setdiff(available_variants, selected_variants_by_context[[context]]))
        if (length(selected_variants_by_context[[context]])< min_var_selection){
          remaining_var_num <- min_var_selection-length(selected_variants_by_context[[context]])
          selected_variants_by_context[[context]] <- c(selected_variants_by_context[[context]], 
                      setdiff(names(context_pip[order(-context_pip)]), selected_variants_by_context[[context]])[1:remaining_var_num])
        }
      }
    }
    return(selected_variants_by_context)
  }
  get_ctwas_weights <- function(post_qc_twas_data, selected_variant_list, twas_weights_data) {
    chrom <- unique(find_data(post_qc_twas_data, c(2, "chrom")))
    if (length(chrom) != 1) stop("Data provided contains more than one chromosome. ")
    genes <- names(post_qc_twas_data)
    # Compile weights_list
    weights_list <- list()
    for (gene in genes) {
      contexts <- names(post_qc_twas_data[[gene]][["weights_qced"]])
      for (context in contexts) {
        data_type <- post_qc_twas_data[[gene]][["data_type"]][[context]]
        postqc_scaled_weight <- post_qc_twas_data[[gene]][["weights_qced"]][[context]][["scaled_weights"]][selected_variant_list[[context]], ,drop=FALSE]
        # adjust susie weights for selected number of weight variants 
        if (colnames(post_qc_twas_data[[gene]][["weights_qced"]][[context]][["scaled_weights"]]) == "susie_weights"){
          gene_data <- find_data(twas_weights_data, c(2, gene))
          ld_variants <- colnames(post_qc_twas_data[[gene]][["LD"]])
          pre_qc_variants <- gene_data[[context]][["variant_names"]]
          qced_original_variants <- allele_qc(pre_qc_variants, ld_variants, pre_qc_variants, match_min_prop = 0)
          idx_in_original_variants <- match(selected_variant_list[[context]], 
                paste0("chr", qced_original_variants$qc_summary$variants_id_qced))
          gene_data[[context]]$susie_weights_intermediate$mu <- gene_data[[context]]$susie_weights_intermediate$mu[, idx_in_original_variants, drop=FALSE]
          gene_data[[context]]$susie_weights_intermediate$lbf_variable <- gene_data[[context]]$susie_weights_intermediate$lbf_variable[,idx_in_original_variants, drop=FALSE]
          gene_data[[context]]$susie_weights_intermediate$X_column_scale_factors <- gene_data[[context]]$susie_weights_intermediate$X_column_scale_factors[idx_in_original_variants]
          gene_data[[context]][["variant_names"]] <- selected_variant_list[[context]]
          gene_data[[context]][["model_weights"]] <- postqc_scaled_weight
          adjusted_susie_weights <- adjust_susie_weights(gene_data[[context]],
              keep_variants = gene_data[[context]][["variant_names"]], allele_qc = TRUE,
              variable_name_obj = c("variant_names"),
              susie_obj = c("susie_weights_intermediate"),
              twas_weights_table = c("model_weights"), ld_variants, match_min_prop = 0.001
          )
          postqc_scaled_weight <- matrix(adjusted_susie_weights$adjusted_susie_weights, ncol=1)
          rownames(postqc_scaled_weight) <- paste0("chr", adjusted_susie_weights$remained_variants_ids)
        }
        colnames(postqc_scaled_weight) <- "weight"
        context_variants <- selected_variant_list[[context]]
        context_range <- sapply(context_variants, function(variant) as.integer(strsplit(variant, "\\:")[[1]][2]))

        weights_list[[paste0(gene, "|", data_type, "_", context)]] <- list(
          chrom = chrom,
          p0 = min(context_range),
          p1 = max(context_range),
          wgt = postqc_scaled_weight,
          #R_wgt = post_qc_twas_data[[gene]]$LD[context_variants, context_variants, drop = FALSE], # ld_list$combined_LD_matrix[context_variants, context_variants],
          molecular_id = gene,
          weight_name = paste0(data_type, "_", context),
          type = data_type,
          context = context,
          n_wgt = length(context_variants)
        )
      }
    }
    # weights_list <- compute_weight_LD_from_ref(weights_list,
    #                                       weight_name=NULL,
    #                                       region_info = region_info,
    #                                       LD_map = LD_info,
    #                                       LD_format="custom", 
    #                                       LD_loader_fun = ctwas_ld_loader,
    #                                       ncore = 1)
    return(weights_list)
  }
  # select variants
  selected_variant_list <- select_variants(post_qc_data, max_var_selection, twas_weights_data, cs_min_cor, min_pip_cutoff, min_var_selection)
  weights <- get_ctwas_weights(post_qc_data, selected_variant_list, twas_weights_data) # reshape weights for all gene-context pairs in the region for cTWAS analysis
  weights <- weights[!sapply(weights, is.null)]

  # load LD variant names
  bim_file_paths <- unique(do.call(c, lapply(1:nrow(region_of_interest), function(region_row){
                        get_regional_ld_meta(LD_meta_file_path, region_of_interest[region_row,,drop=FALSE])$intersections$bim_file_paths})))
  snp_info <- lapply(bim_file_paths, function(file){
                    bimfile <- read.table(file, header = FALSE, sep="\t")[, c(1,2,4:8)]
                    bimfile$V2 <- gsub("chr", "", gsub("_", ":", bimfile$V2))
                    colnames(bimfile) <- c("chrom", "id", "pos", "alt", "ref", "variance", "allele_freq") # A1:alt, A2: ref
                    return(bimfile)})
  names(snp_info)<- do.call(c, lapply(bim_file_paths, function(x) { parts <- strsplit(basename(x), "[_:/.]")[[1]][1:3]
                          gsub("chr", "", paste(parts, collapse="_"))}))
  return(list(weights=weights, snp_info=snp_info))
}