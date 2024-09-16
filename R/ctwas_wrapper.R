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


#' Function to extract information from harmonize_twas output as cTWAS input. For each context within a gene, 
#' we select strong variants based on cs set and scaled weights size for the best TWAS weight method. 
#' 
get_ctwas_input <- function(post_qc_data, twas_weights_data, max_var_selection, LD_meta_file_path){
  # select variants based on cs and size of scaled weights
  select_variants <- function(post_qc_data, max_var_selection, twas_weights_data){
    variant_names_list <- find_data(post_qc_data, c(2, "variant_names"))
    contexts <- names(variant_names_list)
    selected_variants_by_context <- list()
    for (context in contexts){
      scaled_weights <- post_qc_data[[1]][["weights_qced"]][[context]][["scaled_weights"]][,1]
      if (is.null(names(scaled_weights))) names(scaled_weights) <- post_qc_data[[1]][["variant_names"]][[context]]
      if ("cs_variants" %in% names(find_data(twas_weights_data, c(3, context,"susie_weights_intermediate")))){
        cs_variant_list <- setNames(find_data(twas_weights_data, c(3, context,"susie_weights_intermediate", "cs_variants")), NULL)
        LD_ref_variants <- colnames(post_qc_data[[1]][["LD"]])
        cs_variant_list <- lapply(cs_variant_list, function(L){
          cs_variant_qced <- allele_qc(L, LD_ref_variants, L, match_min_prop = 0)
          paste0("chr", cs_variant_qced$target_data_qced$variant_id)
        })
        if (length(unlist(cs_variant_list))<= max_var_selection) {
          selected_variants_by_context[[context]] <- unlist(cs_variant_list)
          remain_var_num <- max_var_selection-length(selected_variants_by_context[[context]])
          if (remain_var_num > 0) {
            # when cs set do not provide enough variants to fill up max_var_selection number, we select remaining variants outside of cs 
            selected_variants_by_context[[context]] <- c(selected_variants_by_context[[context]], 
                names(scaled_weights[order(-abs(scaled_weights))])[!names(scaled_weights[order(-abs(scaled_weights))]) %in% selected_variants_by_context[[context]]][1:remain_var_num])
          }
        } else {
          # when cs set provides greater number of variants, we select top weight variants within cs variants
          top_weight_variants <- names(scaled_weights[order(-abs(scaled_weights))])
          selected_variants_by_context[[context]] <- top_weight_variants[top_weight_variants %in% unlist(cs_variant_list)][1:max_var_selection]
        }
      } else {
        # no cs variants present, we select variants based on scaled weight size 
        selected_variants_by_context[[context]] <- names(scaled_weights[order(-abs(scaled_weights))])[1:max_var_selection]
      }
    }
    return(selected_variants_by_context)
  }

  get_ctwas_weights <- function(post_qc_twas_data, LD_meta_file_path, selected_variant_list) {
    chrom <- unique(find_data(post_qc_twas_data, c(2, "chrom")))
    if (length(chrom) != 1) stop("Data provided contains more than one chromosome. ")
    genes <- names(post_qc_twas_data)
    # Compile weights_list
    weights_list <- list()
    for (gene in genes) {
      contexts <- names(post_qc_twas_data[[gene]][["weights_qced"]])
      for (context in contexts) {
        data_type <- post_qc_twas_data[[gene]][["data_type"]][[context]]
        colnames(post_qc_twas_data[[gene]][["weights_qced"]][[context]][["scaled_weights"]]) <- "weight"
        context_variants <- selected_variant_list[[context]]
        context_range <- sapply(context_variants, function(variant) as.integer(strsplit(variant, "\\:")[[1]][2]))
        weights_list[[paste0(gene, "|", data_type, "_", context)]] <- list(
          chrom = chrom,
          p0 = min(context_range),
          p1 = max(context_range),
          wgt = post_qc_twas_data[[gene]][["weights_qced"]][[context]][["scaled_weights"]][selected_variant_list[[context]], ,drop=FALSE],
          R_wgt = post_qc_twas_data[[gene]]$LD[context_variants, context_variants, drop = FALSE], # ld_list$combined_LD_matrix[context_variants, context_variants],
          molecular_id = gene,
          weight_name = paste0(data_type, "_", context),
          type = data_type,
          context = context,
          n_wgt = length(context_variants)
        )
      }
    }
    return(weights_list)
  }

  selected_variant_list <- select_variants(post_qc_data, max_var_selection, twas_weights_data)
  weights <- get_ctwas_weights(post_qc_data, LD_meta_file_path, selected_variant_list) # reshape weights for all gene-context pairs in the region for cTWAS analysis
  weights <- weights[!sapply(weights, is.null)]
  # load LD variants  
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