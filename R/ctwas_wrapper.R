#' Utility function to load LD in ctwas analyses, to interface with cTWAS package
#' @param ld_matrix_file_path A string of file path to the LD matrix.
#' @export
ctwas_ld_loader <- function(ld_matrix_file_path) {
  ld_loaded <- process_LD_matrix(ld_matrix_file_path, paste0(ld_matrix_file_path, ".bim"))
  ld_loaded <- ld_loaded$LD_matrix
  return(ld_loaded)
}
#' @importFrom data.table fread
#' @export
ctwas_bimfile_loader <- function(bim_file_path) {
  snp_info <- fread(bim_file_path)
  if (ncol(snp_info) == 9) {
    colnames(snp_info) <- c("chrom", "id", "GD", "pos", "alt", "ref", "variance", "allele_freq", "n_nomiss")
  } else {
    colnames(snp_info) <- c("chrom", "id", "GD", "pos", "alt", "ref")
  }
  snp_info$id <- gsub("_", ":", gsub("chr", "", snp_info$id))
  return(snp_info)
}

#' Utility function to format meta data dataframe for cTWAS analyses
#' @importFrom data.table fread
#' @export
get_ctwas_meta_data <- function(ld_meta_data_file, regions_table, xqtl_meta_data = NULL) {
  LD_info <- fread(ld_meta_data_file, header = TRUE, data.table = FALSE)
  colnames(LD_info)[1] <- "chrom"
  LD_info$region_id <- gsub("chr", "", paste(LD_info$chrom, LD_info$start, LD_info$end, sep = "_"))
  LD_info$LD_file <- paste0(dirname(ld_meta_data_file), "/", gsub(",.*$", "", LD_info$path))
  LD_info$SNP_file <- paste0(LD_info$LD_file, ".bim")
  LD_info <- LD_info[, c("region_id", "LD_file", "SNP_file")]
  region_info <- read.table(regions_table, sep = "\t", header = TRUE) # to get exact LD bim file without over-including neighboring LD's info.
  colnames(region_info)[1] <- "chrom"
  region_info$chrom <- as.integer(gsub("chr", "", region_info$chrom))
  region_info$region_id <- paste(region_info$chrom, region_info$start, region_info$stop, sep = "_")
  if (!is.null(xqtl_meta_data)) {
    xqtl_meta_data <- fread(xqtl_meta_data, data.table = FALSE)
    region_info <- do.call(rbind, lapply(1:nrow(region_info), function(x) {
      if (!any(region_info$start[x] < xqtl_meta_data$TSS & region_info$stop[x] > xqtl_meta_data$TSS)) {
        return(NULL)
      } else {
        return(region_info[x, , drop = FALSE])
      }
    }))
  }
  return(list(LD_info = LD_info, region_info = region_info))
}

#' Function to select variants for ctwas weights input
#' @param region_data A list of list containing weights list and snp_info list data for multiple genes/events within a single LD block region.
#' @param export_twas_weight_db A list of list of fine-mapping result data formatted by generate_twas_db function.
#' @param region_block A string for region information for region_weights, consisted of chromosome number, star and end position of LD block conneced with "_".
#' @export
trim_ctwas_variants <- function(region_data, twas_weight_cutoff = 1e-5, cs_min_cor = 0.8,
                                min_pip_cutoff = 0.1, max_num_variants = 1000) {
  # internal functions to select variants for a gene-context pair weight list
  select_variants <- function(group_name, region_data, cs_min_cor, min_pip_cutoff, max_num_variants) {
    weight_list <- region_data$weights[[group_name]]
    context <- weight_list$context
    selected_variants_by_context <- c()
    molecular_id <- gsub("\\|.*", "", group_name)

    if ("cs_variants" %in% names(region_data$susie_weights_intermediate[[molecular_id]][[context]]) & length(region_data$susie_weights_intermediate[[molecular_id]][[context]][["cs_variants"]]) != 0) {
      cs_min_abs_cor <- region_data$susie_weights_intermediate[[molecular_id]][[context]]$cs_purity$min.abs.corr
      for (L in 1:length(region_data$susie_weights_intermediate[[molecular_id]][[context]]$cs_variants)) {
        # we includ all variants in $cs_variant if min_abs_corr > cs_min_cor for the set
        if (cs_min_abs_cor[L] >= cs_min_cor) {
          cs_variants <- gsub("chr", "", region_data$susie_weights_intermediate[[molecular_id]][[context]]$cs_variants[[L]])
          selected_variants_by_context <- cs_variants[cs_variants %in% rownames(weight_list$wgt)]
        }
      }
    }
    context_pip <- region_data$susie_weights_intermediate[[molecular_id]][[context]]$pip
    names(context_pip) <- gsub("chr", "", names(context_pip))
    high_pip_variants <- names(context_pip[context_pip > min_pip_cutoff])[names(context_pip[context_pip > min_pip_cutoff]) %in% rownames(weight_list$wgt)]
    selected_variants_by_context <- unique(c(selected_variants_by_context, high_pip_variants))

    if (length(selected_variants_by_context) < max_num_variants) {
      remaining_var_num <- max_num_variants - length(selected_variants_by_context)
      available_variants <- setdiff(names(context_pip)[names(context_pip) %in% rownames(weight_list$wgt)], selected_variants_by_context)
      context_pip <- context_pip[available_variants]
      selected_variants_by_context <- c(selected_variants_by_context, names(context_pip[order(-context_pip)])[1:remaining_var_num])
    }
    weight_list$wgt <- weight_list$wgt[selected_variants_by_context, , drop = FALSE]
    return(weight_list)
  }

  weights <- setNames(lapply(names(region_data$weights), function(group) {
    region_data$weights[[group]]$wgt <- region_data$weights[[group]]$wgt[abs(region_data$weights[[group]]$wgt[, 1]) > twas_weight_cutoff, , drop = FALSE]
    if (nrow(region_data$weights[[group]]$wgt) < 1) {
      return(NULL)
    }
    if (all(is.na(region_data$weights[[group]]$wgt[, 1])) || all(is.nan(region_data$weights[[group]]$wgt[, 1]))) {
      return(NULL)
    }
    if (nrow(region_data$weights[[group]]$wgt) < max_num_variants) {
      region_data$weights[[group]]$n_wgt <- nrow(region_data$weights[[group]]$wgt)
      return(region_data$weights[[group]])
    }
    region_data$weights[[group]] <- select_variants(group, region_data, cs_min_cor = cs_min_cor, min_pip_cutoff = min_pip_cutoff, max_num_variants = max_num_variants)
    region_data$weights[[group]]$n_wgt <- nrow(region_data$weights[[group]]$wgt)
    return(region_data$weights[[group]])
  }), names(region_data$weights))
  weights <- Filter(Negate(is.null), weights)
  return(weights)
}
