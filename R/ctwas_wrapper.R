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
