#' Utility function to load LD in ctwas analyses, to interface with cTWAS package
#' @param ld_matrix_file_path A string of file path to the LD matrix.
ctwas_ld_loader <- function(ld_matrix_file_path) {
  ld_loaded <- process_LD_matrix(ld_matrix_file_path, paste0(ld_matrix_file_path, ".bim"))
  ld_loaded <- ld_loaded$LD_matrix
  return(ld_loaded)
}
#' @importFrom data.table fread
ctwas_bimfile_loader <- function(ld_matrix_file_path) {
  snp_info <- fread(paste0(ld_matrix_file_path, ".bim"))
  if (ncol(snp_info) == 9) {
    colnames(snp_info) <- c("chrom", "id", "GD", "pos", "alt", "ref", "variance", "allele_freq", "n_nomiss")
  } else {
    colnames(snp_info) <- c("chrom", "id", "GD", "pos", "alt", "ref")
  }
  return(snp_info)
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

#' Function to extract information from harmonize_twas output as cTWAS input. For each imputable context within a gene,
#' we select strong variants based on cs set and pip value.
#' required ctwas package: remotes::install_github("xinhe-lab/ctwas",ref = "multigroup")
get_ctwas_input <- function(post_qc_data, LD_meta_file_path, twas_table, region_block) {
  get_ctwas_weights <- function(post_qc_twas_data) {
    chrom <- unique(find_data(post_qc_twas_data, c(2, "chrom")))
    if (length(chrom) != 1) stop("Data provided contains more than one chromosome. ")
    genes <- names(post_qc_twas_data)
    # Compile weights_list
    weights_list <- list()
    for (gene in genes) {
      contexts <- names(post_qc_twas_data[[gene]][["weights_qced"]])
      for (context in contexts) {
        data_type <- post_qc_twas_data[[gene]][["data_type"]][[context]]
        postqc_scaled_weight <- post_qc_twas_data[[gene]][["weights_qced"]][[context]][["scaled_weights"]]
        colnames(postqc_scaled_weight) <- "weight"
        context_variants <- rownames(post_qc_twas_data[[gene]][["weights_qced"]][[context]][["scaled_weights"]])
        context_range <- sapply(context_variants, function(variant) as.integer(strsplit(variant, "\\:")[[1]][2]))
        weights_list[[paste0(gene, "|", data_type, "_", context)]] <- list(
          chrom = chrom,
          p0 = min(context_range),
          p1 = max(context_range),
          wgt = postqc_scaled_weight,
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
  # format ctwas weights
  weights <- get_ctwas_weights(post_qc_data) # reshape weights for all gene-context pairs in the region for cTWAS analysis
  weights <- weights[!sapply(weights, is.null)]
  # gene_z table
  twas_table <- twas_table[na.omit(twas_table$is_selected_method), , drop = FALSE]
  twas_table$id <- paste0(twas_table$gene, "|", twas_table$type, "_", twas_table$context)
  twas_table$z <- twas_table$twas_z
  twas_table$group <- paste0(twas_table$context, "|", twas_table$type)
  twas_table <- twas_table[, c("id", "z", "type", "context", "group", "gwas_study"), drop = FALSE]
  studies <- unique(twas_table$gwas_study)
  z_gene_list <- list()
  z_snp <- list()
  for (study in studies) {
    z_gene_list[[study]] <- twas_table[twas_table$gwas_study == study, , drop = FALSE]
    z_snp[[study]] <- do.call(rbind, lapply(post_qc_data, function(x) find_data(x, c(1, "gwas_qced", study), docall = rbind)))
    colnames(z_snp[[study]])[which(colnames(z_snp[[study]]) == "variant_id")] <- "id"
    z_snp[[study]] <- z_snp[[study]][, c("id", "A1", "A2", "z")]
    z_snp[[study]] <- z_snp[[study]][!duplicated(z_snp[[study]]$id), ]
  }
  # load LD variant names
  region_of_interest <- region_to_df(region_block)
  bim_file_paths <- unique(do.call(c, lapply(1:nrow(region_of_interest), function(region_row) {
    get_regional_ld_meta(LD_meta_file_path, region_of_interest[region_row, , drop = FALSE])$intersections$bim_file_paths
  })))
  snp_info <- lapply(bim_file_paths, function(file) {
    bimfile <- read.table(file, header = FALSE, sep = "\t")[, c(1, 2, 4:8)]
    bimfile$V2 <- gsub("chr", "", gsub("_", ":", bimfile$V2))
    colnames(bimfile) <- c("chrom", "id", "pos", "alt", "ref", "variance", "allele_freq") # A1:alt, A2: ref
    return(bimfile)
  })
  names(snp_info) <- do.call(c, lapply(bim_file_paths, function(x) {
    parts <- strsplit(basename(x), "[_:/.]")[[1]][1:3]
    gsub("chr", "", paste(parts, collapse = "_"))
  }))
  return(list(weights = weights, snp_info = snp_info, z_gene = z_gene_list, z_snp = z_snp))
}
