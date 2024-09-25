#' Utility function to load LD in ctwas analyses, to interface with cTWAS package
#' @param ld_matrix_file_path A string of file path to the LD matrix.
ctwas_ld_loader <- function(ld_matrix_file_path) {
  ld_loaded <- process_LD_matrix(ld_matrix_file_path, paste0(ld_matrix_file_path, ".bim"))
  ld_loaded <- ld_loaded$LD_matrix
  return(ld_loaded)
}
#' @importFrom data.table fread
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
get_ctwas_meta_data <- function(ld_meta_data_file, regions_table) {
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
  return(list(LD_info = LD_info, region_info = region_info))
}
