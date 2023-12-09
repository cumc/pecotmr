#' Load and Format LD Matrix
#' 
#' @description For each LD LD file, produce a matrix of LD values for the SNPs in the LD file
#' Then merge all the LD matrices into one big matrix
#' @param LD_block_path A path to the LD block files
#' @return A list that contains the LD matrix and the names of the SNPs in the LD matrix
#' @return A list 
#' \describe{
#' \item{names}{The SNPs in the LD matrix}
#' \item{block}{The LD matrix}
#' }
#' @importFrom bedtoolsr bt.intersect
#' @importFrom Matrix bdiag
#' @export

load_LD_matrix <- function(LD_block_path, region) {

    # Load LD matrix
    LD.files <- bt.intersect(a = LD_block_path, b = region)
    LD.files.name <- unique(LD.files$V5)
    LD.list <- list()
    LD.matrix.names <- NULL
    
    for (k in seq_along(LD.files.name)) {
      # Create a temporary file to hold the decompressed data
      temp_file <- tempfile()
      # Decompress the .xz file
      # Requires xz to be installed on system
      system(paste("xzdec", LD.files.name[k], ">", temp_file))
      # Now read the decompressed data
      xz <- readLines(temp_file)
      unlink(temp_file)
      LD.matrix <- xz$f[["arr_0"]]
      #npz <- np$load(LD.files.name[k])
      #LD.matrix <- npz$f[["arr_0"]]
      LD.snps <- str_split(LD.files.name[k], "[.]", simplify = TRUE) %>% 
        .[,-c(length(.), (length(.) - 1))] %>% 
        paste(., collapse = ".") %>% 
        paste0(., ".bim", sep = "") %>% 
        read.table(.)
      
      LD_names <- colnames(LD.matrix) <- rownames(LD.matrix) <- gsub("_", ":", LD.snps$V2)
      snp_merge <- intersect(LD_names, outcome_QC$variant_allele_flip)
      LD.select <- as.matrix(LD.matrix[snp_merge, snp_merge])
      LD.list[[k]] <- LD.select
      LD.matrix.names <- append(LD.matrix.names, snp_merge)
    }
    
    LD.block <- as.matrix(bdiag(LD.list))
    LD.block[upper.tri(LD.block)] <- t(LD.block)[upper.tri(LD.block)]
    #upperTriangle(LD.block, byrow = TRUE) <- lowerTriangle(LD.block)
    colnames(LD.block) <- rownames(LD.block) <- LD.matrix.names
    return list("names" = LD.matrix.names, "block" = LD.block)

}