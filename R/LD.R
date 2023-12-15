# Function to Check if Regions are Consecutive
check_consecutive_regions <- function(df) {
    
  # Remove duplicated rows based on 'chrom' and 'start' columns
  df <- distinct(df, chrom, start, .keep_all = TRUE)

  # Arrange the dataframe by 'chrom' and 'start' columns
  df <- df %>%
    arrange(chrom, start)

  # Check if the input regions are in increasing order based on 'chr' and 'start'
  if (!all(df$chrom[-1] >= df$chrom[-nrow(df)] & df$start[-1] >= df$start[-nrow(df)])) {
    stop("The input list of regions is not in increasing order.")
  }
  return(df)
}

# Function to Find Start and End Rows for Region of Interest
find_start_end_rows <- function(df, region_chrom, region_start, region_end) {
  start_row <- df %>%
    filter(chrom == region_chrom, start <= region_start, end >= region_start) %>%
    slice(1)

  end_row <- df %>%
    filter(chrom == region_chrom, start <= region_end, end >= region_end) %>%
    arrange(desc(end)) %>%
    slice(1)

  if (nrow(start_row) == 0 || nrow(end_row) == 0) {
    stop("Region of interest is not covered by any rows in the data frame.")
  }

  list(start_row = start_row, end_row = end_row)
}

# Function to Validate Selected Region
validate_selected_region <- function(start_row, end_row, region_start, region_end) {
  merged_region_start <- start_row$start
  merged_region_end <- end_row$end

  if (!(merged_region_start <= region_start && merged_region_end >= region_end)) {
    stop("The selected region is not fully covered by the merged region.")
  }
}


# Main Function: bt_intersect
bt_intersect <- function(df, region_strings) {
  # Check if the regions in df and region_strings are in increasing order
  df = check_consecutive_regions(df)
  region_strings = check_consecutive_regions(region_strings)
  # Process each row in region_strings using lapply (if we use variants to merge, the region_strings will be a dataframe)
  result_list <- lapply(seq_len(nrow(region_strings)), function(i) {
    region_chrom <- region_strings$chrom[i]
    region_start <- region_strings$start[i]
    region_end <- region_strings$end[i]

    # Find start and end rows
    rows <- find_start_end_rows(df, region_chrom, region_start, region_end)

    # Validate the selected region
    validate_selected_region(rows$start_row, rows$end_row, region_start, region_end)

    # Get file paths between row_start and row_end
    file_paths <- df$path[which(df$chrom == rows$start_row$chrom & df$start >= rows$start_row$start & df$start <= rows$end_row$start)]

    list(row_start = which(df$chrom == rows$start_row$chrom & df$start == rows$start_row$start),
         row_end = which(df$chrom == rows$end_row$chrom & df$start == rows$end_row$start),
         file_paths = file_paths)
  })
  return(result_list)
}

#' Load and Format LD Matrix
#' 
#' @description For each LD file, produce a matrix of LD values for the SNPs in the LD file
#' Then merge all the LD matrices into one big matrix
#' @param LD_block_path A path to the LD block files
#' @return A list that contains the LD matrix and the names of the SNPs in the LD matrix
#' @return A list 
#' \describe{
#' \item{names}{The SNPs in the LD matrix}
#' \item{block}{The LD matrix}
#' }
#' @importFrom Matrix bdiag
#' @export

load_LD_matrix <- function(LD_block_path, region, variants_select = NULL) {
    
    # list LD block file names (when we use variants to merge with LD block, some variants are located within the same LD block, so we use unique function)
     LD.files.name = NULL
     LD.files <- bt_intersect(LD_block_path, region)
     for (m in seq_along(LD.files))
         {
          LD.files.name = unique(c(LD.files.name, LD.files[[m]]$file_paths))
         }
    LD.list <- list()
    LD.matrix.names <- NULL
    
    for (k in seq_along(LD.files.name)) {
    #load LD matrix
      LD.matrix <- read.table(LD.files.name[k])
      LD_names <- colnames(LD.matrix) <- rownames(LD.matrix) <- gsub("_", ":", rownames(LD.matrix))
    #extract the LD matrix if variants_select are provided
      if (!is.null(variants_select)) {
        snp_merge <- intersect(LD_names, variants_select)
        LD.select <- as.matrix(LD.matrix[snp_merge, snp_merge])
        LD.list[[k]] <- LD.select
        LD.matrix.names <- append(LD.matrix.names, snp_merge)
      } else {
        LD.list[[k]] <- as.matrix(LD.matrix)
        LD.matrix.names <- append(LD.matrix.names, LD_names)
      }
    }
     
    LD.block <- as.matrix(bdiag(LD.list))
    LD.block[upper.tri(LD.block)] <- t(LD.block)[upper.tri(LD.block)]
    colnames(LD.block) <- rownames(LD.block) <- LD.matrix.names
    return(list("names" = LD.matrix.names, "block" = LD.block))

}