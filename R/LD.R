#' Function to Check if Regions are in increasing order and remove duplicated rows
#' @importFrom dplyr arrange
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

#' Function to Find Start and End Rows for Region of Interest
#' @import dplyr
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

#' Function to Validate Selected Region
validate_selected_region <- function(start_row, end_row, region_start, region_end) {
  merged_region_start <- start_row$start
  merged_region_end <- end_row$end

  if (!(merged_region_start <= region_start && merged_region_end >= region_end)) {
    stop("The selected region is not fully covered by the merged region.")
  }
}


# Main Function: bt_intersect
#' Title
#'
#' @param df A data frame with columns "chrom", "start", "end" and "path" which contains the information of LD block.
#' "start" and "end" are the start and end position of LD block, respectively. "path" is the LD block path. 
#' @param region_strings A data frame of the region of interest with columns "chrom", "start" and "end".
#' "start" and "end" are the start and end position of the region, respectively.
#'
#' @return A list 
#' \describe{
#' \item{df}{the regions of df are in incresing order and without duplicated rows.}
#' \item{region_strings}{the regions of input region_strings are in increasing order and without duplicated rows.}
#' \item{result_list}{A list of length equal to the number of rows in `region_strings`. Each element of this list is itself a list containing:
#' \describe{
#' \item{row_start}{The index/row number in df where the intersection starts for each region.}
#' \item{row_end}{The index/row number in df where the intersection ends.}
#' \item{file_paths}{A vector of file paths from df corresponding to the intersected regions.}}}
#' }
#' @export
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
  return(list(result_list = result_list, df = df, region_strings= region_strings))
}


#' Function to load and process LD matrix
#'
#' @param LD_meta_file A data frame specifying the information of LD blocks with the olumns "chrom", "start", "end" and "path".
#' "start" and "end" are the start and end positions of LD blocks, respectively. "path" is the path of each LD block. 
#' @param region A data frame specifying the regions of interest with the columns "chrom", "start" and "end".
#' @param extract_variants A vector of variants to be extracted (optional).
#'
#' @return A list containing:
#' \describe{
#' \item{variants}{merge the variants within each LD block to a vector.}
#' \item{LD}{The linkage disequilibrium block matrix, the row and column names are identical to `variants_id_all`}.
#' }
#' @importFrom Matrix bdiag
#' @importFrom data.table fread
#' @importFrom stringr str_split
#' @import dplyr
#' 
#' @export
load_LD_matrix <- function(LD_meta_file, region, extract_variants = NULL) {
    
    # Intersect LD metadata file with the specified region
    region.LD.files <- bt_intersect(LD_meta_file, region)

    # Initialize variables
    region_merge <- list()
    variants_id_all <- list()

    # Process each region in LD.files
    region_merge <- lapply(seq_along(region.LD.files$result_list), function(i) {
        # Extract file paths for each LD block
        region_LD_file_path <- region.LD.files$result_list[[i]]$file_paths

        # Combine region strings with file paths into a data frame
        data.frame(region.LD.files$region_strings[rep(i, each = length(region_LD_file_path)), ],
                   region_LD_file_path = region_LD_file_path)
    })

    # Combine all regions into a single data frame
    region_merge <- do.call(rbind, region_merge)

    # Create an index for each unique region
    region_index <- mutate(region_merge, index = match(region_LD_file_path, unique(region_LD_file_path)))

    # Process each unique region
    LD.list <- lapply(unique(region_index$index), function(j) {
        # Filter region index for the current region
        LD.index <- filter(region_index, index == j)

        # Read the LD matrix from the file
        LD.matrix <- as.matrix(fread(cmd = paste("xzcat", unique(LD.index$region_LD_file_path)), 
                                     header = TRUE, sep = "\t")[, -1])

        # Process variant names from file paths
        LD.variants <- str_split(unique(LD.index$region_LD_file_path), "\\.", simplify = TRUE) %>%
            .[, -c(length(.), (length(.) - 1), (length(.) - 2))] %>%
            paste(., collapse = ".") %>%
            paste0(., ".bim", sep = "") %>%
            read.table(.) %>%
            mutate(V2 = gsub("_", ":", V2))

        # Set column and row names of the LD matrix
        colnames(LD.matrix) <- rownames(LD.matrix) <- LD.variants$V2

        # Select SNP indices within the specified genomic range
        LD.variants.index.selected <- lapply(seq_len(nrow(LD.index)), function(m) {
            filter(LD.variants, V4 >= LD.index$start[m], V4 <= LD.index$end[m])
        })

        # Combine SNP index metadata into a single data frame
        LD.variants.index.selected <- do.call(rbind, LD.variants.index.selected)

        # Select variants if specified and store their IDs
        if (!is.null(extract_variants)) {
            variants_selected <- intersect(LD.variants.index.selected$V2, extract_variants)
            #variants_id_all[[j]] <- variants_selected
            return(list(LD.block = LD.matrix[variants_selected, variants_selected],variants_id_all = variants_selected))
        } else {
            #variants_id_all[[j]] <- LD.variants.index.selected$V2
            return(list(LD.block = LD.matrix[LD.variants.index.selected$V2, LD.variants.index.selected$V2],variants_id_all = LD.variants.index.selected$V2))
        }
    })

    # Extract and combine LD matrices and variant IDs from the results
    LD.blocks <- lapply(LD.list, function(x) x$LD.block)
    variants_id_all <- unlist(lapply(LD.list, function(x) x$variants_id_all))

    # Create a block matrix with correct names
    LD.block <- as.matrix(bdiag(LD.blocks))
    LD.block[upper.tri(LD.block)] <- t(LD.block)[upper.tri(LD.block)]                              
    rownames(LD.block) <- colnames(LD.block) <- variants_id_all

    # Return a list containing the unique variant IDs and the LD block matrix
    return(list("variants" = variants_id_all, "LD" = LD.block))
}