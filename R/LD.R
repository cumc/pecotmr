#' Function to Check if Regions are in increasing order and remove duplicated rows
#' @importFrom dplyr arrange
check_consecutive_regions <- function(df) {
    
  # Ensure that 'chrom' values are integers
  df$chrom <- ifelse(grepl("^chr", df$chrom), 
                     as.integer(sub("^chr", "", df$chrom)), # Remove 'chr' and convert to integer
                     as.integer(df$chrom)) # Convert to integer if not already
  # Remove duplicated rows based on 'chrom' and 'start' columns
  df <- distinct(df, chrom, start, .keep_all = TRUE)

  # Arrange the dataframe by 'chrom' and 'start' columns
  df <- df %>%
    arrange(chrom, start)

  # Check if the input regions are in increasing order based on 'chr' and 'start'
  if (!all(df$start[-1] >= df$start[-nrow(df)])) {
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
#' \item{file_paths}{A vector of file paths from df corresponding to the intersected regions.}
#' \item{bim_file_paths}{A vector of bim.file paths from df corresponding to the intersected 
#' regions.}}}
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
    # check if there exits bim_path, if the path is present, we use directly, else is NULL
    #if(!is.null(df$bim_path))
    if(any(names(df) == 'bim_path'))
    {
       bim_file_paths <- df$bim_path[which(df$chrom == rows$start_row$chrom & df$start >= rows$start_row$start & df$start <= rows$end_row$start)]

    }else{
       bim_file_paths <-NULL
    }
    list(row_start = which(df$chrom == rows$start_row$chrom & df$start == rows$start_row$start),
         row_end = which(df$chrom == rows$end_row$chrom & df$start == rows$end_row$start),
         file_paths = file_paths,
         bim_file_paths = bim_file_paths)
  })
  return(list(result_list = result_list, df = df, region_strings= region_strings))
}


#' Function to load and process LD matrix
#'
#' @param LD_meta_file A data frame specifying the information of LD blocks with the olumns "chrom",
#' "start", "end" and "path", "bim_path"(optional). "start" and "end" are the start and end positions #' of LD blocks, respectively. "path" is the path of each LD block. "bim_path" is the bim.file path of #' each LD matrix.
#' @param region A data frame specifying the regions of interest with the columns "chrom", "start" and #' "end".
#' @param extract_coordinate A data frame with the columns "chrom","pos" to be extracted (optional).
#'
#' @return A list containing:
#' \describe{
#' \item{variants_df}{merge the variants_selected_df within each LD block to a data frame in the 
#' format of bim file with the columns "chrom", "variants", 
#' "GD", "pos", "A1" and "A2".}
#' \item{LD}{The linkage disequilibrium block matrix, the row and column names are identical to `variants_id_all`}.
#' }
#' @importFrom Matrix bdiag
#' @importFrom data.table fread
#' @importFrom stringr str_split
#' @import dplyr
#' 
#' @export
# Function to load and process LD matrix
load_LD_matrix <- function(LD_meta_file, region, extract_coordinate = NULL) {
    # Intersect LD metadata file with the specified regions
    region_LD_files <- bt_intersect(LD_meta_file, region)

    # Process each region in LD.files
    LD_all_list <- lapply(seq_along(region_LD_files$result_list), function(i) {
        # Extract file paths for each LD block
        region_LD_file_paths <- region_LD_files$result_list[[i]]$file_paths
        region_bim_file_paths <- region_LD_files$result_list[[i]]$bim_file_paths

        # Process each region
        LD_list <- lapply(seq_along(region_LD_file_paths), function(j) {

           # Read the LD matrix from the file
           LD_matrix <- as.matrix(read_delim(region_LD_file_paths[j], col_names = FALSE, delim = " "))
 
           # Determine the .bim file to use
            bim_file_name <- if (!is.null(region_bim_file_paths[j])) {
                region_bim_file_paths[j]
            } else {
                # Construct .bim filename if bim_file_paths are not available
                sub("\\D*$", "", region_LD_file_paths[j]) %>%
                paste(., collapse = ".") %>%
                paste0(., ".bim")
            }
           # Process variant names from file paths
           LD_variants <- bim_file_name%>%
              read.table(.) %>%
              setNames(c("chrom","variants","GD","pos","A2","A1"))%>%
              mutate(chrom = ifelse(grepl("^chr[0-9]+", chrom), sub("^chr", "", chrom), chrom)) %>%
              mutate(variants = gsub("[_]", ":", variants)) %>%
              mutate(variants = ifelse(grepl("^chr[0-9]+:", variants), gsub("^chr", "", variants), variants))

            # Set column and row names of the LD matrix
           colnames(LD_matrix) <- rownames(LD_matrix) <- LD_variants$variants

            # Extract variants within the specified genomic range
            region_data <- region_LD_files$region_strings[i,]
            LD_variants_region_selected <- filter(LD_variants, pos >= region_data$start, pos <= region_data$end)
    

            # Extract variants within each LD matrix if extract_coordinate is specified and store their variants table in the format of bim file
            if (!is.null(extract_coordinate)) {
               # Preprocess 'extract_coordinate' to ensure 'chrom' is numeric and without 'chr'
               extract_coordinate <- extract_coordinate %>%
                          mutate(chrom = ifelse(grepl("^chr", chrom), as.integer(sub("^chr", "", chrom)), chrom))

               # Now merge with 'LD_variants_region_selected'
               variants_selected_df <- LD_variants_region_selected %>%
                      # Ensure that 'chrom' values in 'LD_variants_region_selected' are numeric, remove 'chr' if present
                       mutate(chrom = ifelse(grepl("^chr", chrom), as.integer(sub("^chr", "", chrom)), as.integer(chrom))) %>%
                      # Merge with 'extract_coordinate' after 'chrom' adjustment
                       merge(extract_coordinate, by = c("chrom", "pos")) %>%
                      # Select the desired columns, assuming 'variants' column is equivalent to the 'variants' in 'LD_variants_region_selected'
                       select(chrom, variants, pos, GD, A1, A2)
              LD_matrix_list = LD_matrix[variants_selected_df$variants, variants_selected_df$variants]
            } else {
               LD_matrix_list = LD_matrix[LD_variants_region_selected$variants, LD_variants_region_selected$variants]
               variants_selected_df = LD_variants_region_selected
            }

           return(list(LD_matrix_list = LD_matrix_list,variants_selected_df = variants_selected_df))

       })

           # combine LD matrices and variant_selected_df from the LD_list
           LD_matrices <- lapply(LD_list, function(x) x$LD_matrix_list)
           variants_all_df <- do.call(rbind, lapply(LD_list, function(x) x$variants_selected_df))
           
           # Create a block matrix with correct names
           LD_block <- as.matrix(bdiag(LD_matrices))

           # Check if the matrix is upper diagonal
           # We assume a matrix is upper diagonal if all elements below the main diagonal are zero
           is_upper_diagonal <- all(LD_block[lower.tri(LD_block)] == 0)

           if (is_upper_diagonal) {
           # If the matrix is upper diagonal, transpose the upper triangle to the lower triangle
           LD_block[lower.tri(LD_block)] <- t(LD_block)[lower.tri(LD_block)]
           } else {
           # If the matrix is lower diagonal, transpose the lower triangle to the upper triangle
           LD_block[upper.tri(LD_block)] <- t(LD_block)[upper.tri(LD_block)]
           }                             
           rownames(LD_block) <- colnames(LD_block) <- variants_all_df$variants
           # Return a list containing the variants_all_df and the LD block matrix
           return(list("variants_df" = variants_all_df, "LD" = LD_block))
    })
    return(LD_all_list = LD_all_list)
}