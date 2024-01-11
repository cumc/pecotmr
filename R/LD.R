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

#' Intersect Genomic Data with Regions of Interest
#'
#' @param genomic_data A data frame with columns "chrom", "start", "end", and "path" representing genomic regions.
#' "chrom" is the chromosome, "start" and "end" are the positions of the LD block, and "path" is the file path for the LD block.
#' @param regions_of_interest A data frame with columns "chrom", "start", and "end" specifying regions of interest.
#' "start" and "end" are the positions of these regions.
#'
#' @return A list containing processed genomic data, regions of interest, and a detailed result list.
#' The result list contains, for each region:
#' - The start and end row indices in the genomic data
#' - File paths from the genomic data corresponding to intersected regions
#' - Optionally, bim file paths if available
#' @importFrom stringr str_split
#' @export
intersect_genomic_regions <- function(genomic_data, regions_of_interest) {

  # Set column names
  names(genomic_data) <- c("chrom", "start", "end", "path")
  names(regions_of_interest) <- c("chrom", "start", "end")

  # Order and deduplicate regions
  genomic_data <- check_consecutive_regions(genomic_data)
  regions_of_interest <- check_consecutive_regions(regions_of_interest)

  # Process file paths
  file_path_data <- genomic_data$path %>%
                    str_split(",", simplify = TRUE) %>%
                    data.frame() %>%
                    `colnames<-`(if(ncol(.) == 2) c("LD_file_path", "bim_file_path") else c("LD_file_path"))

  genomic_data <- cbind(genomic_data, file_path_data) %>% select(-path)

  # Process each region
  intersection_results <- lapply(seq_len(nrow(regions_of_interest)), function(i) {
    region <- regions_of_interest[i, ]

    # Find intersection rows
    intersection_rows <- find_intersection_rows(genomic_data, region$chrom, region$start, region$end)

    # Validate region
    validate_region(intersection_rows$start_row, intersection_rows$end_row, region$start, region$end)

    # Extract file paths
    LD_paths <- extract_file_paths(genomic_data, intersection_rows, "LD_file_path")
    bim_paths <- if ("bim_file_path" %in% names(genomic_data)) {
                 extract_file_paths(genomic_data, intersection_rows, "bim_file_path")
               } else {
                 NULL
               }

    list(start_index = intersection_rows$start_index,
         end_index = intersection_rows$end_index,
         LD_file_paths = LD_paths,
         bim_file_paths = bim_paths)
  })

  return(list(intersections = intersection_results, genomic_data = genomic_data, regions_of_interest = regions_of_interest))
}

# Process an LD matrix from a file path
process_LD_matrix <- function(LD_file_path, bim_file_path) {
    # Read the LD matrix
    LD_file_con <- xzfile(LD_file_path)
    LD_matrix <- scan(LD_file_con, quiet = TRUE)
    close(LD_file_con)
    LD_matrix <- matrix(LD_matrix, ncol = sqrt(length(LD_matrix)), byrow = TRUE)

    # Process the bim file to extract variant information
    bim_file_name <- if (!is.null(bim_file_path)) {
        bim_file_path
    } else {
        paste0(LD_file_path, ".bim", sep = "")
    }
    variants_info <- read.table(bim_file_name, header = FALSE, stringsAsFactors = FALSE)
    names(variants_info) <- c("chrom", "variants", "GD", "pos", "A1", "A2")

    list(LD_matrix = LD_matrix, variants = variants_info)
}

# Extract variants for a specific region
extract_variants_for_region <- function(variants_info, region, extract_coordinates) {
    # Filter variants based on region
    region_variants <- subset(variants_info, chrom == region$chrom & pos >= region$start & pos <= region$end)

    if (!is.null(extract_coordinates)) {
        # Filter based on additional coordinates
        extract_coordinates$chrom <- as.integer(grepl("^chr", extract_coordinates$chrom))
        region_variants <- merge(region_variants, extract_coordinates, by = c("chrom", "pos"))
    }

    region_variants
}

# Subset an LD matrix based on selected variants
subset_LD_matrix <- function(LD_matrix, variants) {
    variant_positions <- match(variants$variants, colnames(LD_matrix))
    LD_matrix[variant_positions, variant_positions]
}

# Create a combined LD matrix from multiple matrices
create_combined_LD_matrix <- function(LD_matrices, variants) {
    combined_LD_matrix <- do.call(bdiag, LD_matrices)
    rownames(combined_LD_matrix) <- colnames(combined_LD_matrix) <- variants$variants
    combined_LD_matrix
}

#' Load and Process Linkage Disequilibrium (LD) Matrix
#'
#' @param LD_metadata A data frame specifying LD blocks with columns "chrom", "start", "end", and "path".
#' "start" and "end" denote the positions of LD blocks. "path" is the path of each LD block, optionally including bim file paths.
#' @param regions A data frame specifying regions of interest with columns "chrom", "start", and "end".
#' @param extract_coordinates Optional data frame with columns "chrom" and "pos" for specific coordinates extraction.
#'
#' @return A list of processed LD matrices and associated variant data frames for each region.
#' Each element of the list contains:
#' \describe{
#' \item{variants_df}{A data frame merging selected variants within each LD block in bim file format with columns "chrom", "variants", 
#' "GD", "pos", "A1", and "A2".}
#' \item{LD_matrix}{The LD matrix for each region, with row and column names matching variant identifiers.}
#' }
#' @importFrom Matrix bdiag
#' @import dplyr
#' @export
load_LD_matrices <- function(LD_metadata, regions, extract_coordinates = NULL) {
    # Intersect LD metadata with specified regions using updated function
    intersected_LD_files <- intersect_genomic_regions(LD_metadata, regions)

    # Process each region
    LD_matrices_processed <- lapply(seq_along(intersected_LD_files$intersections), function(i) {
        # Extract file paths for LD and bim files
        LD_file_paths <- intersected_LD_files$intersections[[i]]$LD_file_paths
        bim_file_paths <- intersected_LD_files$intersections[[i]]$bim_file_paths

        # Process each file path
        LD_matrices_list <- lapply(seq_along(LD_file_paths), function(j) {
            # Read and process the LD matrix
            LD_matrix_processed <- process_LD_matrix(LD_file_paths[j], bim_file_paths[j])

            # Extract variants for the region
            region_info <- intersected_LD_files$regions_of_interest[i,]
            selected_variants_df <- extract_variants_for_region(LD_matrix_processed$variants, region_info, extract_coordinates)

            # Subset LD matrix with selected variants
            LD_subset <- subset_LD_matrix(LD_matrix_processed$LD_matrix, selected_variants_df)

            list(LD_matrix_subset = LD_subset, variants_df = selected_variants_df)
        })

        # Combine LD matrices and variants data frames
        combined_LD_matrices <- lapply(LD_matrices_list, function(x) x$LD_matrix_subset)
        combined_variants_dfs <- do.call(rbind, lapply(LD_matrices_list, function(x) x$variants_df))

        # Create a combined LD matrix
        combined_LD_matrix <- create_combined_LD_matrix(combined_LD_matrices, combined_variants_dfs)

        list(variants_df = combined_variants_dfs, LD_matrix = combined_LD_matrix)
    })

    return(LD_matrices_processed)
}
