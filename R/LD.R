#' Function to Check if Regions are in increasing order and remove duplicated rows
#' @importFrom dplyr arrange
check_consecutive_regions <- function(df) {
  # Ensure that 'chrom' values are integers, df can be genomic_data or regions_of_interest
  df$chrom <- ifelse(grepl("^chr", df$chrom), 
                     as.integer(sub("^chr", "", df$chrom)), # Remove 'chr' and convert to integer
                     as.integer(df$chrom)) # Convert to integer if not already
  # Remove duplicated rows based on 'chrom' and 'start' columns
  df <- distinct(df, chrom, start, .keep_all = TRUE)

  # Arrange the genomic regions by 'chrom' and 'start' columns
  df <- df %>%
    arrange(chrom, start)
    
  # Group by chromosome and check if start positions are in ascending order
  start_ordered_check <- df %>% 
    group_by(chrom) %>% 
    mutate(start_order = start >= lag(start, default = first(start))) %>% 
    ungroup()

  if (any(!start_ordered_check$start_order, na.rm = TRUE)) {
    stop("The input list of regions is not in increasing order within each chromosome.")
  }
  return(df)
}


#' Function to Find Start and End Rows of Genomic Data for Region of Interest
#' @import dplyr
#' @noRd
find_intersection_rows <- function(genomic_data, region_chrom, region_start, region_end) {
  start_row <- genomic_data %>%
    filter(chrom == region_chrom, start <= region_start, end >= region_start) %>%
    slice(1)

  end_row <- genomic_data %>%
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
  if (!(start_row$start <= region_start && end_row$end >= region_end)) {
    stop("The selected region is not fully covered by the merged region.")
  }
}

#' Extract File Paths Based on Intersection Criteria
#' @noRd
extract_file_paths <- function(genomic_data, intersection_rows, column_to_extract) {
  # Ensure the file_path_column exists in genomic_data
  if (!column_to_extract %in% names(genomic_data)) {
    stop(paste("Column", column_to_extract, "not found in genomic data"))
  }
  # Extract rows based on intersection criteria
  extracted_paths <- genomic_data[[column_to_extract]][which(
    genomic_data$chrom == intersection_rows$start_row$chrom &
    genomic_data$start >= intersection_rows$start_row$start &
    genomic_data$start <= intersection_rows$end_row$start
  )]

  # Return the extracted paths
  return(extracted_paths)
}

#' Intersect LD reference with Regions of Interest
#'
#' @param ld_reference_meta_file A file of data frame with columns "chrom", "start", "end", and "path" representing genomic regions.
#' "chrom" is the chromosome, "start" and "end" are the positions of the LD block, and "path" is the file path for the LD block.
#' @param region A data frame with columns "chrom", "start", and "end" specifying regions of interest.
#' "start" and "end" are the positions of these regions. Or it can take the form of `chr:start-end`
#'
#' @return A list containing processed genomic data, region of interest, and a detailed result list.
#' The result list contains, for each region:
#' - The start and end row indices in the genomic data
#' - File paths from the genomic data corresponding to intersected regions
#' - Optionally, bim file paths if available
#' @importFrom stringr str_split
#' @noRd
get_regional_ld_meta <- function(ld_reference_meta_file, region) {
  genomic_data <- fread(ld_reference_meta_file, header = "auto")
  region <- parse_region(region) 
  # Set column names
  names(genomic_data) <- c("chrom", "start", "end", "path")
  names(region) <- c("chrom", "start", "end")

  # Order and deduplicate regions
  genomic_data <- check_consecutive_regions(genomic_data)
  region <- check_consecutive_regions(region)

  # Process file paths
  file_path <- genomic_data$path %>%
               str_split(",", simplify = TRUE) %>%
               data.frame() %>%
               `colnames<-`(if(ncol(.) == 2) c("LD_file_path", "bim_file_path") else c("LD_file_path"))

  genomic_data <- cbind(genomic_data, file_path) %>% select(-path)

  # Find intersection rows
  intersection_rows <- find_intersection_rows(genomic_data, region$chrom, region$start, region$end)

  # Validate region
  validate_selected_region(intersection_rows$start_row, intersection_rows$end_row, region$start, region$end)

  # Extract file paths
  LD_paths <- find_valid_file_paths(ld_reference_meta_file, extract_file_paths(genomic_data, intersection_rows, "LD_file_path"))
  bim_paths <- if ("bim_file_path" %in% names(genomic_data)) {
               find_valid_file_paths(ld_reference_meta_file, extract_file_paths(genomic_data, intersection_rows, "bim_file_path"))
             } else {
               NULL
             }

  return(list(intersections = list(start_index = intersection_rows$start_index,
                             end_index = intersection_rows$end_index,
                             LD_file_paths = LD_paths,
                             bim_file_paths = bim_paths), 
              ld_meta_data = genomic_data, 
              region = region))
}

# Process an LD matrix from a file path
process_LD_matrix <- function(LD_meta_file_path, LD_file_path, bim_file_path) {
    # Check whether the LD_file_path and bim_file_path exist or not
    LD_file_path <- find_valid_file_path(LD_meta_file_path, LD_file_path)
    bim_file_path <- find_valid_file_path(LD_meta_file_path, bim_file_path)
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

    # Process variant names from file paths
    LD_variants <- bim_file_name %>%
              read.table(.) %>%
              setNames(c("chrom","variants","GD","pos","A2","A1"))%>%
              mutate(chrom = ifelse(grepl("^chr[0-9]+", chrom), sub("^chr", "", chrom), chrom)) %>%
              mutate(variants = gsub("[_]", ":", variants)) %>%
              mutate(variants = ifelse(grepl("^chr[0-9]+:", variants), gsub("^chr", "", variants), variants))
    # Set column and row names of the LD matrix
    colnames(LD_matrix) <- rownames(LD_matrix) <- LD_variants$variants
    list(LD_matrix = LD_matrix, LD_variants = LD_variants)
}


# Extract LD matrix and variants for a specific region
extract_LD_for_region <- function(LD_matrix, variants, region, extract_coordinates) {
    # Filter variants based on region
    extracted_LD_variants <- subset(variants, chrom == region$chrom & pos >= region$start & pos <= region$end)
    if (!is.null(extract_coordinates)) {
        # Preprocess 'extract_coordinate' to ensure 'chrom' is numeric and without 'chr'
         extract_coordinates <- extract_coordinates %>%
                          mutate(chrom = ifelse(grepl("^chr", chrom), as.integer(sub("^chr", "", chrom)), chrom))%>%
                          select(chrom, pos)
        # Now merge with 'LD_variants_region_selected'
         extracted_LD_variants <- extracted_LD_variants %>%
                      # Ensure that 'chrom' values in 'LD_variants_region_selected' are numeric, remove 'chr' if present
                       mutate(chrom = ifelse(grepl("^chr", chrom), as.integer(sub("^chr", "", chrom)), as.integer(chrom))) %>%
                      # Merge with 'extract_coordinate' after 'chrom' adjustment
                       merge(extract_coordinates, by = c("chrom", "pos")) %>%
                      # Select the desired columns, assuming 'variants' column is equivalent to the 'variants' in 'LD_variants_region_selected'
                       select(chrom, variants, pos, GD, A1, A2)
        
    }
    # Extract LD matrix 
    extracted_LD_matrix = LD_matrix[extracted_LD_variants$variants, extracted_LD_variants$variants]
    list(extracted_LD_matrix = extracted_LD_matrix, extracted_LD_variants = extracted_LD_variants)
}

# Create a combined LD matrix from multiple matrices
create_combined_LD_matrix <- function(LD_matrices, variants) {
    # Create a block matrix with correct names
    combined_LD_matrix <- as.matrix(do.call(bdiag, LD_matrices))
    # Check if the matrix is upper diagonal
    # We assume a matrix is upper diagonal if all elements below the main diagonal are zero
    is_upper_diagonal <- all(combined_LD_matrix[lower.tri(combined_LD_matrix)] == 0)
    if (is_upper_diagonal) {
    # If the matrix is upper diagonal, transpose the upper triangle to the lower triangle
    combined_LD_matrix[lower.tri(combined_LD_matrix)] <- t(combined_LD_matrix)[lower.tri(combined_LD_matrix)]
    } else {
    # If the matrix is lower diagonal, transpose the lower triangle to the upper triangle
    combined_LD_matrix[upper.tri(combined_LD_matrix)] <- t(combined_LD_matrix)[upper.tri(combined_LD_matrix)]
    }                             
    rownames(combined_LD_matrix) <- colnames(combined_LD_matrix) <- variants$variants
    combined_LD_matrix
}

#' Load and Process Linkage Disequilibrium (LD) Matrix
#'
#' @param LD_meta_file_path path of LD_metadata, LD_metadata is a data frame specifying LD blocks with 
#' columns "chrom", "start", "end", and "path". "start" and "end" denote the positions of LD blocks. 
#' "path" is the path of each LD block, optionally including bim file paths.
#' @param region A data frame specifying region of interest with columns "chrom", "start", and "end".
#' @param extract_coordinates Optional data frame with columns "chrom" and "pos" for specific coordinates extraction.
#'
#' @return A list of processed LD matrices and associated variant data frames for region of interest.
#' Each element of the list contains:
#' \describe{
#' \item{combined_LD_variants}{A data frame merging selected variants within each LD block in bim file format with columns "chrom", "variants", 
#' "GD", "pos", "A1", and "A2".}
#' \item{combined_LD_matrix}{The LD matrix for each region, with row and column names matching variant identifiers.}
#' }
#' @importFrom Matrix bdiag
#' @import dplyr
#' @export
load_LD_matrix <- function(LD_meta_file_path, region, extract_coordinates = NULL) {
    # Intersect LD metadata with specified regions using updated function
    intersected_LD_files <- get_regional_ld_meta(LD_meta_file_path, region)

    # Extract file paths for LD and bim files
    LD_file_paths <- intersected_LD_files$intersections$LD_file_paths
    bim_file_paths <- intersected_LD_files$intersections$bim_file_paths

    LD_matrices_list <- lapply(seq_along(LD_file_paths), function(j) {
    # Read and process the LD matrix
    LD_matrix_processed <- process_LD_matrix(LD_meta_file_path, LD_file_paths[j], bim_file_paths[j])

    # Extract LD for the region
    region_info <- intersected_LD_files$region
    extracted_LD_list <- extract_LD_for_region(LD_matrix_processed$LD_matrix,LD_matrix_processed$LD_variants, region_info, extract_coordinates)

    list(extracted_LD_matrix = extracted_LD_list$extracted_LD_matrix, extracted_LD_variants = extracted_LD_list$extracted_LD_variants)
    })

    # Combine LD matrices and variants data frames
    combined_LD_matrices <- lapply(LD_matrices_list, function(x) x$extracted_LD_matrix)
    combined_LD_variants <- do.call(rbind, lapply(LD_matrices_list, function(x) x$extracted_LD_variants))

    # Create a combined LD matrix
    combined_LD_matrix <- create_combined_LD_matrix(combined_LD_matrices, combined_LD_variants)

    # LD list for region
    combined_LD_list <- list(combined_LD_variants = combined_LD_variants, combined_LD_matrix = combined_LD_matrix)

    return(combined_LD_list)
}

#' Filter Genotype Matrix by LD Reference
#'
#' @param X A genotype matrix with col names as variant names in the format chr:pos_ref_alt or chr:pos:ref:alt.
#' @param ld_reference_meta_file A data frame similar to 'genomic_data' in get_regional_ld_meta function.
#' @return A subset of the genotype matrix X, filtered based on LD reference data.
#' @importFrom stringr str_split
#' @importFrom dplyr select
#' @export
filter_genotype_by_ld_reference <- function(X, ld_reference_meta_file) {
  # Step 1: Process variant IDs into a data frame and filter out non-standard nucleotides
  variant_ids <- colnames(X)
  variants_df <- data.frame(
    chrom = gsub("^(chr[^:]+):.*", "\\1", variant_ids),
    pos = as.integer(gsub("^chr[^:]+:(\\d+).*", "\\1", variant_ids)),
    ref = gsub("^chr[^:]+:\\d+[:_](.)[:_].*", "\\1", variant_ids),
    alt = gsub("^chr[^:]+:\\d+[:_].[:_](.)", "\\1", variant_ids),
    stringsAsFactors = FALSE
  )
  variants_df$chrom <- ifelse(grepl("^chr", variants_df$chrom), 
                     as.integer(sub("^chr", "", variants_df$chrom)), # Remove 'chr' and convert to integer
                     as.integer(variants_df$chrom))

  valid_nucleotides <- c("A", "T", "C", "G")
  variants_df <- variants_df[variants_df$ref %in% valid_nucleotides & variants_df$alt %in% valid_nucleotides,]

  # Step 2: Derive region information from the variants data frame
  region_df <- variants_df %>%
               group_by(chrom) %>%
               summarise(start = min(pos), end = max(pos))
  # Step 3: Call get_regional_ld_meta to get bim_file_paths
  bim_file_paths <- get_regional_ld_meta(ld_reference_meta_file, region_df)$intersections$bim_file_paths

  # Step 4: Load bim files and consolidate into a single data frame
  bim_data <- lapply(bim_file_paths, function(path) {
               bim_df <- read.table(path, header = FALSE, stringsAsFactors = FALSE)
               data.frame(chrom = bim_df$V1, pos = bim_df$V4, stringsAsFactors = FALSE)
             }) %>%
             do.call("rbind", .)

  # Step 5: Overlap the variants data frame with bim_data
  keep_indices <- which(paste(variants_df$chrom, variants_df$pos) %in% paste(bim_data$chrom, bim_data$pos))
  X_filtered <- X[,keep_indices,drop=F]

  message("Number of variants dropped: ", ncol(X) - length(keep_indices), 
          " out of ", ncol(X), " total columns.")

  return(X_filtered)
}