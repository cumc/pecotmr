#' Function to Check if Regions are in increasing order and remove duplicated rows
#' @importFrom dplyr distinct arrange group_by mutate ungroup
#' @importFrom magrittr %>%
#' @importFrom stats lag
check_consecutive_regions <- function(df) {
  # Ensure that 'chrom' values are integers, df can be genomic_data or regions_of_interest
  df$chrom <- ifelse(grepl("^chr", df$chrom),
    as.integer(sub("^chr", "", df$chrom)), # Remove 'chr' and convert to integer
    as.integer(df$chrom)
  ) # Convert to integer if not already
  # Remove duplicated rows based on 'chrom' and 'start' columns
  df <- distinct(df, chrom, start, .keep_all = TRUE)

  # Arrange the genomic regions by 'chrom' and 'start' columns
  df <- df %>%
    arrange(chrom, start)

  # Group by chromosome and check if start positions are in ascending order
  # use stats::lag here explicitly to avoid warning message
  # due to a dplyr issue: https://github.com/tidyverse/dplyr/issues/2195
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
#' @importFrom dplyr filter arrange slice
#' @noRd
find_intersection_rows <- function(genomic_data, region_chrom, region_start, region_end) {
  # Adjusting region_start if it's smaller than the smallest start position in genomic_data
  min_start <- min(genomic_data %>% filter(chrom == region_chrom) %>% pull(start))
  if (!is.na(min_start) && region_start < min_start) {
    region_start <- min_start
  }

  # Adjusting region_end if it's larger than the largest end position in genomic_data
  max_end <- max(genomic_data %>% filter(chrom == region_chrom) %>% pull(end))
  if (!is.na(max_end) && region_end > max_end) {
    region_end <- max_end
  }

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
#' @noRd
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
#' @importFrom dplyr select
#' @importFrom data.table fread
#' @noRd
get_regional_ld_meta <- function(ld_reference_meta_file, region, complete_coverage_required = FALSE) {
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
    `colnames<-`(if (ncol(.) == 2) c("LD_file_path", "bim_file_path") else c("LD_file_path"))

  genomic_data <- cbind(genomic_data, file_path) %>% select(-path)

  # Find intersection rows
  intersection_rows <- find_intersection_rows(genomic_data, region$chrom, region$start, region$end)

  # Validate region
  if (complete_coverage_required) {
    validate_selected_region(intersection_rows$start_row, intersection_rows$end_row, region$start, region$end)
  }

  # Extract file paths
  LD_paths <- find_valid_file_paths(ld_reference_meta_file, extract_file_paths(genomic_data, intersection_rows, "LD_file_path"))
  bim_paths <- if ("bim_file_path" %in% names(genomic_data)) {
    find_valid_file_paths(ld_reference_meta_file, extract_file_paths(genomic_data, intersection_rows, "bim_file_path"))
  } else {
    NULL
  }

  return(list(
    intersections = list(
      start_index = intersection_rows$start_row,
      end_index = intersection_rows$end_row,
      LD_file_paths = LD_paths,
      bim_file_paths = bim_paths
    ),
    ld_meta_data = genomic_data,
    region = region
  ))
}

#' @importFrom dplyr mutate
#' @importFrom utils read.table
#' @importFrom stats setNames
# Process an LD matrix from a file path
process_LD_matrix <- function(LD_file_path, bim_file_path) {
  # Order the variants by position
  order_variants_by_position <- function(strings) {
    # Apply the function to each variant to get a vector of positions
    positions <- sapply(strings, function(variant) as.integer(strsplit(variant, ":")[[1]][2]))
    # Check whether the merged variants is orderd
    # when diff() returns 0, it is multiallelic at the position (same position but >1 variations)
    if (!all(diff(positions[order(positions)]) >= 0)) {
      stop("The positions are not in non-decreasing order")
    }
    # Order the variants by position
    strings_ordered <- strings[order(positions)]
    return(strings_ordered)
  }
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
  LD_variants <- read.table(bim_file_name)
  if (ncol(LD_variants) == 9) {
    LD_variants <- LD_variants %>%
      setNames(c("chrom", "variants", "GD", "pos", "A1", "A2", "variance", "allele_freq", "n_nomiss")) %>%
      mutate(chrom = ifelse(grepl("^chr[0-9]+", chrom), sub("^chr", "", chrom), chrom)) %>%
      mutate(variants = format_variant_id(variants)) %>%
      mutate(variants = ifelse(grepl("^chr[0-9]+:", variants), gsub("^chr", "", variants), variants))
  } else if (ncol(LD_variants) == 6) {
    LD_variants <- LD_variants %>%
      setNames(c("chrom", "variants", "GD", "pos", "A1", "A2")) %>%
      mutate(chrom = ifelse(grepl("^chr[0-9]+", chrom), sub("^chr", "", chrom), chrom)) %>%
      mutate(variants = format_variant_id(variants)) %>%
      mutate(variants = ifelse(grepl("^chr[0-9]+:", variants), gsub("^chr", "", variants), variants))
  } else {
    stop("Unexpected number of columns in the input file.")
  }


  # Set column and row names of the LD matrix
  colnames(LD_matrix) <- rownames(LD_matrix) <- LD_variants$variants
  # Check if the matrix is upper diagonal
  # We assume a matrix is upper diagonal if all elements below the main diagonal are zero
  is_upper_diagonal <- all(LD_matrix[lower.tri(LD_matrix)] == 0)
  if (is_upper_diagonal) {
    # If the matrix is upper diagonal, transpose the upper triangle to the lower triangle
    LD_matrix[lower.tri(LD_matrix)] <- t(LD_matrix)[lower.tri(LD_matrix)]
  } else {
    # If the matrix is lower diagonal, transpose the lower triangle to the upper triangle
    LD_matrix[upper.tri(LD_matrix)] <- t(LD_matrix)[upper.tri(LD_matrix)]
  }
  LD_variants_ordered <- LD_variants[match(order_variants_by_position(LD_variants$variants), LD_variants$variants), ]
  LD_matrix <- LD_matrix[match(LD_variants_ordered$variants, rownames(LD_matrix)), match(LD_variants_ordered$variants, rownames(LD_matrix))]
  list(LD_matrix = LD_matrix, LD_variants = LD_variants_ordered)
}

#' Extract LD matrix and variants for a specific region
#' @importFrom dplyr mutate select
#' @importFrom magrittr %>%
#' @importFrom utils tail
extract_LD_for_region <- function(LD_matrix, variants, region, extract_coordinates) {
  # Filter variants based on region
  extracted_LD_variants <- subset(variants, chrom == region$chrom & pos >= region$start & pos <= region$end)
  if (!is.null(extract_coordinates)) {
    # Preprocess 'extract_coordinate' to ensure 'chrom' is numeric and without 'chr'
    extract_coordinates <- extract_coordinates %>%
      mutate(chrom = ifelse(grepl("^chr", chrom), as.integer(sub("^chr", "", chrom)), chrom)) %>%
      select(chrom, pos)
    # Now merge with 'LD_variants_region_selected'
    extracted_LD_variants <- extracted_LD_variants %>%
      # Ensure that 'chrom' values in 'LD_variants_region_selected' are numeric, remove 'chr' if present
      mutate(chrom = ifelse(grepl("^chr", chrom), as.integer(sub("^chr", "", chrom)), as.integer(chrom))) %>%
      # Merge with 'extract_coordinate' after 'chrom' adjustment
      merge(extract_coordinates, by = c("chrom", "pos"))
    # Select the desired columns, assuming 'variants' column is equivalent to the 'variants' in 'LD_variants_region_selected'
    # Select columns dynamically based on the presence of 'variance'
    cols_to_select <- c("chrom", "variants", "pos", "GD", "A1", "A2") # select(chrom, variants, pos, GD, A1, A2)
    if ("variance" %in% names(extracted_LD_variants)) {
      cols_to_select <- c(cols_to_select, "variance")
    }
    extracted_LD_variants <- select(extracted_LD_variants, all_of(cols_to_select))
  }
  # Extract LD matrix
  extracted_LD_matrix <- LD_matrix[extracted_LD_variants$variants, extracted_LD_variants$variants, drop = FALSE]
  list(extracted_LD_matrix = extracted_LD_matrix, extracted_LD_variants = extracted_LD_variants)
}

# Create a combined LD matrix from multiple matrices
create_combined_LD_matrix <- function(LD_matrices, variants) {
  # Extract unique variant names from the list of variants
  mergeVariants <- function(LD_variants_list) {
    # Initialize an empty vector to store the merged variants
    mergedVariants <- character(0)

    # Loop over the list of LD matrices using sapply
    sapply(LD_variants_list, function(LD_variants) {
      # Extract the variants from the current LD matrix
      currentVariants <- get_nested_element(LD_variants, "variants")
      if (length(currentVariants) == 0) {
        return(NULL)
      }

      # Merge variants with the previously merged variants vector
      # Checking if the last variant is the same as the first of the current, if so, skip the first
      if (length(mergedVariants) > 0 && tail(mergedVariants, 1) == currentVariants[1]) {
        mergedVariants <<- c(mergedVariants, currentVariants[-1])
      } else {
        mergedVariants <<- c(mergedVariants, currentVariants)
      }
    })

    # Return the merged vector of variants
    return(mergedVariants)
  }
  unique_variants <- mergeVariants(variants)
  # Initialize an empty combined LD matrix with the unique variants
  combined_LD_matrix <- matrix(0, nrow = length(unique_variants), ncol = length(unique_variants))
  rownames(combined_LD_matrix) <- unique_variants
  colnames(combined_LD_matrix) <- unique_variants
  # Define a function to align the values from each LD matrix to the combined matrix
  align_matrix <- function(ld_matrix, combined_matrix, variant_names) {
    # Find the indices of the variant names in the combined matrix
    indices <- match(variant_names, rownames(combined_matrix))
    # Fill in the values for both rows and columns
    combined_matrix[indices, indices] <- ld_matrix
    return(combined_matrix)
  }
  # Apply the fill_matrix function to each LD matrix and accumulate the results
  combined_LD_matrix <- Reduce(
    function(x, y) align_matrix(y[[1]], x, y[[2]]),
    Map(list, LD_matrices, lapply(LD_matrices, rownames)),
    combined_LD_matrix
  )
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
#' @export
load_LD_matrix <- function(LD_meta_file_path, region, extract_coordinates = NULL) {
  # Intersect LD metadata with specified regions using updated function
  intersected_LD_files <- get_regional_ld_meta(LD_meta_file_path, region)

  # Extract file paths for LD and bim files
  LD_file_paths <- intersected_LD_files$intersections$LD_file_paths
  bim_file_paths <- intersected_LD_files$intersections$bim_file_paths

  # Using a for loop here to allow for rm() in each loop to save memory
  extracted_LD_matrices_list <- list()
  extracted_LD_variants_list <- list()

  # Process each LD block individually
  for (j in seq_along(LD_file_paths)) {
    LD_matrix_processed <- process_LD_matrix(LD_file_paths[j], bim_file_paths[j])
    extracted_LD_list <- extract_LD_for_region(
      LD_matrix = LD_matrix_processed$LD_matrix,
      variants = LD_matrix_processed$LD_variants,
      region = intersected_LD_files$region,
      extract_coordinates = extract_coordinates
    )
    extracted_LD_matrices_list[[j]] <- extracted_LD_list$extracted_LD_matrix
    extracted_LD_variants_list[[j]] <- extracted_LD_list$extracted_LD_variants
    # Remove large objects to free memory
    rm(LD_matrix_processed, extracted_LD_list)
  }
  combined_LD_matrix <- create_combined_LD_matrix(
    LD_matrices = extracted_LD_matrices_list,
    variants = extracted_LD_variants_list
  )
  # Remove large objects to free memory
  rm(extracted_LD_matrices_list)

  ref_panel <- do.call(rbind, lapply(strsplit(rownames(combined_LD_matrix), ":"), function(x) {
    data.frame(chrom = x[1], pos = as.integer(x[2]), A2 = x[3], A1 = x[4])
  }))
  merged_variant_list <- do.call(rbind, extracted_LD_variants_list)
  ref_panel$variant_id <- rownames(combined_LD_matrix)
  if ("variance" %in% colnames(merged_variant_list)) ref_panel$variance <- merged_variant_list$variance[match(rownames(combined_LD_matrix), merged_variant_list$variants)]

  # LD list for region
  combined_LD_list <- list(combined_LD_variants = rownames(combined_LD_matrix), combined_LD_matrix = combined_LD_matrix, ref_panel = ref_panel)

  return(combined_LD_list)
}

#' Filter variants by LD Reference
#'
#' @param variant_ids variant names in the format chr:pos_ref_alt or chr:pos:ref:alt.
#' @param ld_reference_meta_file A data frame similar to 'genomic_data' in get_regional_ld_meta function.
#' @return A subset of variants, filtered based on LD reference data.
#' @importFrom stringr str_split
#' @importFrom dplyr select group_by summarise
#' @importFrom data.table fread
#' @export
filter_variants_by_ld_reference <- function(variant_ids, ld_reference_meta_file, keep_indel = TRUE) {
  # Step 1: Process variant IDs into a data frame and filter out non-standard nucleotides
  variants_df <- do.call(rbind, lapply(strsplit(variant_ids, ":"), function(x) {
    data.frame(chrom = x[1], pos = as.integer(x[2]), ref = x[3], alt = x[4])
  }))

  variants_df$chrom <- ifelse(grepl("^chr", variants_df$chrom),
    as.integer(sub("^chr", "", variants_df$chrom)), # Remove 'chr' and convert to integer
    as.integer(variants_df$chrom)
  )
  # Step 2: Derive region information from the variants data frame
  region_df <- variants_df %>%
    group_by(chrom) %>%
    summarise(start = min(pos), end = max(pos))
  # Step 3: Call get_regional_ld_meta to get bim_file_paths
  bim_file_paths <- get_regional_ld_meta(ld_reference_meta_file, region_df)$intersections$bim_file_paths

  # Step 4: Load bim files and consolidate into a single data frame
  bim_data <- lapply(bim_file_paths, function(path) {
    bim_df <- fread(path, header = FALSE, stringsAsFactors = FALSE)
    data.frame(chrom = bim_df$V1, pos = bim_df$V4, stringsAsFactors = FALSE)
  }) %>%
    do.call("rbind", .)

  # Step 5: Overlap the variants data frame with bim_data
  keep_indices <- which(paste(variants_df$chrom, variants_df$pos) %in% paste(bim_data$chrom, bim_data$pos))
  if (!keep_indel) {
    valid_nucleotides <- c("A", "T", "C", "G")
    snp_idx <- which((variants_df$ref %in% valid_nucleotides) & (variants_df$alt %in% valid_nucleotides))
    keep_indices <- intersect(keep_indices, snp_idx)
  }
  variants_filtered <- variant_ids[keep_indices]

  message(length(variant_ids) - length(keep_indices), " out of ", length(variant_ids), " total variants dropped due to absence on the reference LD panel.")

  return(list(data = variants_filtered, idx = keep_indices))
}
