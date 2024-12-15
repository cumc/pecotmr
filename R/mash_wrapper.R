#' @export
filter_invalid_summary_stat <- function(dat_list, z = TRUE, sig_p_cutoff = 1E-6, filter_by_missing_rate = 0.2) {
  replace_values <- function(df, replace_with) {
    df %>%
      mutate(across(everything(), as.numeric)) %>%
      mutate(across(everything(), ~ replace(., is.nan(.) | is.infinite(.) | is.na(.), replace_with)))
  }

  # Function to process z-scores
  process_z <- function(z_data) {
    z_data <- as.matrix(replace_values(z_data, 0))
    
    if (!is.null(filter_by_missing_rate)) {
      proportion_nonzero <- apply(z_data, 1, function(row) mean(row != 0))
      z_data <- z_data[proportion_nonzero >= filter_by_missing_rate, , drop = FALSE]
    }
    
    return(z_data)
  }

  # Process each component if it exists
  if (!is.null(dat_list$strong) && !is.null(dat_list$strong$z)) {
    dat_list$strong$z <- process_z(dat_list$strong$z)
  }
  
  if (!is.null(dat_list$random) && !is.null(dat_list$random$z)) {
    dat_list$random$z <- process_z(dat_list$random$z)
  }
  
  if (!is.null(dat_list$null) && !is.null(dat_list$null$z)) {
    dat_list$null$z <- process_z(dat_list$null$z)
  }

  # Apply significance cutoff to strong signals if applicable
  if (!is.null(dat_list$strong) && !is.null(dat_list$strong$z) && !is.null(sig_p_cutoff)) {
    chi_square_stat <- qchisq(sig_p_cutoff, df = 1, lower.tail = FALSE)
    z_score <- sqrt(chi_square_stat)
    keep_index <- which(apply(dat_list$strong$z, 1, function(row) any(abs(row) >= z_score)))
    dat_list$strong$z <- dat_list$strong$z[keep_index, , drop = FALSE]
  }

  return(dat_list)
}

#' @export
filter_mixture_components <- function(conditions_to_keep, U, w = NULL, w_cutoff = 1e-04) {
  # Identify conditions not to keep (to be removed)
  conditions_to_filter <- setdiff(colnames(U[[1]]), conditions_to_keep)
  sum_w <- sum(w) # Original total sum of weights

  # Filter U by removing unwanted phenotypes (conditions)
  U <- lapply(U, function(mat, to_keep) {
    missing_conditions <- setdiff(to_keep, colnames(mat))
    if (length(missing_conditions) > 0) {
      stop(paste("Condition(s)", paste(missing_conditions,
        collapse = ", "
      ), "not found in matrix"))
    }
    mat[to_keep, to_keep] # Keep only relevant conditions
  }, conditions_to_keep)

  # Remove matrices where all values are zero or weight is below cutoff
  for (mat_name in names(U)) {
    if (all(U[[mat_name]] == 0)) {
      U[[mat_name]] <- NULL
      if (!is.null(w)) {
        w <- w[!names(w) %in% mat_name]
      }
      next
    }
    if (!is.null(w)) {
      if (w[mat_name] < w_cutoff) {
        w <- w[!names(w) %in% mat_name]
        U[[mat_name]] <- NULL
      }
    }
  }

  # Note: Matrices in U may contain very small values on the diagonal
  # even when contexts are not present, due to EM algorithm adjustments.
  # This makes the matrix not exactly zero, so it won't be removed even though it may not
  # have strong context relevance. This behavior arises because the algorithm attempts to
  # ensure matrices are full-rank, slightly changing initial values.

  # We cannot simply remove diagonal matrices as signals on the diagonal can be strong and relevant.
  # So we manually remove the U components that are driven by non-relevant contexts.
  U[conditions_to_filter] <- NULL
  w <- w[!names(w) %in% conditions_to_filter]

  # Recalculate the sum of remaining weights
  sum_w_new <- sum(w)

  # Adjust weights to maintain the original sum_w
  w <- (w / sum_w_new) * sum_w

  message(paste(length(U), "components of matrices remained after filtering."))

  return(list(U = U, w = w))
}
# This function extracts tensorQTL results for given region for multiple
# summary statistics files
#' @export
load_multitrait_tensorqtl_sumstat <- function(
    sumstats_paths, region, gene = NULL,
    trait_names = NULL, top_loci = FALSE, filter_file = NULL, remove_any_missing = FALSE,
    max_rows_selected = 300, nan_remove = TRUE) {
  if (!is.vector(sumstats_paths) || !all(file.exists(sumstats_paths))) {
    stop("sumstats_paths must be a vector of existing file paths.")
  }
  if (!is.character(region) || length(region) != 1) {
    stop("region must be a single character string.")
  }
  if (!is.character(trait_names)) {
    stop("trait_names must be a vector of character strings.")
  }


  extract_tensorqtl_data <- function(path, region) {
    tabix_region(path, region) %>%
      # first four columns are 'chrom','pos','alt','ref'
      mutate(variants = paste(.[[1]], .[[2]], .[[3]], .[[4]], sep = ":")) %>%
      distinct(variants, molecular_trait_id, .keep_all = TRUE)
  }


  merge_matrices <- function(matrix_list, value_column, id_column = "variants",
                             remove_any_missing = FALSE) {
    # Convert matrices to data frames
    df_list <- lapply(seq_along(matrix_list), function(i) {
      df <- as.data.frame(matrix_list[[i]])
      df2 <- df[, c(id_column, value_column)]
      # Rename columns to avoid duplication
      colnames(df2) <- c(id_column, paste0(value_column, "_", i))
      return(df2)
    })

    # Iteratively merge the data frames
    merged_df <- Reduce(
      function(x, y) merge(x, y, by = id_column, all = TRUE),
      df_list
    )

    # Optionally, remove rows with any missing values
    if (remove_any_missing) {
      merged_df <- merged_df[complete.cases(merged_df), ]
    }
    return(merged_df)
  }

  split_variants_and_match <- function(variant, filter_file, max_rows_selected) {
    if (!file.exists(filter_file)) {
      stop("Filter file does not exist.")
    }

    # Split the variant vector into components
    variant_split <- strsplit(variant, ":")
    variant_df <- data.frame(chr = sapply(variant_split, `[`, 1), pos = sapply(
      variant_split,
      `[`, 2
    ), stringsAsFactors = FALSE)
    variant_df$pos <- as.numeric(variant_df$pos)

    # get the region of interest
    min_pos <- min(variant_df$pos)
    max_pos <- max(variant_df$pos)
    chrom <- unique(variant_df$chr)
    if (length(chrom) != 1) {
      stop("Variants are from multiple chromosomes. Cannot create a single range string.")
    }
    region <- paste0(chrom, ":", min_pos, "-", max_pos)
    ref_table <- tabix_region(filter_file, region)
    if (is.null(ref_table)) {
      stop("No variants in the region.")
    }
    colnames(ref_table)[1:2] <- c("#CHROM", "POS")
    if (!all(c("#CHROM", "POS") %in% colnames(ref_table))) {
      stop("Filter file must contain columns: #CHROM, POS.")
    }
    matched_indices <- which(variant_df$chr %in% ref_table$`#CHROM` & variant_df$pos %in%
      ref_table$POS)
    if (!is.null(max_rows_selected) && max_rows_selected > 0 && max_rows_selected <
      length(matched_indices)) {
      selected_rows <- sample(length(matched_indices), max_rows_selected)
      matched_indices <- matched_indices[selected_rows]
    }
    return(matched_indices)
  }

  Y <- lapply(sumstats_paths, function(x, region, gene) {
    out <- extract_tensorqtl_data(x, region)
    if (!is.null(gene)) {
      out <- out[which(out$molecular_trait_id %in% gene), ]
    }
    sorted <- out[order(-abs(out$beta / out$se)), c("variants", "molecular_trait_id")]
    top_v <- apply(sorted[1:2, ], 1, function(row) paste(row, collapse = "_")) # paste the variant (chr:pos:alt:ref) with gene_id with '_'
    out <- as.list(out)
    out$top_variants <- top_v
    return(out)
  }, region = region, gene = gene)

  ## Y is list of data frames where colnames(Y[[1]])=
  ## c('chrom','pos','alt','ref','variant_id','molecular_trait_id','start_distance',
  ## 'end_distance','af','ma_samples','ma_count','pvalue','beta','se','molecular_trait_object_id',
  ## 'n','variant''top_variants')

  ### The step below assigns condition names to the Y; in case the filename
  ### itself does not contain any condition names, users can input the
  ### condition names via assigning trait_names if trait_names left blank, then
  ### the file names will be assigned as trait_names, where if the text before
  ### '.' (extension) can differentiate the conditions, we will use the shorter
  ### names as trait_name
  if (is.null(trait_names)) {
    trait_names <- gsub("\\..*", "", basename(sumstats_paths)) # extract condition name that is listed before the first appearance of '.'

    if (length(trait_names[duplicated(trait_names)]) >= 1) {
      trait_names <- basename(sumstats_paths)
    }
  }
  names(Y) <- trait_names

  bhat <- merge_matrices(Y,
    value_column = "beta", id_column = c("variants", "molecular_trait_id"),
    remove_any_missing
  )
  sbhat <- merge_matrices(Y,
    value_column = "se", id_column = c("variants", "molecular_trait_id"),
    remove_any_missing
  )
  out <- list(bhat = bhat, sbhat = sbhat)

  # Check if variants are the same in both bhat and sbhat
  if (!identical(out$bhat$variants, out$sbhat$variants)) {
    stop("Error: Variants in bhat and sbhat are not the same.")
  }

  var_idx <- 1:nrow(out$bhat)

  # match with filter_file
  if (!is.null(filter_file)) {
    if (!file.exists(filter_file)) {
      stop("Filter file does not exist.")
    }
    variants <- paste0(out$bhat$variants, "_", out$bhat$molecular_trait_id)
    var_idx <- split_variants_and_match(out$bhat$variants, filter_file, max_rows_selected)
  }

  if (top_loci) {
    union_top_loci <- unique(unlist(lapply(Y, function(item) item$top_variants)))
    var_idx <- which(variants %in% union_top_loci) # var_idx may end up empty if max_rows_selected number too small
  }

  # Extract only subset of data
  variants <- paste0(out$bhat$variants[var_idx], "_", out$bhat$molecular_trait_id[var_idx])
  out$bhat <- out$bhat[var_idx, ]
  out$sbhat <- out$sbhat[var_idx, ]

  if (nan_remove) {
    out <- filter_invalid_summary_stat(out, bhat = "beta", sbhat = "se", nan_remove)
  }

  rownames(out$bhat) <- rownames(out$sbhat) <- variants
  colnames(out$bhat)[which(startsWith(colnames(out$bhat), "beta"))] <- colnames(out$sbhat)[which(startsWith(
    colnames(out$sbhat),
    "se"
  ))] <- trait_names
  out$region <- region
  out$top_variants <- lapply(Y, function(x) x$top_variants)
  return(out)
}

merge_susie_cs <- function(susie_fit, coverage = "cs_coverage_0.95", complementary = FALSE) {
  # Initialize an empty list for the combined_sets
  combined_sets <- list()

  # Identify variant IDs that are associated with more than one credible set
  identify_overlap_sets <- function(variants_sets_and_pips_list) {
    overlap_sets <- list()
    for (variant_id in names(variants_sets_and_pips_list)) {
      sets <- variants_sets_and_pips_list[[variant_id]][["sets"]]
      if (length(sets) > 1) {
        overlap_sets[[variant_id]] <- sets
      }
    }
    return(overlap_sets)
  }
  # Merge overlapping credible sets and update the credible sets in
  # variants_list
  merge_and_update_overlap_sets <- function(variants_sets_and_pips_list, overlap_sets) {
    # Combine and identify unique combined sets
    combined_sets <- unique(unlist(lapply(overlap_sets, function(x) {
      paste(sort(x),
        collapse = ","
      )
    })))
    unique_combined_sets <- unique(combined_sets)

    # Split each combined set into individual sets
    split_sets <- lapply(unique_combined_sets, function(x) strsplit(x, ",")[[1]])

    # Identify and merge overlapping credible sets
    if (length(split_sets) != 1) {
      for (i in 1:(length(split_sets) - 1)) {
        for (j in (i + 1):length(split_sets)) {
          if (!is.null(split_sets[[i]]) && !is.null(split_sets[[j]])) {
            # Check both sets exist
            if (length(intersect(split_sets[[i]], split_sets[[j]])) > 0) {
              # Merge overlapping sets Update both i-th and j-th elements with
              # the merged set
              split_sets[[i]] <- unique(c(split_sets[[i]], split_sets[[j]]))
              split_sets[[j]] <- unique(c(split_sets[[i]], split_sets[[j]]))
            }
          }
        }
      }
    }
    # Eliminate duplicates from the list of sets Convert each set into a string
    # to facilitate comparison
    set_strings <- sapply(split_sets, function(set) paste(sort(set), collapse = ","))
    # Identify unique sets based on their string representation
    unique_set_strings <- unique(set_strings)
    # Retain only the largest combined sets
    final_combined_sets <- character()
    for (set in unique_set_strings) {
      if (length(unique_set_strings) == 1) {
        final_combined_sets <- unique_set_strings[[1]]
      } else {
        included_in_other_set <- FALSE
        set_elements <- unlist(strsplit(set, ","))
        for (other_set in unique_set_strings) {
          if (set != other_set && all(set_elements %in% unlist(strsplit(
            other_set,
            ","
          )))) {
            included_in_other_set <- TRUE
            break
          }
        }
        if (!included_in_other_set) {
          final_combined_sets <- c(final_combined_sets, set)
        }
      }
    }

    # Create a mapping from original set names to combined set names
    set_name_map <- list()
    for (combined_set in final_combined_sets) {
      original_sets <- unlist(strsplit(combined_set, ","))
      for (set in original_sets) {
        set_name_map[[set]] <- combined_set
      }
    }

    # Update the credible_set_names in variants_sets_and_pips_list for each
    # variant_id
    updated_credible_sets <- list()
    for (variant_id in names(variants_sets_and_pips_list)) {
      current_sets <- variants_sets_and_pips_list[[variant_id]][["sets"]] # All credible sets for the current variant_id
      combined_set_found <- FALSE
      # Check if any of the current variant_id's credible sets exist in the
      # set_name_map
      for (set_name in current_sets) {
        if (set_name %in% names(set_name_map)) {
          # If at least one credible set corresponds to a combined credible
          # set, update accordingly
          updated_credible_sets[[variant_id]] <- set_name_map[[set_name]]
          combined_set_found <- TRUE
          break # Exit the loop once the combined credible set is found
        }
      }
      # If no combined credible set is found, keep the original credible sets
      # unchanged
      if (!combined_set_found) {
        updated_credible_sets[[variant_id]] <- current_sets
      }
    }
    return(updated_credible_sets)
  }
  # Loop through each condition and their credible sets
  extract_top_loci <- function(susie_fit, complementary, coverage) {
    results <- list()
    for (i in 1:length(names(susie_fit[[1]]))) {
      if (!is.null(susie_fit[[1]][[i]][["top_loci"]]) && nrow(susie_fit[[1]][[i]][["top_loci"]]) !=
        0) {
        if (!complementary) {
          set_num <- unique(get_nested_element(susie_fit[[1]][[i]], c(
            "top_loci",
            coverage
          )))
          set_num <- set_num[set_num != 0]
        } else {
          set_num <- 0
        }
        num_cs <- length(set_num)
        if (num_cs > 0) {
          for (j in 1:num_cs) {
            variants_df <- get_nested_element(susie_fit[[1]][[i]], c("top_loci")) %>%
              filter(!!sym(coverage) == set_num[j]) %>%
              select(variant_id, pip)
            # Iterate through the rows of the variants_df
            if (dim(variants_df)[1] != 0) {
              for (row in 1:nrow(variants_df)) {
                variant_id <- variants_df$variant_id[row]
                variant_pip <- variants_df$pip[row]

                # Prepare the set name
                set_name <- paste0("cs_", i, "_", set_num[j])

                # If the variant_id is not in the results list, add it with the
                # current set_name and pip
                if (!variant_id %in% names(results)) {
                  results[[variant_id]] <- list(sets = set_name, pips = variant_pip)
                } else {
                  # If the variant_id is already in the results, append the current
                  # set_name and pip
                  results[[variant_id]]$sets <- c(
                    results[[variant_id]]$sets,
                    set_name
                  )
                  results[[variant_id]]$pips <- c(
                    results[[variant_id]]$pips,
                    variant_pip
                  )
                }
              }
            }
          }
        }
      }
    }
    return(results)
  }

  combine_top_loci <- function(extracted_result) {
    top_loci_df <- do.call(rbind, lapply(names(extracted_result), function(variant_id) {
      max_pip <- max(unlist(extracted_result[[variant_id]]$pips))
      median_pip <- median(unlist(extracted_result[[variant_id]]$pips))
      if (length(identify_overlap_sets(extracted_result)) != 0) {
        credible_set_names <- merge_and_update_overlap_sets(extracted_result,
          overlap_sets = identify_overlap_sets(extracted_result)
        )[[variant_id]]
      } else {
        credible_set_names <- paste(sort(unique(unlist(extracted_result[[variant_id]]$sets))),
          collapse = ","
        )
      }
      data.frame(
        variant_id = variant_id, credible_set_names = credible_set_names,
        max_pip = max_pip, median_pip = median_pip, stringsAsFactors = FALSE # Avoid factors for strings
      )
    }))
    return(top_loci_df)
  }

  extracted_top_loci <- extract_top_loci(susie_fit, complementary, coverage = coverage)
  combined_top_loci_df <- combine_top_loci(extracted_top_loci)
  # Clean up row names and make sure variant_id is unique
  combined_top_loci_df <- combined_top_loci_df[!duplicated(combined_top_loci_df$variant_id), ]
  rownames(combined_top_loci_df) <- NULL # Clean up row names
  return(combined_top_loci_df)
}

#' @importFrom data.table as.data.table setnames
#' @export
load_multitrait_R_sumstat <- function(susie_fit, sumstats_db, coverage = NULL, top_loci = FALSE, filter_file = NULL, exclude_condition = NULL, ld_meta_file = NULL, remove_any_missing = TRUE, max_rows_selected = 300, nan_remove = FALSE, condition_filter = FALSE) {
    # Internal recursive filtering function
  filter_nested_list <- function(input_list, valid_conditions) {
    if (is.null(valid_conditions) || length(valid_conditions) == 0) {
      return(input_list)
    }
    
    filtered_list <- list()
    
    # Recursively process list
    for (name in names(input_list)) {
      # Check if current name matches any valid condition
      if (name %in% valid_conditions) {
        filtered_list[[name]] <- input_list[[name]]
      } else {
        # Recursively check nested lists
        if (is.list(input_list[[name]])) {
          nested_result <- filter_nested_list(input_list[[name]], valid_conditions)
          if (length(nested_result) > 0) {
            filtered_list[[name]] <- nested_result
          }
        }
      }
    }
    
    return(filtered_list)
  }

  # Apply condition filtering before processing
  if (!is.null(condition_filter)) {
    # Convert condition_filter to character vector if it's not already
    if (!is.character(condition_filter)) {
      condition_filter <- as.character(condition_filter)
    }
    
    # Split if comma-separated string
    if (length(condition_filter) == 1 && grepl(",", condition_filter)) {
      condition_filter <- strsplit(condition_filter, ",")[[1]]
    }
    
    # Trim whitespace
    condition_filter <- trimws(condition_filter)
    
    # Filter susie_fit and sumstats_db
    susie_fit <- filter_nested_list(susie_fit, condition_filter)
    sumstats_db <- filter_nested_list(sumstats_db, condition_filter)
  }
    
    extract_data <- function(sumstats_db) {
      process_element <- function(element, path = character()) {
        if (is.list(element)) {
          if (all(c("variant_names", "sumstats") %in% names(element))) {
            variant_names <- element$variant_names
            sumstats <- element$sumstats
            
            # Calculate z_scores
            if (all(c("betahat", "sebetahat") %in% names(sumstats))) {
              z_scores <- sumstats$betahat / sumstats$sebetahat
            } else if ("z" %in% names(sumstats)) {
              z_scores <- sumstats$z
            } else {
              return(NULL)  # If we can't find z-scores or calculate them, return NULL
            }
            
            return(data.frame(
              variants = variant_names,
              z = z_scores
            ))
          } else {
            results <- list()
            for (name in names(element)) {
              result <- process_element(element[[name]], c(path, name))
              if (!is.null(result)) {
                # Use the original name as the list name
                results[[name]] <- result
              }
            }
            return(results)
          }
        }
        return(NULL)
      }
      
      results <- process_element(sumstats_db)
      
      if (!is.null(results) && length(results) > 0) {
        return(results)
      } else {
        return(list())
      }
    }
      
  merge_matrices <- function(matrix_list, value_column, ld_meta_file, id_column = "variants",
                               remove_any_missing = FALSE) {
    # Input validation
    if (!is.list(matrix_list) || length(matrix_list) == 0) {
      stop("matrix_list must be a non-empty list")
    }
    if (!is.character(value_column) || length(value_column) != 1) {
      stop("value_column must be a single string")
    }
    if (!is.character(id_column) || length(id_column) != 1) {
      stop("id_column must be a single string")
    }
    
    df_list <- lapply(seq_along(matrix_list), function(i) {
      tryCatch({
        # Step 1: Convert matrix to data frame and extract relevant columns
        df <- as.data.frame(matrix_list[[i]])
        if (!(id_column %in% colnames(df)) || !(value_column %in% colnames(df))) {
          stop(paste("Required columns", id_column, "or", value_column, "not found in dataset", i))
        }
        df2 <- df[, c(id_column, value_column)]
   
        if (!is.null(ld_meta_file)) {
            # Step 2: Split 'variants' to extract chromosomal info
            cohort_variants_df <- parse_variant_id(df2[, c(id_column)])
            # Step 3: Combine extracted chromosomal info with value column
            cohort_df <- cbind(cohort_variants_df, bhat = df2[, value_column, drop = FALSE])
    
            # Step 4: Merge with LD reference and filter
            variants_ld_block_match <- merge(cohort_df, ld_meta_file, by = "chrom", allow.cartesian = TRUE) %>%
              filter(pos > start & pos < end) %>%
              select(-path)
    
            # Function to process each group
            process_group <- function(data) {
              # Construct file path
              bim_file_path <- unique(data$bim_path)
              ld_bim_file <- fread(bim_file_path)
    
              # Perform allele quality control
              flipped_data <- allele_qc(data[, 1:4], ld_bim_file$V2, data,
                col_to_flip = c(value_column),
                match_min_prop = 0, remove_dups = FALSE,
                remove_indels = FALSE, remove_strand_ambiguous = FALSE,
                flip_strand = FALSE, remove_unmatched = TRUE, target_gwas = FALSE
              )$target_data_qced
              return(flipped_data)
            }
    
            final_df <- variants_ld_block_match %>%
              group_by(start, end) %>%
              group_map(~ process_group(.x)) %>%
              bind_rows() %>%
              mutate(variant_id = paste0("chr", variant_id)) %>%
              select(c("variant_id", value_column)) %>%
              rename("variants" = "variant_id")
            # Rename columns to avoid duplication
            colnames(final_df) <- c(id_column, paste0(value_column, "_", i))
        } else {
          final_df <- df2
          colnames(final_df) <- c(id_column, paste0(value_column, "_", i))
        }
        return(final_df)
      }, error = function(e) {
        message(paste("Error processing dataset", i, ":", e$message))
        return(NULL)
      })
    })
    
    # Remove any NULL results from errors
    df_list <- df_list[!sapply(df_list, is.null)]
    
    if (length(df_list) == 0) {
      stop("No valid datasets after processing")
    }
    
    # Iteratively merge the data frames
    merged_df <- Reduce(
      function(x, y) merge(x, y, by = id_column, all = TRUE),
      df_list
    )
    
    # Optionally, remove rows with any missing values
    if (remove_any_missing) {
      merged_df <- merged_df[complete.cases(merged_df), ]
    }
    return(merged_df)
  }

  split_variants_and_match <- function(variant, filter_file, max_rows_selected) {
    if (!file.exists(filter_file)) {
      stop("Filter file does not exist.")
    }

    # Split the variant vector into components
    variant_df <- parse_variant_id(variant)

    # get the region of interest
    min_pos <- min(variant_df$pos)
    max_pos <- max(variant_df$pos)
    chrom <- unique(variant_df$chrom)
    if (length(chrom) != 1) {
      stop("Variants are from multiple chromosomes. Cannot create a single range string.")
    }
    region <- paste0(chrom, ":", min_pos, "-", max_pos)
    ref_table <- tabix_region(filter_file, region)
    if (is.null(ref_table)) {
      stop("No variants in the region.")
    }
    colnames(ref_table)[1:2] <- c("#CHROM", "POS")
    if (!all(c("#CHROM", "POS") %in% colnames(ref_table))) {
      stop("Filter file must contain columns: #CHROM, POS.")
    }
    matched_indices <- which(variant_df$chrom %in% ref_table$`#CHROM` & variant_df$pos %in% ref_table$POS)
    if (!is.null(max_rows_selected) && max_rows_selected > 0 && max_rows_selected < length(matched_indices)) {
      selected_rows <- sample(length(matched_indices), max_rows_selected)
      matched_indices <- matched_indices[selected_rows]
    }
    return(matched_indices)
  }

  results <- lapply(sumstats_db[[1]], function(data) extract_data(data))
  trait_names <- names(results)
  z_scores <- merge_matrices(results, value_column = "z", ld_meta_file, id_column = "variants", remove_any_missing)
  out <- list(z = z_scores)
  
  var_idx <- 1:nrow(out[[1]])
  if (!is.null(filter_file)) {
    variants <- out[[1]]$variants
    print(paste("variants:", variants))
    var_idx <- split_variants_and_match(variants, filter_file, max_rows_selected)
  }

  if (top_loci) {
    union_top_loci <- merge_susie_cs(susie_fit, coverage)
    if (!is.null(union_top_loci)) {
      strong_signal_df <- union_top_loci %>%
        group_by(credible_set_names) %>%
        filter(median_pip == max(median_pip)) %>%
        slice(1) %>%
        ungroup()
      var_idx <- which(out[[1]]$variants %in% strong_signal_df$variant_id)
    } else {
      var_idx <- NULL
    }
  }

  # Extract only subset of data
  variants <- out[[1]]$variants[var_idx]
  for (key in names(out)) {
    out[[key]] <- out[[key]][var_idx, , drop = FALSE]
    rownames(out[[key]]) <- variants
    colnames(out[[key]])[2:ncol(out[[key]])] <- trait_names
    out[[key]] <- out[[key]][, -which(names(out[[key]]) == "variants"), drop = FALSE]
  }
  out$region <- names(susie_fit)

  if (!is.null(exclude_condition) && length(exclude_condition) != 0) {
    for (key in setdiff(names(out), "region")) {
      if (all(exclude_condition %in% colnames(out[[key]]))) {
        out[[key]] <- out[[key]][, -exclude_condition]
      } else {
        stop(paste("Error: exclude_condition are not present in", out$region))
      }
    }
  }

  return(out)
}
                    
#' @export
mash_rand_null_sample <- function(dat, n_random, n_null, exclude_condition, seed = NULL) {
  # Function to extract one data set
  extract_one_data <- function(dat, n_random, n_null) {
    if (is.null(dat)) {
      return(NULL)
    }
    
    abs_z <- abs(dat$z)
    sample_idx <- 1:nrow(abs_z)
    random_idx <- sample(sample_idx, min(n_random, length(sample_idx)), replace = FALSE)
    random <- list(z = dat$z[random_idx, , drop = FALSE])
    
    null.id <- which(apply(abs_z, 1, max) < 2)
    if (length(null.id) == 0) {
      warning(paste("no variants are included in the null dataset because abs_z > 2 for all variants in", dat$region))
      null <- list()
    } else {
      if (length(null.id) < ncol(dat$z)) {
        warning(paste("not enough null data to estimate null correlation in", dat$region))
        null <- list()
      } else {
        null_idx <- sample(null.id, min(n_null, length(null.id)), replace = FALSE)
        null <- list(z = dat$z[null_idx, , drop = FALSE])
      }
    }
    dat <- list(random = random, null = null)
    return(dat)
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  if (length(exclude_condition) > 0) {
    if (all(exclude_condition %in% colnames(dat$z))) {
      dat$z <- dat$z[, -exclude_condition, drop = FALSE]
    } else {
      stop(paste("Error: exclude_condition are not present in", dat$region))
    }
  }
  
  result <- extract_one_data(dat, n_random, n_null)
  return(result)
}

#' @export
merge_mash_data <- function(res_data, one_data) {
  combined_data <- list()
  if (length(res_data) == 0 | is.null(res_data)) {
    return(one_data)
  } else if (length(one_data) == 0 | is.null(one_data)) {
    return(res_data)
  } else {
    for (d in names(one_data)) {
      if (length(one_data[[d]]) == 0 | is.null(one_data[[d]])) {
        combined_data[[d]] <- res_data[[d]] # Keep res_data[[d]] when one_data[[d]] is NULL or empty
        next
      } else {
        # Check if the res_data is NULL
        if (!is.null(res_data[[d]]) | length(res_data[[d]]) != 0) {
          # Check if the number of columns matches
          if (!identical(colnames(res_data[[d]]), colnames(one_data[[d]]))) {
            # Get all column names from both data frames
            all_cols <- union(colnames(res_data[[d]]), colnames(one_data[[d]]))

            # Align res[[d]]
            res_aligned <- setNames(as.data.frame(matrix(NaN,
              nrow = nrow(res_data[[d]]),
              ncol = length(all_cols)
            )), all_cols)
            rownames(res_aligned) <- rownames(res_data[[d]])
            common_cols_res <- intersect(colnames(res_data[[d]]), all_cols)
            res_aligned[common_cols_res] <- res_data[[d]][common_cols_res]

            # Align one_data[[d]]
            one_data_aligned <- setNames(as.data.frame(matrix(NaN,
              nrow = nrow(one_data[[d]]),
              ncol = length(all_cols)
            )), all_cols)
            rownames(one_data_aligned) <- rownames(one_data[[d]])
            common_cols_one_data <- intersect(colnames(one_data[[d]]), all_cols)
            one_data_aligned[common_cols_one_data] <- one_data[[d]][common_cols_one_data]

            # Now both have the same columns, we can rbind them
            combined_data[[d]] <- rbind(res_aligned, one_data_aligned)
          } else {
            # If they already have the same number of columns, just rbind
            combined_data[[d]] <- rbind(as.data.frame(res_data[[d]]), as.data.frame(one_data[[d]]))
          }
        } else {
          combined_data[[d]] <- one_data[[d]]
        }
      }
    }
    return(combined_data)
  }
}

#' @importFrom udr ud_init ud_fit
#' @importFrom mashr mash_set_data cov_canonical estimate_null_correlation_simple
#' @export
mash_pipeline <- function(mash_input, alpha, residual_correlation = NULL, unconstrained.update = "ted", set_seed = 999) {
  set.seed(set_seed)
  if (length(mash_input$null.b) == 0 && length(mash_input$null.s) == 0) {
    if (!is.null(residual_correlation)) {
      vhat <- residual_correlation
    } else {
      condition_num <- ncol(mash_input$random.b)
      vhat <- diag(rep(1, condition_num))
    }
  } else {
    vhat <- estimate_null_correlation_simple(mash_set_data(mash_input$null.b,
      Shat = mash_input$null.s,
      alpha, zero_Bhat_Shat_reset = 1000
    ))
  }

  # mash data Fit mixture model using udr package
  mash_data <- mash_set_data(mash_input$strong.b,
    Shat = mash_input$strong.s, V = vhat,
    alpha, zero_Bhat_Shat_reset = 1000
  )
  # Canonical matrices
  U.can <- cov_canonical(mash_data)

  # Penalty strength
  lambda <- ncol(mash_input$strong.z)
  # Initialize udr
  fit0 <- ud_init(mash_data, n_unconstrained = 50, U_scaled = U.can)
  # Fit udr and use penalty as default as suggested by Yunqi penalty is
  # necessary in small sample size case, and there won't be a difference in
  # large sample size
  fit2 <- ud_fit(fit0, control = list(unconstrained.update,
    scaled.update = "fa",
    resid.update = "none", lambda = lambda, penalty.type = "iw", maxiter = 1000,
    tol = 0.01, tol.lik = 0.01
  ))

  # extract data-driven covariance from udr model. (A list of covariance
  # matrices)
  U.ud <- lapply(fit2$U, function(e) e[["mat"]])
  return(mixture_prior = list(U = U.ud, w = fit2$w, loglik = fit2$loglik))
}