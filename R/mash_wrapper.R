matxMax <- function(mtx) {
  return(arrayInd(which.max(mtx), dim(mtx)))
}

handle_invalid_summary_stat = function(x, bhat="beta", sbhat="se") { 
  x$bhat <- x$bhat %>%
    mutate(across(starts_with(bhat), as.numeric)) %>%
    mutate(across(starts_with(bhat), ~replace(., is.nan(.), 0)))
  x$sbhat <- x$sbhat %>%
    mutate(across(starts_with(sbhat), as.numeric)) %>%
    mutate(across(starts_with(sbhat), ~replace(., is.nan(.) | is.infinite(.), 1E3)))
  return(x)
}

# This function extracts tensorQTL results for given region for multiple summary statistics files
#' @import dplyr
#' @importFrom data.table fread
#' @export 
load_multitrait_tensorqtl_sumstat <- function(sumstats_paths, 
                                              region, 
                                              gene=NULL,
                                              trait_names=NULL, 
                                              top_loci=FALSE, 
                                              filter_file = NULL, 
                                              remove_any_missing = TRUE, 
                                              max_rows_selected = 300, 
                                              nan_remove=TRUE) {
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
      mutate(variants = paste(.[[1]] , .[[2]], .[[3]],.[[4]], sep = ":")) %>% #first four columns are "chrom","pos","alt","ref"
      distinct(variants, molecular_trait_id, .keep_all = TRUE) 
  }
  
  
  merge_matrices <- function(matrix_list, value_column, id_column = "variants", remove_any_missing = FALSE) {
    # Convert matrices to data frames
    df_list <- lapply(seq_along(matrix_list), function(i) {
      df <- as.data.frame(matrix_list[[i]])
      df2 <- df[, c(id_column, value_column)]
      # Rename columns to avoid duplication
      colnames(df2) <- c(id_column, paste0(value_column, "_", i))
      return(df2)
    })
    
    # Iteratively merge the data frames
    merged_df <- Reduce(function(x, y) merge(x, y, by = id_column, all = TRUE), df_list)
    
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
    variant_df <- data.frame(
      chr = sapply(variant_split, `[`, 1),
      pos = sapply(variant_split, `[`, 2),
      stringsAsFactors = FALSE
    )
    variant_df$pos <- as.numeric(variant_df$pos)
    
    # get the region of interest
    min_pos <- min(variant_df$pos)
    max_pos <- max(variant_df$pos)
    chrom <- unique(variant_df$chr)
    if (length(chrom) != 1) {
      stop("Variants are from multiple chromosomes. Cannot create a single range string.")
    }
    region = paste0(chrom, ":", min_pos, "-", max_pos)
    ref_table <- tabix_region(filter_file, region)
    if (is.null(ref_table)){
      stop("No variants in the region.")
    }
    colnames(ref_table)[1:2] = c("#CHROM", "POS")
    if (!all(c("#CHROM", "POS") %in% colnames(ref_table))) {
      stop("Filter file must contain columns: #CHROM, POS.")
    }
    matched_indices <- which(variant_df$chr %in% ref_table$`#CHROM` & variant_df$pos %in% ref_table$POS)
    if (!is.null(max_rows_selected) && max_rows_selected > 0 && max_rows_selected < length(matched_indices)) {
      selected_rows <- sample(length(matched_indices), max_rows_selected)
      matched_indices <- matched_indices[selected_rows]
    }
    return(matched_indices)
  }
  
  Y <- lapply(sumstats_paths, function(x, region, gene){
    out <- extract_tensorqtl_data(x, region)
    if (!is.null(gene)){
        out <- out[which(out$molecular_trait_id %in% gene), ]
    }
    sorted <- out[order(-abs(out$beta/out$se)),c("variants", "molecular_trait_id")]
    top_v <- apply(sorted[1:2,],1,function(row) paste(row, collapse = "_")) #paste the variant (chr:pos:alt:ref) with gene_id with "_"
    out <- as.list(out)
    out$top_variants <- top_v
    return(out)}, region=region, gene=gene)
  
  ## Y is list of data frames where
  ## colnames(Y[[1]])= c('chrom','pos','alt','ref','variant_id','molecular_trait_id','start_distance',
  # 'end_distance','af','ma_samples','ma_count','pvalue','beta','se','molecular_trait_object_id',
  # 'n','variant''top_variants')
  
  ### The step below assigns condition names to the Y; 
  # in case the filename itself does not contain any condition names, users can input the condition names via assigning trait_names 
  # if trait_names left blank, then the file names will be assigned as trait_names, 
  # where if the text before "." (extension) can differentiate the conditions, we will use the shorter names as trait_name
  if (is.null(trait_names)){
    trait_names <- gsub("\\..*", "", basename(sumstats_paths)) # extract condition name that is listed before the first appearance of "."
    
    if (length(trait_names[duplicated(trait_names)])>=1){
      trait_names <- basename(sumstats_paths)
    }
  }
  names(Y) <- trait_names 

  bhat <- merge_matrices(Y, value_column="beta",  id_column = c("variants", "molecular_trait_id"), remove_any_missing)
  sbhat <- merge_matrices(Y, value_column="se",  id_column = c("variants", "molecular_trait_id"), remove_any_missing)
  out <- list(bhat = bhat, sbhat = sbhat)
  
  # Check if variants are the same in both bhat and sbhat
  if (!identical(out$bhat$variants, out$sbhat$variants)) {
    stop("Error: Variants in bhat and sbhat are not the same.")
  }
  
  var_idx = 1:nrow(out$bhat)
  
  # match with filter_file
  if (!is.null(filter_file)) {
    if (!file.exists(filter_file)) {
      stop("Filter file does not exist.")
    }
    variants = paste0(out$bhat$variants, "_", out$bhat$molecular_trait_id)
    var_idx = split_variants_and_match(out$bhat$variants, filter_file, max_rows_selected)
  }
  
  
  if (top_loci) {
    union_top_loci <- unique(unlist(lapply(Y, function(item) item$top_variants)))
    var_idx <- which(variants %in% union_top_loci) #var_idx may end up empty if max_rows_selected number too small
  }

  # Extract only subset of data
  variants <- paste0(out$bhat$variants[var_idx], "_", out$bhat$molecular_trait_id[var_idx])
  out$bhat <- out$bhat[var_idx,]
  out$sbhat <- out$sbhat[var_idx,]

  if (nan_remove) out <- handle_invalid_summary_stat(out, bhat="beta", sbhat="se")

  rownames(out$bhat) <- rownames(out$sbhat) <- variants
  colnames(out$bhat)[which(startsWith(colnames(out$bhat), "beta"))] <- colnames(out$sbhat)[which(startsWith(colnames(out$sbhat),  "se"))] <- trait_names
  out$region = region
  out$top_variants <- lapply(Y, function(x) x$top_variants)
  return(out)
}
                             
#' @export
merge_susie_cs <- function(susie_fit, coverage = "cs_coverage_0.95"){
# Initialize an empty list for the results
results <- list()
combined_sets <- list()

#Identify variant IDs that are associated with more than one credible set
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
#Merge overlapping credible sets and update the credible sets in variants_list
merge_and_update_overlap_sets <- function(variants_sets_and_pips_list,overlap_sets) {
  # Combine and identify unique combined sets
  combined_sets <- unique(unlist(lapply(overlap_sets, function(x) paste(sort(x), collapse = ","))))
  unique_combined_sets <- unique(combined_sets)

  # Split each combined set into individual sets
  split_sets <- lapply(unique_combined_sets, function(x) strsplit(x, ",")[[1]])

  #Identify and merge overlapping credible sets
 if(length(split_sets)!=1){
  for (i in 1:(length(split_sets) - 1)) {
   for (j in (i + 1):length(split_sets)) {
    if (!is.null(split_sets[[i]]) && !is.null(split_sets[[j]])) { # Check both sets exist
      if (length(intersect(split_sets[[i]], split_sets[[j]])) > 0) {
        # Merge overlapping sets
        # Update both i-th and j-th elements with the merged set
        split_sets[[i]] <- unique(c(split_sets[[i]], split_sets[[j]]))
        split_sets[[j]] <- unique(c(split_sets[[i]], split_sets[[j]]))
      }
    }
  }
 }
}
# Eliminate duplicates from the list of sets
# Convert each set into a string to facilitate comparison
set_strings <- sapply(split_sets, function(set) paste(sort(set), collapse = ","))
# Identify unique sets based on their string representation
unique_set_strings <- unique(set_strings)
# Retain only the largest combined sets
  final_combined_sets <- character()
  for (set in unique_set_strings) {
    if(length(unique_set_strings)==1){
      final_combined_sets <- unique_set_strings[[1]]
    } else {
      included_in_other_set <- FALSE
      set_elements <- unlist(strsplit(set, ","))
      for (other_set in unique_set_strings) {
        if (set != other_set && all(set_elements %in% unlist(strsplit(other_set, ",")))) {
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

# Update the credible_set_names in variants_sets_and_pips_list for each variant_id                                      
updated_credible_sets <- list()
for (variant_id in names(variants_sets_and_pips_list)) {
    current_sets <- variants_sets_and_pips_list[[variant_id]][["sets"]] # All credible sets for the current variant_id
    combined_set_found <- FALSE
    # Check if any of the current variant_id's credible sets exist in the set_name_map
    for (set_name in current_sets) {
        if (set_name %in% names(set_name_map)) {
            # If at least one credible set corresponds to a combined credible set, update accordingly
            updated_credible_sets[[variant_id]] <- set_name_map[[set_name]]
            combined_set_found <- TRUE
            break # Exit the loop once the combined credible set is found
        }
    }
    # If no combined credible set is found, keep the original credible sets unchanged
    if (!combined_set_found) {
        updated_credible_sets[[variant_id]] <- current_sets
    }
  }
 return(updated_credible_sets)
}                       
# Loop through each condition and their credible sets
for (i in 1:length(names(susie_fit[[1]]))) {
  if (!is.null(susie_fit[[1]][[i]][["top_loci"]])){
  num_cs <- length(which(unique(get_nested_element(susie_fit[[1]][[i]],c("top_loci",coverage)))>0))

  if (num_cs > 0) {
    for (j in 1:num_cs) {
      variants_df <- get_nested_element(susie_fit[[1]][[i]],c("top_loci")) %>%
                     filter(!!sym(coverage) == j) %>%
                     select(variant_id, pip)

      # Iterate through the rows of the variants_df
      for (row in 1:nrow(variants_df)) {
        variant_id <- variants_df$variant_id[row]
        variant_pip <- variants_df$pip[row]

        # Prepare the set name
        set_name <- paste0("cs_", i, "_", j)

        # If the variant_id is not in the results list, add it with the current set_name and pip
        if (!variant_id %in% names(results)) {
          results[[variant_id]] <- list(sets = set_name, pips = variant_pip)
        } else {
          # If the variant_id is already in the results, append the current set_name and pip
          results[[variant_id]]$sets <- c(results[[variant_id]]$sets, set_name)
          results[[variant_id]]$pips <- c(results[[variant_id]]$pips, variant_pip)
          
        }
      }
    }
  }
 }
}
top_loci_df <- do.call(rbind, lapply(names(results), function(variant_id) {
    max_pip <- max(unlist(results[[variant_id]]$pips))
    median_pip <- median(unlist(results[[variant_id]]$pips))
    if(length(identify_overlap_sets(results))!=0){                                   
    credible_set_names <- merge_and_update_overlap_sets(results,overlap_sets = identify_overlap_sets(results))[[variant_id]]
    } else {
    credible_set_names <- paste(sort(unique(unlist(results[[variant_id]]$sets))), collapse = ",")
    }
    data.frame(
      variant_id = variant_id,
      credible_set_names = credible_set_names,
      max_pip = max_pip,
      median_pip = median_pip,
      stringsAsFactors = FALSE  # Avoid factors for strings
    )
}))

# Clean up row names and make sure variant_id is unique
top_loci_df <- top_loci_df[!duplicated(top_loci_df$variant_id), ]
rownames(top_loci_df) <- NULL  # Clean up row names
return(top_loci_df)
}



                             
#' @export
load_multitrait_R_sumstat <- function(susie_fit, sumstats_db, coverage = "cs_coverage_0.95", top_loci = FALSE, filter_file = NULL, remove_any_missing = TRUE, max_rows_selected = 300, nan_remove=FALSE) {
  
   extract_data <- function(sumstats_db) {
      bhat <- as.data.table(cbind(sumstats_db$variant_names, sumstats_db$sumstats$betahat))
      sbhat <- as.data.table(cbind(sumstats_db$variant_names, sumstats_db$sumstats$sebetahat))
      setnames(bhat, c("variants", "bhat"))
      setnames(sbhat, c("variants", "sbhat"))
      bhat[, bhat := as.numeric(bhat)]
      sbhat[, sbhat := as.numeric(sbhat)]
      list(
        bhat = bhat$bhat,
        sbhat = sbhat$sbhat,
        variants = bhat$variants
      )
    }  
  split_variants_and_match <- function(variant, filter_file, max_rows_selected) {
    if (!file.exists(filter_file)) {
      stop("Filter file does not exist.")
    }
    
    # Split the variant vector into components
    variant_split <- strsplit(variant, ":")
    variant_df <- data.frame(
      chr = sapply(variant_split, `[`, 1),
      pos = sapply(variant_split, `[`, 2),
      stringsAsFactors = FALSE
    )
    variant_df$pos <- as.numeric(variant_df$pos)
    
    # get the region of interest
    min_pos <- min(variant_df$pos)
    max_pos <- max(variant_df$pos)
    chrom <- unique(variant_df$chr)
    if (length(chrom) != 1) {
      stop("Variants are from multiple chromosomes. Cannot create a single range string.")
    }
    region = paste0(chrom, ":", min_pos, "-", max_pos)
    ref_table <- tabix_region(filter_file, region)
    if (is.null(ref_table)){
      stop("No variants in the region.")
    }
    colnames(ref_table)[1:2] = c("#CHROM", "POS")
    if (!all(c("#CHROM", "POS") %in% colnames(ref_table))) {
      stop("Filter file must contain columns: #CHROM, POS.")
    }
    matched_indices <- which(variant_df$chr %in% ref_table$`#CHROM` & variant_df$pos %in% ref_table$POS)
    if (!is.null(max_rows_selected) && max_rows_selected > 0 && max_rows_selected < length(matched_indices)) {
      selected_rows <- sample(length(matched_indices), max_rows_selected)
      matched_indices <- matched_indices[selected_rows]
    }
    return(matched_indices)
  }
  
  merge_matrices <- function(matrix_list, value_column, id_column = "variants", remove_any_missing = FALSE) {
    # Convert matrices to data frames
    df_list <- lapply(seq_along(matrix_list), function(i) {
      df <- as.data.frame(matrix_list[[i]])
      df2 <- df[, c(id_column, value_column)]
      # Rename columns to avoid duplication
      colnames(df2) <- c(id_column, paste0(value_column, "_", i))
      return(df2)
    })
    
    # Iteratively merge the data frames
    merged_df <- Reduce(function(x, y) merge(x, y, by = id_column, all = TRUE), df_list)
    
    # Optionally, remove rows with any missing values
    if (remove_any_missing) {
      merged_df <- merged_df[complete.cases(merged_df), ]
    }
    return(merged_df)
  }
  
  results <- lapply(sumstats_db_list[[1]], function(data) extract_data(data))
  # results <- do.call("c", results)
  trait_names <- names(results)
  
  bhat = merge_matrices(results, value_column="bhat",  id_column = "variants", remove_any_missing)
  sbhat = merge_matrices(results, value_column="sbhat", id_column = "variants", remove_any_missing)
  out <- list(bhat = bhat, sbhat = sbhat)
  
  # Check if variants are the same in both bhat and sbhat
  if (!identical(out$bhat$variants, out$sbhat$variants)) {
    stop("Error: Variants in bhat and sbhat are not the same.")
  }
  var_idx = 1:nrow(out$bhat)
  if (!is.null(filter_file)) {
    variants = out$bhat$variants
    var_idx = split_variants_and_match(variants, filter_file, max_rows_selected)
  }
  
  if (top_loci) {
    union_top_loci <- merge_susie_cs(susie_fit,coverage)
    strong_signal_df <- union_top_loci %>% 
                 group_by(credible_set_names) %>%
                 filter(median_pip == max(median_pip)) %>%
                 slice(1) %>% ungroup()
    var_idx <- which(out$bhat$variants %in% strong_signal_df$variant_id)
  }

  
  # Extract only subset of data
  variants <- out$bhat$variants[var_idx]
  out$bhat <- out$bhat[var_idx,]
  out$sbhat <- out$sbhat[var_idx,]
  
   if (nan_remove) out <- handle_invalid_summary_stat(out)
  #out$z <- out$bhat / out$sbhat
  #rownames(out$bhat) <- rownames(out$sbhat) <- rownames(out$z) <- variants
  #colnames(out$bhat) <- colnames(out$sbhat) <- colnames(out$z) <- trait_names
   rownames(out$bhat) <- rownames(out$sbhat) <- variants
   colnames(out$bhat)[2:ncol(out$bhat)] <- colnames(out$sbhat)[2:ncol(out$bhat)] <- trait_names
   #out$region = rds_files
   return(out)
 }

#' @export
mash_ran_null_sample <- function(dat, n_random, n_null, expected_ncondition, exclude_condition, z_only = FALSE, seed=NULL) {
  # Function to extract one data set
  extract_one_data <- function(dat, n_random, n_null) {
    if (is.null(dat)) return(NULL)
    abs_z <- abs(dat$bhat / dat$sbhat)
    sample_idx <- 1:nrow(abs_z)
    random_idx <- sample(sample_idx, min(n_random, length(sample_idx)), replace = FALSE)
    random <- list(bhat = dat$bhat[random_idx,,drop=FALSE], sbhat = dat$sbhat[random_idx,,drop=FALSE])
    null.id <- which(apply(abs_z, 1, max) < 2)
    if (length(null.id) == 0) {
      warning(paste("Null data is empty for input", dat$region))
      null <- list()
    } else {
      null_idx <- sample(null.id, min(n_null, length(null.id)), replace = FALSE)
      null <- list(bhat = dat$bhat[null_idx,,drop=FALSE], sbhat = dat$sbhat[null_idx,,drop=FALSE])
    }
    dat <- list(random = random, null = null)
    return(dat)
  }

  # Function to reformat data
  reformat_data <- function(dat, z_only) {
    res <- list(random.z = dat$random$bhat / dat$random$sbhat, 
                null.z = dat$null$bhat / dat$null$sbhat)
    if (!z_only) {
      res <- c(res, list(random.b = dat$random$bhat,
                         null.b = dat$null$bhat,
                         null.s = dat$null$sbhat,
                         random.s = dat$random$sbhat))
    }
    return(res)
  }

  if (!is.null(seed)) set.seed(seed)

  # Main processing logic
  if (expected_ncondition > 0 && (ncol(dat$bhat) != expected_ncondition || ncol(dat$sbhat) != expected_ncondition)) {
    stop(paste("Input data has", ncol(dat$bhat), "columns different from required", expected_ncondition))
  }

  if (length(exclude_condition) > 0) {
    dat$bhat <- dat$bhat[,-exclude_condition]
    dat$sbhat <- dat$sbhat[,-exclude_condition]
  }

  result <- reformat_data(extract_one_data(dat, n_random, n_null), z_only)
  return(result)
}
