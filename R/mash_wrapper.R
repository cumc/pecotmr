matxMax <- function(mtx) {
  return(arrayInd(which.max(mtx), dim(mtx)))
}

handle_invalid_summary_stat = function(x) {
    x$bhat[which(is.nan(x$bhat))] = 0
    x$sbhat[which(is.nan(x$sbhat) | is.infinite(x$sbhat))] = 1E3
    return(x)
}

# This function extracts tensorQTL results for given region for multiple summary statistics files
#' @import dplyr
#' @importFrom data.table fread
#' @export 
load_multitrait_tensorqtl_sumstat <- function(sumstats_paths, region, trait_names, filter_file = NULL, remove_any_missing = TRUE, max_rows_selected = 300, nan_remove=TRUE) {
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
            mutate(variant = paste(`#CHROM`, POS, REF, ALT, sep = ":")) %>%
            select(-c(3, 6:9)) %>%
            distinct(variant, .keep_all = TRUE) %>%
            as.matrix
    }

    extract_component <- function(df, component_index) {
        df %>%
        select(6:ncol(df)) %>%
        mutate(across(everything(), ~as.numeric(strsplit(as.character(.), ":")[[1]][component_index]))) %>%
        as.matrix
    }
    
    Y <- lapply(sumstats_paths, extract_tensorqtl_data, region)

    combined_matrix <- Reduce(function(x, y) merge(x, y, by = c("variant", "#CHROM", "POS", "REF", "ALT")), Y) %>%
        distinct(variant, .keep_all = TRUE)

    if (remove_any_missing) {
        combined_matrix <- combined_matrix[complete.cases(combined_matrix), ]
    }

    if (!is.null(filter_file)) {
        if (!file.exists(filter_file)) {
            stop("Filter file does not exist.")
        }
        filter_df <- tabix_region(filter_file, region)
        if (!all(c("#CHROM", "POS") %in% colnames(filter_df))) {
            stop("Filter file must contain columns: #CHROM, POS.")
        }
        combined_matrix <- inner_join(combined_matrix, filter_df %>% select(`#CHROM`, POS), by = c("#CHROM", "POS"))
    }

    if (!is.null(max_rows_selected) && max_rows_selected > 0 && max_rows_selected < nrow(combined_matrix)) {
        selected_rows <- sample(nrow(combined_matrix), max_rows_selected)
        combined_matrix <- combined_matrix[selected_rows, ]
    }

    dat <- list(
        bhat = extract_component(combined_matrix, 1),
        sbhat = extract_component(combined_matrix, 2)
    )
    if (nan_remove) dat <- handle_invalid_summary_stat(dat)
    #dat$z <- dat$bhat / dat$sbhat
    #rownames(dat$bhat) <- rownames(dat$sbhat) <- rownames(dat$z) <- combined_matrix$variant
    #colnames(dat$bhat) <- colnames(dat$sbhat) <- colnames(dat$z) <- trait_names
    rownames(dat$bhat) <- rownames(dat$sbhat) <- combined_matrix$variant
    colnames(dat$bhat) <- colnames(dat$sbhat) <- trait_names
    dat$region = region
    return(dat)
}

#' @export
load_multitrait_R_sumstat <- function(rds_files, top_loci = FALSE, filter_file = NULL, remove_any_missing = TRUE, max_rows_selected = 300, nan_remove=FALSE) {
  
  read_and_extract <- function(rds_file) {
    dat <- readRDS(rds_file)[[1]]
    
    extract_data <- function(item) {
      bhat <- as.data.table(cbind(item$variant_names, item$sumstats$betahat))
      sbhat <- as.data.table(cbind(item$variant_names, item$sumstats$sebetahat))
      setnames(bhat, c("variants", "bhat"))
      setnames(sbhat, c("variants", "sbhat"))
      bhat[, bhat := as.numeric(bhat)]
      sbhat[, sbhat := as.numeric(sbhat)]
      list(
        bhat = bhat$bhat,
        sbhat = sbhat$sbhat,
        variants = bhat$variants,
        top_variants = item$top_loci[, "variant_id"]
      )
    }
    lapply(dat, extract_data)
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
  
  results <- lapply(rds_files, function(file) read_and_extract(file))
  results <- do.call("c", results)
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
    union_top_loci <- unique(unlist(lapply(results, function(item) item$top_variants)))
    var_idx <- which(variants %in% union_top_loci)
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
  colnames(out$bhat) <- colnames(out$sbhat) <- trait_names
  out$region = rds_files
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
