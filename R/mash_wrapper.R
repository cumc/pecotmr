

matxMax <- function(mtx) {
  return(arrayInd(which.max(mtx), dim(mtx)))
}

handle_nan_etc = function(x) {
    x$bhat[which(is.nan(x$bhat))] = 0
    x$sbhat[which(is.nan(x$sbhat) | is.infinite(x$sbhat))] = 1E3
    return(x)
}

# This function extracts tensorQTL results for given region for multiple summary statistics files
#' @import dplyr
#' @importFrom data.table fread
#' @export 
load_multitrait_tensorqtl_sumstat <- function(sumstats_paths, region, trait_names, filter_file = NULL, remove_any_missing = TRUE, max_rows_selected = 300) {
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

    if (remove_any_missing) {
        combined_matrix <- combined_matrix[complete.cases(combined_matrix), ]
    }

    if (!is.null(max_rows_selected) && max_rows_selected > 0 && max_rows_selected < nrow(combined_matrix)) {
        selected_rows <- sample(nrow(combined_matrix), max_rows_selected)
        combined_matrix <- combined_matrix[selected_rows, ]
    }

    dat <- list(
        bhat = extract_component(combined_matrix, 1),
        sbhat = extract_component(combined_matrix, 2)
    )
    dat <- handle_invalid_summary_stat(dat)
    #dat$z <- dat$bhat / dat$sbhat
    #rownames(dat$bhat) <- rownames(dat$sbhat) <- rownames(dat$z) <- combined_matrix$variant
    #colnames(dat$bhat) <- colnames(dat$sbhat) <- colnames(dat$z) <- trait_names
    rownames(dat$bhat) <- rownames(dat$sbhat) <- combined_matrix$variant
    colnames(dat$bhat) <- colnames(dat$sbhat) <- trait_names
    dat$region = region
    return(dat)
}

#' @export
load_multitrait_R_sumstat <- function(rds_files, trait_names, top_loci = FALSE, filter_file = NULL, remove_any_missing = TRUE, max_rows_selected = 300) {
  read_and_extract <- function(rds_file, top_loci) {
      dat <- readRDS(rds_file)
      if (top_loci) {
        bhats <- lapply(dat, function(x) as.data.table(x$top_loci)[, .(variants, bhat)])
        sbhats <- lapply(dat, function(x) as.data.table(x$top_loci)[, .(variants, sbhat)])
      } else {
        bhats <- lapply(dat, function(x) as.data.table(cbind(x$variant_names, x$univariate_sumstats$bhat)))
        sbhats <- lapply(dat, function(x) as.data.table(cbind(x$variant_names, x$univariate_sumstats$sbhat)))
      }
      list(
          bhat = do.call(cbind, lapply(bhats, `[[`, "bhat")),
          sbhat = do.call(cbind, lapply(sbhats, `[[`, "sbhat")),
          variants = bhats[[1]]$variants
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
    if (!all(c("#CHROM", "POS") %in% colnames(filter_df))) {
        stop("Filter file must contain columns: #CHROM, POS.")
    }
    matched_indices <- which(variant_df$chr %in% ref_table$`#CHROM` & variant_df$pos %in% ref_table$POS)
    if (!is.null(max_rows_selected) && max_rows_selected > 0 && max_rows_selected < length(matched_indices)) {
        selected_rows <- sample(length(matched_indices), max_rows_selected)
        matched_indices <- matched_indices[selected_rows]
    }
    return(matched_indices)
  }

  results <- lapply(rds_files, read_and_extract)

  out <- list(
      bhat = rbindlist(lapply(results, function(x) data.table(variants = x$variants, x$bhat)), fill = TRUE),
      sbhat = rbindlist(lapply(results, function(x) data.table(variants = x$variants, x$sbhat)), fill = TRUE)
  )

  variants <- out$bhat$variants
  if (!is.null(filter_file)) {
    var_idx = split_variants_and_match(variants, filter_file, max_rows_selected)
    variants = variants[var_idx]
  }

  out <- handle_invalid_summary_stat(out)
  #out$z <- out$bhat / out$sbhat
  #rownames(out$bhat) <- rownames(out$sbhat) <- rownames(out$z) <- variants
  #colnames(out$bhat) <- colnames(out$sbhat) <- colnames(out$z) <- trait_names
  rownames(out$bhat) <- rownames(out$sbhat) <- variants
  colnames(out$bhat) <- colnames(out$sbhat) <- trait_names
  dat$region = rds_files
  return(out)
}

#' @export
mash_preprocessing <- function(dat, n_random, n_null, expected_ncondition, exclude_condition, z_only = FALSE, seed=NULL) {
  # Function to remove row names from a list of data frames
  remove_rownames <- function(x) {
    for (name in names(x)) rownames(x[[name]]) <- NULL
    return(x)
  }

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