#' @export
wald_test_pval <- function(beta, se, n) {
  # Calculate the t statistic
  t_value <- beta / se
  # Degrees of freedom
  df <- n - 2
  # Calculate two-tailed p-value
  p_value <- 2 * pt(-abs(t_value), df = df, lower.tail = TRUE)

  return(p_value)
}

pval_acat <- function(pvals) {
  if (length(pvals) == 1) {
    return(pvals[1])
  }
  stat <- 0.00
  pval_min <- 1.00

  stat <- sum(qcauchy(pvals))
  pval_min <- min(pval_min, min(qcauchy(pvals)))

  return(pcauchy(stat / length(pvals), lower.tail = FALSE))
}

#' @importFrom harmonicmeanp pLandau
pval_hmp <- function(pvals) {
  # https://search.r-project.org/CRAN/refmans/harmonicmeanp/html/pLandau.html
  pvalues <- unique(pvals)
  L <- length(pvalues)
  HMP <- L / sum(pvalues^-1)

  LOC_L1 <- 0.874367040387922
  SCALE <- 1.5707963267949

  return(pLandau(1 / HMP, mu = log(L) + LOC_L1, sigma = SCALE, lower.tail = FALSE))
}

pval_global <- function(pvals, comb_method = "HMP", naive = FALSE) {
  # assuming sstats has tissues as columns and rows as pvals
  min_pval <- min(pvals)
  n_total_tests <- pvals %>%
    unique() %>%
    length() # There should be one unique pval per tissue
  global_pval <- if (comb_method == "HMP") pval_hmp(pvals) else pval_acat(pvals) # pval vector
  naive_pval <- min(n_total_tests * min_pval, 1.0)
  return(if (naive) naive_pval else global_pval) # global_pval and naive_pval
}

matxMax <- function(mtx) {
  return(arrayInd(which.max(mtx), dim(mtx)))
}

compute_maf <- function(geno) {
  f <- mean(geno, na.rm = TRUE) / 2
  return(min(f, 1 - f))
}

compute_missing <- function(geno) {
  miss <- sum(is.na(geno)) / length(geno)
  return(miss)
}

compute_non_missing_y <- function(y) {
  nonmiss <- sum(!is.na(y))
  return(nonmiss)
}

compute_all_missing_y <- function(y) {
  allmiss <- all(is.na(y))
  return(allmiss)
}

mean_impute <- function(geno) {
  f <- apply(geno, 2, function(x) mean(x, na.rm = TRUE))
  for (i in 1:length(f)) geno[, i][which(is.na(geno[, i]))] <- f[i]
  return(geno)
}

is_zero_variance <- function(x) {
  if (length(unique(x)) == 1) {
    return(T)
  } else {
    return(F)
  }
}

compute_LD <- function(X) {
  if (is.null(X)) {
    stop("X must be provided.")
  }

  # Mean impute X
  genotype_data_imputed <- apply(X, 2, function(x) {
    pos <- which(is.na(x))
    if (length(pos) != 0) {
      x[pos] <- mean(x, na.rm = TRUE)
    }
    return(x)
  })

  # Check if Rfast package is installed
  if (requireNamespace("Rfast", quietly = TRUE)) {
    # Use Rfast::cora for faster correlation calculation
    R <- Rfast::cora(genotype_data_imputed, large = TRUE)
  } else {
    # Use base R cor function if Rfast is not installed
    R <- cor(genotype_data_imputed)
  }

  colnames(R) <- rownames(R) <- colnames(genotype_data_imputed)
  R
}

#' @importFrom matrixStats colVars
filter_X <- function(X, missing_rate_thresh, maf_thresh, var_thresh = 0) {
  rm_col <- which(apply(X, 2, compute_missing) > missing_rate_thresh)
  if (length(rm_col)) X <- X[, -rm_col]
  rm_col <- which(apply(X, 2, compute_maf) <= maf_thresh)
  if (length(rm_col)) X <- X[, -rm_col]
  rm_col <- which(apply(X, 2, is_zero_variance))
  if (length(rm_col)) X <- X[, -rm_col]
  X <- mean_impute(X)
  if (var_thresh > 0) {
    rm_col <- which(matrixStats::colVars(X) < var_thresh)
    if (length(rm_col)) X <- X[, -rm_col]
  }
  return(X)
}

filter_Y <- function(Y, n_nonmiss) {
  rm_col <- which(apply(Y, 2, compute_non_missing_y) < n_nonmiss)
  if (length(rm_col)) Y <- Y[, -rm_col]
  rm_rows <- NULL
  if (is.matrix(Y)) {
    rm_rows <- which(apply(Y, 1, compute_all_missing_y))
    if (length(rm_rows)) Y <- Y[-rm_rows, ]
  } else {
    Y <- Y[which(!is.na(Y))]
  }
  return(list(Y = Y, rm_rows = rm_rows))
}

format_variant_id <- function(names_vector) {
  gsub("_", ":", names_vector)
}

#' Converted  Variant ID into a properly structured data frame
#' @param variant_id A data frame or character vector representing variant IDs.
#'   Expected formats are a data frame with columns "chrom", "pos", "A1", "A2",
#'   or a character vector in "chr:pos:A2:A1" or "chr:pos_A2_A1" format.
#' @return A data frame with columns "chrom", "pos", "A1", "A2", where 'chrom'
#'   and 'pos' are integers, and 'A1' and 'A2' are allele identifiers.
#' @noRd
variant_id_to_df <- function(variant_id) {
  # Check if target_variants is already a data.frame with the required columns
  if (is.data.frame(variant_id)) {
    if (!all(c("chrom", "pos", "A1", "A2") %in% names(variant_id))) {
      names(variant_id) <- c("chrom", "pos", "A2", "A1")
    }
    # Ensure that 'chrom' values are integers
    variant_id$chrom <- ifelse(grepl("^chr", variant_id$chrom),
      as.integer(sub("^chr", "", variant_id$chrom)), # Remove 'chr' and convert to integer
      as.integer(variant_id$chrom)
    ) # Convert to integer if not already
    variant_id$pos <- as.integer(variant_id$pos)
    return(variant_id)
  }
  # Function to split a string and create a data.frame
  create_dataframe <- function(string) {
    string <- gsub("_", ":", string)
    parts <- strsplit(string, ":", fixed = TRUE)
    data <- data.frame(do.call(rbind, parts), stringsAsFactors = FALSE)
    colnames(data) <- c("chrom", "pos", "A2", "A1")
    # Ensure that 'chrom' values are integers
    data$chrom <- ifelse(grepl("^chr", data$chrom),
      as.integer(sub("^chr", "", data$chrom)), # Remove 'chr' and convert to integer
      as.integer(data$chrom)
    ) # Convert to integer if not already
    data$pos <- as.integer(data$pos)
    return(data)
  }
  return(create_dataframe(variant_id))
}

#' @importFrom stringr str_split
#' @export
parse_region <- function(region) {
  if (!is.character(region) || length(region) != 1) {
    return(region)
  }

  if (!grepl("^chr[0-9XY]+:[0-9]+-[0-9]+$", region)) {
    stop("Input string format must be 'chr:start-end'.")
  }
  parts <- str_split(region, "[:-]")[[1]]
  df <- data.frame(
    chrom = gsub("^chr", "", parts[1]),
    start = as.integer(parts[2]),
    end = as.integer(parts[3])
  )

  return(df)
}

# Retrieve a nested element from a list structure
#' @export
get_nested_element <- function(nested_list, name_vector) {
  if (is.null(name_vector)) {
    return(NULL)
  }
  current_element <- nested_list
  for (name in name_vector) {
    if (is.null(current_element[[name]])) {
      stop("Element not found in the list")
    }
    current_element <- current_element[[name]]
  }
  return(current_element)
}

thisFile <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmdArgs)
  if (length(match) > 0) {
    ## Rscript
    path <- cmdArgs[match]
    path <- gsub("\\~\\+\\~", " ", path)
    return(normalizePath(sub(needle, "", path)))
  } else {
    ## 'source'd via R console
    return(sys.frames()[[1]]$ofile)
  }
}

load_script <- function() {
  fileName <- thisFile()
  return(ifelse(!is.null(fileName) && file.exists(fileName),
    readChar(fileName, file.info(fileName)$size), ""
  ))
}

#' Find Valid File Path
find_valid_file_path <- function(reference_file_path, target_file_path) {
  # Check if the reference file path exits
  try_reference <- function() {
    if (file.exists(reference_file_path)) {
      return(reference_file_path)
    } else {
      return(NULL)
    }
  }
  # Check if the target file path exists
  try_target <- function() {
    if (file.exists(target_file_path)) {
      return(target_file_path)
    } else {
      # If not, construct a new target path by combining the directory of the reference file path with the target file path
      target_full_path <- file.path(dirname(reference_file_path), target_file_path)
      if (file.exists(target_full_path)) {
        return(target_full_path)
      } else {
        return(NULL)
      }
    }
  }

  target_result <- try_target()
  if (!is.null(target_result)) {
    return(target_result)
  }

  reference_result <- try_reference()
  if (!is.null(reference_result)) {
    return(reference_result)
  }

  stop(sprintf(
    "Both reference and target file paths do not work. Tried paths: '%s' and '%s'",
    reference_file_path, file.path(dirname(reference_file_path), target_file_path)
  ))
}

find_valid_file_paths <- function(reference_file_path, target_file_paths) sapply(target_file_paths, function(x) find_valid_file_path(reference_file_path, x))

#' Filter a vector based on a correlation matrix
#'
#' This function filters a vector `z` based on a correlation matrix `LD` and a correlation threshold `rThreshold`.
#' It keeps only one element among those having an absolute correlation value greater than the threshold.
#'
#' @param z A numeric vector to be filtered.
#' @param LD A square correlation matrix with dimensions equal to the length of `z`.
#' @param rThreshold The correlation threshold for filtering.
#'
#' @return A list containing the following elements:
#'   \describe{
#'     \item{filteredZ}{The filtered vector `z` based on the correlation threshold.}
#'     \item{filteredLD}{The filtered matrix `LD` based on the correlation threshold.}
#'     \item{dupBearer}{A vector indicating the duplicate status of each element in `z`.}
#'     \item{corABS}{A vector storing the absolute correlation values of duplicates.}
#'     \item{sign}{A vector storing the sign of the correlation values (-1 for negative, 1 for positive).}
#'     \item{minValue}{The minimum absolute correlation value encountered.}
#'   }
#'
#' @examples
#' z <- c(1, 2, 3, 4, 5)
#' LD <- matrix(c(
#'   1.0, 0.8, 0.2, 0.1, 0.3,
#'   0.8, 1.0, 0.4, 0.2, 0.5,
#'   0.2, 0.4, 1.0, 0.6, 0.1,
#'   0.1, 0.2, 0.6, 1.0, 0.3,
#'   0.3, 0.5, 0.1, 0.3, 1.0
#' ), nrow = 5, ncol = 5)
#' rThreshold <- 0.5
#'
#' result <- find_duplicate_variants(z, LD, rThreshold)
#' print(result)
#'
#' @export
find_duplicate_variants <- function(z, LD, rThreshold) {
  p <- length(z)
  dupBearer <- rep(-1, p)
  corABS <- rep(0, p)
  sign <- rep(1, p)
  count <- 1
  minValue <- 1

  for (i in 1:(p - 1)) {
    if (dupBearer[i] != -1) next

    idx <- (i + 1):p
    corVec <- abs(LD[i, idx])
    dupIdx <- which(dupBearer[idx] == -1 & corVec > rThreshold)

    if (length(dupIdx) > 0) {
      j <- idx[dupIdx]
      sign[j] <- ifelse(LD[i, j] < 0, -1, sign[j])
      corABS[j] <- corVec[dupIdx]
      dupBearer[j] <- count
    }

    minValue <- min(minValue, min(corVec))
    count <- count + 1
  }

  # Filter z based on dupBearer
  filteredZ <- z[dupBearer == -1]
  filteredLD <- LD[dupBearer == -1, dupBearer == -1, drop = F]

  return(list(filteredZ = filteredZ, filteredLD = filteredLD, dupBearer = dupBearer, corABS = corABS, sign = sign, minValue = minValue))
}
