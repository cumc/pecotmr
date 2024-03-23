# This needs pgenlibr package
# devtools::install_github("chrchang/plink-ng", subdir="/2.0/pgenlibr")

# read PLINK files
#' @importFrom dplyr rename
#' @importFrom data.table fread
#' @importFrom tools file_path_sans_ext
read_pvar <- function(pgen) {
  pvarf <- paste0(file_path_sans_ext(pgen), ".pvar")
  pvardt <- fread(pvarf, skip = "#CHROM")
  pvardt <- rename(pvardt,
    "chrom" = "#CHROM", "pos" = "POS",
    "alt" = "ALT", "ref" = "REF", "id" = "ID"
  )
  pvardt <- pvardt[, c("chrom", "id", "pos", "alt", "ref")]
  return(pvardt)
}

#' @importFrom data.table fread
#' @importFrom tools file_path_sans_ext
read_bim <- function(bed) {
  bimf <- paste0(file_path_sans_ext(bed), ".bim")
  bim <- fread(bimf)
  colnames(bim) <- c("chrom", "id", "gpos", "pos", "a1", "a0")
  return(bim)
}

#' @importFrom data.table fread
#' @importFrom tools file_path_sans_ext
read_psam <- function(pgen) {
  psamf <- paste0(file_path_sans_ext(pgen), ".psam")
  psam <- fread(psamf, header = T)
  colnames(psam)[1:2] <- c("FID", "IID")
  return(psam)
}

#' @importFrom data.table fread
#' @importFrom tools file_path_sans_ext
read_fam <- function(bed) {
  famf <- paste0(file_path_sans_ext(bed), ".fam")
  return(fread(famf, header = F))
}

# open pgen/pvar PLINK 2 data format
#' @importFrom pgenlibr NewPgen
open_pgen <- function(pgenf) {
  return(NewPgen(pgenf))
}

# open bed/bim/fam: A PLINK 1 .bed is a valid .pgen
#' @importFrom pgenlibr NewPgen
open_bed <- function(bed) {
  raw_s_ct <- nrow(read_fam(bed))
  return(NewPgen(bed, raw_sample_ct = raw_s_ct))
}

#' @importFrom pgenlibr GetVariantCt ReadList
read_pgen <- function(pgen, variantidx = NULL, meanimpute = F) {
  if (is.null(variantidx)) {
    variantidx <- 1:GetVariantCt(pgen)
  }

  ReadList(pgen,
    variant_subset = variantidx,
    meanimpute = meanimpute
  )
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

#' @importFrom data.table fread
#' @importFrom dplyr as_tibble mutate
#' @importFrom tibble tibble
#' @importFrom magrittr %>%
tabix_region <- function(file, region, tabix_header = "auto") {
  # Execute tabix command and capture the output
  cmd_output <- tryCatch(
    {
      fread(cmd = paste0("tabix -h ", file, " ", region), sep = "auto", header = tabix_header)
    },
    error = function(e) NULL
  )

  # Check if the output is empty and return an empty tibble if so
  if (is.null(cmd_output) || nrow(cmd_output) == 0) {
    return(tibble())
  }
  cmd_output %>%
    as_tibble() %>%
    mutate(
      !!names(.)[1] := as.character(.[[1]]),
      !!names(.)[2] := as.numeric(.[[2]])
    )
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

#' Load summary statistic data
#'
#' This function formats the input summary statistics dataframe with uniform column names
#' to fit into the SuSiE pipeline. The mapping is performed through the specified column file.
#' Additionally, it extracts sample size, case number, control number, and variance of Y.
#' Missing values in n_sample, n_case, and n_control are backfilled with median values.
#'
#' @param sumstat_path File path to the summary statistics.
#' @param column_file_path File path to the column file for mapping.
#' @param n_sample User-specified sample size. If unknown, set as 0 to retrieve from the sumstat file.
#' @param n_case User-specified number of cases.
#' @param n_control User-specified number of controls.
#'
#' @return A list of rss_input, including the column-name-formatted summary statistics,
#' sample size (n), and var_y.
#'
#' @importFrom data.table fread
#' @export
load_rss_data <- function(sumstat_path, column_file_path, n_sample = 0, n_case = 0, n_control = 0) {
  var_y <- NULL
  sumstats <- fread(sumstat_path)
  column_data <- read.table(column_file_path, header = FALSE, sep = ":", stringsAsFactors = FALSE)
  colnames(column_data) <- c("standard", "original")

  # Map column names using the column file
  for (name in colnames(sumstats)) {
    if (name %in% column_data$original) {
      index <- which(column_data$original == name)
      colnames(sumstats)[colnames(sumstats) == name] <- column_data$standard[index]
    }
  }

  # Calculate z-score if missing
  if (!"z" %in% colnames(sumstats) && all(c("beta", "se") %in% colnames(sumstats))) {
    sumstats$z <- sumstats$beta / sumstats$se
  }

  # Set beta and se if missing
  if (!"beta" %in% colnames(sumstats) && "z" %in% colnames(sumstats)) {
    sumstats$beta <- sumstats$z
    sumstats$se <- 1
  }

  # Backfill missing values with median for n_sample, n_case, and n_control
  for (col in c("n_sample", "n_case", "n_control")) {
    if (col %in% colnames(sumstats)) {
      sumstats[[col]][is.na(sumstats[[col]])] <- median(sumstats[[col]], na.rm = TRUE)
    }
  }

  # Validate and calculate sample size and variance of Y
  if (n_sample != 0 && (n_case + n_control) != 0) {
    stop("Please provide sample size, or case number with control number, but not both")
  } else if (n_sample != 0) {
    n <- n_sample
  } else if ((n_case + n_control) != 0) {
    n <- n_case + n_control
    phi <- n_case / n
    var_y <- 1 / (phi * (1 - phi))
  } else {
    if ("n_sample" %in% colnames(sumstats)) {
      n <- median(sumstats$n_sample)
    } else if (all(c("n_case", "n_control") %in% colnames(sumstats))) {
      n <- median(sumstats$n_case + sumstats$n_control)
      phi <- median(sumstats$n_case / n)
      var_y <- 1 / (phi * (1 - phi))
    } else {
      warning("Sample size and variance of Y could not be determined from the summary statistics.")
      n <- NULL
    }
  }

  return(list(sumstats = sumstats, n = n, var_y = var_y))
}
