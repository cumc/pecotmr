# This needs pgenlibr package
# devtools::install_github("chrchang/plink-ng", subdir="/2.0/pgenlibr")
# and PLINK2R
# Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS"=TRUE); remotes::install_github('gabraham/plink2R', subdir='plink2R', ref='d74be015e8f54d662b96c6c2a52a614746f9030d')

# read PLINK files
#' @importFrom dplyr rename
#' @importFrom data.table fread
#' @importFrom tools file_path_sans_ext
read_pvar <- function(pgen){
  pvarf <- paste0(file_path_sans_ext(pgen), ".pvar")
  pvardt <- fread(pvarf, skip = "#CHROM")
  pvardt <- rename(pvardt, "chrom" = "#CHROM", "pos" = "POS",
                "alt" = "ALT", "ref" = "REF", "id" = "ID")
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
  psam = fread(psamf, header=T)
  colnames(psam)[1:2] = c("FID", "IID")
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
open_pgen <- function(pgenf){
    return(NewPgen(pgenf))
} 

# open bed/bim/fam: A PLINK 1 .bed is a valid .pgen
#' @importFrom pgenlibr NewPgen
open_bed <- function(bed){
    raw_s_ct <- nrow(read_fam(bed))
    return(NewPgen(bed, raw_sample_ct = raw_s_ct))
}
#' @importFrom pgenlibr GetVariantCt ReadList
read_pgen <- function(pgen, variantidx = NULL, meanimpute = F ) {
  if (is.null(variantidx)){
    variantidx <- 1:GetVariantCt(pgen)}

    ReadList(pgen, variant_subset = variantidx,
                     meanimpute = meanimpute)
}

thisFile <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  needle  <- "--file="
  match   <- grep(needle,cmdArgs)
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
                readChar(fileName,file.info(fileName)$size),""))
}

#' @importFrom data.table fread
#' @import dplyr
tabix_region <- function(file, region, tabix_header = "auto"){
  # Execute tabix command and capture the output
  cmd_output <- tryCatch(
    {
      fread(cmd = paste0("tabix -h ", file, " ", region), sep="auto", header = tabix_header)
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
      #If not, construct a new target path by combining the directory of the reference file path with the target file path
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
    
  stop(sprintf("Both reference and target file paths do not work. Tried paths: '%s' and '%s'", 
               reference_file_path, file.path(dirname(reference_file_path), target_file_path)))
}

find_valid_file_paths <- function(reference_file_path, target_file_paths) sapply(target_file_paths, function(x) find_valid_file_path(reference_file_path, x))
