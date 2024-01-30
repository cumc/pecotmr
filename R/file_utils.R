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
tabix_region <- function(file, region, tabix_header = FALSE ){
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
find_valid_file_path <- function(primary_file_path, fallback_file_path) {
  # Check if the primary file path exits
  try_primary <- function() {
    if (file.exists(primary_file_path)) {
      return(primary_file_path)
    } else {
      return(NULL)
    }
  }
 # Check if the fallback file path exists
  try_fallback <- function() {
      if (file.exists(fallback_file_path)) {
      return(fallback_file_path)
    } else {
      #If not, construct a new fallback path by combining the directory of the primary file path with the fallback file path
     fallback_full_path <- file.path(dirname(primary_file_path), fallback_file_path)
      if (file.exists(fallback_full_path)) {
      return(fallback_full_path)
      } else {
      return(NULL)
      }
    }
  }
    
  fallback_result <- try_fallback()
  if (!is.null(fallback_result)) {
    return(fallback_result)
  }

  primary_result <- try_primary()
  if (!is.null(primary_result)) {
    return(primary_result)
  }
    
  stop(sprintf("Both primary and fallback file paths do not work. Tried paths: '%s' and '%s'", 
               primary_file_path, file.path(dirname(primary_file_path), fallback_file_path)))
}