#' Converted  Variant ID into a properly structured data frame
#' @param variant_id A data frame or character vector representing variant IDs.
#'   Expected formats are a data frame with columns "chrom", "pos", "A1", "A2",
#'   or a character vector in "chr:pos:A1:A2" or "chr:pos_A1_A2" format.
#' @return A data frame with columns "chrom", "pos", "A1", "A2", where 'chrom' 
#'   and 'pos' are integers, and 'A1' and 'A2' are allele identifiers.
#' @import stats
#' @noRd
convert_to_dataframe <- function(variant_id) {
  # Check if target_variants is already a data.frame with the required columns
  if (is.data.frame(variant_id)){
    if (!all(c("chrom", "pos", "A1", "A2") %in% names(variant_id))) {
        names(variant_id) <- c("chrom", "pos", "A1", "A2")
      }
    # Ensure that 'chrom' values are integers
    variant_id$chrom <- ifelse(grepl("^chr", variant_id$chrom), 
                       as.integer(sub("^chr", "", variant_id$chrom)), # Remove 'chr' and convert to integer
                       as.integer(variant_id$chrom))# Convert to integer if not already
    variant_id$pos <- as.integer(variant_id$pos)
    return(variant_id)
  }
  # Function to split a string and create a data.frame
  create_dataframe <- function(string, pattern) {
    # If the pattern is for "chr:pos_ref_at", replace '_' with ':'
    if(pattern == "colon_underscore") {
        string <- gsub("_", ":", string)
    }
    parts <- strsplit(string, ":", fixed = TRUE)
    data <- data.frame(do.call(rbind, parts),stringsAsFactors = FALSE)
    colnames(data) <- c("chrom", "pos", "A1", "A2")
    #Ensure that 'chrom' values are integers
    data$chrom <- ifelse(grepl("^chr", data$chrom), 
                       as.integer(sub("^chr", "", data$chrom)), # Remove 'chr' and convert to integer
                       as.integer(data$chrom))# Convert to integer if not already 
    data$pos <- as.integer(data$pos)
    return(data)
  }

  # Check if id1 is in the first vector format
  if (any(grepl(":", variant_id[1])) && any(grepl("_", variant_id[1]))) {
    return(create_dataframe(variant_id, "colon_underscore"))
  }

  # Check if id1 is in the second vector format
  if (all(grepl(":", variant_id[1]))) {
    return(create_dataframe(variant_id, ":"))
  }
  # If none of the conditions are met, stop and print an error
  stop("Input does not match any expected format. Please provide a valid data frame or a character vector in the specified formats.")
}

#' Match alleles between target_variants and ref_variants
#'
#' Match by ("chrom", "A1", "A2" and "pos"), accounting for possible
#' strand flips and major/minor allele flips (opposite effects and zscores).
#'
#' @param target_variants A data frame with columns "chrom", "pos", "A1", "A2" or strings in the format of "chr:pos:A1:A2"/"chr:pos_A1_A2".
#' @param ref_variants A data frame with columns "chrom", "pos", "A1", "A2" or strings in the format of "chr:pos:A1:A2"/"chr:pos_A1_A2".
#' @param target_data A data frame on which QC procedures will be applied..
#' @param col_to_flip The name of the column in target_data where flips are to be applied.
#' @param match_min_prop Minimum proportion of variants in the smallest data
#'   to be matched, otherwise stops with an error. Default is 20%.
#' @param remove_dups Whether to remove duplicates, default is TRUE.
#' @param flip Whether the alleles must be flipped: A <--> T & C <--> G, in which case 
#'   corresponding `col_to_flip` are multiplied by -1. Default is `TRUE`.
#' @param remove_strand_ambiguous Whether to remove strand SNPs (if any). Default is `TRUE`.
#' @return A single data frame with matched variants.
#' @import dplyr
#' @import tidyr
#' @importFrom vctrs vec_duplicate_detect
#' @export
allele_qc <- function(target_variants, ref_variants, target_data, col_to_flip, 
                      match.min.prop = 0.2, remove_dups = TRUE, flip = TRUE, 
                      remove_strand_ambiguous = TRUE) {
  
  target_variants <- convert_to_dataframe(target_variants)
  ref_variants <- convert_to_dataframe(ref_variants)
  
  matched <- merge(target_variants, ref_variants, by = c("chrom", "pos"), all = FALSE, suffixes = c(".target", ".ref")) %>%
    as.data.frame() %>%
    mutate(variants_id_qced = paste(chrom, pos, A1.ref, A2.ref, sep = ":"))
  
  target_data_qced <- target_data[target_variants$chrom %in% matched$chrom & target_variants$pos %in% matched$pos, , drop = FALSE] %>%
    as.data.frame()
  
  # Align the rownames to the target_data_qced
  target_data_qced$variant = matched$variants_id_qced
  
  a1 <- toupper(matched$A1.target)
  a2 <- toupper(matched$A2.target)
  ref1 <- toupper(matched$A1.ref)
  ref2 <- toupper(matched$A2.ref)
  
  # Strand flip function
  strand_flip <- function(ref) {
    flip <- ref
    flip[ref == "A"] <- "T"
    flip[ref == "T"] <- "A"
    flip[ref == "G"] <- "C"
    flip[ref == "C"] <- "G"
    return(flip)
  }
  
  flip1 <- strand_flip(ref1)
  flip2 <- strand_flip(ref2)
  
  snp <- list()
  if (remove_strand_ambiguous) {
    strand_unambiguous <- !((a1 == "A" & a2 == "T") | (a1 == "T" & a2 == "A") | 
                     (a1 == "C" & a2 == "G") | (a1 == "G" & a2 == "C"))
  } else {
    strand_unambiguous <- rep(TRUE, length(a1))
  }
  snp[["keep"]] <- strand_unambiguous
  # Remove non-ATCG coding
  non_ATCG <- !(a1 %in% c("A", "T", "G", "C") & a2 %in% c("A", "T", "G", "C"))
  snp[["keep"]][non_ATCG] <- FALSE
  snp[["sign_flip"]] <- ((a1 == ref2 & a2 == ref1) | (a1 == flip2 & a2 == flip1)) & (a1 != ref1 & a2 != ref2)
  snp[["strand_flip"]] <- ((a1 == flip1 & a2 == flip2) | (a1 == flip2 & a2 == flip1)) & (a1 != ref1 & a2 != ref2)  
  exact_match <- (a1 == ref1 & a2 == ref2)
  snp[["keep"]][!(exact_match | snp[["sign_flip"]] | snp[["strand_flip"]])] <- FALSE

  if (!any(snp[["strand_flip"]][which(strand_unambiguous)])) {
    # we conclude that strand flip does not exists in the data at all
    # so we can bring back those previous marked to drop because of strand ambiguous
    snp[["keep"]][which(!strand_unambiguous)] <- TRUE
  }

  qc_summary <- matched %>%
    mutate(keep = snp[["keep"]],
           sign_flip = snp[["sign_flip"]],
           strand_flip = snp[["strand_flip"]])
  
  # Apply allele flip if required
  if (flip) {
    if (!is.null(target_data_qced[, col_to_flip])) {
      target_data_qced[qc_summary$sign_flip, col_to_flip] <- -1 * target_data_qced[qc_summary$sign_flip, col_to_flip]
    } else {
      stop("Column '", col_to_flip, "' not found in target_data.")
    }
  }
 
  # Keep SNPs based on the 'keep' flag
  target_data_qced <- target_data_qced[qc_summary$keep, ,drop=FALSE]
  
  # Remove duplicates if specified
  if (remove_dups) {
    dups <- vec_duplicate_detect(qc_summary[, c("chrom", "pos", "A1.target", "A2.target")])
    if (any(dups)) {
      target_data_qced <- target_data_qced[!dups, ,drop=FALSE]
      message("Some duplicates were removed.")
    }
  }
  
  # Check if the minimum proportion of variants is matched
  min_match <- match.min.prop * min(nrow(target_variants), nrow(ref_variants))
  if (nrow(target_data_qced) < min_match) {
    stop("Not enough variants have been matched.")
  }
    # change A1 and A2 so that it can fit the reference, and rearrange the columns so that the four are at the very front
    target_data_qced = target_data_qced %>% tidyr::separate(variant, into = c("chrom", "pos", "A1", "A2"), sep = ":", remove = FALSE) %>% select(chrom, pos, A1, A2, everything()) %>% mutate(chrom = as.integer(chrom), pos = as.integer(pos))
  return(list(target_data_qced = target_data_qced, qc_summary = qc_summary))
}