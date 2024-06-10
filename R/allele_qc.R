#' Match alleles between target_variants and ref_variants
#'
#' Match by ("chrom", "A1", "A2" and "pos"), accounting for possible
#' strand flips and major/minor allele flips (opposite effects and zscores).
#'
#' @param target_variants A data frame with columns "chrom", "pos", "A1", "A2" or strings in the format of "chr:pos:A2:A1"/"chr:pos_A2_A1".
#' @param ref_variants A data frame with columns "chrom", "pos", "A1", "A2" or strings in the format of "chr:pos:A2:A1"/"chr:pos_A2_A1".
#' @param target_data A data frame on which QC procedures will be applied.
#' @param col_to_flip The name of the column in target_data where flips are to be applied.
#' @param match_min_prop Minimum proportion of variants in the smallest data
#'   to be matched, otherwise stops with an error. Default is 20%.
#' @param remove_dups Whether to remove duplicates, default is TRUE.
#' @param remove_indels Whether to remove INDELs, default is FALSE.
#' @param flip Whether the alleles must be flipped: A <--> T & C <--> G, in which case
#'   corresponding `col_to_flip` are multiplied by -1. Default is `TRUE`.
#' @param remove_strand_ambiguous Whether to remove strand SNPs (if any). Default is `TRUE`.
#' @param flip_strand Whether to output the variants after strand flip. Default is `FALSE`.
#' @param remove_unmatched Whether to remove unmatched variants. Default is `TRUE`.
#' @return A single data frame with matched variants.
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate inner_join filter pull select everything row_number
#' @importFrom vctrs vec_duplicate_detect
#' @importFrom tidyr separate
#' @export
allele_qc <- function(target_variants, ref_variants, target_data, col_to_flip = NULL,
                      match_min_prop = 0.2, remove_dups = TRUE,
                      remove_indels = FALSE, remove_strand_ambiguous = TRUE,
                      flip_strand = FALSE, remove_unmatched = TRUE, target_gwas = TRUE) {
  target_variants <- variant_id_to_df(target_variants)
  ref_variants <- variant_id_to_df(ref_variants)
  if (isTRUE(target_gwas)) target_data <- variant_id_to_df(target_data)

  matched <- merge(target_variants, ref_variants, by = c("chrom", "pos"), all = FALSE, suffixes = c(".target", ".ref")) %>%
    as.data.frame() %>%
    mutate(variants_id_qced = paste(chrom, pos, A2.ref, A1.ref, sep = ":"))

  matched_indices <- target_variants %>%
    mutate(index = row_number()) %>%
    inner_join(matched, by = c("chrom" = "chrom", "pos" = "pos", "A1" = "A1.target", "A2" = "A2.target")) %>%
    filter(duplicated(.) | !duplicated(.)) # this will keep all rows, duplicates and non-duplicates alike
  variants_id_qced <- matched_indices %>% pull(variants_id_qced)
  target_data_qced <- target_data[matched_indices$index, , drop = FALSE] %>%
    as.data.frame() %>%
    mutate(variant_id = variants_id_qced)

  a1 <- toupper(matched$A1.target)
  a2 <- toupper(matched$A2.target)
  ref1 <- toupper(matched$A1.ref)
  ref2 <- toupper(matched$A2.ref)

  strand_flip <- function(ref) {
    # Define a mapping for complementary bases
    base_mapping <- c("A" = "T", "T" = "A", "G" = "C", "C" = "G")

    # Function to complement a single base
    complement_base <- function(base) {
      complement <- base_mapping[base]
      return(complement)
    }

    # Function to reverse and complement a DNA sequence
    reverse_complement <- function(sequence) {
      reversed <- rev(strsplit(sequence, NULL)[[1]])
      complemented <- sapply(reversed, complement_base, USE.NAMES = FALSE)
      return(paste(complemented, collapse = ""))
    }

    complemented_sequence <- sapply(ref, reverse_complement, USE.NAMES = FALSE)

    return(complemented_sequence)
  }

  flip1 <- strand_flip(ref1)
  flip2 <- strand_flip(ref2)
  remove_strand_ambiguous <- TRUE
  snp <- list()

  if (remove_strand_ambiguous) {
    strand_unambiguous <- !((a1 == "A" & a2 == "T") | (a1 == "T" & a2 == "A") |
      (a1 == "C" & a2 == "G") | (a1 == "G" & a2 == "C"))
  } else {
    strand_unambiguous <- rep(TRUE, length(a1))
  }

  snp[["keep"]] <- strand_unambiguous

  # Remove non-ATCG coding
  check_ATCG <- function(vec) {
    pattern <- "^[ATCGDI]+$"

    # Function to check if a single element matches the pattern
    check_element <- function(element) {
      grepl(pattern, element)
    }

    result <- sapply(vec, check_element, USE.NAMES = FALSE)
    return(result)
  }

  non_ATCG <- !(check_ATCG(a1) & check_ATCG(a2))


  snp[["keep"]][non_ATCG] <- FALSE
  snp[["sign_flip"]] <- ((a1 == ref2 & a2 == ref1) | (a1 == flip2 & a2 == flip1)) & (a1 != ref1 & a2 != ref2)
  snp[["strand_flip"]] <- ((a1 == flip1 & a2 == flip2) | (a1 == flip2 & a2 == flip1)) & (a1 != ref1 & a2 != ref2)
  snp[["INDEL"]] <- (a2 == "I" | a2 == "D" | nchar(a2) > 1 | nchar(a1) > 1)
  exact_match <- (a1 == ref1 & a2 == ref2)
  ID_match <- ((a2 == "D" | a2 == "I") & (nchar(ref1) > 1 | nchar(ref2) > 1))

  snp[["keep"]][!(exact_match | snp[["sign_flip"]] | snp[["strand_flip"]] | ID_match)] <- FALSE

  if (!any(snp[["strand_flip"]][which(strand_unambiguous)])) {
    # we conclude that strand flip does not exists in the data at all
    # so we can bring back those previous marked to drop because of strand ambiguous
    snp[["keep"]][which(!strand_unambiguous)] <- TRUE
    snp[["keep"]][!(exact_match | snp[["sign_flip"]] | ID_match)] <- FALSE
  }

  if (remove_indels) {
    snp[["keep"]][snp[["INDEL"]]] <- FALSE
  }

  qc_summary <- matched %>%
    mutate(
      keep = snp[["keep"]],
      sign_flip = snp[["sign_flip"]],
      strand_flip = snp[["strand_flip"]],
      indel = snp[["INDEL"]]
    )

  # Apply allele flip if required
  if (!is.null(col_to_flip)) {
    if (!is.null(target_data_qced[, col_to_flip])) {
      target_data_qced[qc_summary$sign_flip, col_to_flip] <- -1 * target_data_qced[qc_summary$sign_flip, col_to_flip]
    } else {
      stop("Column '", col_to_flip, "' not found in target_data.")
    }
  }

  # Keep SNPs based on the 'keep' flag
  target_data_qced_matched <- target_data_qced[qc_summary$keep, , drop = FALSE]

  # Output the variants after strand flip if specified
  if (flip_strand) {
    strand_flipped_indices <- which(qc_summary$strand_flip)
    target_data_qced_matched[strand_flipped_indices, "A1"] <- strand_flip(target_data_qced_matched[strand_flipped_indices, "A1"])
    target_data_qced_matched[strand_flipped_indices, "A2"] <- strand_flip(target_data_qced_matched[strand_flipped_indices, "A2"])
  }

  # Remove duplicates if specified
  if (remove_dups) {
    dups <- vec_duplicate_detect(qc_summary[, c("chrom", "pos", "A1.target", "A2.target")])
    if (any(dups)) {
      target_data_qced_matched <- target_data_qced_matched[!dups, , drop = FALSE]
      message("Some duplicates were removed.")
    }
  }

  # Check if the minimum proportion of variants is matched
  min_match <- match_min_prop * min(nrow(target_variants), nrow(ref_variants))
  if (nrow(target_data_qced_matched) < min_match) {
    stop("Not enough variants have been matched.")
  }

  # Keep unmatched variants if specified
  if (!remove_unmatched) {
    unmatched_indices <- setdiff(seq_len(nrow(target_variants)), matched_indices$index)
    if (length(unmatched_indices) > 0) {
      target_data_qced_unmatched <- target_data[unmatched_indices, , drop = FALSE] %>%
        as.data.frame() %>%
        mutate(variant_id = paste(target_variants$chrom[unmatched_indices],
          target_variants$pos[unmatched_indices],
          target_variants$A2[unmatched_indices],
          target_variants$A1[unmatched_indices],
          sep = ":"
        ))
      target_data_qced <- rbind(target_data_qced_matched, target_data_qced_unmatched)
    } else {
      target_data_qced <- target_data_qced_matched
    }
  } else {
    target_data_qced <- target_data_qced_matched
  }

  # Change A1 and A2 so that it can fit the reference, and rearrange the columns so that the four are at the very front
  target_data_qced <- target_data_qced %>%
    separate(variant_id, into = c("chrom", "pos", "A2", "A1"), sep = ":", remove = FALSE) %>%
    select(chrom, pos, A1, A2, everything()) %>%
    mutate(chrom = as.integer(chrom), pos = as.integer(pos))

  if (remove_dups & any(duplicated(target_data_qced$variant_id))) {
    stop("In the input, there are duplicated varaints with different z scores. Please check the data and determine which to keep.")
  }

  return(list(target_data_qced = target_data_qced, qc_summary = qc_summary))
}

#' Align Variant Names
#'
#' This function aligns variant names from two strings containing variant names in the format of
#' "chr:pos:A1:A2" or "chr:pos_A1_A2". The first string should be the "source" and the second
#' should be the "reference".
#'
#' @param source A character vector of variant names in the format "chr:pos:A2:A1" or "chr:pos_A2_A1".
#' @param reference A character vector of variant names in the format "chr:pos:A2:A1" or "chr:pos_A2_A1".
#'
#' @return A list with two elements:
#' - aligned_variants: A character vector of aligned variant names.
#' - unmatched_indices: A vector of indices for the variants in the source that could not be matched.
#'
#' @examples
#' source <- c("1:123:A:C", "2:456:G:T", "3:789:C:A")
#' reference <- c("1:123:A:C", "2:456:T:G", "4:101:G:C")
#' align_variant_names(source, reference)
#'
#' @export
align_variant_names <- function(source, reference) {
  # Check if source and reference follow the expected pattern
  source_pattern <- grepl("^(chr)?[0-9]+:[0-9]+:[ATCG*]+:[ATCG*]+$|^(chr)?[0-9]+:[0-9]+_[ATCG*]+_[ATCG*]+$", source)
  reference_pattern <- grepl("^(chr)?[0-9]+:[0-9]+:[ATCG*]+:[ATCG*]+$|^(chr)?[0-9]+:[0-9]+_[ATCG*]+_[ATCG*]+$", reference)

  if (!all(source_pattern) && !all(reference_pattern)) {
    # Both source and reference do not follow the expected pattern
    warning("Cannot unify variant names because they do not follow the expected variant naming convention chr:pos:A2:A1 or chr:pos_A2_A1.")
    return(list(aligned_variants = source, unmatched_indices = integer(0)))
  }

  if ((!all(source_pattern) && all(reference_pattern)) || (all(source_pattern) && !all(reference_pattern))) {
    # One of source or reference follows the pattern, while the other does not
    stop("Source and reference have different variant naming conventions. They cannot be aligned.")
  }
  source_has_chr_prefix <- grepl("^chr", source[1])
  reference_has_chr_prefix <- grepl("^chr", reference[1])

  source_df <- variant_id_to_df(source)
  reference_df <- variant_id_to_df(reference)

  qc_result <- allele_qc(
    target_variants = source_df,
    ref_variants = reference_df,
    target_data = source_df,
    col_to_flip = NULL,
    match_min_prop = 0,
    remove_dups = FALSE,
    flip_strand = TRUE,
    remove_indels = FALSE,
    remove_strand_ambiguous = FALSE,
    remove_unmatched = FALSE
  )

  aligned_df <- qc_result$target_data_qced

  # Determine the output format based on the reference convention
  if (!grepl("_", reference[1])) {
    # Reference follows "chr:pos:A2:A1" convention
    output_format <- "colon"
  } else {
    # Reference follows "chr:pos_A2_A1" convention
    output_format <- "colon_underscore"
  }

  # Format the aligned variants based on the reference convention
  if (output_format == "colon") {
    aligned_variants <- apply(aligned_df[, c("chrom", "pos", "A2", "A1")], 1, paste, collapse = ":")
  } else {
    aligned_variants <- apply(aligned_df[, c("chrom", "pos", "A2", "A1")], 1, paste, collapse = ":")
    aligned_variants <- gsub(":(\\d+):(\\w):(\\w)", ":\\1_\\2_\\3", aligned_variants)
  }

  names(aligned_variants) <- NULL

  # Adjust the chr prefix in aligned_variants based on the reference convention
  aligned_variants_has_chr_prefix <- grepl("^chr", aligned_variants[1])
  if (reference_has_chr_prefix && !aligned_variants_has_chr_prefix) {
    aligned_variants <- paste0("chr", aligned_variants)
  } else if (!reference_has_chr_prefix && aligned_variants_has_chr_prefix) {
    aligned_variants <- sub("^chr", "", aligned_variants)
  }

  unmatched_indices <- which(match(aligned_variants, reference, nomatch = 0) == 0)

  list(
    aligned_variants = aligned_variants,
    unmatched_indices = unmatched_indices
  )
}
