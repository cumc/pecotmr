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
#' @importFrom dplyr mutate inner_join filter pull select everything row_number if_else
#' @importFrom vctrs vec_duplicate_detect
#' @importFrom dolyr dplyr::if_else
#' @importFrom tidyr separate
#' @export
allele_qc <- function(target_variants, ref_variants, target_data, col_to_flip = NULL,
                      match_min_prop = 0.2, remove_dups = TRUE,
                      remove_indels = FALSE, remove_strand_ambiguous = TRUE,
                      flip_strand = FALSE, remove_unmatched = TRUE, remove_same_vars = FALSE, target_gwas = TRUE) {
    
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

    # check if the pattern is ATCG/DI
    check_ATCG <- function(vec) {
        pattern <- "^[ATCGDI]+$"

        # Function to check if a single element matches the pattern
        check_element <- function(element) {
            grepl(pattern, element)
        }

        result <- sapply(vec, check_element, USE.NAMES = FALSE)
        return(result)
    }
    
    # transform all inputs to dataframe
    target_variants <- variant_id_to_df(target_variants)
    ref_variants <- variant_id_to_df(ref_variants)
    if (isTRUE(target_gwas)) target_data <- variant_id_to_df(target_data)
    
    # remove the variant_id column in target data to avoid conflict
    if ("variant_id" %in% colnames(target_data)) {
        target_data <- select(target_data, -variant_id)
    }

    match_result <- merge(target_data, ref_variants, by = c("chrom", "pos"), all = FALSE, suffixes = c(".target", ".ref")) %>%
        # match target & ref by chrom and position
        as.data.frame() %>%
        mutate(variants_id_original = paste(chrom, pos, A2.target, A1.target, sep = ":")) %>%
        mutate(variants_id_qced = paste(chrom, pos, A2.ref, A1.ref, sep = ":")) %>%
        # filter out totally same rows.
        filter(duplicated(.) | !duplicated(.)) %>%
        # upper case target/reference A1 A2
        mutate(across(c(A1.target, A2.target, A1.ref, A2.ref), toupper)) %>%
        mutate(flip1.ref = strand_flip(A1.ref), flip2.ref = strand_flip(A2.ref)) %>%
        # these pairings are ambiguous: because we cannot tell it's an sign flip / strand flip
        mutate(strand_unambiguous = dplyr::if_else((A1.target == "A" & A2.target == "T") | (A1.target == "T" & A2.target == "A") |
                                            (A1.target == "C" & A2.target == "G") | (A1.target == "G" & A2.target == "C"), FALSE, TRUE)) %>%
        # filter out non-ATCG coded alleles
        mutate(non_ATCG = !(check_ATCG(A1.target) & check_ATCG(A2.target))) %>%
        # exact match should be kept all the time
        mutate(exact_match = A1.target == A1.ref & A2.target == A2.ref) %>%
        mutate(sign_flip = ((A1.target == A2.ref & A2.target == A1.ref) | (A1.target == flip2.ref & A2.target == flip1.ref)) & (A1.target != A1.ref & A2.target != A2.ref)) %>%
        mutate(strand_flip = ((A1.target == flip1.ref & A2.target == flip2.ref) | (A1.target == flip2.ref & A2.target == flip1.ref)) & (A1.target != A1.ref & A2.target != A2.ref)) %>%
        mutate(INDEL = (A2.target == "I" | A2.target == "D" | nchar(A2.target) > 1 | nchar(A1.target) > 1)) %>%
        mutate(ID_match = ((A2.target == "D" | A2.target == "I") & (nchar(A1.ref) > 1 | nchar(A2.ref) > 1)))
    
    # if not remove, then this should'nt be a condition to filter out any variants
    if (!remove_strand_ambiguous) {
        match_result <- match_result %>% mutate(strand_unambiguous = TRUE)
    }
    # if all strand flip is un-ambigous, then we know ambigous cases are indeed a strand flip
    # not a sign flip, then we infer there is no ambigous in the whole dataset, and keep those ambigous ones
    if (nrow(match_result %>% filter(strand_flip == TRUE) %>% filter(strand_unambiguous == TRUE)) == 0) {
        match_result <- match_result %>% mutate(strand_unambiguous = TRUE)
    }
    
    # To keep variants: if it's a strand flip, we will keep those unambigous (because if ambigous, cannot know it's trand / sign flip, so discard all)
    # or exact match or indel match (ID_match)
    # If not a strand flip, then we will keep those that are exact match / those are sign flip / INDEL matched
    match_result <- match_result %>% mutate(keep = dplyr::if_else(strand_flip, true = (strand_unambiguous | exact_match | ID_match), false = 
                                                           (exact_match | sign_flip | ID_match)))

    if (remove_indels) {
        match_result <- match_result %>% mutate(keep = dplyr::if_else(INDEL == FALSE, FALSE, TRUE))
    }
    
    # flip the signs of the column col_to_flip if there is a sign flip
    if (!is.null(col_to_flip)) {
        if (!is.null(match_result[, col_to_flip])) {
            match_result[match_result$sign_flip, col_to_flip] <- -1 * match_result[match_result$sign_flip, col_to_flip]
        } else {
            stop("Column '", col_to_flip, "' not found in target_data.")
        }
    }
    # flip the strands if there is a strand flip
    if (flip_strand) {
        strand_flipped_indices <- which(match_result$strand_flip)
        match_result[strand_flipped_indices, "A1.target"] <- strand_flip(match_result[strand_flipped_indices, "A1.target"])
        match_result[strand_flipped_indices, "A2.target"] <- strand_flip(match_result[strand_flipped_indices, "A2.target"])
    }
    
    # FIXME: I think this parameter is confusing. I inheritated directly from our function, whose default setting is TRUE.
    # It is removing all multi-allelic alleles which is unnecessary. I suggest remove this parameter directly. 
    # What we are trying to avoid is the SAME allele having diferent z score. I defined one parameter remove_same_vars later, but I can re-use this
    # remove_dup name
    if (remove_dups) {
        dups <- vec_duplicate_detect(match_result[, c("chrom", "pos", "variants_id_qced")])
        if (any(dups)) {
            match_result <- match_result[!dups, , drop = FALSE]
            message("Some duplicates were removed.")
        }
    }

    # Remove all unnecessary columns used to determine qc status
    # Finally keep those variants with FLAG keep = TRUE
    result <- match_result[match_result$keep, , drop = FALSE]
    
    
    result <- result %>% select(-(flip1.ref:keep)) %>% select(-variants_id_original, -A1.target, -A2.target) %>%
        rename(A1 = A1.ref, A2 = A2.ref, variant_id = variants_id_qced)
    
    # default FALSE, but if want to remove same variants having different z score, then set as TRUE
    if (remove_same_vars) {
        same_vars <- vec_duplicate_detect(result[, c("chrom", "pos", "variant_id")])
        if (any(same_vars)) {
            result <- result[!same_vars, , drop = FALSE]
            message("Same variants with different z scores are removed.")
        }
    }

    if (!remove_unmatched) {
        match_variant <- match_result %>% pull(variants_id_original)
        match_result <- select(match_result, -(flip1.ref:keep)) %>% select(-variants_id_original, -A1.target, -A2.target) %>%
            rename(A1 = A1.ref, A2 = A2.ref, variant_id = variants_id_qced)
        target_data <- target_data %>% mutate(variant_id = paste(chrom, pos, A2, A1, sep = ":"))
        if (length(setdiff(target_data %>% pull(variant_id), match_variant)) > 0) {
            unmatch_data <- target_data %>% filter(!variant_id %in% match_variant)
            result <- rbind(result, unmatch_data)
        }
    }
    
    min_match <- match_min_prop *  nrow(ref_variants)
        if (nrow(result) < min_match) {
        stop("Not enough variants have been matched.")
    }

    # throw an error if there are any variant_id that are duplicated (meaning that same variant having different other infos for example z score)
    if (!remove_same_vars & any(duplicated(result$variant_id))) {
        stop("In the input, there are duplicated variants with different z scores. Please check the data and determine which to keep.")
    }

    return(list(target_data_qced = result, qc_summary = match_result))
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
    remove_dups = TRUE,
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
