#' Converted  Variant ID into a properly structured data frame
#' @param variant_id A data frame or character vector representing variant IDs.
#'   Expected formats are a data frame with columns "chrom", "pos", "A1", "A2",
#'   or a character vector in "chr:pos:A1:A2" or "chr:pos_A1_A2" format.
#' @return A data frame with columns "chrom", "pos", "A1", "A2", where 'chrom' 
#'   and 'pos' are integers, and 'A1' and 'A2' are allele identifiers.
#' @import stats
#' @export
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
#' @param remove Whether to remove strand SNPs (if any). Default is `TRUE`.
#' @return A single data frame with matched variants.
#' @importFrom dplyr mutate filter
#' @importFrom vctrs vec_duplicate_detect
#' @export
allele_qc = function(target_variants, ref_variants, target_data, col_to_flip, match.min.prop, remove_dups,flip,remove){

    target_variants = convert_to_dataframe(target_variants)
    ref_variants = convert_to_dataframe(ref_variants)
    
    matched <- as.data.frame(merge(target_variants, ref_variants,
                   by = c("chrom","pos"), all = FALSE, suffixes = c(".target", ".ref")))
    target_data_merged <- as.data.frame(target_data)[target_variants$chrom %in% matched$chrom & target_variants$pos %in% matched$pos,]

    a1 = toupper(matched$A1.target)
    a2 = toupper(matched$A2.target)
    ref1 = toupper(matched$A1.ref)
    ref2 = toupper(matched$A2.ref)
    # Strand flip, to change the allele representation in the 2nd data-set
	strand_flip = function(ref) {
		flip = ref
		flip[ref == "A"] = "T"
		flip[ref == "T"] = "A"
		flip[ref == "G"] = "C"
		flip[ref == "C"] = "G"
		return(flip)
	}
    flip1 = strand_flip(ref1)
	flip2 = strand_flip(ref2)
	snp = list()
    
	# Remove strand ambiguous SNPs (scenario 3)
 	snp[["keep"]] = !((a1=="A" & a2=="T") | (a1=="T" & a2=="A") | (a1=="C" & a2=="G") | (a1=="G" & a2=="C"))
 	# Remove non-ATCG coding
 	snp[["keep"]][ a1 != "A" & a1 != "T" & a1 != "G" & a1 != "C" ] = F
 	snp[["keep"]][ a2 != "A" & a2 != "T" & a2 != "G" & a2 != "C" ] = F
    
 	# as long as scenario 1 is involved, sign_flip will return TRUE
 	snp[["sign_flip"]] = (a1 == ref2 & a2 == ref1) | (a1 == flip2 & a2 == flip1)
    
 	# as long as scenario 2 is involved, strand_flip will return TRUE
 	snp[["strand_flip"]] = (a1 == flip1 & a2 == flip2) | (a1 == flip2 & a2 == flip1)
    
 	# remove other cases, eg, tri-allelic, one dataset is A C, the other is A G, for example.
 	exact_match = (a1 == ref1 & a2 == ref2) 
 	snp[["keep"]][!(exact_match | snp[["sign_flip"]] | snp[["strand_flip"]])] = F

    # Add the 'keep', 'sign_flip', and 'strand_flip' flags to the matched data
    matched = matched%>%
            mutate(keep = snp[["keep"]])%>%
            mutate(sign_flip= snp[["sign_flip"]])%>%
            mutate(strand_flip=snp[["strand_flip"]])
    
    # Apply allele flip if required
    if(flip) {
      if(!is.null(target_data_merged[,col_to_flip])){
         target_data_merged[matched$sign_flip,col_to_flip] = -1*target_data_merged[matched$sign_flip,col_to_flip]
      } else {
         stop("Column '", col_to_flip, "' not found in target_data.")
      }
    }
    
    # Keep SNPs based on the 'keep' flag
    if(remove) {
       target_data_merged = target_data_merged[matched$keep,]
    }

    # Remove duplicates if specified
    if (remove_dups) {
       dups <- vec_duplicate_detect(matched[, c("chrom", "pos","A1.target","A2.target")])
       if (any(dups)) {
          target_data_merged <- target_data_merged[!dups, ]
          message2("Some duplicates were removed.")
        }
     }

    # Check if the minimum proportion of variants is matched
    min_match <- match.min.prop * min(nrow(target_variants), nrow(ref_variants))
    if (nrow(target_data_merged) < min_match){
       stop("Not enough variants have been matched.")
     }
     return(target_data_merged)
}