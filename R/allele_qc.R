#' Match alleles between summary statistics and SNP information.
#' 
#' Match by ("chr", "A1", "A2" and "pos"), accounting for possible
#' strand flips and major/minor allele flips (opposite effects and zscores).
#' 
#' @param sumstats A data frame with columns "chr", "pos", "A1", "A2", "beta" and "z". 
#' @param info_snp A data frame with columns "chr", "pos", "A1" and "A2".
#' @param match.min.prop  Minimum proportion of variants in the smallest data
#' to be matched, otherwise stops with an error. Default is 20\%.
#' @param remove_dups Whether to remove duplicates, default is True.
#' @param flip Whether the alleles must be flipped: A <--> T & C <--> G, in which case corresponding `$beta` and `$z` are multiplied by -1
#' Default is `TRUE`.
#' @param remove Whether to remove strand SNPs (if any)
#'
#' @return A single data frame with matched variants. Values in columns "beta" and "z"
#'   are multiplied by -1 for variants with alleles reversed
#' @importFrom vctrs vec_duplicate_detect
#' @importFrom Rlab message2
#' @export
#'
allele_qc = function(sumstats,info_snp,match.min.prop,remove_dups,flip,remove) {
  sumstats  =  as.data.frame(sumstats)
  info_snp  =  as.data.frame(info_snp)

  sumstats <- sumstats %>% mutate(chrom = ifelse(grepl("^chr", chrom), as.integer(sub("^chr", "", chrom)), chrom))
  
  matched <- merge(as.data.table(sumstats), as.data.table(info_snp),
                   by = c("chrom","pos"), all = FALSE, suffixes = c(".sumstats", ".ref"))

  a1 = toupper(matched$A1.sumstats)
  a2 = toupper(matched$A2.sumstats)
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
   
  matched = matched%>%
            mutate(keep = snp[["keep"]])%>%
            mutate(sign_flip= snp[["sign_flip"]])%>%
            mutate(strand_flip=snp[["strand_flip"]])
  if(flip) {
    if(!is.null(matched$beta)) {
    matched$beta[matched$sign_flip] = -1 * matched$beta[matched$sign_flip]
    }
    matched$beta[matched$sign_flip] = -1 * matched$beta[matched$sign_flip]
    matched$z[matched$sign_flip] = -1 * matched$z[matched$sign_flip]
    matched$A1.sumstats[matched$sign_flip] = matched$A1.ref[matched$sign_flip]
    matched$A2.sumstats[matched$sign_flip] = matched$A2.ref[matched$sign_flip]
  }  
  if(remove) {
    matched = matched[matched$keep,]
  }
  if (remove_dups) {
    dups <- vec_duplicate_detect(matched[, c("chrom", "pos","A1.sumstats","A2.sumstats")])
    if (any(dups)) {
      matched <- matched[!dups, ]
      message2("Some duplicates were removed.")
    }
  }
  #match.min.prop Minimum proportion of variants in the smallest data to be matched, otherwise stops with an error. Default is `20%`.
  min_match <- match.min.prop * min(nrow(sumstats), nrow(info_snp))
  if (nrow(matched) < min_match)
  stop("Not enough variants have been matched.")
  return(matched)
}