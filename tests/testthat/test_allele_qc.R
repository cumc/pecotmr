context("allele_qc")
library(data.table)

create_allele_data <- function(seed, n=100, match_min_prop=0.8, ambiguous=FALSE, non_actg=FALSE, edge_cases=FALSE) {
  set.seed(seed)
  num_pass <- n*match_min_prop
  sumstat_A1 <- sample(c("A", "T", "G", "C"), num_pass, replace = TRUE)
  sumstat_A2 <- unlist(lapply(sumstat_A1, function(x) {
    if (x == "A") {
      return(sample(c("G", "C"), 1))
    } else if (x == "T") {
      return(sample(c("G", "C"), 1))
    } else if (x == "G") {
      return(sample(c("A", "T"), 1))
    } else if (x == "C") {
      return(sample(c("A", "T"), 1))
    }
  }))

  if (ambiguous) {
    # Strand Ambiguous SNPs 
    sumstat_A1 <- c(sumstat_A1, sample(c("A", "T", "G", "C"), n-num_pass, replace = TRUE))
    sumstat_A2 <- unlist(c(sumstat_A2, lapply(sumstat_A1[(num_pass+1):length(sumstat_A1)], function(x) {
      if (x == "A") {
        return("T")
      } else if (x == "T") {
        return("A")
      } else if (x == "G") {
        return("C")
      } else if (x == "C") {
        return("G")
      }
    })))
  } else if (non_actg) {
    # Non-ATCG coding SNPs
    sumstat_A1 <- c(sumstat_A1, sample(c("ATG", "TAC", "GACA", "CTAA"), n-num_pass, replace = TRUE))
    sumstat_A2 <- unlist(c(sumstat_A2, lapply(sumstat_A1[(num_pass+1):length(sumstat_A1)], function(x) {
      if (x == "ATG") {
        return("TAC")
      } else if (x == "TAC") {
        return("ATG")
      } else if (x == "GACA") {
        return("CTGT")
      } else if (x == "CTAA") {
        return("GATT")
      }
    })))
  }

  # Info SNPs
  info_A1 <- lapply(sumstat_A1[1:num_pass], function(x) {
    if(runif(1) < 0.2) {
      # flip a small proportion of the alleles
      if (x == "A") {
        return("T")
      } else if (x == "T") {
        return("A")
      } else if (x == "G") {
        return("C")
      } else if (x == "C") {
        return("G")
      }
    } else {
      return(x)
    }
  })
  info_A2 <- sumstat_A2[1:num_pass]
  # Handle random flips
  info_A2[info_A1 != sumstat_A1[1:num_pass]] <- unlist(lapply(info_A2[info_A1 != sumstat_A1[1:num_pass]], function(x) {
    if (x == "A") {
      return("T")
    } else if (x == "T") {
      return("A")
    } else if (x == "G") {
      return("C")
    } else if (x == "C") {
      return("G")
    }
  }))

  # Create the rest of the alleles
  info_A1 <- unlist(c(info_A1, sample(c("A", "T", "G", "C"), n-num_pass, replace = TRUE)))
  info_A2 <- unlist(c(info_A2, lapply(info_A1[(num_pass+1):length(info_A1)], function(x) {
    if (x == "A") {
      return(sample(c("G", "C"), 1))
    } else if (x == "T") {
      return(sample(c("G", "C"), 1))
    } else if (x == "G") {
      return(sample(c("A", "T"), 1))
    } else if (x == "C") {
      return(sample(c("A", "T"), 1))
    }
  })))

  chromosome <- unlist(rep(sample(1:20, 1), n))
  snp_positions <- sample(1:1000000, n)
  target_variants <- data.frame(
    chrom = chromosome,
    pos = snp_positions,
    A1 = sumstat_A1,
    A2 = sumstat_A2
  )
  ref_variants <- data.frame(
    chrom = chromosome,
    pos = snp_positions,
    A1 = info_A1,
    A2 = info_A2
  )
  target_data <- data.frame(
    chrom = chromosome,
    pos = snp_positions,
    A1 = sumstat_A1,
    A2 = sumstat_A2,
    beta = rnorm(n),
    z = rnorm(n)
  )

  return(list(target_data = target_data, target_variants = target_variants, ref_variants = ref_variants))
}

test_that("Check that we correctly remove stand ambiguous SNPs",{
  res <- create_allele_data(1, n=100, match_min_prop=0.8, ambiguous=TRUE)
  output <- allele_qc(
    res$target_variants, res$ref_variants, res$target_data, "beta", match.min.prop = 0.2,
    TRUE, TRUE, FALSE, TRUE)
  expect_equal(nrow(output$target_data_qced), 80)
})

test_that("Check that we correctly remove non-ACTG coding SNPs",{
  res <- create_allele_data(1, n=100, match_min_prop=0.4, non_actg=TRUE)
  output <- allele_qc(
    res$target_variants, res$ref_variants, res$target_data, "beta", match.min.prop = 0.2,
    TRUE, TRUE, FALSE, TRUE)
  expect_equal(nrow(output$target_data_qced), 40)
})

test_that("Check that execution stops if not enough variants are matched",{
  res <- create_allele_data(1, n=100, match_min_prop=0.1, ambiguous=TRUE)
  expect_error(allele_qc(
    res$target_variants, res$ref_variants, res$target_data, "beta", match.min.prop = 0.2,
    TRUE, TRUE, FALSE, TRUE), "Not enough variants have been matched.")
})