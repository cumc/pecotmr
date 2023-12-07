context("allele_qc")

create_allele_data <- function(seed, n=100, match_min_prop=0.8, ambiguous=FALSE, non_actg=FALSE, edge_cases=FALSE) {
  set.seed(seed)
  num_pass <- n*match_min_prop
  sumstat_A1 <- sample(c("A", "T", "G", "C"), num_pass, replace = TRUE)
  sumstat_A2 <- lapply(sumstat_A1, function(x) {
    if (x == "A") {
      return("T")
    } else if (x == "T") {
      return("A")
    } else if (x == "G") {
      return("C")
    } else if (x == "C") {
      return("G")
    }
  })

  if (ambiguous) {
    # Strand Ambiguous SNPs 
    sumstat_A1 <- c(sumstat_A1, sample(c("A", "T", "G", "C"), n-num_pass, replace = TRUE))
    sumstat_A2 <- c(sumstat_A2, lapply(sumstat_A1[n-num_pass+1:], function(x) {
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
  } else if (non_actg) {
    # Non-ATCG coding SNPs
    sumstat_A1 <- c(sumstat_A1, sample(c("ATG", "TAC", "GACA", "CTAA"), n-num_pass, replace = TRUE))
    sumstat_A2 <- c(sumstat_A2, lapply(sumstat_A1[n-num_pass+1:], function(x) {
      if (x == "ATG") {
        return("TAC")
      } else if (x == "TAC") {
        return("ATG")
      } else if (x == "GACA") {
        return("CTGT")
      } else if (x == "CTAA") {
        return("GATT")
      }
    }))
  }

  # Info SNPs
  info_A1 <- lapply(sumstat_A1[:num_pass], function(x) {
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
  info_A2 <- lapply(info_A1, function(x) {
    if (x == "A") {
      return("T")
    } else if (x == "T") {
      return("A")
    } else if (x == "G") {
      return("C")
    } else if (x == "C") {
      return("G")
    }
  })

  # Create the rest of the alleles
  info_A1 <- c(info_A1, sample(c("A", "T", "G", "C"), n-num_pass, replace = TRUE))
  info_A2 <- c(info_A2, lapply(info_A1[n-num_pass+1:], function(x) {
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

  chromosome <- rep(sample(1:20, 1), n)
  snp_positions <- sample(1:1000000, n)

  sumstats <- matrix(
    data = c(
      # chr, pos, A1, A2, beta, z
      chromosome,
      snp_positions,
      sumstat_A1,
      sumstat_A2,
      rnorm(n),
      rnorm(n)
    ),
    ncol = 6,
    byrow = TRUE
  )
  info_snp <- matrix(
    data = c(
      # chr, pos, A1, A2
      chromosome,
      snp_positions,
      info_A1,
      info_A2
    ),
    ncol = 4,
    byrow = TRUE
  )

  return(list(sumstats = sumstats, info_snp = info_snp))
}

test_that("Check that we correctly remove stand ambiguous SNPs",{
  res <- create_allele_data(1, n=100, match_min_prop=0.8, ambiguous=True)
  output <- allele_qc(res[1],res[2],0.2,TRUE,TRUE,TRUE) 
  expect_equal(nrow(output), 80)
})

test_that("Check that we correctly remove non-ACTG coding SNPs",{
  res <- create_allele_data(1, n=100, match_min_prop=0.4, non_actg=True)
  output <- allele_qc(res[1],res[2],0.2,TRUE,TRUE,TRUE) 
  expect_equal(nrow(output), 40)
})

test_that("Check that execution stops if not enough variants are matched",{
  res <- create_allele_data(1, n=100, match_min_prop=0.1, ambiguous=True)
  output <- allele_qc(res[1],res[2],0.2,TRUE,TRUE,TRUE) 
  expect_message(output, "Not enough variants have been matched.")
})