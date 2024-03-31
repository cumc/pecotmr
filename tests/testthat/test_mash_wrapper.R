test_that("merge_susie_cs merges credible sets correctly", {
  # Test case 1: No overlapping credible sets
  susie_fit_1 <- list(
    list(
      condition_1 = list(
        top_loci = data.frame(
          variant_id = c("variant1", "variant2"),
          pip = c(0.8, 0.6),
          cs_coverage_0.95 = c(1, 1)
        )
      ),
      condition_2 = list(
        top_loci = data.frame(
          variant_id = c("variant3", "variant4"),
          pip = c(0.9, 0.7),
          cs_coverage_0.95 = c(1, 2)
        )
      )
    )
  )
  
  expected_output_1 <- data.frame(
    variant_id = c("variant1", "variant2", "variant3", "variant4"),
    credible_set_names = c("cs_1_1", "cs_1_1", "cs_2_1", "cs_2_2"),
    max_pip = c(0.8, 0.6, 0.9, 0.7),
    median_pip = c(0.8, 0.6, 0.9, 0.7),
    stringsAsFactors = FALSE
  )
  
  expect_equal(merge_susie_cs(susie_fit_1), expected_output_1)
  
  # Test case 2: Overlapping credible sets
  susie_fit_2 <- list(
    list(
      condition_1 = list(
        top_loci = data.frame(
          variant_id = c("variant1", "variant2"),
          pip = c(0.8, 0.6),
          cs_coverage_0.95 = c(1, 1)
        )
      ),
      condition_2 = list(
        top_loci = data.frame(
          variant_id = c("variant2", "variant3"),
          pip = c(0.7, 0.9),
          cs_coverage_0.95 = c(2, 2)
        )
      )
    )
  )
  
  expected_output_2 <- data.frame(
    variant_id = c("variant1", "variant2", "variant3"),
    credible_set_names = c("cs_1_1,cs_2_2", "cs_1_1,cs_2_2", "cs_1_1,cs_2_2"),
    max_pip = c(0.8, 0.7, 0.9),
    median_pip = c(0.8, 0.65, 0.9),
    stringsAsFactors = FALSE
  )

  expect_equal(merge_susie_cs(susie_fit_2), expected_output_2)
  
  # Test case 3: Empty input
  susie_fit_3 <- list(condition_1 = list(top_loci = data.frame(
    variant_id = character(),
    credible_set_names = character(),
    max_pip = numeric(),
    median_pip = numeric(),
    stringsAsFactors = FALSE
  )))
  
  expected_output_3 <- NULL
  
  expect_equal(merge_susie_cs(susie_fit_3), expected_output_3)
  
  # Test case 4: Complementary parameter set to TRUE
  susie_fit_4 <- list(
    list(
      condition_1 = list(
        top_loci = data.frame(
          variant_id = c("variant1", "variant2"),
          pip = c(0.8, 0.6),
          cs_coverage_0.95 = c(1, 1)
        )
      ),
      condition_2 = list(
        top_loci = data.frame(
          variant_id = c("variant3", "variant4"),
          pip = c(0.9, 0.7),
          cs_coverage_0.95 = c(2, 2)
        )
      )
    )
  )
  
  expected_output_4 <- NULL
  
  expect_equal(merge_susie_cs(susie_fit_4, complementary = TRUE), expected_output_4)
  
  # Test case 5: Different coverage parameter
  susie_fit_5 <- list(
    list(
      condition_1 = list(
        top_loci = data.frame(
          variant_id = c("variant1", "variant2"),
          pip = c(0.8, 0.6),
          cs_coverage_0.90 = c(1, 1)
        )
      ),
      condition_2 = list(
        top_loci = data.frame(
          variant_id = c("variant3", "variant4"),
          pip = c(0.9, 0.7),
          cs_coverage_0.90 = c(2, 2)
        )
      )
    )
  )
  
  expected_output_5 <- data.frame(
    variant_id = c("variant1", "variant2", "variant3", "variant4"),
    credible_set_names = c("cs_1_1", "cs_1_1", "cs_2_2", "cs_2_2"),
    max_pip = c(0.8, 0.6, 0.9, 0.7),
    median_pip = c(0.8, 0.6, 0.9, 0.7),
    stringsAsFactors = FALSE
  )
  
  expect_equal(merge_susie_cs(susie_fit_5, coverage = "cs_coverage_0.90"), expected_output_5)
  
  # Test case 6: Multiple top_loci tables with mixed coverage indices
  susie_fit_6 <- list(
    list(
      condition_1 = list(
        top_loci = data.frame(
          variant_id = c("variant1", "variant2", "variant3"),
          pip = c(0.8, 0.6, 0.7),
          cs_coverage_0.95 = c(1, 1, 2)
        )
      ),
      condition_2 = list(
        top_loci = data.frame(
          variant_id = c("variant4", "variant5"),
          pip = c(0.9, 0.7),
          cs_coverage_0.95 = c(2, 3)
        )
      ),
      condition_3 = list(
        top_loci = data.frame(
          variant_id = c("variant6", "variant7", "variant8"),
          pip = c(0.85, 0.75, 0.8),
          cs_coverage_0.95 = c(1, 3, 2)
        )
      )
    )
  )
  
  expected_output_6 <- data.frame(
    variant_id = c("variant1", "variant2", "variant3", "variant4", "variant5", "variant6", "variant7", "variant8"),
    credible_set_names = c("cs_1_1", "cs_1_1", "cs_1_2", "cs_2_2", "cs_2_3", "cs_3_1", "cs_3_3", "cs_3_2"),
    max_pip = c(0.8, 0.6, 0.7, 0.9, 0.7, 0.85, 0.75, 0.8),
    median_pip = c(0.8, 0.6, 0.7, 0.9, 0.7, 0.85, 0.75, 0.8),
    stringsAsFactors = FALSE
  )
  
  expect_equal(merge_susie_cs(susie_fit_6), expected_output_6)
  
  # Test case 7: Multiple top_loci tables with overlapping sets and mixed coverage indices
  susie_fit_7 <- list(
    list(
      condition_1 = list(
        top_loci = data.frame(
          variant_id = c("variant1", "variant2", "variant3"),
          pip = c(0.8, 0.6, 0.7),
          cs_coverage_0.95 = c(1, 1, 2)
        )
      ),
      condition_2 = list(
        top_loci = data.frame(
          variant_id = c("variant2", "variant3", "variant4"),
          pip = c(0.7, 0.9, 0.85),
          cs_coverage_0.95 = c(2, 2, 1)
        )
      ),
      condition_3 = list(
        top_loci = data.frame(
          variant_id = c("variant4", "variant5"),
          pip = c(0.75, 0.8),
          cs_coverage_0.95 = c(3, 2)
        )
      )
    )
  )
  
  expected_output_7 <- data.frame(
    variant_id = c("variant1", "variant2", "variant3", "variant4", "variant5"),
    credible_set_names = c("cs_1_1,cs_1_2,cs_2_2", "cs_1_1,cs_1_2,cs_2_2", "cs_1_1,cs_1_2,cs_2_2","cs_2_1,cs_3_3", "cs_3_2"),
    max_pip = c(0.8, 0.7, 0.9, 0.85, 0.8),
    median_pip = c(0.8, 0.65, 0.8, 0.8, 0.8),
    stringsAsFactors = FALSE
  )
  
  expect_equal(merge_susie_cs(susie_fit_7), expected_output_7)
  
  # Test case 8: Multiple top_loci tables with different coverage indices and no overlapping sets
  susie_fit_8 <- list(
    list(
      condition_1 = list(
        top_loci = data.frame(
          variant_id = c("variant1", "variant2", "variant3"),
          pip = c(0.8, 0.6, 0.7),
          cs_coverage_0.95 = c(1, 2, 3)
        )
      ),
      condition_2 = list(
        top_loci = data.frame(
          variant_id = c("variant4", "variant5"),
          pip = c(0.9, 0.7),
          cs_coverage_0.95 = c(3, 1)
        )
      ),
      condition_3 = list(
        top_loci = data.frame(
          variant_id = c("variant6", "variant7", "variant8"),
          pip = c(0.85, 0.75, 0.8),
          cs_coverage_0.95 = c(2, 3, 1)
        )
      )
    )
  )
  
  expected_output_8 <- data.frame(
    variant_id = c("variant1", "variant2", "variant3", "variant4", "variant5", "variant6", "variant7", "variant8"),
    credible_set_names = c("cs_1_1", "cs_1_2", "cs_1_3", "cs_2_3", "cs_2_1", "cs_3_2", "cs_3_3", "cs_3_1"),
    max_pip = c(0.8, 0.6, 0.7, 0.9, 0.7, 0.85, 0.75, 0.8),
    median_pip = c(0.8, 0.6, 0.7, 0.9, 0.7, 0.85, 0.75, 0.8),
    stringsAsFactors = FALSE
  )
  
  expect_equal(merge_susie_cs(susie_fit_8), expected_output_8)
  
  # Test case 9: Single top_loci table with mixed coverage indices
  susie_fit_9 <- list(
    list(
      condition_1 = list(
        top_loci = data.frame(
          variant_id = c("variant1", "variant2", "variant3", "variant4", "variant5"),
          pip = c(0.8, 0.6, 0.7, 0.9, 0.85),
          cs_coverage_0.95 = c(1, 1, 2, 3, 2)
        )
      )
    )
  )
  
  expected_output_9 <- data.frame(
    variant_id = c("variant1", "variant2", "variant3", "variant5", "variant4"),
    credible_set_names = c("cs_1_1", "cs_1_1", "cs_1_2", "cs_1_2", "cs_1_3"),
    max_pip = c(0.8, 0.6, 0.7, 0.85, 0.9),
    median_pip = c(0.8, 0.6, 0.7, 0.85, 0.9),
    stringsAsFactors = FALSE
  )
  
  expect_equal(merge_susie_cs(susie_fit_9), expected_output_9)
  
  # Test case 10: Multiple top_loci tables with mixed coverage indices and overlapping sets
  susie_fit_10 <- list(
    list(
      condition_1 = list(
        top_loci = data.frame(
          variant_id = c("variant1", "variant2", "variant3"),
          pip = c(0.8, 0.6, 0.7),
          cs_coverage_0.95 = c(1, 2, 1)
        )
      ),
      condition_2 = list(
        top_loci = data.frame(
          variant_id = c("variant2", "variant4", "variant5"),
          pip = c(0.75, 0.9, 0.85),
          cs_coverage_0.95 = c(2, 1, 3)
        )
      ),
      condition_3 = list(
        top_loci = data.frame(
          variant_id = c("variant3", "variant5", "variant6"),
          pip = c(0.65, 0.8, 0.7),
          cs_coverage_0.95 = c(3, 2, 1)
        )
      )
    )
  )
  
  expected_output_10 <- data.frame(
    variant_id = c("variant1", "variant3", "variant2", "variant4", "variant5", "variant6"),
    credible_set_names = c("cs_1_1,cs_3_3", "cs_1_1,cs_3_3", "cs_1_2,cs_2_2", "cs_2_1", "cs_2_3,cs_3_2", "cs_3_1"),
    max_pip = c(0.8, 0.7, 0.75, 0.9, 0.85, 0.7),
    median_pip = c(0.8, 0.675, 0.675, 0.9, 0.825, 0.7),
    stringsAsFactors = FALSE
  )
  
  expect_equal(merge_susie_cs(susie_fit_10), expected_output_10)
})