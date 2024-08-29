#' Perform udr Analysis on Multigene Data
#'
#' @param combined_susie_list A list containing the combined SuSiE and summary statistics results.
#' @param coverage A numeric vector representing the coverage values.
#' @param independent_variant_list A list of independent variants to be used as a filter.
#' @param n_random An integer specifying the number of random samples to generate.
#' @param n_null An integer specifying the number of null samples to generate.
#' @param seed An integer specifying the seed for random number generation.
#' @param exclude_condition A character vector specifying conditions to be excluded from the analysis. Defaults to NULL.
#'
#' @return A list containing the results of the prior, or NULL if conditions are not met.
#' @import dplyr
#' @export
multigene_udr <- function(combined_susie_list, coverage, independent_variant_list, n_random, n_null, seed, exclude_condition = NULL) {
  # Default to an empty vector if exclude_condition is NULL
  if (is.null(exclude_condition)) {
    exclude_condition <- c()
  }
  reformat_data <- function(dat) {
    res <- list(
      strong.b = dat$strong$bhat,
      random.b = dat$random$bhat,
      null.b = dat$null$bhat,
      strong.s = dat$strong$sbhat,
      null.s = dat$null$sbhat,
      random.s = dat$random$sbhat
    )
    return(res)
  }
  # Load strong and random null summary statistics
  strong_file <- load_multitrait_R_sumstat(
    combined_susie_list$extracted_regional_window_combined_susie_result,
    combined_susie_list$extracted_regional_window_combined_sumstats_result,
    coverage,
    top_loci = TRUE,
    exclude_condition = exclude_condition
  )

  ran_null_file <- load_multitrait_R_sumstat(
    combined_susie_list$extracted_regional_window_combined_susie_result,
    combined_susie_list$extracted_regional_window_combined_sumstats_result,
    filter_file = independent_variant_list,
    exclude_condition = exclude_condition
  )

  # Generate random null samples
  ran_null <- mash_rand_null_sample(
    ran_null_file,
    n_random = n_random,
    n_null = n_null,
    exclude_condition = exclude_condition,
    seed = seed
  )


  # Prepare the strong summary statistics
  strong <- list(strong = list(bhat = strong_file$bhat, sbhat = strong_file$sbhat))

  # Combine strong and random null samples
  res <- c(strong, ran_null)

  # Reformat data for MASH analysis
  mash_input <- reformat_data(res)

  # Filter invalid summary statistics for each condition
  conditions <- c("strong", "random", "null")
  for (cond in conditions) {
    mash_input <- filter_invalid_summary_stat(
      mash_input,
      bhat = paste0(cond, ".b"),
      sbhat = paste0(cond, ".s"),
      z = TRUE,
      filter_by_missing_rate = NULL,
      sig_p_cutoff = NULL
    )
  }

  # Calculate ZtZ matrix
  mash_input$ZtZ <- t(as.matrix(mash_input$strong.z)) %*% as.matrix(mash_input$strong.z) / nrow(mash_input$strong.z)

  # Perform MASH analysis if conditions are met
  dd_prior <- if (nrow(mash_input$strong.b) < 2 || ncol(mash_input$strong.b) < 2) {
    NULL
  } else {
    mash_pipeline(mash_input, alpha = 1)
  }
  return(dd_prior)
}
