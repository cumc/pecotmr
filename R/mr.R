#' heterogeneity:  calculate I2 statistics based on the Cochran's Q statistic
#' @noRd 
calc_I2 <- function(Q, Est) {
  Q <- Q[[1]]
  Est <- length(unique(Est))
  I2 <- if (Q > 1e-3) (Q - Est + 1) / Q else 0
  return(if (I2 < 0) 0 else I2)
}

#' MR Format Function
#'
#' Description of what the function does.
#'
#' @param susie_result A list containing the results of SuSiE analysis. This list should include nested elements such as 'susie_results', 'susie_result_trimmed', and 'top_loci', containing details about the statistical analysis of genetic variants.
#' @param condition A character string specifying the conditions. This is used to select the corresponding subset of results within 'susie_result'.
#' @param gwas_sumstats_db A data frame containing summary statistics from GWAS studies. It should include columns for variant id and their associated statistics such as beta coefficients and standard errors.
#' @param sets A character string indicating the method used to define sets of genetic variants. Defaults to "sets". This parameter is used to specify the type of sets to extract from the 'susie_result' object.
#' @param coverage A character string specifying the coverage threshold for credible sets, used when 'sets' is not "sets". Defaults to "coverage_0.95", indicating a 95% coverage credible set.
#' @param allele_qc Optional. A logical value indicating whether allele qc should be performed on the variants. When TRUE, allele qc are applied to the variants based on the GWAS summary statistics database ('gwas_sumstats_db').
#' @return A data frame formatted for MR analysis or NULL if cs_list is empty.
#' @export
mr_format <- function(susie_result, condition, gwas_sumstats_db, coverage = "cs_coverage_0.95", allele_qc = TRUE) {
  # Create null mr_format_input
  create_null_mr_input <- function(gene_name) {
    mr_format_input <- data.frame(
      gene_name = gene_name,
      variant_id = as.character(rep(NA, length(gene_name))),
      bhat_x = as.numeric(rep(NA, length(gene_name))),
      sbhat_x = as.numeric(rep(NA, length(gene_name))),
      cs = as.numeric(rep(NA, length(gene_name))),
      pip = as.numeric(rep(NA, length(gene_name))),
      bhat_y = as.numeric(rep(NA, length(gene_name))),
      sbhat_y = as.numeric(rep(NA, length(gene_name))),
      stringsAsFactors = FALSE # Optional, to prevent factors
    )
  }
  gene_name <- unique(get_nested_element(susie_result, c("susie_results", condition, "region_info", "region_name")))
  # Attempt to retrieve top_loci; if not found, return NULL
  top_loci <- tryCatch(
    {
      get_nested_element(susie_result, c("susie_results", condition, "top_loci"))
    },
    error = function() { # No parameter here
      message("top_loci does not exist for the specified condition in susie_result.")
      return(NULL)
    }
  )
  if (is.data.frame(top_loci)) {
    if (all(unique(get_nested_element(top_loci, coverage)) != 0)) {
      susie_cs_result_formatted <- top_loci %>%
        mutate(gene_name = gene_name) %>%
        filter(coverage >= 1) %>%
        mutate(variant = ifelse(grepl("^chr[0-9]+:", variant_id), gsub("^chr", "", variant_id), variant_id)) %>%
        select(gene_name, variant, betahat, sebetahat, all_of(coverage), pip) %>%
        rename("bhat_x" = "betahat", "sbhat_x" = "sebetahat", "cs" = coverage)
      if (allele_qc == TRUE) {
        susie_cs_result_formatted <- allele_qc(susie_cs_result_formatted$variant, gwas_sumstats_db$variant_id, susie_cs_result_formatted, c("bhat_x", "sbhat_x"))$target_data_qced[, c("gene_name", "variant_id", "bhat_x", "sbhat_x", "cs", "pip")]
      }

      susie_cs_gwas_variants_merge <- intersect(susie_cs_result_formatted$variant, gwas_sumstats_db$variant_id)

      mr_format_input <- susie_cs_result_formatted[match(susie_cs_gwas_variants_merge, susie_cs_result_formatted$variant), ] %>%
        cbind(., gwas_sumstats_db[match(susie_cs_gwas_variants_merge, gwas_sumstats_db$variant_id), ] %>%
          select(beta, se) %>%
          rename("bhat_y" = "beta", "sbhat_y" = "se"))
    } else {
      mr_format_input <- create_null_mr_input(gene_name)
    }
  } else {
    mr_format_input <- create_null_mr_input(gene_name)
  }
  return(mr_format_input)
}

#' Mendelian Randomization (MR)
#'
#' @param mr_formatted_input the output of twas_mr_format_input function
#' @param cpip_cutoff the threshold of cumulative posterior inclusion probability, default is 0.5
#' @return A single data frame of output with columns "gene_name", "num_CS", "num_IV",
#' "meta_eff", "se_meta_eff", "meta_pval", "Q", "Q_pval" and "I2". "gene_name" is ensemble ID. "num_CS" is the number of credible sets
#' contained in each gene, "num_IV" is the number of variants contained in each gene. "meta_eff", "se_meta_eff" and "meta_pval" are the MR estimate, standard error and pvalue.
#' "Q" is Cochran’s Q statistic, "I2" quantifies the heterogeneity, range from 0 to 1.
#' @importFrom dplyr mutate group_by filter ungroup distinct arrange select
#' @importFrom magrittr %>%
#' @importFrom stats pnorm pchisq
#' @export
mr_analysis <- function(mr_formatted_input, cpip_cutoff = 0.5) {
  create_null_output <- function(gene_name) {
    data.frame(
      gene_name = gene_name,
      num_CS = as.integer(rep(NA, length(gene_name))),
      num_IV = as.integer(rep(NA, length(gene_name))),
      cpip = as.numeric(rep(NA, length(gene_name))),
      meta_eff = as.numeric(rep(NA, length(gene_name))),
      se_meta_eff = as.numeric(rep(NA, length(gene_name))),
      meta_pval = as.numeric(rep(NA, length(gene_name))),
      Q = as.numeric(rep(NA, length(gene_name))),
      Q_pval = as.numeric(rep(NA, length(gene_name))),
      I2 = as.numeric(rep(NA, length(gene_name))),
      stringsAsFactors = FALSE
    )
  }
  if (all(is.na(mr_formatted_input[, -1]))) {
    return(create_null_output(unique(mr_formatted_input$gene_name)))
  }
  output <- mr_formatted_input %>%
    mutate(
      bhat_x = bhat_x / sbhat_x,
      sbhat_x = 1
    ) %>%
    group_by(gene_name, cs) %>%
    mutate(cpip = sum(pip)) %>%
    filter(cpip >= cpip_cutoff) # Cumulative pip greater than a user defined cumulative pip threshold

  if (dim(output)[1] != 0) {
    output <- output %>%
      group_by(gene_name, cs) %>%
      mutate(
        beta_yx = bhat_y / bhat_x,
        se_yx = sqrt((sbhat_y^2 / bhat_x^2) + ((bhat_y^2 * sbhat_x^2) / bhat_x^4)),
        composite_bhat = sum((beta_yx * pip) / cpip),
        composite_sbhat = sum((beta_yx^2 + se_yx^2) * pip / cpip)
      ) %>%
      mutate(
        composite_sbhat = sqrt(composite_sbhat - composite_bhat^2),
        wv = composite_sbhat^-2
      ) %>%
      ungroup() %>%
      mutate(
        meta_eff = sum(unique(wv) * unique(composite_bhat)),
        sum_w = sum(unique(wv)),
        se_meta_eff = sqrt(sum_w^-1),
        num_CS = length(unique(cs))
      ) %>%
      mutate(
        num_IV = length(variant_id),
        meta_eff = meta_eff / sum_w,
        meta_pval = 2 * pnorm(abs(meta_eff) / se_meta_eff, lower.tail = FALSE),
        Q = sum(unique(wv) * (unique(composite_bhat) - unique(meta_eff))^2),
        I2 = calc_I2(Q, composite_bhat),
        Q_pval = pchisq(Q, df = length(unique(composite_bhat)) - 1, lower = F)
      ) %>%
      ungroup() %>%
      distinct(gene_name, .keep_all = TRUE) %>%
      mutate(
        cpip = round(cpip, 3),
        meta_pval = round(meta_pval, 3),
        meta_eff = round(meta_eff, 3),
        se_meta_eff = round(se_meta_eff, 3),
        Q = round(Q, 3),
        Q_pval = round(Q_pval, 3),
        I2 = round(I2, 3)
      ) %>%
      arrange(meta_pval) %>%
      select(gene_name, num_CS, num_IV, cpip, meta_eff, se_meta_eff, meta_pval, Q, Q_pval, I2)
  } else {
    return(create_null_output(unique(mr_formatted_input$gene_name)))
  }
  output
}
