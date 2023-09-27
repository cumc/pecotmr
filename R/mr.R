######heterogeneity:  calculate I2 statistics based on the Cochran's Q statistic
calc_I2 = function(Q, Est) {
    Q = Q[[1]]
    Est = length(unique(Est))
    I2 = if (Q > 1e-3) (Q - Est + 1)/Q else 0
    return(if (I2 < 0) 0 else I2)
}

#####ptwas_est
fine_mr <- function (formatted_input, cpip_cutoff=0.5) {
output = formatted_input %>%
    mutate(
        beta_eQTL = beta_eQTL/se_eQTL,
        se_eQTL = 1) %>%
    group_by(gene_id, cluster) %>%
    mutate(spip = sum(pip)) %>%
    filter(spip >= cpip_cutoff) %>% # Cumulative cluster pip greater than a user defined cumulative pip threshold
    group_by(gene_id, cluster) %>%
    mutate(
        beta_yx = beta_GWAS/beta_eQTL,
        se_yx = sqrt(
            (se_GWAS^2/beta_eQTL^2) + ((beta_GWAS^2*se_eQTL^2)/beta_eQTL^4)),
        grp_beta = sum((beta_yx*pip)/spip),
        grp_se = sum((beta_yx^2 + se_yx^2)*pip/spip)) %>%
    mutate(
        grp_se = sqrt(grp_se - grp_beta^2),
        wv = grp_se^-2) %>%
    ungroup() %>%
    group_by(gene_id) %>%
    mutate(
        meta = sum(unique(wv) * unique(grp_beta)),
        sum_w = sum(unique(wv)),
        se_meta = sqrt(sum_w^-1),
        num_cluster = length(unique(cluster))) %>%
    mutate(
        num_instruments = length(snp),
        meta = meta/sum_w,
        ######sum(unique(wv)*(unique(grp_beta) - unique(meta))^2)
        Q = sum(unique(wv)*(unique(grp_beta) - unique(meta))^2),
        I2 = calc_I2(Q, grp_beta)) %>%
        ungroup() %>%
    distinct(gene_id, .keep_all = TRUE) %>%
    mutate(
        spip = round(spip, 3),
        grp_beta = round(grp_beta, 3),
        meta = round(meta, 3),
        se_meta = round(se_meta, 3),
        Q = round(Q, 3),
        I2 = round(I2, 3)) %>%
    select(gene_id, num_cluster, num_instruments, spip, grp_beta, grp_se, meta, se_meta, Q, I2)
    return(output)
}