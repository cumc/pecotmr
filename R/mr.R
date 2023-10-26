# heterogeneity:  calculate I2 statistics based on the Cochran's Q statistic
calc_I2 = function(Q, Est) {
    Q = Q[[1]]
    Est = length(unique(Est))
    I2 = if (Q > 1e-3) (Q - Est + 1)/Q else 0
    return(if (I2 < 0) 0 else I2)
}

#' @export
fine_mr <- function(formatted_input, cpip_cutoff=0.5) {
output = formatted_input %>%
    mutate(
        bhat_x = bhat_x/sbhat_x,
        sbhat_x = 1) %>%
    group_by(X_ID, cs) %>%
    mutate(cpip = sum(pip)) %>%
    filter(cpip >= cpip_cutoff) %>% # Cumulative pip greater than a user defined cumulative pip threshold
    group_by(X_ID, cs) %>%
    mutate(
        beta_yx = bhat_y/bhat_x,
        se_yx = sqrt(
            (sbhat_y^2/bhat_x^2) + ((bhat_y^2*sbhat_x^2)/bhat_x^4)),
        composite_bhat = sum((beta_yx*pip)/cpip),
        composite_sbhat = sum((beta_yx^2 + se_yx^2)*pip/cpip)) %>%
    mutate(
        composite_sbhat = sqrt(composite_sbhat - composite_bhat^2),
        wv = composite_sbhat^-2) %>%
    ungroup() %>%
    group_by(X_ID) %>%
    mutate(
        meta_eff = sum(unique(wv) * unique(composite_bhat)),
        sum_w = sum(unique(wv)),
        se_meta_eff = sqrt(sum_w^-1),
        num_CS = length(unique(cs))) %>%
    mutate(
        num_IV = length(snp),
        meta_eff = meta_eff/sum_w,
        ######sum(unique(wv)*(unique(composite_bhat) - unique(meta_eff))^2)
        Q = sum(unique(wv)*(unique(composite_bhat) - unique(meta_eff))^2),
        I2 = calc_I2(Q, composite_bhat)) %>%
        ungroup() %>%
    distinct(X_ID, .keep_all = TRUE) %>%
    mutate(
        cpip = round(cpip, 3),
        composite_bhat = round(composite_bhat, 3),
        meta_eff = round(meta_eff, 3),
        se_meta_eff = round(se_meta_eff, 3),
        Q = round(Q, 3),
        I2 = round(I2, 3)) %>%
    select(X_ID, num_CS, num_IV, cpip, composite_bhat, composite_sbhat, meta_eff, se_meta_eff, Q, I2)
    return(output)
}