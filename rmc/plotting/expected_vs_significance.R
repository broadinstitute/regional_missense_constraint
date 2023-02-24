library(ggplot2)

# Start with the following scenario:
# O/E = 0 on side A, O/E = 1 on side B
# Side A and side B have the same number of expected variants

# Calculate chi-square value in the above scenario with a specified
# expected number of variants in side A (equal to the expected in side B)
get_chisq_single_break = function(obs1, exp1, obs2, exp2) {
    return(
        -2 *
            (
                log(dpois(obs1, exp1 * (obs1 + obs2) / (exp1 + exp2))) +
                    log(dpois(obs2, exp2 * (obs1 + obs2) / (exp1 + exp2))) -
                    log(dpois(obs1, obs1)) -
                    log(dpois(obs2, obs2))
            )
    )
}

# Plot p-values of LRT in the above scenario, varying over
# expected number of variants in side A (equal to the expected in side B)
df_lrt = data.frame(exp = seq(1, 50))
df_lrt$chisq = sapply(df_lrt$exp, function(x) {
    get_chisq_single_break(0, x, x, x)
})
df_lrt$pval = sapply(df_lrt$chisq, function(x) 1 - pchisq(x, df = 1))
df_lrt$nlog10_pval = -log10(df_lrt$pval)

p = ggplot(df_lrt, aes(x = exp, y = nlog10_pval)) +
    geom_smooth(se = FALSE, method = "gam", formula = y ~ s(log(x))) +
    geom_point(size = 0.5) +
    ggtitle("Single break: O/E=0 vs. O/E=1\nEqual # of expected vars per section") +
    xlab("Expected number of variants per section") +
    ylab("-log10(p)") +
    scale_y_continuous(breaks = seq(0, round(max(df_lrt$nlog10_pval), 1), 2)) +
    theme_bw()
