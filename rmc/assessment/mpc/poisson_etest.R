## Code used to create rate ratio table of MPC scores
## in de novo variants from
## neurodevelopmental disorders (NDD) cases vs controls
## Designed to be run in RStudio
## (And manually adjust input/output paths)
# Load libraries
library(dplyr)
library(formattable)
library(purrr)

# Read input table
df <-read.table('/Users/kchao/Desktop/all_ndd_mpc_rate.tsv.bgz', header=T)
head(df)
# This should print something like this:
#   mpc_bin n_case n_control total_cases total_controls rate_per_case rate_per_control
# 1       0  17603      1007      37488          2179      0.469560        0.4621400
# 2       1   7295       350      37488          2179      0.194600        0.1606200
# 3       2   3719        77      37488          2179      0.099205        0.0353370
# 4       3   1418         8      37488          2179      0.037825        0.0036714

# Run two-sided Poisson exact test
df2 <- df %>%
  mutate(
    ptest=map2(
        c(case_count, control_count),
        c(total_case, total_control),
        poisson.test
      )
    )

# Extract rate ratio estimate and p-value from Poisson exact test
df <- df2 %>% mutate(p=map_dbl(ptest, ~ .x$p.value), estimate=map_dbl(ptest, ~ .x$estimate))
# Drop field with full Poisson test results
df <- subset(df, select=-ptest)
write.table(df, '/Users/kchao/Desktop/all_ndd_mpc_rate_ptest.tsv')

# Display slightly prettier table
formattable(df, align="l")
