# -----------------------------------------------------------------------------
# Scenario-specific validation: E1
# -----------------------------------------------------------------------------
# This script reproduces the benchmark calculation used for the E1 scenario.
# In E1, the prompt imposes:
# - a fixed target effect size of Cohen's d = 2.5
# - a selectively favourable variance input based on dataset D2 at day 14
#
# The purpose of this script is to make the scenario-specific benchmark
# transparent and reproducible. It is kept separate from the main reference
# pipeline because it validates one prompt scenario rather than the general
# dataset-day reference space.
#
# Note: this script assumes that the function get_comprehensive_n() has already
# been loaded from the main reference pipeline.
# -----------------------------------------------------------------------------

library(dplyr)

# Define the target effect size and the scenario-specific variance inputs
target_d <- 2.5
v1 <- 11.3
v2 <- 19.0

# Compute the pooled standard deviation implied by the two variances
s_p <- sqrt((v1 + v2) / 2)

# Convert Cohen's d into an absolute mean difference for the calculator
delta_for_test <- target_d * s_p

# Run the reference calculator under the E1 assumptions
res <- get_comprehensive_n(
  m_diff = delta_for_test,
  v1 = v1,
  v2 = v2,
  A_val = 0.96,  # approximately pnorm(2.5 / sqrt(2))
  alpha = 0.05,
  power = 0.8
)

# Print the rounded reference values used for the scenario validation
print(res %>% mutate(across(starts_with("n_"), ceiling)))