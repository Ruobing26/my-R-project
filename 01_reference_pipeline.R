# -----------------------------------------------------------------------------
# Reference pipeline for Chapters 3 and 4
# -----------------------------------------------------------------------------
# This script prepares the anonymised pilot datasets used in the thesis,
# computes the reference sample size table for the benchmark analyses, and
# builds an extended reproducibility space for alternative model outputs.
#
# The script contains three parts:
# 1. Main reference calculations based on absolute tumour area by dataset and day
# 2. Auxiliary pairwise day-to-day contrasts prepared in advance for broader
#    reproducibility checks
# 3. Additional merged-dataset calculations to benchmark model-side pooling
# -----------------------------------------------------------------------------

library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(pwr)
library(MKpower)


# Import anonymised pilot datasets
ds1 <- read_csv("public_D1.csv")
ds2 <- read_csv("public_D2.csv")
ds3 <- read_csv("public_D3.csv")
ds2_3 <- read_csv("public_merged_sorted.csv")

# Define the string patterns used to identify the two treatment groups
TREAT_KEY <- "Therapie_A"
CTRL_KEY  <- "Therapie_B"

# -----------------------------------------------------------------------------
# Dataset import and harmonisation
# -----------------------------------------------------------------------------
# Convert each raw pilot dataset from wide to long format and create a common
# structure with dataset, day, mouse ID, treatment group, and tumour area.
# -----------------------------------------------------------------------------

clean_pilot_dataset <- function(raw_data, dataset_id) {
  
  # Harmonise column names and rename the first column to "day"
  names(raw_data) <- make.unique(names(raw_data), sep = "_")
  raw_data <- raw_data %>% rename(day = 1)
  
  # Reshape the dataset from wide to long format
  raw_data %>%
    pivot_longer(
      cols = -day,
      names_to = "original_col",
      values_to = "tumor"
    ) %>%
    mutate(
      dataset = dataset_id,
      
      # Assign each original column to treatment group A or B
      group = case_when(
        str_detect(original_col, TREAT_KEY) ~ "A",
        str_detect(original_col, CTRL_KEY)  ~ "B",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(group), !is.na(tumor)) %>%
    
    # Within each dataset, day, and group, assign running mouse IDs
    group_by(dataset, day, group) %>%
    mutate(mouse_id = row_number()) %>% 
    ungroup() %>%
    select(dataset, day, mouse_id, group, tumor) %>%
    arrange(dataset, day, group, mouse_id)
}

# Apply the harmonisation function to the three original pilot datasets
ds1_clean <- clean_pilot_dataset(ds1, "D1")
ds2_clean <- clean_pilot_dataset(ds2, "D2")
ds3_clean <- clean_pilot_dataset(ds3, "D3")
pilot_clean <- bind_rows(ds1_clean, ds2_clean, ds3_clean)

# -----------------------------------------------------------------------------
# Complete-day restriction
# -----------------------------------------------------------------------------
# For the prompt input and the main reference calculations, only dataset-day
# combinations with complete 1:1 pilot observations are retained.
# This means that a day is kept only if all 5 animals in group A and all 5
# animals in group B are still observed on that day.
# -----------------------------------------------------------------------------
pilot_clean <- pilot_clean %>%
  group_by(dataset, day) %>%
  filter(
    sum(group == "A") == 5 & 
      sum(group == "B") == 5
  ) %>%
  ungroup()

# -----------------------------------------------------------------------------
# Probabilistic effect size for the Wilcoxon-Mann-Whitney reference calculation
# -----------------------------------------------------------------------------
# This function computes the pilot-based probabilistic effect used in the
# Noether approximation for the WMW sample size calculation.
# -----------------------------------------------------------------------------
estimate_A <- function(x_A, x_B) {
  x_A <- x_A[is.finite(x_A)]
  x_B <- x_B[is.finite(x_B)]
  if (length(x_A) == 0 || length(x_B) == 0) return(NA_real_)
  
  # Compare all pairwise combinations between group A and group B
  cmp <- outer(x_A, x_B, FUN = "-")
  (sum(cmp > 0) + 0.5 * sum(cmp == 0)) / (length(x_A) * length(x_B))
}

# -----------------------------------------------------------------------------
# Auxiliary reference space: pairwise day-to-day contrasts
# -----------------------------------------------------------------------------
# This section constructs all possible day-pair contrasts. These calculations
# were prepared as part of the broader reproducibility framework, so that model
# outputs based on derived contrasts can also be checked against a pre-specified
# reference space if needed.
# -----------------------------------------------------------------------------
all_days <- sort(unique(pilot_clean$day))
day_combos <- combn(all_days, 2)

cfb_all_pairs <- list()
for(i in 1:ncol(day_combos)) {
  d_start <- day_combos[1, i]
  d_end   <- day_combos[2, i]
  
  cfb_all_pairs[[i]] <- pilot_clean %>%
    filter(day %in% c(d_start, d_end)) %>%
    group_by(dataset, mouse_id, group) %>%
    filter(n() == 2) %>%  # retain only mice observed at both time points
    summarise(diff = tumor[day == d_end] - tumor[day == d_start], .groups = "drop") %>%
    mutate(time_pair = paste0("D", d_end, "-D", d_start))
}

cfb_matrix <- bind_rows(cfb_all_pairs)

# -----------------------------------------------------------------------------
# Core sample size calculator for benchmark comparisons
# -----------------------------------------------------------------------------
# For a given effect estimate and variance input, this function computes
# reference sample sizes for:
# - pooled t
# - Welch t
# - Wilcoxon-Mann-Whitney (Noether approximation)
# - pooled z
# - Welch z
# Each calculation is returned for both two-sided and one-sided settings.
# -----------------------------------------------------------------------------
get_comprehensive_n <- function(m_diff, v1, v2, A_val, alpha = 0.05, power = 0.8) {
  delta_abs <- abs(m_diff)
  
  # Pooled variance is simplified because the pilot datasets satisfy n1 = n2
  s_pooled  <- sqrt((v1 + v2) / 2)
  
  calc_by_side <- function(label) {
    
    # Map the side specification to the syntax expected by pwr and MKpower
    alt_pwr <- ifelse(label == "Two-Sided", "two.sided", "greater")
    alt_mk  <- ifelse(label == "Two-Sided", "two.sided", "one.sided")
    
    alpha_t <- alpha
    
    tryCatch({
      # Pooled t reference calculation
      n_c <- pwr.t.test(
        d = delta_abs / s_pooled,
        sig.level = alpha_t,
        power = power,
        type = "two.sample",
        alternative = alt_pwr
      )$n
      
      # Welch t reference calculation
      n_w <- power.welch.t.test(
        delta = delta_abs,
        sd1 = sqrt(v1), sd2 = sqrt(v2),
        sig.level = alpha_t,
        power = power,
        alternative = alt_mk
      )$n
      
      # WMW reference calculation using the Noether approximation
      sided <- ifelse(label == "Two-Sided", 2, 1)
      n_m <- (qnorm(1 - alpha / sided) + qnorm(power))^2 / (6 * (A_val - 0.5)^2)
      
      # Welch z normal approximation
      z_alpha <- if (label == "Two-Sided") qnorm(1 - alpha / 2) else qnorm(1 - alpha)
      z_beta  <- qnorm(power)  # since power = 1 - beta
      n_z <- (v1 + v2) * (z_alpha + z_beta)^2 / (delta_abs^2)
      
      # Pooled z normal approximation
      n_pooled_z <- 2 * (s_pooled^2) * (z_alpha + z_beta)^2 / (delta_abs^2)
      
      data.frame(
        Side = label,
        n_pooled_t = n_c,
        n_Welch_t = n_w,
        n_WMW = n_m,
        n_welch_z = n_z,
        n_pooled_z = n_pooled_z
      )
    }, error = function(e) {
      data.frame(
        Side = label,
        n_pooled_t = NA,
        n_Welch_t = NA,
        n_WMW = NA,
        n_welch_z = NA,
        n_pooled_z = NA
      )
    })
  }
  
  bind_rows(calc_by_side("Two-Sided"), calc_by_side("One-Sided"))
}

# -----------------------------------------------------------------------------
# Main reference calculations based on absolute tumour area
# -----------------------------------------------------------------------------
# For each dataset-day combination, compute group means, group variances,
# the mean difference, the pilot-based WMW effect, and the corresponding
# reference sample sizes.
# -----------------------------------------------------------------------------
final_results <- pilot_clean %>%
  group_by(dataset, day) %>%
  summarise(
    m_A = mean(tumor[group == "A"], na.rm = TRUE),
    m_B = mean(tumor[group == "B"], na.rm = TRUE),
    v_A = var(tumor[group == "A"], na.rm = TRUE),
    v_B = var(tumor[group == "B"], na.rm = TRUE),
    delta = m_A - m_B,
    A_val = estimate_A(tumor[group == "A"], tumor[group == "B"]),
    .groups = "drop"
  ) %>%
  filter(
    is.finite(delta), is.finite(v_A), is.finite(v_B),
    abs(delta) > 0, v_A > 0, v_B > 0,
    is.finite(A_val),
    abs(A_val - 0.5) > 1e-6
  ) %>%
  rowwise() %>%
  do({
    res <- get_comprehensive_n(.$delta, .$v_A, .$v_B, .$A_val)
    data.frame(., res)
  }) %>%
  ungroup()

# Apply the same reference calculations to the auxiliary day-to-day contrasts
final_cfb_analysis <- cfb_matrix %>%
  group_by(dataset, time_pair) %>%
  summarise(
    m_A = mean(diff[group == "A"], na.rm = TRUE),
    m_B = mean(diff[group == "B"], na.rm = TRUE),
    v_A = var(diff[group == "A"], na.rm = TRUE),
    v_B = var(diff[group == "B"], na.rm = TRUE),
    delta = m_A - m_B,
    A_val = estimate_A(diff[group == "A"], diff[group == "B"]),
    .groups = "drop"
  ) %>%
  filter(
    is.finite(delta), is.finite(v_A), is.finite(v_B),
    abs(delta) > 0, v_A > 0, v_B > 0,
    is.finite(A_val),
    abs(A_val - 0.5) > 1e-6
  ) %>%
  rowwise() %>%
  do({
    res <- get_comprehensive_n(.$delta, .$v_A, .$v_B, .$A_val)
    data.frame(., res)
  }) %>%
  ungroup() %>%
  rename(time_point = time_pair)

# -----------------------------------------------------------------------------
# Additional merged dataset for model-side pooling checks
# -----------------------------------------------------------------------------
# The merged D2+D3 dataset is not part of the standard prompt input.
# It is included here only to benchmark model-side aggregation decisions
# against a transparent reference calculation.
# -----------------------------------------------------------------------------
ds_D2D3_clean <- clean_pilot_dataset(ds2_3, "D2+D3")
ds_D2D3_clean <- ds_D2D3_clean %>%
  group_by(dataset, day) %>%
  filter(
    sum(group == "A") == 10 & 
      sum(group == "B") == 10
  ) %>%
  ungroup()

# Combine the original pilot datasets with the merged D2+D3 dataset
pilot_plus <- bind_rows(pilot_clean, ds_D2D3_clean)
pilot_plus <- pilot_plus %>%
  mutate(dataset = factor(dataset, levels = c("D1", "D2", "D3", "D2+D3")))

# Compute the reference table for the merged D2+D3 dataset
B <- ds_D2D3_clean %>%
  group_by(dataset, day) %>%
  summarise(
    m_A = mean(tumor[group == "A"], na.rm = TRUE),
    m_B = mean(tumor[group == "B"], na.rm = TRUE),
    v_A = var(tumor[group == "A"], na.rm = TRUE),
    v_B = var(tumor[group == "B"], na.rm = TRUE),
    delta = m_A - m_B,
    A_val = estimate_A(tumor[group == "A"], tumor[group == "B"]),
    .groups = "drop"
  ) %>%
  filter(
    is.finite(delta), is.finite(v_A), is.finite(v_B),
    abs(delta) > 0, v_A > 0, v_B > 0,
    is.finite(A_val),
    abs(A_val - 0.5) > 1e-6
  ) %>%
  rowwise() %>%
  do({
    res <- get_comprehensive_n(.$delta, .$v_A, .$v_B, .$A_val)
    data.frame(., res)
  }) %>%
  ungroup()

# -----------------------------------------------------------------------------
# Final reference planning table
# -----------------------------------------------------------------------------
# Combine the main absolute-endpoint reference table, the merged D2+D3 table,
# and the auxiliary day-to-day contrast table into one master object.
# All sample sizes are rounded up to the next integer.
# -----------------------------------------------------------------------------
master_planning_table <- bind_rows(
  final_results %>% rename(time_point = day) %>% mutate(time_point = as.character(time_point)),
  B %>% rename(time_point = day) %>% mutate(time_point = as.character(time_point)),
  final_cfb_analysis
) %>%
  mutate(across(starts_with("n_"), ceiling)) %>%
  select(dataset, time_point, Side, m_A, m_B, v_A, v_B, A_val, everything())

# Export the final reference planning table
write.csv(master_planning_table, file = "Referentable.csv", row.names = FALSE)