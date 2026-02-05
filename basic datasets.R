library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(pwr)
library(MKpower)

ds1 <- read_csv("public_D1.csv")
ds2 <- read_csv("public_D2.csv")
ds3 <- read_csv("public_D3.csv")
ds2_3 <- read_csv("public_merged_sorted.csv")

TREAT_KEY <- "Therapie_A"
CTRL_KEY  <- "Therapie_B"

clean_pilot_dataset <- function(raw_data, dataset_id) {
  
  names(raw_data) <- make.unique(names(raw_data), sep = "_")
  raw_data <- raw_data %>% rename(day = 1)
  
  # Convert to a long table
  raw_data %>%
    pivot_longer(
      cols = -day,
      names_to = "original_col",
      values_to = "tumor"
    ) %>%
    mutate(
      dataset = dataset_id,
      # Identify the group 
      group = case_when(
        str_detect(original_col, TREAT_KEY) ~ "A",
        str_detect(original_col, CTRL_KEY)  ~ "B",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(group), !is.na(tumor)) %>%
    # In each group and on each day, number 1 to 5 was assigned to the mice
    group_by(dataset, day, group) %>%
    mutate(mouse_id = row_number()) %>% # 确保每个组都从 1 开始
    ungroup() %>%
    select(dataset, day, mouse_id, group, tumor) %>%
    arrange(dataset, day, group, mouse_id)
}


ds1_clean <- clean_pilot_dataset(ds1, "D1")
ds2_clean <- clean_pilot_dataset(ds2, "D2")
ds3_clean <- clean_pilot_dataset(ds3, "D3")
pilot_clean <- bind_rows(ds1_clean, ds2_clean, ds3_clean)

# --- Ensure that there are 5 mice in Group A/B every day ---
pilot_clean <- pilot_clean %>%
  group_by(dataset, day) %>%
  filter(
    sum(group == "A") == 5 & 
      sum(group == "B") == 5
  ) %>%
  ungroup()

## Visualize the distribution of source data (as the basis for screening data)

P1 <- ggplot(pilot_clean, aes(x = factor(day), y = tumor, fill = group)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  # Superimpose scattered points to clearly identify isolated mouse individuals
  geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
  # Divide by dataset
  facet_wrap(~dataset, scales = "free_x") + 
  theme_bw() +
  labs(title = "Growth Dynamics by Dataset",
       subtitle = "Recommended layout for quality control and data filtering",
       x = "Days Post-Implantation",
       y = "Tumor Volume (mm³)")

compute_params_long <- function(data_long) {
  # Calculate the basic statistics for each day and each group
  # Group by group and the number of days
  stats_summary <- data_long %>%
    group_by(day, group) %>%
    summarise(
      n    = n(),
      mean = mean(tumor, na.rm = TRUE), 
      var  = var(tumor, na.rm = TRUE),
      .groups = "drop"
    )
  
  # The long table was split into the A and the B
  s_a <- stats_summary %>% filter(group == "A")
  s_b <- stats_summary %>% filter(group == "B")
  
  # Security Check: Ensure that the number of days is exactly matched
  common_days <- intersect(s_a$day, s_b$day)
  s_a <- s_a %>% filter(day %in% common_days)
  s_b <- s_b %>% filter(day %in% common_days)
  
  # Build a wide table and calculate Delta and combined variance
  tibble(
    day        = common_days,
    n_A        = s_a$n,
    n_B        = s_b$n,
    mean_A     = s_a$mean,
    mean_B     = s_b$mean,
    delta      = s_b$mean - s_a$mean, 
    var_A      = s_a$var,
    var_B      = s_b$var,
    pooled_var = ((s_a$n - 1) * s_a$var + (s_b$n - 1) * s_b$var) / 
      (s_a$n + s_b$n - 2)
  ) %>%
    arrange(day)
}

res_D1 <- pilot_clean %>%
  filter(dataset == "D1") %>%
  compute_params_long() %>%
  mutate(datasets = "D1", .before = 1)

res_D2 <- pilot_clean %>%
  filter(dataset == "D2") %>%
  compute_params_long() %>%
  mutate(datasets = "D2", .before = 1)

res_D3 <- pilot_clean %>%
  filter(dataset == "D3") %>%
  compute_params_long() %>%
  mutate(datasets = "D3", .before = 1)

# ---------------------------------------
# Combine datasets by row-binding
# ---------------------------------------

results_all <- bind_rows(
  res_D1,
  res_D2,
  res_D3
) %>%
  mutate(
    datasets = factor(
      datasets,
      levels = c("D1", "D2", "D3")
    )
  ) %>%
  arrange(datasets, day)


# -------------------------------------------------------
# Parameter-level combination of D2 and D3 (no raw pooling)
# -------------------------------------------------------
res_D2D3_param <- results_all %>%
  filter(datasets %in% c("D2", "D3")) %>%
  group_by(day) %>%
  # Ensure that only the days jointly owned by D2 and D3 are merged
  filter(n() == 2) %>% 
  summarise(
    # Calculate the degrees of freedom (df = n_A + n_B -2)
    df_D2 = (n_A[datasets == "D2"] + n_B[datasets == "D2"] - 2),
    df_D3 = (n_A[datasets == "D3"] + n_B[datasets == "D3"] - 2),
    
    # Combined variance (weighted average based on degrees of freedom
    pooled_var = (df_D2 * pooled_var[datasets == "D2"] + 
                    df_D3 * pooled_var[datasets == "D3"]) / (df_D2 + df_D3),
    
    # Merge other key statistics
    delta     = mean(delta),
    n_A       = mean(n_A),
    n_B       = mean(n_B),
    mean_A    = mean(mean_A),
    mean_B    = mean(mean_B),
    var_A     = mean(var_A),
    var_B     = mean(var_B),
    .groups = "drop"
  ) %>%
  mutate(datasets = "D2 + D3") %>%
  # Adjust the column order to be exactly the same as the results_all structure
  select(datasets, everything(), -df_D2, -df_D3)

# -------------------------------------------------------
# Merge back into the main table
# -------------------------------------------------------
results_all <- bind_rows(results_all, res_D2D3_param) %>%
  mutate(datasets = factor(datasets, levels = c("D1", "D2", "D3", "D2 + D3"))) %>%
  arrange(datasets, day)

# Visualize the evolution of Delta
P2 <- ggplot(results_all, aes(x = day, y = abs(delta), group = datasets, color = datasets)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = c(5, 10, 15), linetype = "dashed", color = "red") +
  annotate("text", x = 5, y = 12, label = "Candidate Thresholds", color = "red") +
  labs(title = "Trend of Mean Difference (Delta) Over Time", 
       y = "Absolute Delta", x = "Days") +
  theme_light()


# -----------------------------------------------------------
# Calculate the standard t, welch (pwr, MKpower)
# -----------------------------------------------------------
results_ss <- results_all %>%
  filter(abs(delta) >= 1, pooled_var > 0) %>%
  
  rowwise() %>%
  mutate(
    # Standard t-test sample size
    n_standard = tryCatch({
      power.t.test(delta = abs(delta), sd = sqrt(pooled_var), 
                   sig.level = 0.05, power = 0.8)$n
    }, error = function(e) NA_real_),
    
    # Welch tests the sample size 
    n_welch = {
      res <- tryCatch({
        MKpower::power.welch.t.test(delta = abs(delta), 
                                    sd1 = sqrt(var_A), sd2 = sqrt(var_B),
                                    sig.level = 0.05, power = 0.8)
      }, error = function(e) NULL)
      
      if (is.null(res)) {
        NA_real_
      } else if (!is.null(res[["n1"]])) {
        as.numeric(res[["n1"]])
      } else if (!is.null(res[["n"]])) {
        as.numeric(res[["n"]])
      } else {
        NA_real_
      }
    }
  ) %>%
  ungroup()

# --- Final screening: Delta >= 10 determined based on visualization ---
final_planning_table <- results_ss %>%
  filter(abs(delta) >= 10) 



# -------------------------------------------------------
# Compute sample size per valid time point
# -------------------------------------------------------

plot_data <- results_ss %>%
  select(datasets, day, n_standard, n_welch) %>%
  pivot_longer(cols = starts_with("n_"), 
               names_to = "test_type", 
               values_to = "n_per_group")


# Calculate the true median sample size for each dataset
real_ref_lines <- results_ss %>%
  group_by(datasets) %>%
  summarise(
    real_target_n = ceiling(median(n_welch, na.rm = TRUE)) 
  )





# -------------------------------------------------------
# Robust WMW Sample Size (Noether's Formula)
# -------------------------------------------------------
# 1 Non-parametric effect size calculation function.

result_ss_full <- pilot_clean %>%
  # 1. Group by dataset and days (handle absolute volume)
  group_by(dataset, day) %>%
  summarise(
    n_A = sum(group == "A"),
    n_B = sum(group == "B"),
    mean_A = mean(tumor[group == "A"], na.rm = TRUE),
    mean_B = mean(tumor[group == "B"], na.rm = TRUE),
    sd_A   = sd(tumor[group == "A"], na.rm = TRUE),
    sd_B   = sd(tumor[group == "B"], na.rm = TRUE),
    # Calculate the non-parametric effect size A_val
    A_val  = estimate_A(tumor[group == "A"], tumor[group == "B"]),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    # Call the integrated computing function to obtain 8 variants
    # This function should return data frames containing both single-tail and double-tail versions
    res = list(get_all_n_variations(mean_A, mean_B, sd_A, sd_B, A_val))
  ) %>%
  # Expand the column of calculation results
  unnest(cols = res)



## Calculate the difference parameters for different days
df <- result_ss_full %>%
  mutate(
    delta_mean = abs(mean_B - mean_A),
    sd_ratio   = sd_B / sd_A,
    var_ratio  = (sd_B^2) / (sd_A^2)
  )

# Pull the n_* column into a long table and parse the method/tail
n_long <- df %>%
  pivot_longer(
    cols = starts_with("n_"),
    names_to = "n_type",
    values_to = "n"
  ) %>%
  # n_type: n_cohen_2t / n_wmw_1t
  extract(n_type, into = c("method","tail"), regex = "^n_([^_]+)_([^_]+)$") %>%
  mutate(
    tail = recode(tail, "1t" = "one-sided", "2t" = "two-sided"),
    method = factor(method, levels = c("cohen","pooled","welch","wmw")),
    method_label = recode(method,
                          "cohen"  = "t-test (Cohen d)",
                          "pooled" = "t-test (pooled var)",
                          "welch"  = "Welch t-test",
                          "wmw"    = "Wilcoxon/MW"
    )
  )


## The mean trajectory, how does the mean difference change over time, and is the SD large
mean_long <- df %>%
  select(dataset, day, mean_A, mean_B, sd_A, sd_B) %>%
  pivot_longer(
    cols = c(mean_A, mean_B, sd_A, sd_B),
    names_to = c(".value","group"),
    names_pattern = "(mean|sd)_(A|B)"
  ) %>%
  mutate(group = recode(group, "A"="Group A", "B"="Group B"))


## The sample size varies with day
P3 <- ggplot(n_long, aes(x = day, y = n, group = interaction(method_label, tail))) +
  geom_line(aes(linetype = tail)) +
  geom_point(size = 1.6) +
  facet_wrap(~dataset, scales = "free_x") +
  scale_y_log10() +
  labs(
    x = "Day", y = "Required sample size per group (log10 scale)",
    linetype = "Tail",
    title = "Sample size sensitivity across methods and tail choices"
  ) +
  theme_bw()

P4 <- ggplot(mean_long, aes(x = day, y = mean, group = group, color = group)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.3) +
  facet_wrap(~dataset, scales = "free_y") +
  labs(
    x = "Day", y = "Mean ± SD",
    color = "Group",
    title = "Tumor volume trajectories (mean ± SD)"
  ) +
  theme_bw()


