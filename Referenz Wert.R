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


# 1 非参数效应量计算函数
estimate_A <- function(x_A, x_B) {
  x_A <- x_A[is.finite(x_A)]
  x_B <- x_B[is.finite(x_B)]
  if (length(x_A) == 0 || length(x_B) == 0) return(NA_real_)
  cmp <- outer(x_A, x_B, FUN = "-")
  (sum(cmp > 0) + 0.5 * sum(cmp == 0)) / (length(x_A) * length(x_B))
}

# 自动生成所有可能的天数组合
all_days <- sort(unique(pilot_clean$day))
day_combos <- combn(all_days, 2)

cfb_all_pairs <- list()
for(i in 1:ncol(day_combos)) {
  d_start <- day_combos[1, i]; d_end <- day_combos[2, i]
  
  cfb_all_pairs[[i]] <- pilot_clean %>%
    filter(day %in% c(d_start, d_end)) %>%
    group_by(dataset, mouse_id, group) %>%
    filter(n() == 2) %>% # 确保两点都有数据
    summarise(diff = tumor[day == d_end] - tumor[day == d_start], .groups = "drop") %>%
    mutate(time_pair = paste0("D", d_end, "-D", d_start))
}

cfb_matrix <- bind_rows(cfb_all_pairs)

# -------------------------------------------------------
# Calculator for LLM Benchmarking
# -------------------------------------------------------

get_comprehensive_n <- function(m_diff, v1, v2, A_val, alpha=0.05, power=0.8) {
  delta_abs <- abs(m_diff)
  # pooled variance simplified because all pilot datasets satisfy n1 = n2
  s_pooled  <- sqrt((v1 + v2) / 2)
  
  calc_by_side <- function(label) {
    
    alt_pwr <- ifelse(label == "Two-Sided", "two.sided", "greater")
    alt_mk  <- ifelse(label == "Two-Sided", "two.sided", "one.sided")
    
    # pwr / MKpower
    alpha_t <- alpha
    
    tryCatch({
      # pooled t, t-based numeric via pwr.t.test
      n_c <- pwr.t.test(
        d = delta_abs / s_pooled,
        sig.level = alpha_t,
        power = power,
        type = "two.sample",
        alternative = alt_pwr
      )$n
      # Welch t, t-based numeric
      n_w <- power.welch.t.test(
        delta = delta_abs,
        sd1 = sqrt(v1), sd2 = sqrt(v2),
        sig.level = alpha_t,
        power = power,
        alternative = alt_mk
      )$n
      
      # WMW (Noether) — n per group
      sided <- ifelse(label == "Two-Sided", 2, 1)
      n_m <- (qnorm(1 - alpha / sided) + qnorm(power))^2 / (6 * (A_val - 0.5)^2)
      
      # Welch z approximation, unpooled variance term
      z_alpha <- if (label == "Two-Sided") qnorm(1 - alpha/2) else qnorm(1 - alpha)
      z_beta  <- qnorm(power)  # since power = 1 - beta
      n_z <- (v1 + v2) * (z_alpha + z_beta)^2 / (delta_abs^2)
      
      # pooled z approximation
      n_pooled_z <- 2 * (s_pooled^2) * (z_alpha + z_beta)^2 / (delta_abs^2)
      
      data.frame(Side = label, n_pooled_t = n_c, n_Welch_t = n_w, n_WMW = n_m, n_welch_z = n_z, n_pooled_z = n_pooled_z)
    }, error = function(e) {
      data.frame(Side = label, n_pooled_t = NA, n_Welch_t = NA, n_WMW = NA, n_welch_z = NA, n_pooled_z = NA)
    })
  }
  
  bind_rows(calc_by_side("Two-Sided"), calc_by_side("One-Sided"))
}

# 处理绝对体积数据
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
    res <- get_comprehensive_n(.$delta, .$v_A, .$v_B, .$A_val)   # <-- 这里是正确调用
    data.frame(., res)
  }) %>%
  ungroup() %>%
  rename(time_point = time_pair)


# 合并数据(结果层面)
ds_D2D3_clean <- clean_pilot_dataset(ds2_3, "D2+D3")
ds_D2D3_clean <- ds_D2D3_clean %>%
  group_by(dataset, day) %>%
  filter(
    sum(group == "A") == 10 & 
      sum(group == "B") == 10
  ) %>%
  ungroup()


pilot_plus <- bind_rows(pilot_clean, ds_D2D3_clean)
pilot_plus <- pilot_plus %>%
  mutate(dataset = factor(dataset, levels = c("D1","D2","D3","D2+D3")))


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

# 合并结果
master_planning_table <- bind_rows(
  final_results %>% rename(time_point = day) %>% mutate(time_point = as.character(time_point)),
  B %>% rename(time_point = day) %>% mutate(time_point = as.character(time_point)),
  final_cfb_analysis
) %>%
  mutate(across(starts_with("n_"), ceiling)) %>%
  select(dataset, time_point, Side, m_A, m_B, v_A, v_B, A_val, everything())



write.csv(master_planning_table, file = "Referentable.csv", row.names = FALSE)
