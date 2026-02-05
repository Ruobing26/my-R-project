# -------------------------------------------------------
# Robust WMW Sample Size (Noether's Formula)
# -------------------------------------------------------
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

get_comprehensive_n <- function(m_diff, s1, s2, A_val, alpha=0.05, power=0.8) {
  delta_abs <- abs(m_diff)
  s_pooled  <- sqrt((s1^2 + s2^2) / 2)
  
  calc_by_side <- function(label) {
    
    alt_pwr <- ifelse(label == "Two-Sided", "two.sided", "greater")
    alt_mk  <- ifelse(label == "Two-Sided", "two.sided", "one.sided")
    
    # pwr / MKpower: 直接用 alpha（two.sided 会由函数内部处理）
    alpha_t <- alpha
    
    tryCatch({
      n_c <- pwr.t.test(
        d = delta_abs / s_pooled,
        sig.level = alpha_t,
        power = power,
        type = "two.sample",
        alternative = alt_pwr
      )$n
      
      n_w <- power.welch.t.test(
        delta = delta_abs,
        sd1 = s1, sd2 = s2,
        sig.level = alpha_t,
        power = power,
        alternative = alt_mk
      )$n
      
      # WMW (Noether) — n per group
      sided <- ifelse(label == "Two-Sided", 2, 1)
      n_m <- (qnorm(1 - alpha / sided) + qnorm(power))^2 / (12 * (A_val - 0.5)^2)
      
      # z-approx for two-sample mean difference (per group)
      z_alpha <- if (label == "Two-Sided") qnorm(1 - alpha/2) else qnorm(1 - alpha)
      z_beta  <- qnorm(power)  # since power = 1 - beta
      n_z <- (2 * (s_pooled^2) * (z_alpha + z_beta)^2) / (delta_abs^2)
      
      
      data.frame(Side = label, n_CohenD = n_c, n_Welch = n_w, n_WMW = n_m, n_z_approx = n_z)
    }, error = function(e) {
      data.frame(Side = label, n_CohenD = NA, n_Welch = NA, n_WMW = NA, n_z_approx = NA)
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
    s_A = sd(tumor[group == "A"], na.rm = TRUE),
    s_B = sd(tumor[group == "B"], na.rm = TRUE),
    delta = m_A - m_B,
    A_val = estimate_A(tumor[group == "A"], tumor[group == "B"]),
    .groups = "drop"
  ) %>%
  filter(
    is.finite(delta), is.finite(s_A), is.finite(s_B),
    abs(delta) > 0, s_A > 0, s_B > 0,
    is.finite(A_val),
    abs(A_val - 0.5) > 1e-6
  ) %>%
  rowwise() %>%
  do({
    res <- get_comprehensive_n(.$delta, .$s_A, .$s_B, .$A_val)
    data.frame(., res)
  }) %>%
  ungroup() %>%
  mutate(across(starts_with("n_"), ceiling),
         analysis_type = "Baseline")


final_cfb_analysis <- cfb_matrix %>%
  group_by(dataset, time_pair) %>%
  summarise(
    m_A = mean(diff[group == "A"], na.rm = TRUE),
    m_B = mean(diff[group == "B"], na.rm = TRUE),
    s_A = sd(diff[group == "A"], na.rm = TRUE),
    s_B = sd(diff[group == "B"], na.rm = TRUE),
    delta = m_A - m_B,
    A_val = estimate_A(diff[group == "A"], diff[group == "B"]),
    .groups = "drop"
  ) %>%
  filter(
    is.finite(delta), is.finite(s_A), is.finite(s_B),
    abs(delta) > 0, s_A > 0, s_B > 0,
    is.finite(A_val),
    abs(A_val - 0.5) > 1e-6
  ) %>%
  rowwise() %>%
  do({
    res <- get_comprehensive_n(.$delta, .$s_A, .$s_B, .$A_val)   # <-- 这里是正确调用
    data.frame(., res)
  }) %>%
  ungroup() %>%
  rename(time_point = time_pair) %>%
  mutate(
    across(starts_with("n_"), ceiling),
    analysis_type = "Change_From_Baseline"
  )



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


P5 <- ggplot(pilot_plus, aes(x = factor(day), y = tumor, fill = group)) +
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


B <- pilot_plus %>%
  group_by(dataset, day) %>%
  summarise(
    m_A = mean(tumor[group == "A"], na.rm = TRUE),
    m_B = mean(tumor[group == "B"], na.rm = TRUE),
    s_A = sd(tumor[group == "A"], na.rm = TRUE),
    s_B = sd(tumor[group == "B"], na.rm = TRUE),
    delta = m_A - m_B,
    A_val = estimate_A(tumor[group == "A"], tumor[group == "B"]),
    .groups = "drop"
  ) %>%
  filter(
    is.finite(delta), is.finite(s_A), is.finite(s_B),
    abs(delta) > 0, s_A > 0, s_B > 0,
    is.finite(A_val),
    abs(A_val - 0.5) > 1e-6
  ) %>%
  rowwise() %>%
  do({
    res <- get_comprehensive_n(.$delta, .$s_A, .$s_B, .$A_val)
    data.frame(., res)
  }) %>%
  ungroup() %>%
  mutate(across(starts_with("n_"), ceiling),
         analysis_type = "Baseline")

# 合并结果
master_planning_table <- bind_rows(
  final_results %>% rename(time_point = day) %>% mutate(time_point = as.character(time_point)),
  B %>% rename(time_point = day) %>% mutate(time_point = as.character(time_point)),
  final_cfb_analysis
) %>%
  mutate(across(starts_with("n_"), ceiling)) %>%
  select(dataset, analysis_type, time_point, Side, m_A, m_B, s_A, s_B, A_val, everything())

## 筛选数据
B %>%
  filter(abs(delta) > 10) %>%
  select(1:7)
