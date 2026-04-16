###############################################################################
# manual_reference_checks.R
# Manual reference calculations aligned with the thesis reference table
###############################################################################

# ---------------------------------------------------------------------------
# 0) Helper functions
# ---------------------------------------------------------------------------

normalize_side <- function(side) {
  side <- trimws(as.character(side))
  if (side %in% c("Two-Sided", "two-sided", "two_sided", "two sided")) return("Two-Sided")
  if (side %in% c("One-Sided", "one-sided", "one_sided", "one sided")) return("One-Sided")
  stop("Unknown Side value: ", side)
}

z_alpha <- function(alpha, side) {
  side <- normalize_side(side)
  if (side == "Two-Sided") qnorm(1 - alpha / 2) else qnorm(1 - alpha)
}

z_beta <- function(power) {
  qnorm(power)
}

pooled_variance <- function(v_A, v_B, pilot_n_A = 5, pilot_n_B = 5) {
  ((pilot_n_A - 1) * v_A + (pilot_n_B - 1) * v_B) / (pilot_n_A + pilot_n_B - 2)
}

power_z <- function(ncp, alpha, side) {
  side <- normalize_side(side)
  zcrit <- z_alpha(alpha, side)
  if (side == "Two-Sided") {
    pnorm(-zcrit, mean = ncp, sd = 1) + (1 - pnorm(zcrit, mean = ncp, sd = 1))
  } else {
    1 - pnorm(zcrit, mean = ncp, sd = 1)
  }
}

power_t <- function(ncp, df, alpha, side) {
  side <- normalize_side(side)
  if (side == "Two-Sided") {
    tcrit <- qt(1 - alpha / 2, df = df)
    pt(-tcrit, df = df, ncp = ncp) + (1 - pt(tcrit, df = df, ncp = ncp))
  } else {
    tcrit <- qt(1 - alpha, df = df)
    1 - pt(tcrit, df = df, ncp = ncp)
  }
}

# ---------------------------------------------------------------------------
# 1) pooled z
# Formula class: z-based approximation under equal variances
# Returns per-group sample size
# ---------------------------------------------------------------------------

calc_n_pooled_z <- function(delta, v_A, v_B, alpha = 0.05, power = 0.80,
                            side = "Two-Sided", pilot_n_A = 5, pilot_n_B = 5) {
  delta <- abs(delta)
  sp2 <- pooled_variance(v_A, v_B, pilot_n_A, pilot_n_B)
  z_a <- z_alpha(alpha, side)
  z_b <- z_beta(power)
  
  n_raw <- 2 * sp2 * (z_a + z_b)^2 / delta^2
  n <- ceiling(n_raw)
  
  # achieved power after rounding
  se <- sqrt(2 * sp2 / n)
  ncp <- delta / se
  actual_power <- power_z(ncp, alpha, side)
  
  list(n = n, actual_power = actual_power)
}

# ---------------------------------------------------------------------------
# 2) Welch z
# Formula class: z-based approximation under unequal variances
# Returns per-group sample size
# ---------------------------------------------------------------------------

calc_n_welch_z <- function(delta, v_A, v_B, alpha = 0.05, power = 0.80,
                           side = "Two-Sided") {
  delta <- abs(delta)
  z_a <- z_alpha(alpha, side)
  z_b <- z_beta(power)
  
  n_raw <- (v_A + v_B) * (z_a + z_b)^2 / delta^2
  n <- ceiling(n_raw)
  
  se <- sqrt(v_A / n + v_B / n)
  ncp <- delta / se
  actual_power <- power_z(ncp, alpha, side)
  
  list(n = n, actual_power = actual_power)
}

# ---------------------------------------------------------------------------
# 3) WMW / Noether approximation
# A_val is the probabilistic effect from your reference table
# Returns per-group sample size for balanced design
# ---------------------------------------------------------------------------

calc_n_wmw <- function(A_val, alpha = 0.05, power = 0.80, side = "Two-Sided") {
  effect_abs <- abs(A_val - 0.5)
  if (effect_abs <= 0) stop("A_val must differ from 0.5.")
  
  z_a <- z_alpha(alpha, side)
  z_b <- z_beta(power)
  
  # total N for r = 1
  N_raw <- 4 * (z_a + z_b)^2 / (12 * effect_abs^2)
  n <- ceiling(N_raw / 2)
  
  N <- 2 * n
  actual_info <- sqrt(12 * effect_abs^2 * N / 4)
  actual_power <- pnorm(actual_info - z_a)
  
  list(n = n, actual_power = actual_power)
}

# ---------------------------------------------------------------------------
# 4) pooled t
# Manual integer search using noncentral t power
# Returns per-group sample size
# ---------------------------------------------------------------------------

calc_n_pooled_t <- function(delta, v_A, v_B, alpha = 0.05, power = 0.80,
                            side = "Two-Sided", pilot_n_A = 5, pilot_n_B = 5,
                            max_n = 10000) {
  delta <- abs(delta)
  sp2 <- pooled_variance(v_A, v_B, pilot_n_A, pilot_n_B)
  
  for (n in 2:max_n) {
    df <- 2 * n - 2
    se <- sqrt(2 * sp2 / n)
    ncp <- delta / se
    pwr <- power_t(ncp, df, alpha, side)
    
    if (pwr >= power) {
      return(list(n = n, actual_power = pwr))
    }
  }
  
  stop("No pooled t solution found up to max_n.")
}

# ---------------------------------------------------------------------------
# 5) Welch t
# Manual integer search using Welch SE + Satterthwaite df + noncentral t power
# Returns per-group sample size
# ---------------------------------------------------------------------------

calc_n_welch_t <- function(delta, v_A, v_B, alpha = 0.05, power = 0.80,
                           side = "Two-Sided", max_n = 10000) {
  delta <- abs(delta)
  
  for (n in 2:max_n) {
    se2 <- v_A / n + v_B / n
    se <- sqrt(se2)
    
    df_num <- se2^2
    df_den <- (v_A / n)^2 / (n - 1) + (v_B / n)^2 / (n - 1)
    df <- df_num / df_den
    
    ncp <- delta / se
    pwr <- power_t(ncp, df, alpha, side)
    
    if (pwr >= power) {
      return(list(n = n, actual_power = pwr, df = df))
    }
  }
  
  stop("No Welch t solution found up to max_n.")
}

# ---------------------------------------------------------------------------
# 6) One-row calculation aligned with your table
# ---------------------------------------------------------------------------

calc_reference_row <- function(dataset, time_point, Side,
                               m_A, m_B, v_A, v_B, A_val, delta,
                               alpha = 0.05, power = 0.80,
                               pilot_n_A = 5, pilot_n_B = 5) {
  Side <- normalize_side(Side)
  
  pooled_t <- calc_n_pooled_t(delta, v_A, v_B, alpha, power, Side, pilot_n_A, pilot_n_B)
  welch_t  <- calc_n_welch_t(delta, v_A, v_B, alpha, power, Side)
  wmw      <- calc_n_wmw(A_val, alpha, power, Side)
  welch_z  <- calc_n_welch_z(delta, v_A, v_B, alpha, power, Side)
  pooled_z <- calc_n_pooled_z(delta, v_A, v_B, alpha, power, Side, pilot_n_A, pilot_n_B)
  
  data.frame(
    dataset = dataset,
    time_point = time_point,
    Side = Side,
    m_A = m_A,
    m_B = m_B,
    v_A = v_A,
    v_B = v_B,
    A_val = A_val,
    delta = delta,
    n_pooled_t = pooled_t$n,
    n_Welch_t = welch_t$n,
    n_WMW = wmw$n,
    n_welch_z = welch_z$n,
    n_pooled_z = pooled_z$n,
    stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------------------
# 7) Build full reference table from an aggregated input table
# Input must contain:
# dataset, time_point, Side, m_A, m_B, v_A, v_B, A_val, delta
# ---------------------------------------------------------------------------

build_reference_table <- function(df,
                                  alpha = 0.05,
                                  power = 0.80,
                                  pilot_n_A = 5,
                                  pilot_n_B = 5) {
  required_cols <- c("dataset", "time_point", "Side", "m_A", "m_B",
                     "v_A", "v_B", "A_val", "delta")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop("Missing columns: ", paste(missing_cols, collapse = ", "))
  }
  
  out <- vector("list", nrow(df))
  
  for (i in seq_len(nrow(df))) {
    out[[i]] <- calc_reference_row(
      dataset = df$dataset[i],
      time_point = df$time_point[i],
      Side = df$Side[i],
      m_A = df$m_A[i],
      m_B = df$m_B[i],
      v_A = df$v_A[i],
      v_B = df$v_B[i],
      A_val = df$A_val[i],
      delta = df$delta[i],
      alpha = alpha,
      power = power,
      pilot_n_A = pilot_n_A,
      pilot_n_B = pilot_n_B
    )
  }
  
  do.call(rbind, out)
}

# ---------------------------------------------------------------------------
# 8) Compare manual results with an existing reference table
# ---------------------------------------------------------------------------

compare_reference_tables <- function(original_df, manual_df) {
  key_cols <- c("dataset", "time_point", "Side")
  
  merged <- merge(
    original_df, manual_df,
    by = key_cols,
    suffixes = c("_orig", "_manual")
  )
  
  compare_cols <- c("n_pooled_t", "n_Welch_t", "n_WMW", "n_welch_z", "n_pooled_z")
  
  for (col in compare_cols) {
    merged[[paste0(col, "_diff")]] <- merged[[paste0(col, "_manual")]] - merged[[paste0(col, "_orig")]]
  }
  
  merged
}

# ---------------------------------------------------------------------------
# 9) Example usage
# ---------------------------------------------------------------------------

ref_input <- read.csv("Referentable.csv", stringsAsFactors = FALSE)

manual_ref <- build_reference_table(ref_input)

comparison <- compare_reference_tables(ref_input, manual_ref)

compare_cols <- c("n_pooled_t_diff", "n_Welch_t_diff", "n_WMW_diff", "n_welch_z_diff", "n_pooled_z_diff")
comparison[apply(comparison[, compare_cols] != 0, 1, any), ]

nrow(comparison)
colSums(comparison[, compare_cols] != 0)


# Save
write.csv(manual_ref, "Referenztabelle_manual.csv", row.names = FALSE)
write.csv(comparison, "Referenztabelle_comparison.csv", row.names = FALSE)
