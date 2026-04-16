# -----------------------------------------------------------------------------
# Exploratory figure for Chapter 3
# -----------------------------------------------------------------------------
# This script generates the exploratory figure used in Chapter 3 to show the
# tumour trajectories in the three original pilot datasets.
#
# The figure highlights:
# - the earliest day on which a clear group separation becomes visible
# - the later day from which the prompt-input time window begins
#
# The plot is exported as both PDF and PNG for use in the thesis.
#
# Note: this script assumes that the object `pilot_plus` has already been
# created in the main reference pipeline.
# -----------------------------------------------------------------------------

library(ggplot2)

# Earliest day with visually recognisable group separation in each dataset
sep_df <- data.frame(
  dataset = c("D1", "D2", "D3"),
  sep_day = c(8, 12, 12)
)

# First day included in the prompt-input window for each dataset
input_df <- data.frame(
  dataset   = c("D1", "D2", "D3"),
  input_day = c(13, 14, 14)
)

# Compute mean tumour area and standard error by dataset, day, and group
plot_df <- pilot_plus %>%
  filter(dataset %in% c("D1", "D2", "D3")) %>%
  group_by(dataset, day, group) %>%
  summarise(
    mean_tumor = mean(tumor, na.rm = TRUE),
    sd_tumor   = sd(tumor, na.rm = TRUE),
    n          = sum(!is.na(tumor)),
    se_tumor   = sd_tumor / sqrt(n),
    .groups    = "drop"
  ) %>%
  mutate(
    dataset = factor(dataset, levels = c("D1", "D2", "D3")),
    group   = factor(group, levels = c("A", "B"))
  )

# Build the exploratory line plot with error bars and reference lines
P_line <- ggplot(
  plot_df,
  aes(
    x = day,
    y = mean_tumor,
    group = group,
    color = group,
    linetype = group,
    shape = group
  )
) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.2, stroke = 0.7) +
  geom_errorbar(
    aes(ymin = mean_tumor - se_tumor, ymax = mean_tumor + se_tumor),
    width = 0.22,
    linewidth = 0.45,
    show.legend = FALSE
  ) +
  geom_vline(
    data = sep_df,
    aes(xintercept = sep_day),
    linetype = "dotted",
    linewidth = 0.5,
    colour = "grey30",
    inherit.aes = FALSE
  ) +
  geom_vline(
    data = input_df,
    aes(xintercept = input_day),
    linetype = "dashed",
    linewidth = 0.6,
    colour = "grey10",
    inherit.aes = FALSE
  ) +
  facet_wrap(~ dataset, scales = "free_x", nrow = 1) +
  scale_color_manual(
    values = c("A" = "#0072B2", "B" = "#D55E00"),
    name = "Gruppe"
  ) +
  scale_linetype_manual(
    values = c("A" = "solid", "B" = "longdash"),
    name = "Gruppe"
  ) +
  scale_shape_manual(
    values = c("A" = 16, "B" = 17),
    name = "Gruppe"
  ) +
  labs(
    x = "Tag nach Tumorzellinjektion",
    y = expression("Mittlere Tumorfläche (mm"^2*")")
  ) +
  guides(
    color = guide_legend(
      nrow = 1,
      byrow = TRUE,
      override.aes = list(linewidth = 0.9, size = 2.4)
    ),
    linetype = "none",
    shape = "none"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey95", colour = "black", linewidth = 0.4),
    strip.text = element_text(face = "bold", size = 11),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 10, colour = "black"),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box.margin = margin(t = 2),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(colour = "grey85", linewidth = 0.3),
    panel.border = element_rect(linewidth = 0.5),
    plot.margin = margin(5.5, 8, 5.5, 5.5)
  )


# Print the plot in the current R session
P_line

# Export the figure in high quality for the thesis
ggsave(
  filename = "tumorverlaeufe_prompt_input.pdf",
  plot = P_line,
  width = 16,
  height = 7,
  units = "cm",
  device = cairo_pdf
)

ggsave(
  filename = "tumorverlaeufe_prompt_input.png",
  plot = P_line,
  width = 16,
  height = 7,
  units = "cm",
  dpi = 600
)



