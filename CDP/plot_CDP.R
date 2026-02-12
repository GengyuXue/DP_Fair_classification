# --- required packages ---
library(ggplot2)
library(ggh4x)
library(dplyr)

source("oracle.R")

# ============================================================
# Toggle selection (subset vs full)
# ============================================================
USE_SELECTION <- TRUE           
n_keep   <- c(5000, 7000, 9000) # selected n_train
eps_keep <- c(0.75, 1, 2, 3, 4) # selected epsilon


stopifnot(all(c("n_train","alpha","epsilon","risk","disparity") %in% names(all_df)))

# ============================================================
# summarise over runs
# ============================================================
df_use <- all_df %>%
  { if (USE_SELECTION) dplyr::filter(., n_train %in% n_keep, epsilon %in% eps_keep) else . }

summ_df <- df_use %>%
  group_by(n_train, alpha, epsilon) %>%
  summarise(
    error_rate_mean = mean(risk, na.rm = TRUE),
    error_rate_sd   = sd(risk,   na.rm = TRUE),
    disparity_mean  = mean(disparity, na.rm = TRUE),
    disparity_sd    = sd(disparity,   na.rm = TRUE),
    mc = n(), .groups = "drop"
  ) %>%
  mutate(
    error_rate_se   = error_rate_sd / sqrt(mc),
    disparity_se    = disparity_sd / sqrt(mc),
    error_lower     = pmax(0, error_rate_mean - 1.96 * error_rate_se),
    error_upper     = pmin(1, error_rate_mean + 1.96 * error_rate_se),
    disparity_lower = pmax(0, disparity_mean - 1.96 * disparity_se),
    disparity_upper = pmin(1, disparity_mean + 1.96 * disparity_se)
  )

# ============================================================
# Long format: stack Error / Disparity rows
# ============================================================
df_err_long <- summ_df %>%
  transmute(
    n_train, alpha, epsilon,
    metric = factor("Error", levels = c("Error","Disparity")),
    value  = error_rate_mean,
    lower  = error_lower,
    upper  = error_upper
  )

df_dd_long <- summ_df %>%
  transmute(
    n_train, alpha, epsilon,
    metric = factor("Disparity", levels = c("Error","Disparity")),
    value  = disparity_mean,
    lower  = disparity_lower,
    upper  = disparity_upper
  )

df_long <- bind_rows(df_err_long, df_dd_long)

# ============================================================
# color
# ============================================================
pal_color_all <- c(
  "0.75" = "#748EC2",
  "1"    = "#7262AC",
  "2"    = "#748B7E",
  "3"    = "#CB6030",
  "4"    = "#E0A34A"
)
eps_levels <- intersect(names(pal_color_all), as.character(sort(unique(df_long$epsilon))))
pal_color  <- pal_color_all[eps_levels]
pal_fill   <- paste0(pal_color, "33")

df_long$epsilon <- factor(df_long$epsilon, levels = eps_levels)

ntrain_levels <- sort(unique(df_long$n_train))
df_abline <- data.frame(
  metric   = factor("Disparity", levels = c("Error","Disparity")),
  n_train  = ntrain_levels,
  slope    = 1,
  intercept= 0
)

# ============================================================
# Oracle
# ============================================================


if (!exists("oracle_dp_montecarlo")) {
  stop("oracle_dp_montecarlo() not found.")
}
alpha_grid <- sort(unique(df_long$alpha))
df_oracle <- oracle_dp_montecarlo(
  alpha_list = alpha_grid,
  n_pop      = 200000,
  d          = 2,
  p_a        = 0.3,
  alpha1 = 4, beta1 = 2,
  alpha0 = 4.5, beta0 = 2,
  w = c(1,1), b = 0, gamma = -0.3, k = 12, eps_clip = 1e-3
)

if (!("error" %in% names(df_oracle))) {
  stop("Expected column 'error' in oracle output.")
}
df_oracle_all <- merge(data.frame(n_train = ntrain_levels), df_oracle, by = NULL)
df_oracle_all$metric <- factor("Error", levels = c("Error","Disparity"))


# ============================================================
# Plot
# ============================================================
p_combined <- ggplot(
  df_long,
  aes(x = alpha, y = value, color = epsilon, group = epsilon)
) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = epsilon),
              alpha = 0.25, colour = NA, show.legend = FALSE) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.2, shape = 16) +
  geom_line(
    data = df_oracle_all,
    aes(x = alpha, y = error, linetype = "Oracle"),
    inherit.aes = FALSE,
    color = "#58708A", linewidth = 0.8
  ) +
  scale_linetype_manual(name = NULL, values = c("Oracle" = "dashed")) +
  scale_colour_manual(values = pal_color, name = "\u03F5",
                      limits = eps_levels, breaks = eps_levels) +
  scale_fill_manual(values = pal_fill, guide = "none",
                    limits = eps_levels) +
  facet_grid(
    rows = vars(metric),
    cols = vars(n_train),
    scales = "free_y",
    switch = "y",
    labeller = labeller(n_train = function(x) paste0("N = ", x))
  ) +
  scale_x_continuous(
    breaks = seq(0.2, 1, 0.2),
    limits = c(min(df_long$alpha, na.rm = TRUE), max(df_long$alpha, na.rm = TRUE))
  ) +
  ggh4x::facetted_pos_scales(
    y = list(
      metric == "Error" ~ scale_y_continuous(
        limits = range(df_err_long$value, df_oracle_all$error, na.rm = TRUE) + c(-0.02, 0.02),
        expand = expansion(mult = c(0.02, 0.06))
      ),
      metric == "Disparity" ~ scale_y_continuous(
        limits = c(0, max(0.02, max(df_dd_long$upper, na.rm = TRUE))),
        breaks = scales::pretty_breaks(n = 6)
      )
    )
  ) +
  labs(x = expression(alpha), y = NULL, color = expression(epsilon)) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "right",
    strip.background = element_rect(fill = "white", colour = NA),
    strip.placement = "outside",
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    strip.text.y.left = element_text(angle = 90, face = "bold"),
    strip.text.x = element_text(face = "bold", color = "black")
  ) +
  geom_abline(
    data = df_abline,
    aes(slope = slope, intercept = intercept),
    linetype = "dashed", color = "gray60",
    inherit.aes = FALSE
  )

p_combined <- p_combined +
  theme(
    strip.text.x = element_text(size = 16, face = "bold", color = "black"),
    strip.text.y.left = element_text(size = 14,face = "plain", angle = 90)
  )+
  theme(
    axis.title.x = element_text(size = 16),
    legend.key.width = unit(0.8, "cm"),   
    legend.spacing.x = unit(0.25, "cm"),    
    legend.title = element_text(size = 14, face = "bold"),
    legend.text  = element_text(size = 13),
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13)
  )
print(p_combined)
