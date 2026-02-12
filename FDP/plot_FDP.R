# --- required packages ---
library(ggplot2)
library(ggh4x)
library(dplyr)

source("oracle.R")


# ------------------------------------------------------------
# 1. Filter the data
# ------------------------------------------------------------

S_keep   <- c(4,5,6)
eps_keep <- c(0.75, 1, 2, 3, 4, 10)

df_use <- all_df %>%
  filter(S %in% S_keep, epsilon %in% eps_keep)

# ------------------------------------------------------------
# 2. Summary
# ------------------------------------------------------------
summ_df <- df_use %>%
  group_by(S, alpha, epsilon) %>%
  summarise(
    error_rate_mean = mean(err,     na.rm = TRUE),
    error_rate_sd   = sd(err,       na.rm = TRUE),
    disparity_mean  = mean(abs_DD,  na.rm = TRUE),
    disparity_sd    = sd(abs_DD,    na.rm = TRUE),
    mc = n(), .groups = "drop"
  ) %>%
  mutate(
    error_rate_se   = error_rate_sd / sqrt(mc),
    disparity_se    = disparity_sd  / sqrt(mc),
    error_lower     = pmax(0, error_rate_mean - 1.96 * error_rate_se),
    error_upper     = pmin(1, error_rate_mean + 1.96 * error_rate_se),
    disparity_lower = pmax(0, disparity_mean  - 1.96 * disparity_se),
    disparity_upper = pmin(1, disparity_mean  + 1.96 * disparity_se)
  )

# ------------------------------------------------------------
# 3. Long format
# ------------------------------------------------------------
df_err_long <- summ_df %>%
  transmute(
    S, alpha, epsilon,
    metric = factor("Error", levels = c("Error","Disparity")),
    value  = error_rate_mean,
    lower  = error_lower,
    upper  = error_upper
  )

df_dd_long <- summ_df %>%
  transmute(
    S, alpha, epsilon,
    metric = factor("Disparity", levels = c("Error","Disparity")),
    value  = disparity_mean,
    lower  = disparity_lower,
    upper  = disparity_upper
  )

df_long <- bind_rows(df_err_long, df_dd_long)

# ------------------------------------------------------------
# 4. Colour and label
# ------------------------------------------------------------
pal_color <- c(
  "0.75" = "#748EC2",
  "1"    = "#7262AC",
  "2"    = "#748B7E",  
  "3"    = "#CB6030",
  "4"    = "#E0A34A",
  "10"   = "#9F6E45"
)
pal_fill <- paste0(pal_color, "33")


df_long$epsilon <- factor(df_long$epsilon, levels = names(pal_color))

fmt_eps <- function(x) paste0("\u03F5 = ", x)

# ------------------------------------------------------------
# 5. Reference line y=x
# ------------------------------------------------------------
S_levels <- sort(unique(df_long$S))

df_abline <- data.frame(
  metric = factor("Disparity", levels = c("Error","Disparity")),
  S = S_levels,
  slope = 1, intercept = 0
)

# ------------------------------------------------------------
# 6. Base plot
# ------------------------------------------------------------
p_combined <- ggplot(
  df_long,
  aes(x = alpha, y = value, color = epsilon, group = epsilon)
) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = epsilon),
              alpha = 0.25, colour = NA, show.legend = FALSE) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.2, shape = 16) +
  scale_colour_manual(
    values = pal_color,
    name   = "\u03F5",         
    breaks = levels(df_long$epsilon),
    limits = levels(df_long$epsilon)
  ) +
  scale_fill_manual(values = pal_fill, guide = "none") +
  facet_grid(
    rows = vars(metric),
    cols = vars(S),
    scales = "free_y",
    switch = "y",
    labeller = labeller(S = function(x) paste0("S = ", x))
  ) +
  scale_x_continuous(
    breaks = seq(0.2, 1, 0.2),
    limits = c(min(df_long$alpha, na.rm = TRUE), max(df_long$alpha, na.rm = TRUE))
  ) +
  facetted_pos_scales(
    y = list(
      metric == "Error" ~ scale_y_continuous(
        limits = range(df_err_long$value, na.rm = TRUE) + c(-0.02, 0.02),
        expand = expansion(mult = c(0.02, 0.06))
      ),
      metric == "Disparity" ~ scale_y_continuous(
        limits = c(0, max(0.02, max(df_dd_long$upper, na.rm = TRUE))),
        breaks = scales::pretty_breaks(n = 6)
      )
    )
  ) +
  labs(x = expression(alpha), y = NULL, color = "\u03F5") +  
  theme_bw(base_size = 12) +
  theme(
    legend.position = "right",
    strip.background = element_rect(fill = "white", colour = NA),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 90, face = "bold"),
    strip.text.x = element_text(face = "bold", color = "black"),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.minor = element_blank()
  ) +
  geom_abline(
    data = df_abline,
    aes(slope = slope, intercept = intercept),
    color = "gray60",
    linetype = "dashed",
    inherit.aes = FALSE
  )

# ------------------------------------------------------------
# 7. Oracle curve 
# ------------------------------------------------------------
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

df_oracle_all <- merge(data.frame(S = S_levels), df_oracle, by = NULL)
df_oracle_all$metric <- factor("Error", levels = c("Error","Disparity"))

p_combined <- p_combined +
  geom_line(
    data = df_oracle_all,
    aes(x = alpha, y = error, linetype = "Oracle"),
    inherit.aes = FALSE,
    color = "#58708A", linewidth = 0.7
  ) +
  scale_linetype_manual(name = NULL, values = c("Oracle" = "dashed"))

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

# ------------------------------------------------------------
print(p_combined)
