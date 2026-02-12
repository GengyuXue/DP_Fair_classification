# --- pkgs ---
library(ggplot2)
library(dplyr)
library(patchwork)


df_base <- all_df

pal_color <- c(
  "0.75"="#748EC2","1"="#7262AC","2"="#5F7D68",
  "3"="#924144","4"="#CB6030","5"="#E0A34A"
)

eps_order   <- names(pal_color)
eps_present <- intersect(eps_order, as.character(sort(unique(df_base$epsilon))))
df_base <- df_base %>% filter(as.character(epsilon) %in% eps_present)

# --- summarise over runs ---
summ_df <- df_base %>%
  group_by(Cw, epsilon) %>%
  summarise(
    err_mean = mean(err,    na.rm = TRUE),
    err_sd   = sd(err,      na.rm = TRUE),
    dd_mean  = mean(abs_DD, na.rm = TRUE),
    dd_sd    = sd(abs_DD,   na.rm = TRUE),
    mc = dplyr::n(), .groups = "drop"
  ) %>%
  mutate(
    err_se = ifelse(mc > 1, err_sd / sqrt(mc), 0),
    dd_se  = ifelse(mc > 1, dd_sd  / sqrt(mc), 0),
    err_lo = pmax(0, err_mean - 1.96 * err_se),
    err_hi = pmin(1, err_mean + 1.96 * err_se),
    dd_lo  = pmax(0, dd_mean  - 1.96 * dd_se),
    dd_hi  = pmax(0, dd_mean  + 1.96 * dd_se)
  )

fmt_eps <- function(x) sub("(?<=\\d)\\.0+$", "", as.character(x), perl = TRUE)


base_cols  <- scale_colour_manual(
  values = pal_color[eps_present],
  name   = "\u03F5",                                 
  limits = eps_present,
  breaks = eps_present,
)

base_fills <- scale_fill_manual(
  values = pal_color[eps_present],
  guide  = "none",
  limits = eps_present
)

base_theme <- theme_bw(base_size = 12) +
  theme(
    legend.position = "right",
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.minor = element_blank()
  )

# ---------------- Error plot ----------------
err_lim <- range(c(summ_df$err_lo, summ_df$err_hi), na.rm = TRUE)

p_err <- ggplot(
  summ_df,
  aes(x = Cw, y = err_mean, color = as.factor(epsilon), group = epsilon)
) +
  geom_ribbon(
    aes(ymin = err_lo, ymax = err_hi, fill = as.factor(epsilon)),
    alpha = 0.22, colour = NA, show.legend = FALSE
  ) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.6, shape = 16) +
  base_cols + base_fills +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) +
  scale_y_continuous(
    limits = err_lim + c(-0.02, 0.02),
    expand = expansion(mult = c(0.02, 0.06))
  ) +
  labs(x = expression(C[omega]), y = "Error")

# ---------------- Disparity plot  ----------------
dd_lim <- range(c(summ_df$dd_lo, summ_df$dd_hi, 0.3), na.rm = TRUE)
dd_lim[1] <- max(0, dd_lim[1])

p_dd <- ggplot(
  summ_df,
  aes(x = Cw, y = dd_mean, color = as.factor(epsilon), group = epsilon)
) +
  geom_ribbon(
    aes(ymin = dd_lo, ymax = dd_hi, fill = as.factor(epsilon)),
    alpha = 0.22, colour = NA, show.legend = FALSE
  ) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.6, shape = 16) +
  base_cols + base_fills +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) +
  scale_y_continuous(
    limits = c(max(0, dd_lim[1] - 0.02), dd_lim[2] + 0.01),
    breaks = scales::pretty_breaks(n = 6),
    expand = expansion(mult = c(0.02, 0.06))
  ) +
  geom_hline(yintercept = 0.3, linetype = "dashed", color = "grey45", linewidth = 0.6) +
  labs(x = expression(C[omega]), y = "Disparity")

# ---------------- Combine  ----------------
p_combined <- (p_err + p_dd) + plot_layout(guides = "collect") & base_theme
big_theme <- theme(
  axis.title.x = element_text(size = 16, face = "plain"),
  axis.title.y = element_text(size = 16, face = "plain"),
  axis.text.x  = element_text(size = 14),
  axis.text.y  = element_text(size = 14),
  legend.title = element_text(size = 14, face = "bold"),
  legend.text  = element_text(size = 13),
  legend.key.width = unit(0.8, "cm")
)

p_combined <- p_combined & big_theme
print(p_combined)

