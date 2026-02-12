library(ggplot2)
library(dplyr)
library(ggh4x)


df_base <- all_df

if (!"method" %in% names(df_base)) df_base$method <- "FDP"
df_base <- df_base %>% filter((method == "CDP") | (S %in% c(1,2,3,4,5)))

# Keep only these epsilons
eps_keep <- c(0.75, 1, 2)
df_base  <- df_base %>% filter(epsilon %in% eps_keep)


df_base <- df_base %>% mutate(series = ifelse(method == "CDP", "CDP", paste0("S=", S)))

eps_present    <- eps_keep[eps_keep %in% sort(unique(df_base$epsilon))]
stopifnot(length(eps_present) > 0)
S_levels       <- paste0("S=", sort(unique(df_base$S[df_base$method != "CDP"])))
series_levels  <- c(S_levels, if (any(df_base$method == "CDP")) "CDP")

# summary
summ_df <- df_base %>%
  group_by(alpha, epsilon, series) %>%
  summarise(
    err_mean = mean(err,    na.rm = TRUE),
    err_sd   = sd(err,      na.rm = TRUE),
    dd_mean  = mean(abs_DD, na.rm = TRUE),
    dd_sd    = sd(abs_DD,   na.rm = TRUE),
    mc = dplyr::n(),
    .groups = "drop"
  ) %>%
  mutate(
    err_se = ifelse(mc > 1, err_sd / sqrt(mc), 0),
    dd_se  = ifelse(mc > 1, dd_sd / sqrt(mc), 0),
    err_lo = pmax(0, err_mean - 1.96 * err_se),
    err_hi = pmin(1, err_mean + 1.96 * err_se),
    dd_lo  = pmax(0, dd_mean  - 1.96 * dd_se),
    dd_hi  = pmax(0, dd_mean  + 1.96 * dd_se)
  )

summ_df$series  <- factor(summ_df$series,  levels = series_levels)
summ_df$epsilon <- factor(summ_df$epsilon, levels = eps_present)

# Long format
df_err <- summ_df %>%
  transmute(alpha, epsilon, series,
            metric = factor("Error", levels = c("Error","Disparity")),
            value = err_mean, lower = err_lo, upper = err_hi)

df_dd <- summ_df %>%
  transmute(alpha, epsilon, series,
            metric = factor("Disparity", levels = c("Error","Disparity")),
            value = dd_mean, lower = dd_lo, upper = dd_hi)

df_long <- dplyr::bind_rows(df_err, df_dd) %>% arrange(metric, epsilon, alpha)


df_abline <- data.frame(
  metric   = factor("Disparity", levels = c("Error","Disparity")),
  epsilon  = factor(eps_present, levels = eps_present),
  slope    = 1,
  intercept= 0
)


# Colors per S
# base_colors_vec <- c("#748EC2", "#7262AC", "#5F7D68",
#                      "#924144", "#CB6030", "#E0A34A")

base_colors_vec <- c("#748EC2", "#7262AC", "#5F7D68",
                     "#CB6030", "#E0A34A")


col_S <- rep(base_colors_vec, length.out = length(S_levels))
names(col_S) <- S_levels

pal_series <- col_S
if ("CDP" %in% series_levels) pal_series["CDP"] <- "#D65DB1"
pal_series <- pal_series[series_levels] 

base_cols  <- scale_colour_manual(values = pal_series, breaks = series_levels, name = NULL)
base_fills <- scale_fill_manual(  values = pal_series, breaks = series_levels, guide = "none")

#limit
err_lim <- range(df_err$value, df_err$lower, df_err$upper, na.rm = TRUE)
dd_lim  <- range(df_dd$value,  df_dd$lower,  df_dd$upper,  na.rm = TRUE)

alpha_max <- max(df_long$alpha, na.rm = TRUE)
dd_lim[2] <- max(dd_lim[2], alpha_max)  # show y = x fully

# small padding
err_lim <- err_lim + c(-0.02, 0.02)
dd_lim  <- dd_lim  + c(-0.02, 0.02)

# plot
p <- ggplot(
  df_long,
  aes(x = alpha, y = value, color = series, group = series)
) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = series),
              alpha = 0.25, colour = NA, show.legend = FALSE, na.rm = TRUE) +
  geom_line(linewidth = 0.7, na.rm = TRUE) +
  geom_point(size = 1.1, na.rm = TRUE) +
  base_cols + base_fills +
  
  facet_grid(
    rows = vars(metric),
    cols = vars(epsilon),
    switch   = "y",
    labeller = labeller(epsilon = function(x) paste0("\u03F5 = ", x)),
    scales   = "free_y"  
  ) +
  
  ggh4x::facetted_pos_scales(
    y = list(
      metric == "Error" ~ scale_y_continuous(
        limits = c(0.1,0.33),
        breaks = scales::pretty_breaks(n = 6),
        minor_breaks = NULL,
        oob = scales::oob_squish
      ),
      metric == "Disparity" ~ scale_y_continuous(
        limits = c(0,0.6),
        breaks = scales::pretty_breaks(n = 6),
        minor_breaks = NULL,
        oob = scales::oob_squish
      )
    )
  ) +
  
  scale_x_continuous(
    breaks = seq(0.2, 1, 0.2),
    limits = c(0.1, 1),
    minor_breaks = NULL
  ) +
  

  geom_abline(
    data = df_abline,
    aes(slope = slope, intercept = intercept),
    linetype = "dashed", color = "gray60",
    inherit.aes = FALSE
  ) +
  
  labs(x = expression(alpha), y = NULL) +
  theme_bw(base_size = 13) +
  theme(
    legend.position   = "right",
    strip.background  = element_rect(fill = "white", colour = NA),
    strip.placement   = "outside",
    strip.text.y.left = element_text(angle = 90, face = "bold"),
    strip.text.x      = element_text(face = "bold", color = "black", hjust = 0.5),
    panel.grid.major  = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.minor  = element_blank()
  )

p  <- p +
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
print(p)