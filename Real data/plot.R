# --- packages ---
library(ggplot2)
library(ggh4x)
library(dplyr)


eps_list <- c(0.75, 1, 2, 3)   

# --- column names ---
err_col <- if ("err" %in% names(all_df)) "err" else if ("error" %in% names(all_df)) "error" else stop("No error/err column")
dd_col  <- if ("abs_DD" %in% names(all_df)) "abs_DD" else if ("disparity" %in% names(all_df)) "disparity" else stop("No disparity/abs_DD column")


if (is.factor(all_df$epsilon))    all_df$epsilon <- as.numeric(as.character(all_df$epsilon))
if (is.character(all_df$epsilon)) all_df$epsilon <- suppressWarnings(as.numeric(all_df$epsilon))
alpha_num <- all_df$alpha
if (is.factor(alpha_num))    alpha_num <- as.character(alpha_num)
if (is.character(alpha_num)) alpha_num <- gsub("[^0-9.+\\-eE]", "", alpha_num)
all_df$alpha <- suppressWarnings(as.numeric(alpha_num))

# --- choose epsilon ---
if (is.null(eps_list) || !length(eps_list)) eps_list <- sort(unique(na.omit(all_df$epsilon)))
stopifnot(length(eps_list) > 0)


df_base <- all_df %>% filter((epsilon %in% eps_list) | is.na(epsilon))


na_rows <- df_base %>% filter(is.na(epsilon))
if (nrow(na_rows) > 0) {
  na_expanded <- na_rows[rep(seq_len(nrow(na_rows)), each = length(eps_list)), , drop = FALSE]
  na_expanded$epsilon <- rep(eps_list, times = nrow(na_rows))
  df_base <- bind_rows(df_base %>% filter(!is.na(epsilon)), na_expanded)
}

# define series for colour
if ("S" %in% names(df_base)) {
  df_base <- df_base %>% mutate(series = paste0("S=", .data$S))
} else if ("method" %in% names(df_base)) {
  df_base <- df_base %>% mutate(series = as.character(.data$method))
} else {
  df_base$series <- "All"
}

# --- rename methods per mapping ---
df_base <- df_base %>%
  mutate(series = dplyr::case_when(
    series == "FDP(S=1)"              ~ "FDP(S=1) F",
    series == "CDP"                   ~ "CDP F",
    series == "CDP_unfair"            ~ "CDP NF",
    series == "Groupwise_noDP_noFair" ~ "NDP NF",
    series == "FAIR_no_DP"            ~ "NDP F",
    TRUE                              ~ series
  ))

# ---------------- summaries over runs ----------------
summ_df <- df_base %>%
  group_by(alpha, epsilon, series) %>%
  summarise(
    err_mean = mean(.data[[err_col]], na.rm = TRUE),
    err_sd   = sd(  .data[[err_col]], na.rm = TRUE),
    dd_mean  = mean(.data[[dd_col]],  na.rm = TRUE),
    dd_sd    = sd(  .data[[dd_col]],  na.rm = TRUE),
    mc = dplyr::n(),
    .groups = "drop"
  ) %>%
  mutate(
    err_se = ifelse(mc > 1, err_sd / sqrt(mc), 0),
    dd_se  = ifelse(mc > 1, dd_sd  / sqrt(mc), 0),
    floor_err = 0,
    floor_dd  = 0,
    err_lo = pmax(0, err_mean - pmax(1.96 * err_se, floor_err)),
    err_hi = pmin(1, err_mean + pmax(1.96 * err_se, floor_err)),
    dd_lo  = pmax(0, dd_mean  - pmax(1.96 * dd_se,  floor_dd)),
    dd_hi  = pmin(1, dd_mean  + pmax(1.96 * dd_se,  floor_dd))
  )

# split by alpha presence
summ_has_a <- summ_df %>% filter(!is.na(alpha))
summ_no_a  <- summ_df %>% filter(is.na(alpha)) 

# ---------- Long format ----------
df_err <- summ_has_a %>%
  transmute(
    alpha, epsilon, series,
    metric = factor("Error", levels = c("Error","Disparity")),
    value = err_mean, lower = err_lo, upper = err_hi
  )
df_dd <- summ_has_a %>%
  transmute(
    alpha, epsilon, series,
    metric = factor("Disparity", levels = c("Error","Disparity")),
    value = dd_mean, lower = dd_lo, upper = dd_hi
  )
df_long <- bind_rows(df_err, df_dd)
df_long$epsilon <- factor(df_long$epsilon, levels = eps_list)

# ---------- family (DP/NDP), fairness (F/NF), subtype (CDP/FDP/NDP) ----------
get_family <- function(s) ifelse(grepl("^NDP", s), "NDP", "DP")
get_status <- function(s) ifelse(grepl(" NF$", s), "NF", "F")
get_sub    <- function(s) dplyr::case_when(
  grepl("^CDP", s)                  ~ "CDP",
  grepl("^FDP\\(S=1\\)", s)         ~ "FDP(S=1)",
  grepl("^NDP", s)                  ~ "NDP",
  TRUE                              ~ "Other"
)

df_long <- df_long %>%
  mutate(
    family = factor(get_family(series), levels = c("DP","NDP")),
    status = factor(get_status(series), levels = c("F","NF")),
    sub    = factor(get_sub(series), levels = c("CDP","FDP(S=1)","NDP"))
  )


df_noa_lines <- NULL
df_noa_bands <- NULL
if (nrow(summ_no_a) > 0) {
  x_min <- suppressWarnings(min(df_long$alpha, na.rm = TRUE))
  x_max <- suppressWarnings(max(df_long$alpha, na.rm = TRUE))
  df_noa_lines <- bind_rows(
    summ_no_a %>% transmute(
      epsilon = factor(epsilon, levels = eps_list),
      series,
      metric = factor("Error", levels = c("Error","Disparity")),
      y = err_mean,
      y_lo = err_lo, y_hi = err_hi,
      x_start = x_min, x_end = x_max
    ),
    summ_no_a %>% transmute(
      epsilon = factor(epsilon, levels = eps_list),
      series,
      metric = factor("Disparity", levels = c("Error","Disparity")),
      y = dd_mean,
      y_lo = dd_lo, y_hi = dd_hi,
      x_start = x_min, x_end = x_max
    )
  ) %>%
    mutate(
      family = factor(get_family(series), levels = c("DP","NDP")),
      status = factor(get_status(series), levels = c("F","NF")),
      sub    = factor(get_sub(series),    levels = c("CDP","FDP(S=1)","NDP"))
    )
  
  df_noa_bands <- df_noa_lines %>%
    transmute(
      epsilon, series, metric, family, status, sub,
      xmin = x_start, xmax = x_end,
      ymin = pmax(0, y_lo), ymax = pmin(1, y_hi)
    )
}

# ---------- limits ----------
x_min <- suppressWarnings(min(df_long$alpha, na.rm = TRUE))
x_max <- suppressWarnings(max(df_long$alpha, na.rm = TRUE))
if (!is.finite(x_min) || !is.finite(x_max)) { x_min <- 0; x_max <- 1 }

err_vals <- c(df_err$value, if (!is.null(df_noa_lines)) df_noa_lines$y[df_noa_lines$metric=="Error"])
err_lim <- range(err_vals, df_err$lower, df_err$upper, na.rm = TRUE)
if (!all(is.finite(err_lim))) err_lim <- c(0, 1)

dd_vals <- c(df_dd$value, if (!is.null(df_noa_lines)) df_noa_lines$y[df_noa_lines$metric=="Disparity"])
dd_lim  <- range(dd_vals, df_dd$lower, df_dd$upper, na.rm = TRUE)
if (!all(is.finite(dd_lim))) dd_lim <- c(0, 1)
dd_lim <- dd_lim + c(-0.01, 0.01)

# ---------- Palettes + legend labels ----------
pal_family <- c(DP = "#1f77b4", NDP = "#e68a2e")  
family_breaks <- c("DP","NDP")
family_labels <- c("DP","Non-private")

sh_sub     <- c("CDP" = 16, "FDP(S=1)" = 17, "NDP" = 15)
sub_breaks <- c("CDP","FDP(S=1)","NDP")
sub_labels <- c("CDP","FDP(S=1)","Non-private")

lt_vals   <- c("F" = "solid", "NF" = "22")
lt_breaks <- c("F","NF")
lt_labels <- c("Fair","Unfair")

df_long <- df_long %>%
  mutate(
    lw = ifelse(sub %in% c("CDP","FDP(S=1)") & status == "F", 0.6, 1.3)
  )

if (!is.null(df_noa_lines) && nrow(df_noa_lines) > 0) {
  df_noa_lines <- df_noa_lines %>%
    mutate(
      lw = ifelse(sub %in% c("CDP","FDP(S=1)") & status == "F", 0.6, 1.3)
    )
}

# ======================== PLOT ===========================
p <- ggplot(
  df_long,
  aes(x = alpha, y = value,
      color = family, linetype = status, shape = sub,
      group = series)
) +
  geom_ribbon(
    aes(ymin = lower, ymax = upper, fill = family,
        group = interaction(series, epsilon, metric)),
    alpha = 0.25, colour = NA, show.legend = FALSE
  ) +
  { if (!is.null(df_noa_bands) && nrow(df_noa_bands) > 0)
    geom_rect(
      data = df_noa_bands,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = family),
      inherit.aes = FALSE,
      alpha = 0.25,
      colour = NA
    )
    else NULL } +
  geom_line(aes(linewidth = lw), na.rm = TRUE) +
  geom_point(size = 3.3, na.rm = TRUE) +
  { if (!is.null(df_noa_lines) && nrow(df_noa_lines) > 0)
    geom_segment(
      data = df_noa_lines,
      aes(x = x_start, xend = x_end, y = y, yend = y,
          color = family, linetype = status, linewidth = lw),
      inherit.aes = FALSE,
      lineend  = "butt"
    )
    else NULL } +
  geom_abline(
    data = data.frame(metric = factor("Disparity", levels = c("Error","Disparity"))),
    aes(slope = 1, intercept = 0),
    inherit.aes = FALSE,
    linetype = "dashed",
    color = "gray50"
  ) +
  # ---- Scales with updated labels ----
scale_color_manual(values = pal_family,
                   breaks = family_breaks, labels = family_labels,
                   name = "Privacy") +
  scale_fill_manual(values = pal_family, guide = "none") +
  scale_linetype_manual(values = lt_vals,
                        breaks = lt_breaks, labels = lt_labels,
                        name = "Fairness") +
  scale_shape_manual(values = sh_sub,
                     breaks = sub_breaks, labels = sub_labels,
                     name = "Privacy type") +
  scale_linewidth_identity(guide = "none") +  
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
        limits = err_lim, breaks = scales::pretty_breaks(n = 6),
        minor_breaks = NULL, oob = scales::oob_squish
      ),
      metric == "Disparity" ~ scale_y_continuous(
        limits = dd_lim, breaks = scales::pretty_breaks(n = 6),
        minor_breaks = NULL, oob = scales::oob_squish
      )
    )
  ) +
  scale_x_continuous(
    limits = c(0.025, 0.4),
    breaks = seq(0, 0.4, 0.1),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  labs(x = bquote(bold(alpha)), y = NULL) +
  theme_bw(base_size = 13) +
  theme(
    legend.key.width = unit(1.2, "cm"),
    legend.position   = "right",
    strip.background  = element_rect(fill = "white", colour = NA),
    strip.placement   = "outside",
    strip.text.y.left = element_text(angle = 90),
    strip.text.x      = element_text(face = "bold", color = "black"),
    panel.grid.major  = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.minor  = element_blank()
  )

# Make legend samples clear
p <- p +
  guides(
    color    = guide_legend(override.aes = list(linetype = "solid", shape = NA, linewidth = 1)),
    linetype = guide_legend(
      override.aes = list(color = "black", shape = NA, linewidth = 1,
                          linetype = c("solid","22"))
    ),
    shape    = guide_legend(override.aes = list(linetype = "blank"))
  )

p <- p + theme(
  strip.text.x      = element_text(size = 18, color = "black"),
  strip.text.y.left = element_text(size = 16, angle = 90)
)

p <- p +
  theme(
    # ----- x axis title (alpha) -----
    axis.title.x = element_text(size = 16),
    axis.text.x  = element_text(size = 12),   
    axis.text.y  = element_text(size = 12), 
    
    # ----- legend text -----
    legend.title = element_text(size = 16),
    legend.text  = element_text(size = 14),
    
    # ----- legend key sizing -----
    legend.key.size  = unit(0.8, "cm"),
    legend.key.width = unit(1.8, "cm")
  ) +
  guides(
    color    = guide_legend(override.aes = list(linetype = "solid", shape = NA, linewidth = 1.6)),
    linetype = guide_legend(override.aes = list(color = "black", shape = NA, linewidth = 1.6,
                                                linetype = c("solid","22"))),
    shape    = guide_legend(override.aes = list(linetype = "blank"))
  )
print(p)
