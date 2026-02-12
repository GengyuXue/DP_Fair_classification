# -------------------- Environment --------------------
suppressPackageStartupMessages({ library(mvtnorm); library(RSpectra) })
Sys.setenv(
  OMP_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  OPENBLAS_CORETYPE = "NEHALEM" 
)

# -------------------- Settings --------------------
p_a <- 0.3
n_train_list <- c(5000,7000,9000)
n_test <- 4000
alpha_list <- seq(0.05, 1, 0.05)
epsilon_list <- c(0.75, 1, 2, 3, 4)
h_list <- seq(0.1, 0.35, 0.05)
d <- 2
alpha1 <- 4; beta1 <- 2; alpha0 <- 4.5; beta0 <- 2
w <- c(1,1); b <- 0; gamma <- -0.3; k <- 12; C_K <- 1

# -------------------- Array index --------------------
args <- commandArgs(trailingOnly = TRUE)
run_idx <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", NA))
if (is.na(run_idx) && length(args)) run_idx <- as.integer(args[1])
if (is.na(run_idx) || run_idx < 1) stop("run_idx not provided")

# -------------------- Methods --------------------
source("method-function-CDP.R")

# -------------------- Paths --------------------
submit_dir <- Sys.getenv("SLURM_SUBMIT_DIR")
base_dir <- if (nzchar(submit_dir)) submit_dir else getwd()
outdir <- file.path(base_dir, "results")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Prepare per-task output file
run_file <- file.path(outdir, sprintf("run_%03d.txt", run_idx))
header_line <- "run\tn_train\talpha\tepsilon\th_cv\trisk\tdisparity"
writeLines(header_line, con = run_file)


append_row <- function(df_row, file) {
  stopifnot(nrow(df_row) == 1)
  write.table(
    df_row, file,
    sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE
  )
}

# -------------------- One simulation --------------------
set.seed(2000 + run_idx)

test_data <- generate_data(
  n_test, d, p_a, alpha1, beta1, alpha0, beta0,
  eta_fun = eta_arctan, w = w, b = b, gamma = gamma, k = k
)

for (n_train in n_train_list) {
  delta <- (n_train/2)^(-2)
  
  input_data <- generate_data(
    n_train, d, p_a, alpha1, beta1, alpha0, beta0,
    eta_fun = eta_arctan, w = w, b = b, gamma = gamma, k = k
  )
  
  # CV h
  cv_out <- cv_select_h(
    list(X = input_data$X, Y = input_data$Y),
    h_grid = h_list, n_folds = 3, kernel_fun = gaussian_kernel, tol = 1e-12
  )
  h_cv <- cv_out$h_star
  
  for (eps_priv in epsilon_list) {
    for (a_lvl in alpha_list) {
      out_cf <- try(
        CDP_fair_classifier(
          test_data, input_data,
          alpha = a_lvl, h = h_cv, epsilon = eps_priv, delta = delta,
          kernel_fun = gaussian_kernel, tol = 1e-8, C_K = C_K,
          cross_fit = TRUE, pi_mode = "empirical", tau_method = "directional"
        ),
        silent = TRUE
      )
      
      if (!inherits(out_cf, "try-error")) {
        risk <- mean(out_cf$y_hat_out != test_data$Y)
        disparity <- as.numeric(out_cf$DD_hat_out)
        row_now <- data.frame(
          run = run_idx,
          n_train = n_train,
          alpha = a_lvl,
          epsilon = eps_priv,
          h_cv = h_cv,
          risk = risk,
          disparity = disparity
        )
        # append 
        append_row(row_now, run_file)
        
      } else {
        msg <- if (!is.null(attr(out_cf, "condition"))) attr(out_cf, "condition")$message else "unknown"
        cat(sprintf("[warn] n_train=%d eps=%.3g alpha=%.3g failed: %s\n",
                    n_train, eps_priv, a_lvl, msg), file = stderr())
      }
    }
  }
}

cat("[done] wrote incrementally to:", run_file, "\n")
