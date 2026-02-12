#!/usr/bin/env Rscript

# -------------------- Environment --------------------
suppressPackageStartupMessages({ library(mvtnorm); library(RSpectra) })
Sys.setenv(
  OMP_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  OPENBLAS_CORETYPE = "NEHALEM"
)

# -------------------- Array index --------------------
args <- commandArgs(trailingOnly = TRUE)
run_idx <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", NA))
if (is.na(run_idx) && length(args)) run_idx <- as.integer(args[1])
if (is.na(run_idx) || run_idx < 1) stop("run_idx not provided")

# -------------------- Methods --------------------
source("method-function-FDP.R")

# -------------------- Settings --------------------
S_list       <- 3                     
n_s          <- 2000                  
n_test       <- 2000
alpha_list   <- 0.3                
epsilon_list <- c(0.75, 1, 2, 3, 4, 5)
h_list       <- seq(0.1, 0.35, 0.05)

# Sensitivity grid
Cw_list <- seq(0.06, 0.5, 0.02)

# Data-gen parameters
d      <- 2
p_a    <- 0.3
alpha1 <- 4; beta1 <- 2
alpha0 <- 4.5; beta0 <- 2
w      <- c(1,1); b <- 0; gamma <- -0.3; k <- 12

# FDP constants
C_K   <- 1
rho   <- 0.03       
omega <- NULL
eta   <- 1e-2

# Choose tree depth M 
choose_M <- function(n_vec, eps_vec, M_min = 1, M_cap = 6) {
  val <- sum(pmin(n_vec, (n_vec^2) * (eps_vec^2)))
  M <- floor(log2(max(val, 1)))
  M <- max(as.integer(M), as.integer(M_min))
  M <- min(M, as.integer(M_cap))
  M
}

# -------------------- Paths --------------------
submit_dir <- Sys.getenv("SLURM_SUBMIT_DIR")
base_dir <- if (nzchar(submit_dir)) submit_dir else getwd()
outdir <- file.path(base_dir, "results")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Per-task TSV 
out_file <- file.path(outdir, sprintf("cw_sensitivity_run_%03d.txt", run_idx))
header <- "run\tS\tn_s\talpha\tepsilon\th_cv\tM\tCw\terr\tabs_DD"
writeLines(header, con = out_file)

append_row <- function(df_row, file) {
  stopifnot(nrow(df_row) == 1)
  write.table(df_row, file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)
}

# -------------------- One array-task run --------------------
set.seed(3000 + run_idx)

# Fixed test_dat
test_dat <- generate_site_data(
  n = n_test, d = d, p_a = p_a,
  alpha1 = alpha1, beta1 = beta1, alpha0 = alpha0, beta0 = beta0,
  eta_fun = eta_arctan, w = w, b = b, gamma = gamma, k = k
)

for (S in S_list) {
  sites <- generate_sites(
    S = S, n_per_site = n_s, d = d, p_a = p_a, split_ratio = 0.5,
    alpha1 = alpha1, beta1 = beta1, alpha0 = alpha0, beta0 = beta0,
    eta_fun = eta_arctan, w = w, b = b, gamma = gamma, k = k
  )
  
  cv_out <- cv_select_h(data = sites, h_grid = h_list, n_folds = 3,
                        kernel_fun = gaussian_kernel, tol = 1e-12)
  h_cv <- cv_out$h_star
  
  
  delta_vec <- rep((n_s/2)^(-2), S)  
  mu <- rep(1/S, S)
  nu <- rep(1/S, S)
  
  for (Cw_val in Cw_list) {
    for (eps_priv in epsilon_list) {
      eps_vec <- rep(eps_priv, S)
      M_dyn   <- choose_M(n_vec = rep(n_s, S), eps_vec = eps_vec)
      
      for (a_lvl in alpha_list) {
        out_fd <- try(
          FDP_fair_classifier_cross_fit(
            Sites_train = sites$train, Sites_cal = sites$cal,
            alpha = a_lvl, h = h_cv, eps = eps_vec, delta = delta_vec,
            mu = mu, nu = nu,
            C_K = C_K, Cw = Cw_val, rho = rho, M = M_dyn, omega = omega, eta = eta,
            pi_method = "empirical"
          ),
          silent = TRUE
        )
        
        if (!inherits(out_fd, "try-error")) {
          y_hat <- out_fd$predict(test_dat$X, test_dat$A)
          err   <- mean(y_hat != test_dat$Y)
          absDD <- empirical_DD(y_hat, test_dat$A)$abs_DD
          
          append_row(data.frame(
            run = run_idx,
            S = S,
            n_s = n_s,
            alpha = a_lvl,
            epsilon = eps_priv,
            h_cv = h_cv,
            M = M_dyn,
            Cw = Cw_val,
            err = err,
            abs_DD = absDD
          ), out_file)
        } else {
          msg <- if (!is.null(attr(out_fd, "condition"))) attr(out_fd, "condition")$message else "unknown"
          cat(sprintf("[warn] Cw run=%d S=%d eps=%.3g alpha=%.3g Cw=%.3g failed: %s\n",
                      run_idx, S, eps_priv, a_lvl, Cw_val, msg), file = stderr())
        }
      }
    }
  }
}

cat("[done] Cw sensitivity wrote to:", out_file, "\n")
