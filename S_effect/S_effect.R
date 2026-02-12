# -------------------- Environment --------------------
suppressPackageStartupMessages({
  library(mvtnorm)
  suppressWarnings(requireNamespace("RSpectra", quietly = TRUE))
})
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
source("method-function-combined.R")

# -------------------- Settings --------------------
S_list       <- c(1, 2, 3, 4, 5, 6)
n_s          <- 7200
n_test       <- 2500
alpha_list   <- seq(0.05, 1, 0.05)
epsilon_list <- c(0.75, 1, 2, 3, 4, 5)
h_list       <- seq(0.1, 0.35, 0.05)

# Data-gen parameters
d      <- 2
p_a    <- 0.3
alpha1 <- 4; beta1 <- 2
alpha0 <- 4.5; beta0 <- 2
w      <- c(1,1); b <- 0; gamma <- -0.3; k <- 12

# FDP constants
C_K  <- 1
Cw   <- 0.1
rho  <- 0.03
omega <- NULL
eta <- 1e-2

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

# Unified output file (append rows)
out_file <- file.path(outdir, sprintf("run_%03d.txt", run_idx))
header <- "run\tmethod\tS\tn_s\talpha\tepsilon\th_cv\tM\terr\tabs_DD\tchosen_site"
writeLines(header, con = out_file)

append_row <- function(df_row, file) {
  stopifnot(nrow(df_row) == 1)
  write.table(df_row, file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)
}

# -------------------- One simulation (array) --------------------
set.seed(2000 + run_idx)

test_dat <- generate_site_data(
  n = n_test, d = d, p_a = p_a,
  alpha1 = alpha1, beta1 = beta1, alpha0 = alpha0, beta0 = beta0,
  eta_fun = eta_arctan, w = w, b = b, gamma = gamma, k = k
)

for (S in S_list) {
  n_site <- ceiling(n_s/S)
  sites <- generate_sites(
    S = S, n_per_site = n_site, d = d, p_a = p_a, split_ratio = 0.5,
    alpha1 = alpha1, beta1 = beta1, alpha0 = alpha0, beta0 = beta0,
    eta_fun = eta_arctan, w = w, b = b, gamma = gamma, k = k
  )
  
  # CV bandwidth
  cv_out <- cv_select_h(data = sites, h_grid = h_list, n_folds = 3,
                        kernel_fun = gaussian_kernel, tol = 1e-12,
                        site_index = NULL)
  h_cv <- cv_out$h_star
  chosen_site <- cv_out$chosen_site
  
  delta_vec <- rep((n_site/2)^(-2), S)
  mu <- rep(1/S, S)
  nu <- rep(1/S, S)
  
  for (eps_priv in epsilon_list) {
    eps_vec <- rep(eps_priv, S)
    M_dyn <- choose_M(n_vec = rep(n_site, S), eps_vec = eps_vec)
    
    for (a_lvl in alpha_list) {
      out_fd <- try(
        FDP_fair_classifier_cross_fit(
          Sites_train = sites$train, Sites_cal = sites$cal,
          alpha = a_lvl, h = h_cv, eps = eps_vec, delta = delta_vec,
          mu = mu, nu = nu,
          C_K = C_K, Cw = Cw, rho = rho, M = M_dyn, omega = omega, eta = eta,
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
          method = "FDP",
          S = S,
          n_s = n_s,
          alpha = a_lvl,
          epsilon = eps_priv,
          h_cv = h_cv,
          M = M_dyn,
          err = err,
          abs_DD = absDD,
          chosen_site = ifelse(is.na(chosen_site), -1L, chosen_site)
        ), out_file)
        
        # CDP 
        if (S == 1) {
          X_all <- sites$all[[1]]$X
          Y_all <- sites$all[[1]]$Y
          A_all <- sites$all[[1]]$A
          input_cdp <- as_cdp_input(X_all, Y_all, A_all)
          delta_cdp <- (n_site/2)^(-2)
          
          out_cdp <- try(
            CDP_fair_classifier(
              test_data  = list(X = test_dat$X, A = test_dat$A, Y = test_dat$Y),
              input_data = input_cdp,
              alpha = a_lvl, h = h_cv, epsilon = eps_priv, delta = delta_cdp,
              kernel_fun = gaussian_kernel, tol = 1e-8, C_K = C_K,
              cross_fit = TRUE, pi_mode = "empirical",
              tau_method = "directional"
            ),
            silent = TRUE
          )
          if (!inherits(out_cdp, "try-error")) {
            err_cdp <- mean(out_cdp$y_hat_out != test_dat$Y)
            absDD_cdp <- as.numeric(out_cdp$DD_hat_out)
            append_row(data.frame(
              run = run_idx,
              method = "CDP",
              S = 1,
              n_s = n_s,       
              alpha = a_lvl,
              epsilon = eps_priv,
              h_cv = h_cv,
              M = NA_real_,        
              err = err_cdp,
              abs_DD = absDD_cdp,
              chosen_site = -1   
            ), out_file)
          }
        }
        
      } else {
        msg <- if (!is.null(attr(out_fd, "condition"))) attr(out_fd, "condition")$message else "unknown"
        cat(sprintf("[warn] run=%d S=%d eps=%.3g alpha=%.3g failed: %s\n",
                    run_idx, S, eps_priv, a_lvl, msg), file = stderr())
      }
    }
  }
}

cat("[done] wrote to:", out_file, "\n")
