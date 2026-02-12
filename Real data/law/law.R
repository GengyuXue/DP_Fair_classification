#!/usr/bin/env Rscript

# -------------------- Environment --------------------
suppressWarnings(suppressMessages({
}))

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
source("method-fair-no_dp.R")
source("method-dp-no-fair-groupwise.R")
source("method_no_fair_no_dp.R")

# -------------------- Settings --------------------
alpha_list   <- seq(0.025, 0.4, 0.025)
epsilon_list <- c(0.75, 1, 2, 3)

h_grid <- c(0.2, 0.3)
tol    <- 1e-8
C_K    <- 1

# FDP (S=1) tuning
M_fdp    <- 6
Cw_fdp   <- 0.1
eta_iso  <- 1e-2
mu_vec_1 <- 1
nu_vec_1 <- 1
rho      <- 0.02

# -------------------- Paths --------------------
submit_dir <- Sys.getenv("SLURM_SUBMIT_DIR")
base_dir <- if (nzchar(submit_dir)) submit_dir else getwd()
setwd(base_dir)

outdir <- file.path(base_dir, "results")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
run_file <- file.path(outdir, sprintf("law_unfair_run_%03d.txt", run_idx))
header_line <- "run\tmethod\tepsilon\talpha\terror\tdisparity\th_cv"
writeLines(header_line, con = run_file)

append_row <- function(df_row, file) {
  stopifnot(nrow(df_row) == 1)
  write.table(
    df_row, file,
    sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE
  )
}

# -------------------- Data helpers --------------------
read_csv_base <- function(path) {
  utils::read.csv(path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
}

prepare_cdp_full <- function(data, scale_X = FALSE) {
  req <- c("LSAT","GPA","resident","Sensitive","Label")
  missing_cols <- setdiff(req, colnames(data))
  if (length(missing_cols)) stop("Missing columns: ", paste(missing_cols, collapse = ", "))
  df <- data[, req, drop = FALSE]
  df <- df[stats::complete.cases(df), , drop = FALSE]
  X <- as.matrix(df[, c("LSAT","GPA","resident"), drop = FALSE])
  if (scale_X) X <- scale(X)
  list(
    X = X,
    A = as.integer(df[["Sensitive"]]),
    Y = as.integer(df[["Label"]])
  )
}

# -------------------- One run (this array task) --------------------
set.seed(1000 + run_idx)

# --- load data ---
dat_path <- file.path(base_dir, "law_unfair_40000.csv")
if (!file.exists(dat_path)) stop("Cannot find law_unfair_40000.csv in working directory")
law_unfair <- read_csv_base(dat_path)

# one split
n_all   <- nrow(law_unfair)
idx_all <- seq_len(n_all)
idx_tr  <- sample(idx_all, size = ceiling(0.70 * n_all), replace = FALSE)
idx_te  <- setdiff(idx_all, idx_tr)

tr_full   <- prepare_cdp_full(law_unfair[idx_tr, ],  FALSE)
te_full   <- prepare_cdp_full(law_unfair[idx_te, ],  FALSE)
test_data <- list(X = te_full$X, A = te_full$A, Y = te_full$Y)

n_train   <- nrow(tr_full$X)
delta_run <- (n_train/2)^(-2)

# CV for h
cv_out <- cv_select_h(
  data = tr_full,
  h_grid = h_grid,
  n_folds = 3,
  kernel_fun = gaussian_kernel,
  tol = 1e-12
)
h_cv <- as.numeric(cv_out$h_star)
cat(sprintf("[info] run=%d Chosen h: %g (CV err: %.6g)\n", run_idx, h_cv, cv_out$cv_error))

# -------------------- RUNS --------------------
# ---- FAIR_no_DP  ----
for (alp in alpha_list) {
  fair_fit <- try(
    fair_classifier(
      test_data  = test_data,
      input_data = as_cdp_input(tr_full$X, tr_full$Y, tr_full$A),
      alpha      = alp,
      h          = h_cv,
      kernel_fun = gaussian_kernel,
      tol        = tol,
      C_K        = C_K,
      cross_fit  = TRUE
    ),
    silent = TRUE
  )
  if (!inherits(fair_fit, "try-error") && !is.null(fair_fit$y_hat_out)) {
    yhat <- as.integer(fair_fit$y_hat_out)
    row_now <- data.frame(
      run       = run_idx,
      method    = "FAIR_no_DP",
      epsilon   = NA_real_,
      alpha     = alp,
      error     = mean(yhat != test_data$Y),
      disparity = empirical_DD(yhat, test_data$A)$abs_DD,
      h_cv      = h_cv
    )
    append_row(row_now, run_file)
  } else {
    row_now <- data.frame(
      run=run_idx, method="FAIR_no_DP", epsilon=NA_real_, alpha=alp,
      error=NA_real_, disparity=NA_real_, h_cv=h_cv
    )
    append_row(row_now, run_file)
    msg <- if (!is.null(attr(fair_fit, "condition"))) attr(fair_fit, "condition")$message else "unknown"
    cat(sprintf("[warn] FAIR_no_DP failed (alpha=%.3g): %s\n", alp, msg), file = stderr())
  }
}

# ---- Groupwise_noDP_noFair ----
sp <- split_train_cal(list(X=tr_full$X, A=tr_full$A, Y=tr_full$Y), split_ratio = 0.5)
half_train <- sp$train
yhat_g <- try(classifier_xya(half_train, test_data, h = h_cv, kernel_fun = gaussian_kernel), silent = TRUE)
if (!inherits(yhat_g, "try-error")) {
  row_g <- data.frame(
    run       = run_idx,
    method    = "Groupwise_noDP_noFair",
    epsilon   = NA_real_,
    alpha     = NA_real_,
    error     = mean(yhat_g != test_data$Y),
    disparity = empirical_DD(yhat_g, test_data$A)$abs_DD,
    h_cv      = h_cv
  )
  append_row(row_g, run_file)
} else {
  row_g <- data.frame(
    run=run_idx, method="Groupwise_noDP_noFair", epsilon=NA_real_, alpha=NA_real_,
    error=NA_real_, disparity=NA_real_, h_cv=h_cv
  )
  append_row(row_g, run_file)
  emsg_g <- if (!is.null(attr(yhat_g, "condition"))) attr(yhat_g, "condition")$message else "unknown"
  cat(sprintf("[warn] classifier_xya failed: %s\n", emsg_g), file = stderr())
}



for (eps in epsilon_list) {
  
  # CDP_unfair 
  cdp_unfair_fit <- try(
    CDP_unfair_classifier_groupwise(
      test_data  = test_data,
      input_data = as_cdp_input(tr_full$X, tr_full$Y, tr_full$A),
      h = h_cv, epsilon = eps, delta = delta_run,
      kernel_fun = gaussian_kernel, tol = tol, C_K = C_K, pi_mode = "empirical"
    ),
    silent = TRUE
  )
  if (!inherits(cdp_unfair_fit, "try-error") && !is.null(cdp_unfair_fit$y_hat)) {
    yhat_ud <- as.integer(cdp_unfair_fit$y_hat)
    row_now <- data.frame(
      run=run_idx, method="CDP_unfair", epsilon=eps, alpha=NA_real_,
      error=mean(yhat_ud != test_data$Y),
      disparity=empirical_DD(yhat_ud, test_data$A)$abs_DD,
      h_cv=h_cv
    )
    append_row(row_now, run_file)
  } else {
    row_now <- data.frame(
      run=run_idx, method="CDP_unfair", epsilon=eps, alpha=NA_real_,
      error=NA_real_, disparity=NA_real_, h_cv=h_cv
    )
    append_row(row_now, run_file)
    msg <- if (!is.null(attr(cdp_unfair_fit, "condition"))) attr(cdp_unfair_fit, "condition")$message else "unknown"
    cat(sprintf("[warn] CDP_unfair failed (eps=%.3g): %s\n", eps, msg), file = stderr())
  }
  
  for (alp in alpha_list) {
    # CDP
    cdp_fit <- try(
      CDP_fair_classifier(
        test_data  = test_data,
        input_data = list(
          X = tr_full$X, Y = tr_full$Y, A = tr_full$A,
          X_11 = tr_full$X[tr_full$A==1 & tr_full$Y==1,,drop=FALSE],
          X_10 = tr_full$X[tr_full$A==1 & tr_full$Y==0,,drop=FALSE],
          X_01 = tr_full$X[tr_full$A==0 & tr_full$Y==1,,drop=FALSE],
          X_00 = tr_full$X[tr_full$A==0 & tr_full$Y==0,,drop=FALSE]
        ),
        alpha   = alp, h = h_cv, epsilon = eps, delta = delta_run,
        kernel_fun = gaussian_kernel, tol = tol, C_K = C_K, cross_fit = TRUE,
        pi_mode = "empirical", tau_method = "directional"
      ),
      silent = TRUE
    )
    if (!inherits(cdp_fit, "try-error") && !is.null(cdp_fit$y_hat_out)) {
      yhat_cdp <- as.integer(cdp_fit$y_hat_out)
      row_now <- data.frame(
        run=run_idx, method="CDP", epsilon=eps, alpha=alp,
        error=mean(yhat_cdp != test_data$Y),
        disparity=empirical_DD(yhat_cdp, test_data$A)$abs_DD,
        h_cv=h_cv
      )
      append_row(row_now, run_file)
    } else {
      row_now <- data.frame(
        run=run_idx, method="CDP", epsilon=eps, alpha=alp,
        error=NA_real_, disparity=NA_real_, h_cv=h_cv
      )
      append_row(row_now, run_file)
      msg <- if (!is.null(attr(cdp_fit, "condition"))) attr(cdp_fit, "condition")$message else "unknown"
      cat(sprintf("[warn] CDP failed (eps=%.3g, alpha=%.3g): %s\n", eps, alp, msg), file = stderr())
    }
    
    # ----- FDP (S=1)
    split_50 <- split_train_cal(list(X=tr_full$X, A=tr_full$A, Y=tr_full$Y), split_ratio = 0.5)
    fdp_fit <- try(
      FDP_fair_classifier_cross_fit(
        Sites_train = list(split_50$train),
        Sites_cal   = list(split_50$cal),
        alpha = alp, h = h_cv, eps = c(eps), delta = c(delta_run),
        mu = c(mu_vec_1), nu = c(nu_vec_1),
        C_K = C_K, Cw = Cw_fdp, M = M_fdp, rho = rho, omega = NULL, eta = eta_iso,
        pi_method = "empirical"
      ),
      silent = TRUE
    )
    if (!inherits(fdp_fit, "try-error") && !is.null(fdp_fit$predict)) {
      yhat_fdp <- as.integer(fdp_fit$predict(test_data$X, test_data$A))
      row_now <- data.frame(
        run=run_idx, method="FDP(S=1)", epsilon=eps, alpha=alp,
        error=mean(yhat_fdp != test_data$Y),
        disparity=empirical_DD(yhat_fdp, test_data$A)$abs_DD,
        h_cv=h_cv
      )
      append_row(row_now, run_file)
    } else {
      row_now <- data.frame(
        run=run_idx, method="FDP(S=1)", epsilon=eps, alpha=alp,
        error=NA_real_, disparity=NA_real_, h_cv=h_cv
      )
      append_row(row_now, run_file)
      msg <- if (!is.null(attr(fdp_fit, "condition"))) attr(fdp_fit, "condition")$message else "unknown"
      cat(sprintf("[warn] FDP(S=1) failed (eps=%.3g, alpha=%.3g): %s\n", eps, alp, msg), file = stderr())
    }
    
    
  }
}

cat("[done] wrote incrementally to: ", run_file, "\n")
