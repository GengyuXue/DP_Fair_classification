# -------------------- Combine FDP run outputs --------------------
outdir <- getwd()

if (!dir.exists(outdir)) stop("Results folder not found: ", normalizePath(outdir))
cat("Combining from: ", normalizePath(outdir), "\n")

# Find per-run files
files <- list.files(outdir, pattern = "^run_\\d+\\.txt$", full.names = TRUE)
if (!length(files)) {
  files <- list.files(outdir, pattern = "\\.txt$", full.names = TRUE)
}
if (!length(files)) stop("No TXT files found in: ", normalizePath(outdir))

# Required columns
required_cols <- c("run","S","n_s","alpha","epsilon","h_cv","M","err","abs_DD")

dfs <- list()
bad <- 0L
for (f in files) {
  df <- try(read.table(f, header = TRUE, sep = "\t", check.names = FALSE), silent = TRUE)
  if (inherits(df, "try-error") || !nrow(df)) { bad <- bad + 1L; next }
  if (!all(required_cols %in% names(df))) { bad <- bad + 1L; next }
  
  # Coerce types
  df$run     <- as.integer(df$run)
  df$S       <- as.integer(df$S)
  df$n_s     <- as.integer(df$n_s)
  df$alpha   <- as.numeric(df$alpha)
  df$epsilon <- as.numeric(df$epsilon)
  df$h_cv    <- as.numeric(df$h_cv)
  df$M       <- as.integer(df$M)
  df$err     <- as.numeric(df$err)
  df$abs_DD  <- as.numeric(df$abs_DD)
  
  dfs[[length(dfs)+1L]] <- df
}
if (!length(dfs)) stop("No valid result files with required columns in: ", normalizePath(outdir))
if (bad) cat("Skipped", bad, "files that lacked required columns or failed to parse.\n")

all_df <- do.call(rbind, dfs)

# -------------------- Infer grids --------------------
mc           <- max(all_df$run, na.rm = TRUE)
S_list       <- sort(unique(all_df$S))
n_s_list     <- sort(unique(all_df$n_s))  
alpha_list   <- sort(unique(all_df$alpha))
epsilon_list <- sort(unique(all_df$epsilon))

cat("Detected:\n  runs (mc): ", mc,
    "\n  S:         ", paste(S_list, collapse=", "),
    "\n  n_s:       ", paste(n_s_list, collapse=", "),
    "\n  alpha:     ", paste(alpha_list, collapse=", "),
    "\n  epsilon:   ", paste(epsilon_list, collapse=", "), "\n", sep="")

# -------------------- Allocate arrays --------------------
err_arr <- array(NA_real_,
                 dim = c(mc, length(S_list), length(alpha_list), length(epsilon_list)),
                 dimnames = list(paste0("run", 1:mc),
                                 paste0("S=", S_list),
                                 paste0("alpha=", alpha_list),
                                 paste0("eps=", epsilon_list)))
DD_arr  <- array(NA_real_, dim = dim(err_arr), dimnames = dimnames(err_arr))

hsel_mat <- matrix(NA_real_, nrow = mc, ncol = length(S_list),
                   dimnames = list(paste0("run", 1:mc), paste0("S=", S_list)))
M_mat    <- matrix(NA_integer_, nrow = mc, ncol = length(S_list),
                   dimnames = dimnames(hsel_mat))

# -------------------- Fill arrays --------------------
for (i in seq_len(nrow(all_df))) {
  r  <- all_df$run[i]
  if (is.na(r) || r < 1 || r > mc) next
  si <- match(all_df$S[i],       S_list)
  ai <- match(all_df$alpha[i],   alpha_list)
  ei <- match(all_df$epsilon[i], epsilon_list)
  if (any(is.na(c(si, ai, ei)))) next
  
  err_arr[r, si, ai, ei] <- all_df$err[i]
  DD_arr[r,  si, ai, ei] <- all_df$abs_DD[i]
  
  # keep last-seen h and M for that (run,S)
  hsel_mat[r, si] <- all_df$h_cv[i]
  M_mat[r,    si] <- all_df$M[i]
}

# -------------------- Summaries --------------------
mean_err <- apply(err_arr, c(2,3,4), mean, na.rm = TRUE)
sd_err   <- apply(err_arr, c(2,3,4), sd,   na.rm = TRUE)
mean_DD  <- apply(DD_arr,  c(2,3,4), mean, na.rm = TRUE)
sd_DD    <- apply(DD_arr,  c(2,3,4), sd,   na.rm = TRUE)

cat("Summary shapes:\n  mean_err: ", paste(dim(mean_err), collapse=" x"),
    "\n  mean_DD:  ", paste(dim(mean_DD),  collapse=" x"), "\n", sep="")

# -------------------- Save combined artifacts --------------------
timestamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
combined_rds   <- file.path(outdir, paste0("combined_FDP_runs_", timestamp, ".rds"))

saveRDS(all_df, combined_rds)
