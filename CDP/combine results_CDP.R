outdir <- getwd()

if (!dir.exists(outdir)) stop("Results folder not found: ", normalizePath(outdir))

cat("Combining from: ", normalizePath(outdir), "\n")

# Find per-run files 
files <- list.files(outdir, pattern = "^run_\\d+\\.txt$", full.names = TRUE)
if (!length(files)) {
  files <- list.files(outdir, pattern = "\\.txt$", full.names = TRUE)
}
if (!length(files)) stop("No TXT files found in: ", normalizePath(outdir))

# Read and validate columns
required_cols <- c("run","n_train","alpha","epsilon","h_cv","risk","disparity")
dfs <- list()
bad <- 0L
for (f in files) {
  df <- try(read.table(f, header = TRUE, sep = "\t", check.names = FALSE), silent = TRUE)
  if (inherits(df, "try-error") || !nrow(df)) { bad <- bad + 1L; next }
  if (!all(required_cols %in% names(df))) { bad <- bad + 1L; next }
  

  df$run       <- as.integer(df$run)
  df$n_train   <- as.integer(df$n_train)
  df$alpha     <- as.numeric(df$alpha)
  df$epsilon   <- as.numeric(df$epsilon)
  df$h_cv      <- as.numeric(df$h_cv)
  df$risk      <- as.numeric(df$risk)
  df$disparity <- as.numeric(df$disparity)
  
  dfs[[length(dfs)+1L]] <- df
}
if (!length(dfs)) stop("No valid result files with required columns in: ", normalizePath(outdir))
if (bad) cat("Skipped", bad, "files that lacked required columns or failed to parse.\n")

all_df <- do.call(rbind, dfs)

# Infer grids
mc             <- max(all_df$run, na.rm = TRUE)
n_train_list   <- sort(unique(all_df$n_train))
alpha_list     <- sort(unique(all_df$alpha))
epsilon_list   <- sort(unique(all_df$epsilon))

cat("Detected:\n  runs (mc): ", mc,
    "\n  n_train:   ", paste(n_train_list, collapse=", "),
    "\n  alpha:     ", paste(alpha_list, collapse=", "),
    "\n  epsilon:   ", paste(epsilon_list, collapse=", "), "\n", sep="")

# Allocate arrays
errors_arr <- array(NA_real_,
                    dim = c(mc, length(n_train_list), length(alpha_list), length(epsilon_list)),
                    dimnames = list(paste0("run", 1:mc),
                                    paste0("ntr=", n_train_list),
                                    paste0("alpha=", alpha_list),
                                    paste0("eps=", epsilon_list)))
DD_arr <- array(NA_real_, dim = dim(errors_arr), dimnames = dimnames(errors_arr))
hsel_mat <- matrix(NA_real_, nrow = mc, ncol = length(n_train_list),
                   dimnames = list(paste0("run", 1:mc), paste0("ntr=", n_train_list)))

# Fill arrays
for (i in seq_len(nrow(all_df))) {
  r  <- all_df$run[i]
  if (is.na(r) || r < 1 || r > mc) next
  tidx <- match(all_df$n_train[i], n_train_list)
  ai   <- match(all_df$alpha[i], alpha_list)
  ei   <- match(all_df$epsilon[i], epsilon_list)
  if (any(is.na(c(tidx, ai, ei)))) next
  
  errors_arr[r, tidx, ai, ei] <- all_df$risk[i]
  DD_arr[r,     tidx, ai, ei] <- all_df$disparity[i]
  hsel_mat[r, tidx] <- all_df$h_cv[i]
}

# Summaries
mean_error <- apply(errors_arr, c(2,3,4), mean, na.rm = TRUE)
sd_error   <- apply(errors_arr,   c(2,3,4), sd,   na.rm = TRUE)
mean_DD    <- apply(DD_arr,       c(2,3,4), mean, na.rm = TRUE)
sd_DD      <- apply(DD_arr,       c(2,3,4), sd,   na.rm = TRUE)

cat("Summary shapes:\n  mean_error: ", paste(dim(mean_error), collapse=" x"),
    "\n  mean_DD:    ", paste(dim(mean_DD), collapse=" x"), "\n", sep="")

# Save combined artifacts
timestamp <- format(Sys.time(), "%Y%m%d-%H%M%S")

combined_rds <- file.path(outdir, paste0("combined_runs_", timestamp, ".rds"))

saveRDS(all_df, combined_rds)


