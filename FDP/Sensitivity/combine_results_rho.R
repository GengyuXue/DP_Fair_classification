outdir <- getwd()
if (!dir.exists(outdir)) stop("Results folder not found: ", normalizePath(outdir))
cat("Combining from: ", normalizePath(outdir), "\n")

## -------------------- Discover files --------------------
files <- list.files(outdir, pattern = "^rho_sensitivity_run_\\d+\\.txt$", full.names = TRUE)
if (!length(files)) {
  files <- list.files(outdir, pattern = "\\.txt$", full.names = TRUE)
}
if (!length(files)) stop("No TXT files found in: ", normalizePath(outdir))

## -------------------- Read & validate --------------------
required_cols <- c("run","S","n_s","alpha","epsilon","h_cv","M","rho","err","abs_DD")

dfs <- list()
bad <- 0L
for (f in files) {
  df <- try(read.table(f, header = TRUE, sep = "\t", check.names = FALSE), silent = TRUE)
  if (inherits(df, "try-error") || !nrow(df)) { bad <- bad + 1L; next }
  if (!all(required_cols %in% names(df))) { bad <- bad + 1L; next }
  
  # Coerce to proper types
  df$run     <- as.integer(df$run)
  df$S       <- as.integer(df$S)
  df$n_s     <- as.integer(df$n_s)
  df$alpha   <- as.numeric(df$alpha)
  df$epsilon <- as.numeric(df$epsilon)
  df$h_cv    <- as.numeric(df$h_cv)
  df$M       <- as.integer(df$M)
  df$rho     <- as.numeric(df$rho)
  df$err     <- as.numeric(df$err)
  df$abs_DD  <- as.numeric(df$abs_DD)
  
  dfs[[length(dfs)+1L]] <- df
}
if (!length(dfs)) stop("No valid result files with required columns in: ", normalizePath(outdir))
if (bad) cat("Skipped", bad, "files that lacked required columns or failed to parse.\n")

all_df <- do.call(rbind, dfs)

## -------------------- Infer grids --------------------
mc           <- max(all_df$run, na.rm = TRUE)
S_list       <- sort(unique(all_df$S))
n_s_list     <- sort(unique(all_df$n_s))
alpha_list   <- sort(unique(all_df$alpha))
epsilon_list <- sort(unique(all_df$epsilon))
rho_list     <- sort(unique(all_df$rho))

cat("Detected:\n  runs (mc): ", mc,
    "\n  S:          ", paste(S_list, collapse=", "),
    "\n  n_s:        ", paste(n_s_list, collapse=", "),
    "\n  alpha:      ", paste(alpha_list, collapse=", "),
    "\n  epsilon:    ", paste(epsilon_list, collapse=", "),
    "\n  rho:        ", paste(rho_list, collapse=", "), "\n", sep="")

## -------------------- Allocate arrays --------------------
# Dimensions: run × S × n_s × alpha × epsilon × rho
dn <- list(
  paste0("run",    seq_len(mc)),
  paste0("S=",     S_list),
  paste0("n_s=",   n_s_list),
  paste0("alpha=", alpha_list),
  paste0("eps=",   epsilon_list),
  paste0("rho=",   rho_list)
)
err_arr <- array(NA_real_, dim = sapply(dn, length), dimnames = dn)
DD_arr  <- array(NA_real_, dim = sapply(dn, length), dimnames = dn)

hsel_arr <- array(NA_real_,
                  dim = c(mc, length(S_list), length(n_s_list)),
                  dimnames = dn[1:3])

dnM <- dn[1:5]
M_arr <- array(NA_integer_,
               dim = sapply(dnM, length),
               dimnames = dnM)

## -------------------- Fill arrays --------------------
idx <- function(val, vec) match(val, vec)

for (i in seq_len(nrow(all_df))) {
  r  <- all_df$run[i]
  if (is.na(r) || r < 1 || r > mc) next
  
  si  <- idx(all_df$S[i],       S_list)
  ni  <- idx(all_df$n_s[i],     n_s_list)
  ai  <- idx(all_df$alpha[i],   alpha_list)
  ei  <- idx(all_df$epsilon[i], epsilon_list)
  ri  <- idx(all_df$rho[i],     rho_list)
  
  if (any(is.na(c(si,ni,ai,ei,ri)))) next
  
  err_arr[r, si, ni, ai, ei, ri] <- all_df$err[i]
  DD_arr[ r, si, ni, ai, ei, ri] <- all_df$abs_DD[i]

  hsel_arr[r, si, ni]            <- all_df$h_cv[i]
  M_arr[   r, si, ni, ai, ei]    <- all_df$M[i]
}

## -------------------- Summaries over runs --------------------
# Mean/SD for error and disparity
mean_err <- apply(err_arr, c(2,3,4,5,6), mean, na.rm = TRUE)
sd_err   <- apply(err_arr, c(2,3,4,5,6),   sd, na.rm = TRUE)
mean_DD  <- apply(DD_arr,  c(2,3,4,5,6), mean, na.rm = TRUE)
sd_DD    <- apply(DD_arr,  c(2,3,4,5,6),   sd, na.rm = TRUE)

# For h and M, summarize across runs
mean_hcv <- apply(hsel_arr, c(2,3), mean, na.rm = TRUE)  # S × n_s
mode_or_mean_M <- apply(M_arr, c(2,3,4,5), function(v) {
  v <- v[!is.na(v)]
  if (!length(v)) return(NA_integer_)
  if (length(unique(v)) == 1L) unique(v) else as.integer(round(mean(v)))
})

cat("Summary shapes:\n",
    "  mean_err: ", paste(dim(mean_err), collapse=" x"), "\n",
    "  mean_DD:  ", paste(dim(mean_DD), collapse=" x"), "\n",
    "  mean_hcv: ", paste(dim(mean_hcv), collapse=" x"), "\n",
    "  M summary:", paste(dim(mode_or_mean_M), collapse=" x"), "\n", sep="")

## -------------------- Save artifacts --------------------
timestamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
combined_rds  <- file.path(outdir, paste0("rho_combined_rows_", timestamp, ".rds"))

saveRDS(all_df, combined_rds)

