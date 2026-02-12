options(stringsAsFactors = FALSE)

## ------------------ Config ------------------
outdir <- getwd()
if (!dir.exists(outdir)) stop("Working directory not found: ", normalizePath(outdir))
cat("Combining from: ", normalizePath(outdir), "\n", sep = "")

files <- list.files(
  outdir,
  pattern = "adult_unfair_run_\\d+\\.(txt|tsv)$",
  full.names = TRUE,
  ignore.case = TRUE
)
if (!length(files)) {
  files <- list.files(outdir, pattern = "\\.(tsv|txt)$", full.names = TRUE, ignore.case = TRUE)
}
if (!length(files)) stop("No TSV/TXT files found in: ", normalizePath(outdir))

## ------------------ Read helpers ------------------
read_tab_or_ws <- function(path) {
  df <- try(utils::read.table(path, header = TRUE, sep = "\t",
                              check.names = FALSE, stringsAsFactors = FALSE,
                              na.strings = c("NA", "NaN", "", "null")),
            silent = TRUE)
  if (!inherits(df, "try-error") && is.data.frame(df) && ncol(df) > 1) return(df)
  
  utils::read.table(path, header = TRUE, sep = "",
                    check.names = FALSE, stringsAsFactors = FALSE,
                    na.strings = c("NA", "NaN", "", "null"))
}

required_cols <- c("run","method","epsilon","alpha","error","disparity","h_cv")

## ------------------ Read & validate ------------------
dfs <- list()
bad <- 0L

for (f in files) {
  df <- try(read_tab_or_ws(f), silent = TRUE)
  if (inherits(df, "try-error") || !is.data.frame(df) || !nrow(df)) { bad <- bad + 1L; next }
  
  # normalize column names
  names(df) <- trimws(names(df))
  
  if (!all(required_cols %in% names(df))) { bad <- bad + 1L; next }
  
  df$run       <- suppressWarnings(as.integer(df$run))
  df$method    <- as.character(df$method)
  df$epsilon   <- suppressWarnings(as.numeric(df$epsilon))   # FAIR_no_DP => NA ok
  df$alpha     <- suppressWarnings(as.numeric(df$alpha))
  df$error     <- suppressWarnings(as.numeric(df$error))
  df$disparity <- suppressWarnings(as.numeric(df$disparity))
  df$h_cv      <- suppressWarnings(as.numeric(df$h_cv))
  
  # Keep only expected columns
  df <- df[, required_cols]
  
  dfs[[length(dfs)+1L]] <- df
}

if (!length(dfs)) stop("No valid result files with required columns in: ", normalizePath(outdir))
if (bad) cat("Skipped", bad, "file(s) that lacked required columns or failed to parse.\n")

all_df <- do.call(rbind, dfs)

## ------------------ Inspect grids ------------------
mc           <- suppressWarnings(max(all_df$run, na.rm = TRUE))
if (!is.finite(mc) || is.na(mc)) mc <- 0L

method_list  <- sort(unique(all_df$method))
alpha_list   <- sort(unique(all_df$alpha))
epsilon_list <- sort(unique(all_df$epsilon))  # includes NA if present

to_lab <- function(x) ifelse(is.na(x), "NA", as.character(x))

cat("Detected:\n  runs (mc):  ", mc,
    "\n  methods:    ", paste(method_list, collapse = ", "),
    "\n  alpha:      ", paste(to_lab(alpha_list), collapse = ", "),
    "\n  epsilon:    ", paste(to_lab(epsilon_list), collapse = ", "),
    "\n", sep = "")

## ------------------ Per-method arrays ------------------
match_with_na <- function(x, levels) {
  idx <- match(x, levels)
  if (any(is.na(levels))) {
    na_level <- which(is.na(levels))[1L]
    idx[is.na(x)] <- na_level
  }
  idx
}

errors_arr <- list()
DD_arr     <- list()
mean_error <- list()
sd_error   <- list()
mean_DD    <- list()
sd_DD      <- list()

runs_present <- sort(unique(all_df$run))
stopifnot(length(runs_present) > 0)

for (m in method_list) {
  df_m <- all_df[all_df$method == m, , drop = FALSE]
  
  # ensure numeric types
  df_m$alpha   <- suppressWarnings(as.numeric(df_m$alpha))
  df_m$epsilon <- suppressWarnings(as.numeric(df_m$epsilon))
  
  # grids present for this method
  alpha_vals <- sort(unique(df_m$alpha))
  eps_vals   <- unique(df_m$epsilon)
  if (length(eps_vals) == 0L) eps_vals <- NA_real_
  eps_vals   <- eps_vals[order(ifelse(is.na(eps_vals), Inf, eps_vals))]
  
  # labels
  run_labs   <- paste0("run=",   runs_present)
  alpha_labs <- paste0("alpha=", to_lab(alpha_vals))
  eps_labs   <- paste0("eps=",   to_lab(eps_vals))
  
  # allocate
  dims <- c(length(runs_present), length(alpha_vals), length(eps_vals))
  arr_err <- array(NA_real_, dim = dims)
  arr_dd  <- array(NA_real_, dim = dims)
  
  # assign dimnames only for non-zero extents
  dn <- vector("list", 3L)
  if (dims[1] > 0) dn[[1]] <- run_labs
  if (dims[2] > 0) dn[[2]] <- alpha_labs
  if (dims[3] > 0) dn[[3]] <- eps_labs
  dimnames(arr_err) <- dn
  dimnames(arr_dd)  <- dn
  
  # fill
  r_idx <- match(df_m$run,     runs_present)
  a_idx <- match(df_m$alpha,   alpha_vals)
  e_idx <- match_with_na(df_m$epsilon, eps_vals)
  
  keep <- !(is.na(r_idx) | is.na(a_idx) | is.na(e_idx))
  if (any(keep)) {
    sel <- cbind(r_idx[keep], a_idx[keep], e_idx[keep])
    arr_err[sel] <- df_m$error[keep]
    arr_dd[ sel] <- df_m$disparity[keep]
  }
  
  # summaries over runs 
  mean_error[[m]] <- apply(arr_err, c(2,3), mean, na.rm = TRUE)
  sd_error[[m]]   <- apply(arr_err, c(2,3),   sd, na.rm = TRUE)
  mean_DD[[m]]    <- apply(arr_dd,  c(2,3), mean, na.rm = TRUE)
  sd_DD[[m]]      <- apply(arr_dd,  c(2,3),   sd, na.rm = TRUE)
  
  errors_arr[[m]] <- arr_err
  DD_arr[[m]]     <- arr_dd
}


## ------------------ Save artifacts ------------------
timestamp     <- format(Sys.time(), "%Y%m%d-%H%M%S")
combined_rds  <- file.path(outdir, paste0("adult_unfair_combined_rows_", timestamp, ".rds"))

saveRDS(all_df, combined_rds)


