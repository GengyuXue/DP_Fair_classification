clip01 <- function(x, eps = 1e-6) pmin(pmax(x, eps), 1 - eps)
safediv <- function(num, den) if (den > 0) num / den else 0

############# Kernel functions #################

gaussian_kernel <- function(a, b, h) {
  d <- length(b)
  abs_diff <- t(abs(t(a) - b)) / h
  l2_abs_diff <- sqrt(rowSums(abs_diff^2))
  h^(-d) * exp(-(l2_abs_diff^2) / 2)
}

K_cov_gaussian <- function(X, h) {
  scale_X <- X / h
  Gram <- tcrossprod(scale_X)
  gdiag <- diag(Gram)
  Dist <- outer(gdiag, gdiag, "+") - 2 * Gram
  exp(-0.5 * Dist)
}

K_cov_gaussian_block <- function(X, Y, h) {
  Xs <- X / h; Ys <- Y / h
  gx <- rowSums(Xs * Xs); gy <- rowSums(Ys * Ys)
  Kxy <- outer(gx, gy, "+") - 2 * tcrossprod(Xs, Ys)
  exp(-0.5 * Kxy)
}

# eigen / Nyström factorization for GP draws
kernel_eig_factor_trunc <- function(
    X, h, tol = 1e-8, r = 35,
    nystrom_mult = 2, m_min = 50, small_n_threshold = 1800
) {
  n <- nrow(X)
  if (n == 0L) return(list(U = matrix(0, 0, 0), s = numeric(0)))
  if (n <= small_n_threshold) {
    K <- K_cov_gaussian(X, h)
    K <- (K + t(K)) / 2
    diag(K) <- diag(K) + tol
    ev <- eigen(K, symmetric = TRUE)
    lam <- pmax(ev$values, 0)
    r_eff <- min(as.integer(r), sum(lam > tol))
    if (r_eff < 1) return(list(U = matrix(0, n, 0), s = numeric(0)))
    list(U = ev$vectors[, seq_len(r_eff), drop = FALSE],
         s = sqrt(lam[seq_len(r_eff)]))
  } else {
    m <- max(as.integer(nystrom_mult * r), m_min); m <- min(m, n)
    set.seed(1)
    idx <- sample.int(n, m, replace = FALSE)
    X_land <- X[idx, , drop = FALSE]
    W <- K_cov_gaussian(X_land, h); W <- (W + t(W)) / 2; diag(W) <- diag(W) + tol
    C <- K_cov_gaussian_block(X, X_land, h)
    evW <- eigen(W, symmetric = TRUE)
    lamW <- pmax(evW$values, 0)
    keep <- which(lamW > tol)
    if (!length(keep)) return(list(U=matrix(0,n,0), s=numeric(0)))
    keep <- keep[order(lamW[keep], decreasing = TRUE)]
    if (length(keep) > r) keep <- keep[seq_len(r)]
    V <- evW$vectors[, keep, drop = FALSE]; Lam <- lamW[keep]
    inv_sqrt_Lam <- 1 / sqrt(Lam)
    U_tilde <- C %*% (V * rep(inv_sqrt_Lam, each = ncol(V)))
    qrU <- qr(U_tilde); Q <- qr.Q(qrU); R <- qr.R(qrU)
    S_small <- R %*% t(R)
    evS <- eigen((S_small + t(S_small))/2, symmetric = TRUE)
    lamS <- pmax(evS$values, 0); Vsmall <- evS$vectors
    U <- Q %*% Vsmall; s <- sqrt(lamS)
    r_fin <- min(ncol(U), length(s), r)
    if (r_fin < 1) return(list(U = matrix(0, n, 0), s = numeric(0)))
    list(U = U[, seq_len(r_fin), drop = FALSE], s = s[seq_len(r_fin)])
  }
}

GP_draw_from_factor <- function(fac) {
  r <- length(fac$s)
  if (r == 0) return(rep(0, nrow(fac$U)))
  z <- rnorm(r)
  drop(fac$U %*% (fac$s * z))
}

# KDE evaluated at Xeval against Xtrain
kde_vec <- function(Xtrain, Xeval, h, kernel_fun) {
  if (nrow(Xtrain) == 0) return(rep(0, nrow(Xeval)))
  apply(Xeval, 1, function(x) sum(kernel_fun(Xtrain, x, h)))
}

################# Data generation ##################

eta_arctan <- function(x, a, w = c(1.5, 0), b = -0.05, gamma = 0.05, k = 4, eps = 1e-3) {
  s <- b + sum(w * (x - 0.5)) + gamma * (2 * a - 1)
  p <- 0.5 + atan(k * s) / pi
  max(min(p, 1 - eps), eps)
}

generate_data <- function(n, d = 2, p_a,
                          alpha1 = 4, beta1 = 3,
                          alpha0 = 3, beta0 = 4,
                          eta_fun = eta_arctan,
                          w = c(1.5, 0.1), b = 0, gamma = 0.05, k = 4, eps = 1e-3) {
  A <- rbinom(n, 1, p_a)
  X <- matrix(0, n, d)
  X[A == 1, 1] <- rbeta(sum(A == 1), alpha1, beta1)
  X[A == 0, 1] <- rbeta(sum(A == 0), alpha0, beta0)
  if (d >= 2) X[, 2:d] <- matrix(runif(n * (d - 1), 0, 1), n, d - 1)
  Y <- vapply(seq_len(n),
              function(i) rbinom(1, 1, eta_fun(X[i, ], A[i], w, b, gamma, k)),
              numeric(1))
  idx_11 <- which(A == 1 & Y == 1)
  idx_10 <- which(A == 1 & Y == 0)
  idx_01 <- which(A == 0 & Y == 1)
  idx_00 <- which(A == 0 & Y == 0)
  list(
    Y = Y, X = X, A = A,
    X_11 = X[idx_11, , drop = FALSE],
    X_10 = X[idx_10, , drop = FALSE],
    X_01 = X[idx_01, , drop = FALSE],
    X_00 = X[idx_00, , drop = FALSE]
  )
}

# Split by cell into train/cal
sample_split <- function(data, split_ratio = 0.5) {
  X_11 <- data$X_11; X_10 <- data$X_10; X_01 <- data$X_01; X_00 <- data$X_00
  n_00 <- nrow(X_00); n_01 <- nrow(X_01); n_10 <- nrow(X_10); n_11 <- nrow(X_11)
  idx <- sample(1:n_00, ceiling(n_00 * split_ratio), replace = FALSE)
  X_00_train <- X_00[idx, , drop = FALSE]; X_00_cal <- X_00[-idx, , drop = FALSE]
  idx <- sample(1:n_01, ceiling(n_01 * split_ratio), replace = FALSE)
  X_01_train <- X_01[idx, , drop = FALSE]; X_01_cal <- X_01[-idx, , drop = FALSE]
  idx <- sample(1:n_10, ceiling(n_10 * split_ratio), replace = FALSE)
  X_10_train <- X_10[idx, , drop = FALSE]; X_10_cal <- X_10[-idx, , drop = FALSE]
  idx <- sample(1:n_11, ceiling(n_11 * split_ratio), replace = FALSE)
  X_11_train <- X_11[idx, , drop = FALSE]; X_11_cal <- X_11[-idx, , drop = FALSE]
  list(
    X_00_train = X_00_train, X_00_cal = X_00_cal,
    X_01_train = X_01_train, X_01_cal = X_01_cal,
    X_10_train = X_10_train, X_10_cal = X_10_cal,
    X_11_train = X_11_train, X_11_cal = X_11_cal
  )
}

############## π_tilde (empirical or private) ################

compute_pi_tilde <- function(input_data, split = c("train","cal"),
                             epsilon, delta,
                             pi_mode = c("private","empirical")) {
  split <- match.arg(split); pi_mode <- match.arg(pi_mode)
  if (split == "train") {
    n1 <- nrow(input_data$X_10_train) + nrow(input_data$X_11_train)
    n0 <- nrow(input_data$X_00_train) + nrow(input_data$X_01_train)
  } else {
    n1 <- nrow(input_data$X_10_cal) + nrow(input_data$X_11_cal)
    n0 <- nrow(input_data$X_00_cal) + nrow(input_data$X_01_cal)
  }
  n <- n1 + n0
  pi1_emp <- safediv(n1, n); pi0_emp <- safediv(n0, n)
  if (pi_mode == "empirical" || n == 0) {
    return(list(pi_tilde_1 = clip01(pi1_emp), pi_tilde_0 = clip01(pi0_emp)))
  }
  # PRIVATE 
  sd_pi <- if (n > 0) 4 * sqrt(2 * log(5 / delta)) / (n * epsilon) else 0
  noise <- rnorm(2, 0, sd_pi)
  list(pi_tilde_1 = clip01(pi1_emp + noise[1]),
       pi_tilde_0 = clip01(pi0_emp + noise[2]))
}

############## Core components ################

# Evaluate noisy posterior for group a on Xeval
eval_eta_for_a <- function(Xeval, input_data, a, h, epsilon, delta,
                           d = ncol(Xeval),
                           kernel_fun = gaussian_kernel, C_K = 1, tol = 1e-8,
                           p_floor = 1e-12, use_cal = FALSE, fac = NULL,
                           pi_mode = c("private","empirical")) {
  pi_mode <- match.arg(pi_mode)
  
  if (!isTRUE(use_cal)) {
    if (a == 0) { Xa1 <- input_data$X_01_train; Xa0 <- input_data$X_00_train }
    else        { Xa1 <- input_data$X_11_train; Xa0 <- input_data$X_10_train }
  } else {
    if (a == 0) { Xa1 <- input_data$X_01_cal;   Xa0 <- input_data$X_00_cal   }
    else        { Xa1 <- input_data$X_11_cal;   Xa0 <- input_data$X_10_cal   }
  }
  n_a <- nrow(Xa1) + nrow(Xa0); m <- nrow(Xeval)
  if (n_a == 0 || m == 0) return(list(p_hat = rep(NA_real_, m), eta_hat = rep(NA_real_, m)))
  
  S1 <- kde_vec(Xa1, Xeval, h, kernel_fun)
  S0 <- kde_vec(Xa0, Xeval, h, kernel_fun)
  
  use_fac <- (!is.null(fac) && !is.null(fac$U) && nrow(fac$U) == m)
  if (use_fac) { W1 <- GP_draw_from_factor(fac); W2 <- GP_draw_from_factor(fac) }
  else {
    fac_local <- kernel_eig_factor_trunc(Xeval, h, tol = tol, r = min(35, m))
    W1 <- GP_draw_from_factor(fac_local); W2 <- GP_draw_from_factor(fac_local)
  }
  
  # scl <- safediv(8 * sqrt(2 * C_K * log(8 / delta)), n_a * epsilon * h^d)
  scl <- if (pi_mode == "empirical") safediv(4 * sqrt(2 * C_K * log(4 / delta)),
                                             n_a * epsilon * h^d) 
          else safediv(8 * sqrt(2 * C_K * log(8 / delta)), n_a * epsilon * h^d)
  
  p_hat <- (S1 + S0) / ifelse(n_a > 0, n_a, 1) + scl * W1
  bad <- !is.finite(p_hat) | (p_hat <= 0); if (any(bad)) p_hat[bad] <- p_floor
  numer <- S1 / ifelse(n_a > 0, n_a, 1) + scl * W2; numer[!is.finite(numer)] <- 0
  eta_hat <- pmin(pmax(numer / p_hat, 0), 1)
  list(p_hat = p_hat, eta_hat = eta_hat)
}

# ==== τ calibration ====
CDP_calibrate_tau <- function(
    input_data, h, epsilon, delta, alpha,
    kernel_fun = gaussian_kernel, C_K = 1, tol = 1e-8,
    pi_tilde_1 = NULL, pi_tilde_0 = NULL,
    switch = FALSE,
    pi_mode = c("private","empirical"),
    tau_method = c("directional","grid"),
    grid_size = 250,       
    use_range_bounds = TRUE,  
    tol_root = 1e-6
) {
  pi_mode <- match.arg(pi_mode); tau_method <- match.arg(tau_method)
  

  if (!isTRUE(switch)) {
    Xcal_0 <- rbind(input_data$X_00_cal, input_data$X_01_cal)
    Xcal_1 <- rbind(input_data$X_10_cal, input_data$X_11_cal)
    split_used <- "cal"
  } else {
    Xcal_0 <- rbind(input_data$X_00_train, input_data$X_01_train)
    Xcal_1 <- rbind(input_data$X_10_train, input_data$X_11_train)
    split_used <- "train"
  }
  n0 <- nrow(Xcal_0); n1 <- nrow(Xcal_1)
  d  <- if (n0 + n1 > 0) ncol(rbind(Xcal_0, Xcal_1)) else 2
  
  # π-tilde
  if (is.null(pi_tilde_1) || is.null(pi_tilde_0)) {
    pis <- compute_pi_tilde(input_data, split = split_used,
                            epsilon = epsilon, delta = delta, pi_mode = pi_mode)
    pi_tilde_1 <- pis$pi_tilde_1; pi_tilde_0 <- pis$pi_tilde_0
  }
  pi1 <- pi_tilde_1; pi0 <- pi_tilde_0
  
 
  eta0 <- if (n0) eval_eta_for_a(Xcal_0, input_data, 0, h, epsilon, delta, d,
                                 kernel_fun, C_K, tol, use_cal = switch, pi_mode = pi_mode)$eta_hat else numeric(0)
  eta1 <- if (n1) eval_eta_for_a(Xcal_1, input_data, 1, h, epsilon, delta, d,
                                 kernel_fun, C_K, tol, use_cal = switch, pi_mode = pi_mode)$eta_hat else numeric(0)
  
  
  n_cal <- n0 + n1
  sd_DD <- if (n_cal > 0) 2 * sqrt(2 * log(1.25 / delta)) / (n_cal * epsilon) else 0
  w_tilde <- rnorm(1, 0, sd_DD)
  
  eta1s <- sort(eta1); eta0s <- sort(eta0)
  tail_geq <- function(s, thr) {
    if (!length(s)) return(0)
    k <- findInterval(thr, s, left.open = FALSE)
    (length(s) - k) / length(s)
  }
  DD <- function(tau) {
    t1 <- clip01(0.5 + tau/(2*pi1))
    t0 <- clip01(0.5 - tau/(2*pi0))
    tail_geq(eta1s, t1) - tail_geq(eta0s, t0) + w_tilde
  }

  
  L_base <- -min(pi1, pi0); U_base <-  min(pi1, pi0)

  if (use_range_bounds && n1 > 0) { r1 <- range(eta1s); L1 <-  2*pi1*(r1[1]-0.5); U1 <-  2*pi1*(r1[2]-0.5) } else { L1 <- -Inf; U1 <- +Inf }
  if (use_range_bounds && n0 > 0) { r0 <- range(eta0s); L0 <- -2*pi0*(r0[2]-0.5); U0 <- -2*pi0*(r0[1]-0.5) } else { L0 <- -Inf; U0 <- +Inf }
  
  tau_lo <- max(L_base, L1, L0)
  tau_hi <- min(U_base, U1, U0)
  if (!(tau_lo < tau_hi)) {
    warning("No feasible tau interval; returning tau=0.")
    return(list(tauhat = 0, DD = NA_real_, feasible = FALSE,
                feasible_set = c(NA_real_, NA_real_),
                pi_tilde_1 = pi1, pi_tilde_0 = pi0))
  }
  
  # First evaluate at τ = 0
  D0 <- DD(0)
  if (abs(D0) <= alpha) {
    return(list(tauhat = 0, DD = D0, feasible = TRUE,
                feasible_set = c(0,0),
                pi_tilde_1 = pi1, pi_tilde_0 = pi0))
  }
  
  # ---------- directional (bisection on needed side) ----------
  if (tau_method == "directional") {
    solve_root <- function(target, a, b) {
      fa <- DD(a) - target; fb <- DD(b) - target
      if (fa * fb > 0) return(NA_real_)
      for (it in 1:80) {
        m <- 0.5*(a+b); fm <- DD(m) - target
        if (abs(fm) <= 1e-6 || (b-a) <= tol_root) return(m)
        if (fa * fm <= 0) { b <- m; fb <- fm } else { a <- m; fa <- fm }
      }
      0.5*(a+b)
    }
    if (D0 > alpha) {              
      tau_star <- solve_root(alpha, 0, tau_hi)
      feasible <- !is.na(tau_star)
      tau_hat  <- if (feasible) tau_star else tau_hi
      return(list(tauhat = tau_hat, DD = DD(tau_hat), feasible = feasible,
                  feasible_set = if (feasible) c(tau_star, tau_hi) else c(NA_real_, NA_real_),
                  pi_tilde_1 = pi1, pi_tilde_0 = pi0))
    } else {                             
      tau_star <- solve_root(-alpha, tau_lo, 0)
      feasible <- !is.na(tau_star)
      tau_hat  <- if (feasible) tau_star else tau_lo
      return(list(tauhat = tau_hat, DD = DD(tau_hat), feasible = feasible,
                  feasible_set = if (feasible) c(tau_lo, tau_star) else c(NA_real_, NA_real_),
                  pi_tilde_1 = pi1, pi_tilde_0 = pi0))
    }
  }
  
  # ---------- one-sided grid ----------
  grid_size <- as.integer(grid_size)
  if (grid_size < 5L) grid_size <- 5L
  if (grid_size %% 2 == 0L) grid_size <- grid_size + 1L
  
  if (D0 > alpha) {
    tau_grid <- seq(0, tau_hi, length.out = grid_size)
  } else {
    tau_grid <- seq(tau_lo, 0, length.out = grid_size)
  }
  
  tau_grid <- sort(unique(tau_grid))
  
  DD_vals <- vapply(tau_grid, DD, numeric(1))
  feas    <- which(abs(DD_vals) <= alpha)
  
  if (length(feas)) {
    idx <- feas[which.min(abs(tau_grid[feas]))]
    same <- which(abs(tau_grid[feas]) == abs(tau_grid[idx]))
    if (length(same) > 1) {
      cand <- feas[same]
      idx <- cand[which.min(abs(DD_vals[cand]))]
    }
    tauhat <- tau_grid[idx]; feasible <- TRUE
    L <- min(tau_grid[feas]); U <- max(tau_grid[feas])
  } else {
    viol <- abs(abs(DD_vals) - alpha)
    idx  <- which.min(viol)
    tie  <- which(viol == viol[idx])
    if (length(tie) > 1) idx <- tie[which.min(abs(tau_grid[tie]))]
    tauhat <- tau_grid[idx]; feasible <- FALSE
    L <- NA_real_; U <- NA_real_
    warning(sprintf("No τ with |DD| <= α on %s side grid; picked τ minimizing violation.",
                    if (D0 > alpha) "right" else "left"))
  }
  
  list(tauhat = tauhat, DD = DD(tauhat), feasible = feasible,
       feasible_set = c(L, U),
       pi_tilde_1 = pi1, pi_tilde_0 = pi0)
}

DD_train <- function(test_data, y_hat) {
  A <- as.integer(test_data$A)
  p1 <- if (any(A == 1)) mean(y_hat[A == 1] == 1) else 0
  p0 <- if (any(A == 0)) mean(y_hat[A == 0] == 1) else 0
  gap <- p1 - p0
  list(disparity = gap, abs_disparity = abs(gap))
}

# Main classifier 
CDP_fair_classifier <- function(test_data, input_data, alpha, h, epsilon, delta,
                                kernel_fun = gaussian_kernel, tol, C_K = 1,
                                cross_fit = TRUE,
                                pi_mode = c("private","empirical"),
                                tau_method = c("directional","grid")) {
  pi_mode <- match.arg(pi_mode); tau_method <- match.arg(tau_method)
  d <- ncol(input_data$X)
  input_data <- sample_split(input_data, split_ratio = 0.5)
  test_X <- test_data$X
  

  pis_tr <- compute_pi_tilde(input_data, split = "train",
                             epsilon = epsilon, delta = delta, pi_mode = pi_mode)
  pi_tilde_1 <- pis_tr$pi_tilde_1; pi_tilde_0 <- pis_tr$pi_tilde_0
  

  n_1_train <- nrow(input_data$X_10_train) + nrow(input_data$X_11_train)
  n_0_train <- nrow(input_data$X_00_train) + nrow(input_data$X_01_train)
  fac <- kernel_eig_factor_trunc(test_X, h, tol = 1e-8, r = 35)
  W1_0 <- GP_draw_from_factor(fac); W2_0 <- GP_draw_from_factor(fac)
  W1_1 <- GP_draw_from_factor(fac); W2_1 <- GP_draw_from_factor(fac)
  
  scl_const <- if (pi_mode == "empirical") 4 else 8
  scl_log   <- if (pi_mode == "empirical") log(4 / delta) else log(8 / delta)
  
  scl0 <- safediv(scl_const * sqrt(2 * C_K * scl_log), n_0_train * epsilon * h^d)
  scl1 <- safediv(scl_const * sqrt(2 * C_K * scl_log), n_1_train * epsilon * h^d)

  compute_for_a <- function(a) {
    if (a == 0) {
      Xa1 <- input_data$X_01_train; Xa0 <- input_data$X_00_train
      na <- n_0_train ; W1 <- W1_0; W2 <- W2_0; scl <- scl0
    } else {
      Xa1 <- input_data$X_11_train; Xa0 <- input_data$X_10_train
      na <- n_1_train; W1 <- W1_1; W2 <- W2_1; scl <- scl1
    }
    S1 <- kde_vec(Xa1, test_X, h, kernel_fun)
    S0 <- kde_vec(Xa0, test_X, h, kernel_fun)
    p_hat <- pmax((S1 + S0) / ifelse(na > 0, na, 1) + scl * W1, .Machine$double.eps)
    numer <- S1 / ifelse(na > 0, na, 1) + scl * W2
    eta_hat <- pmin(pmax(numer / p_hat, 0), 1)
    list(p_hat = p_hat, eta_hat = eta_hat)
  }
  eta_output_0 <- compute_for_a(0); eta_output_1 <- compute_for_a(1)
  

  tau_data <- CDP_calibrate_tau(
    input_data, h, epsilon, delta, alpha,
    kernel_fun = gaussian_kernel, C_K = 1, tol = 1e-8,
    pi_tilde_1 = pi_tilde_1, pi_tilde_0 = pi_tilde_0,
    switch = FALSE, pi_mode = pi_mode,
    tau_method = tau_method
  )
  tau_hat <- tau_data$tauhat
  
  thr0 <- clip01(0.5 - tau_hat / (2 * pi_tilde_0))
  thr1 <- clip01(0.5 + tau_hat / (2 * pi_tilde_1))
  eta_test <- ifelse(test_data$A == 0, eta_output_0$eta_hat, eta_output_1$eta_hat)
  thr_vec  <- ifelse(test_data$A == 0, thr0, thr1)
  y_hat    <- as.integer(eta_test >= thr_vec)
  DD_hat   <- DD_train(test_data, y_hat)
  
  if (!isTRUE(cross_fit)) {
    return(list(
      pi_tilde_1 = pi_tilde_1, pi_tilde_0 = pi_tilde_0,
      eta_0 = eta_output_0$eta_hat, eta_1 = eta_output_1$eta_hat,
      tau_hat = tau_hat, DD_hat = DD_hat$abs_disparity, y_hat = y_hat,
      tau_feasible = tau_data$feasible, tau_feasible_set = tau_data$feasible_set
    ))
  }
  
  # cross-fit pass 
  pis_cal <- compute_pi_tilde(input_data, split = "cal",
                              epsilon = epsilon, delta = delta, pi_mode = pi_mode)
  pi_tilde_1_cal <- pis_cal$pi_tilde_1; pi_tilde_0_cal <- pis_cal$pi_tilde_0
  n_1_cal <- nrow(input_data$X_10_cal) + nrow(input_data$X_11_cal)
  n_0_cal <- nrow(input_data$X_00_cal) + nrow(input_data$X_01_cal)
  W1_0_cal <- GP_draw_from_factor(fac); W2_0_cal <- GP_draw_from_factor(fac)
  W1_1_cal <- GP_draw_from_factor(fac); W2_1_cal <- GP_draw_from_factor(fac)
  
  
  scl_const <- if (pi_mode == "empirical") 4 else 8
  scl_log   <- if (pi_mode == "empirical") log(4 / delta) else log(8 / delta)
  
  scl0_cal <- safediv(scl_const * sqrt(2 * C_K * scl_log), max(n_0_cal, 1) * epsilon * h^d)
  scl1_cal <- safediv(scl_const * sqrt(2 * C_K * scl_log), max(n_1_cal, 1) * epsilon * h^d)
  
  compute_for_a_cal <- function(a) {
    if (a == 0) {
      Xa1 <- input_data$X_01_cal; Xa0 <- input_data$X_00_cal
      na <- n_0_cal ; W1 <- W1_0_cal; W2 <- W2_0_cal; scl <- scl0_cal
    } else {
      Xa1 <- input_data$X_11_cal; Xa0 <- input_data$X_10_cal
      na <- n_1_cal; W1 <- W1_1_cal; W2 <- W2_1_cal; scl <- scl1_cal
    }
    S1 <- kde_vec(Xa1, test_X, h, kernel_fun)
    S0 <- kde_vec(Xa0, test_X, h, kernel_fun)
    p_hat <- pmax((S1 + S0) / ifelse(na > 0, na, 1) + scl * W1, .Machine$double.eps)
    numer <- S1 / ifelse(na > 0, na, 1) + scl * W2
    eta_hat <- pmin(pmax(numer / p_hat, 0), 1)
    list(p_hat = p_hat, eta_hat = eta_hat)
  }
  eta_output_0_cal <- compute_for_a_cal(0); eta_output_1_cal <- compute_for_a_cal(1)
  
  tau_data_cal <- CDP_calibrate_tau(
    input_data, h, epsilon, delta, alpha,
    kernel_fun = gaussian_kernel, C_K = 1, tol = 1e-8,
    pi_tilde_1 = pi_tilde_1_cal, pi_tilde_0 = pi_tilde_0_cal,
    switch = TRUE, pi_mode = pi_mode,
    tau_method = tau_method
  )
  tau_hat_cal <- tau_data_cal$tauhat
  thr0_cal <- clip01(0.5 - tau_hat_cal / (2 * pi_tilde_0_cal))
  thr1_cal <- clip01(0.5 + tau_hat_cal / (2 * pi_tilde_1_cal))
  eta_test_cal <- ifelse(test_data$A == 0, eta_output_0_cal$eta_hat, eta_output_1_cal$eta_hat)
  thr_vec_cal  <- ifelse(test_data$A == 0, thr0_cal, thr1_cal)
  y_hat_cal    <- as.integer(eta_test_cal >= thr_vec_cal)
  
  # simple blend for out-of-fold style
  y_hat_out <- rbinom(length(y_hat), 1, 0.5 * y_hat + 0.5 * y_hat_cal)
  DD_hat_out <- DD_train(test_data, y_hat_out)
  
  list(
    pi_tilde_1_train = pi_tilde_1, pi_tilde_0_train = pi_tilde_0,
    pi_tilde_1_cal = pi_tilde_1_cal, pi_tilde_0_cal = pi_tilde_0_cal,
    eta_0_train = eta_output_0$eta_hat, eta_1_train = eta_output_1$eta_hat,
    eta_0_cal = eta_output_0_cal$eta_hat, eta_1_cal = eta_output_1_cal$eta_hat,
    tau_hat_train = tau_hat, tau_hat_cal = tau_hat_cal,
    DD_hat_train = DD_hat, DD_hat_cal = tau_data_cal$DD,
    y_hat_train = y_hat, y_hat_cal = y_hat_cal,
    y_hat_out = y_hat_out, DD_hat_out = DD_hat_out$abs_disparity,
    tau_feasible_train = tau_data$feasible, tau_feasible_set_train = tau_data$feasible_set,
    tau_feasible_cal   = tau_data_cal$feasible, tau_feasible_set_cal   = tau_data_cal$feasible_set
  )
}

############## CV bandwidth ################

eta_from_kde <- function(Xtrain, Ytrain, Xeval, h,
                         kernel_fun = gaussian_kernel, tol = 1e-12) {
  X1 <- Xtrain[Ytrain == 1, , drop = FALSE]
  X0 <- Xtrain[Ytrain == 0, , drop = FALSE]
  S1 <- if (nrow(X1)) kde_vec(X1, Xeval, h, kernel_fun) else rep(0, nrow(Xeval))
  S0 <- if (nrow(X0)) kde_vec(X0, Xeval, h, kernel_fun) else rep(0, nrow(Xeval))
  den <- S1 + S0
  prior <- if (length(Ytrain)) mean(Ytrain) else 0.5
  eta <- ifelse(!is.finite(den) | den <= tol, prior, S1 / den)
  pmin(pmax(eta, 0), 1)
}

cv_select_h <- function(data, h_grid, n_folds = 3,
                        kernel_fun = gaussian_kernel, tol = 1e-12) {
  X <- data$X; Y <- as.integer(data$Y); n <- nrow(X)
  fold_ids <- sample(rep(1:n_folds, length.out = n))
  mean_err <- numeric(length(h_grid))
  for (j in seq_along(h_grid)) {
    h <- h_grid[j]; errs <- numeric(n_folds)
    for (k in seq_len(n_folds)) {
      idx_val <- which(fold_ids == k); idx_tr <- which(fold_ids != k)
      eta_vl  <- eta_from_kde(X[idx_tr, , drop = FALSE], Y[idx_tr],
                              X[idx_val, , drop = FALSE],
                              h, kernel_fun, tol)
      yhat_vl <- as.integer(eta_vl >= 0.5)
      errs[k] <- mean(yhat_vl != Y[idx_val])
    }
    mean_err[j] <- mean(errs)
  }
  j_star <- which.min(mean_err)
  list(h_star = h_grid[j_star],
       cv_error = mean_err[j_star],
       errors_by_h = setNames(mean_err, paste0("h=", h_grid)))
}
