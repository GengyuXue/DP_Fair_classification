# ---------------------------- Data Generation -------------------------------
eta_arctan <- function(x, a, w = c(1.5, 0.1), b = 0, gamma = 0.05, k = 4, eps = 1e-3) {
  s <- b + sum(w * (x - 0.5)) + gamma * (2*a - 1)
  p <- 0.5 + atan(k * s) / pi
  max(min(p, 1 - eps), eps)
}


generate_site_data <- function(n, d = 2, p_a = 0.5,
                               alpha1 = 4, beta1 = 3, alpha0 = 3, beta0 = 4,
                               eta_fun = eta_arctan, w = c(1.5,0.1), b = 0, gamma = 0.05, k = 4) {
  A <- rbinom(n, 1, p_a)
  X <- matrix(0, n, d)
  X[A == 1, 1] <- rbeta(sum(A == 1), alpha1, beta1)
  X[A == 0, 1] <- rbeta(sum(A == 0), alpha0, beta0)
  if (d >= 2) X[, 2:d] <- matrix(runif(n * (d - 1), 0, 1), n, d - 1)
  Y <- vapply(seq_len(n), function(i) rbinom(1, 1, eta_fun(X[i,], A[i], w,b,gamma,k)), 1)
  
  list(X = X, A = A, Y = Y)
}

split_train_cal <- function(dat, split_ratio = 0.5) {
  idx <- sample.int(length(dat$Y), size = ceiling(length(dat$Y) * split_ratio))
  list(train = list(X = dat$X[idx,,drop=FALSE], A = dat$A[idx], Y = dat$Y[idx]),
       cal   = list(X = dat$X[-idx,,drop=FALSE], A = dat$A[-idx], Y = dat$Y[-idx]))
}


# distribute data into S servers
generate_sites <- function(S = 3, n_per_site = 400, d = 2, p_a = 0.6,
                           split_ratio = 0.5, alpha1 = 4, beta1 = 3, alpha0 = 3, beta0 = 4,
                           eta_fun = eta_arctan,  w = c(1.5, 0.1), b = 0, gamma = 0.05,
                           k = 4) {
  
  repS <- function(x) if (length(x) == 1) rep(x, S) else x
  
  n_vec   <- repS(n_per_site)
  p_vec   <- repS(p_a)
  split_v <- repS(split_ratio)
  
  sites_all <- vector("list", S)
  sites_tr  <- vector("list", S)
  sites_cal <- vector("list", S)
  
  for (s in seq_len(S)) {
    dat <- generate_site_data(n = n_vec[s], d = d, p_a = p_vec[s], alpha1 = alpha1, beta1 = beta1,
                              alpha0 = alpha0, beta0 = beta0, eta_fun = eta_fun, w = w, b = b,
                              gamma = gamma, k = k)
    parts <- split_train_cal(dat, split_ratio = split_v[s])
    sites_all[[s]] <- dat
    sites_tr[[s]]  <- parts$train
    sites_cal[[s]] <- parts$cal
  }
  
  list(all = sites_all, train = sites_tr, cal = sites_cal)
}




# ---------------------------- Kernel and evaluation -------------------------------
gaussian_kernel <- function(a, b, h) {
  d <- length(b)
  abs_diff <- t(abs(t(a) - b)) / h
  l2_abs_diff <- sqrt(rowSums(abs_diff^2))
  kernel_value <- h^(-d) * exp(-(l2_abs_diff^2) / (2 ))
  return(kernel_value)
}


K_cov_gaussian <- function(X, h) {
  scale_X <- X / h
  Gram  <- tcrossprod(scale_X)
  Gram_diagonal <- diag(Gram)
  Dist <- outer(Gram_diagonal, Gram_diagonal, "+") - 2 * Gram 
  return(exp(-0.5 * Dist))
}



K_cov_gaussian_block <- function(X, Y, h) {
  X <- as.matrix(X); Y <- as.matrix(Y)
  Xs <- X / h; Ys <- Y / h
  gx <- rowSums(Xs * Xs); gy <- rowSums(Ys * Ys)
  Kxy <- outer(gx, gy, "+") - 2 * tcrossprod(Xs, Ys)
  exp(-0.5 * Kxy)
}


kernel_eig_factor_trunc <- function(
    X, h, tol = 1e-8, r = 35,
    nystrom_mult = 2, m_min = 50, small_n_threshold = 1800
) {
  X <- as.matrix(X); n <- nrow(X)
  if (n == 0L) return(list(U = matrix(0, 0, 0), s = numeric(0)))
  
  if (n <= small_n_threshold) {
    K <- K_cov_gaussian(X, h)
    K <- (K + t(K)) / 2; diag(K) <- diag(K) + tol
    ev  <- eigen(K, symmetric = TRUE)
    lam <- pmax(ev$values, 0)
    r_eff <- min(as.integer(r), sum(lam > tol))
    if (r_eff < 1) return(list(U = matrix(0, n, 0), s = numeric(0)))
    return(list(U = ev$vectors[, seq_len(r_eff), drop = FALSE],
                s = sqrt(lam[seq_len(r_eff)])))
  }
  
  m <- max(as.integer(nystrom_mult * r), m_min); m <- min(m, n)
  idx    <- sample.int(n, m, replace = FALSE)
  X_land <- X[idx, , drop = FALSE]
  
  W <- K_cov_gaussian(X_land, h); W <- (W + t(W)) / 2; diag(W) <- diag(W) + tol
  C <- K_cov_gaussian_block(X, X_land, h)
  
  evW  <- eigen(W, symmetric = TRUE)
  lamW <- pmax(evW$values, 0)
  keep <- which(lamW > tol)
  if (!length(keep)) return(list(U = matrix(0, n, 0), s = numeric(0)))
  keep <- keep[order(lamW[keep], decreasing = TRUE)]
  if (length(keep) > r) keep <- keep[seq_len(r)]
  
  V   <- evW$vectors[, keep, drop = FALSE]
  Lam <- lamW[keep]
  
  inv_sqrt_Lam <- 1 / sqrt(Lam)
  U_tilde <- C %*% (V * rep(inv_sqrt_Lam, each = ncol(V)))
  
  qrU <- qr(U_tilde); Q <- qr.Q(qrU); R <- qr.R(qrU)
  S_small <- R %*% t(R)
  evS  <- eigen((S_small + t(S_small))/2, symmetric = TRUE)
  lamS <- pmax(evS$values, 0)
  Vsm  <- evS$vectors
  
  U <- Q %*% Vsm
  s <- sqrt(lamS)
  
  r_fin <- min(ncol(U), length(s), r)
  if (r_fin < 1) return(list(U = matrix(0, n, 0), s = numeric(0)))
  list(U = U[, seq_len(r_fin), drop = FALSE], s = s[seq_len(r_fin)])
}


GP_draw_from_factor <- function(fac) {
  r <- length(fac$s)
  if (r == 0) return(rep(0, nrow(fac$U)))  
  z <- rnorm(r)
  drop(fac$U %*% (fac$s * z))
}


#' Evaluated KDE
kde_vec <- function(Xtrain, Xeval, h, kernel_fun) {
  if (nrow(Xtrain) == 0) return(rep(0, nrow(Xeval)))
  apply(Xeval, 1, function(x) sum(kernel_fun(Xtrain, x, h)))
}




# ---------------------------- Algorithm 2 S1 -------------------------------
site_private_components <- function(train, Xeval, h, eps_s, delta_s,
                                    C_K = 1, d = ncol(Xeval),
                                    kernel_fun = gaussian_kernel,
                                    pi_method = c("dp","empirical")) {
  pi_method <- match.arg(pi_method)
  m   <- nrow(Xeval)
  n   <- length(train$A)
  ns0 <- sum(train$A == 0)
  ns1 <- sum(train$A == 1)
  
  # ---- class priors ----
  if (pi_method == "dp") {
    sigma <- 4 * sqrt(2 * log(5 / delta_s)) / (n * eps_s)
    w <- rnorm(2, 0, sigma)
    pi1 <- min(max(ns1 / n + w[1], 1e-6), 1 - 1e-6)
    pi0 <- 1 - pi1
  } else { # empirical
    pi1 <- min(max(ns1 / n, 1e-6), 1 - 1e-6)
    pi0 <- 1 - pi1
  }
  
  # ---- GP factors for kernel perturbations ----
  fac <- kernel_eig_factor_trunc(Xeval, h, tol = 1e-8, r = min(35, m))
  Wp0 <- GP_draw_from_factor(fac)
  Wn0 <- GP_draw_from_factor(fac)
  Wp1 <- GP_draw_from_factor(fac)
  Wn1 <- GP_draw_from_factor(fac)
  
  # ---- DP scale for kernel terms: conditional on pi_method ----
  scale_mult <- if (pi_method == "empirical") 4 else 8
  
  ns0 <- sum(train$A == 0); ns1 <- sum(train$A == 1)
  scl0 <- scale_mult * sqrt(2 * C_K * log(scale_mult / delta_s)) / (ns0 * eps_s * h^d)
  scl1 <- scale_mult * sqrt(2 * C_K * log(scale_mult / delta_s)) / (ns1 * eps_s * h^d)
  
  comp <- function(a, Wp, Wn, scl) {
    Xa1 <- train$X[train$A == a & train$Y == 1, , drop = FALSE]
    Xa0 <- train$X[train$A == a & train$Y == 0, , drop = FALSE]
    n_a <- nrow(Xa1) + nrow(Xa0)
    
    S1 <- kde_vec(Xa1, Xeval, h, kernel_fun)
    S0 <- kde_vec(Xa0, Xeval, h, kernel_fun)
    
    p   <- (S1 + S0) / n_a + scl * Wp
    p[!is.finite(p) | p <= 0] <- 1e-6
    
    px1 <- S1 / n_a + scl * Wn
    px1[!is.finite(px1)] <- 0
    px1 <- pmax(px1, 0)
    
    list(p = p, px1 = px1)
  }
  
  a0 <- comp(0, Wp0, Wn0, scl0)
  a1 <- comp(1, Wp1, Wn1, scl1)
  
  list(
    p_X_given_A   = list('0' = a0$p,   '1' = a1$p),
    p_XY1_given_A = list('0' = a0$px1, '1' = a1$px1),
    pi_tilde      = c('0' = pi0, '1' = pi1)
  )
}





# ---------------------------- Algorithm 4 -------------------------------
compute_Z <- function(Sites_train, Xcal, Acal, h, eps, delta, nu,
                      C_K = 1, d = ncol(Xcal), pi_method = c("dp","empirical")) {
  pi_method <- match.arg(pi_method)
  S <- length(Sites_train)
  m <- nrow(Xcal)
  
  comps <- lapply(seq_len(S), function(s) {
    site_private_components(Sites_train[[s]], Xcal, h, eps[s], delta[s],
                            C_K, d, kernel_fun = gaussian_kernel,
                            pi_method = pi_method)
  })
  
  pi_agg <- Reduce(`+`, lapply(seq_len(S), function(s) nu[s] * comps[[s]]$pi_tilde))
  
  p_agg <- list(
    '0' = Reduce(`+`, lapply(seq_len(S), function(s) nu[s] * comps[[s]]$p_X_given_A[['0']])),
    '1' = Reduce(`+`, lapply(seq_len(S), function(s) nu[s] * comps[[s]]$p_X_given_A[['1']]))
  )
  px1_agg <- list(
    '0' = Reduce(`+`, lapply(seq_len(S), function(s) nu[s] * comps[[s]]$p_XY1_given_A[['0']])),
    '1' = Reduce(`+`, lapply(seq_len(S), function(s) nu[s] * comps[[s]]$p_XY1_given_A[['1']]))
  )
  eta_agg <- list(
    '0' = pmin(pmax(px1_agg[['0']] / pmax(p_agg[['0']], .Machine$double.eps), 0), 1),
    '1' = pmin(pmax(px1_agg[['1']] / pmax(p_agg[['1']], .Machine$double.eps), 0), 1)
  )
  
  Z <- numeric(m)
  i0 <- which(Acal == 0)
  i1 <- which(Acal == 1)
  if (length(i0)) Z[i0] <- 2 * (-1) * pi_agg[['0']] * (eta_agg[['0']][i0] - 0.5)
  if (length(i1)) Z[i1] <- 2 * (+1) * pi_agg[['1']] * (eta_agg[['1']][i1] - 0.5)
  
  list(Z = Z, pi_agg = pi_agg, eta_agg = eta_agg,
       p_agg = p_agg, px1_agg = px1_agg)
}




FDP_Binary_Tree <- function(Zs, As, eps_s, delta_s, M, tree_noise_sd = NULL) {
  # Uniform grids
  J    <- 2^M + 1
  taus <- seq(-1, 1, length.out = J)
  
  # Histogram counts per a per leaves
  hist_by_a <- function(a) {
    z <- Zs[As == a]
    if (!length(z)) return(rep(0, J - 1))
    z <- pmin(pmax(z, taus[1] + 1e-12), tail(taus, 1) - 1e-12)
    as.integer(tabulate(findInterval(z, taus, rightmost.closed = TRUE), nbins = J) [1:(J-1)])
  }
  counts_M_0 <- hist_by_a(0)
  counts_M_1 <- hist_by_a(1)
  
  # Build trees
  build_tree <- function(counts_M) {
    counts <- vector("list", M)
    counts[[M]] <- counts_M
    if (M > 1) {
      for (ell in (M-1):1) {
        child  <- counts[[ell+1]]
        parent <- as.integer(tapply(child, rep(1:(2^ell), each = 2), sum))
        counts[[ell]] <- parent
      }
    }
    counts
  }
  counts0 <- build_tree(counts_M_0)
  counts1 <- build_tree(counts_M_1)
  
  # Add Gaussian noise
  if (is.null(tree_noise_sd)) {
    var_leaf <- (4 * log(1/delta_s) / eps_s + 2) / (eps_s / M)
    tree_noise_sd <- sqrt(var_leaf)
  }
  
  N0 <- vector("list", M)
  N1 <- vector("list", M)
  for (ell in 1:M) {
    K <- 2^ell
    N0[[ell]] <- counts0[[ell]] + rnorm(K, 0, tree_noise_sd)
    N1[[ell]] <- counts1[[ell]] + rnorm(K, 0, tree_noise_sd)
  }
  
  list(N0 = N0, N1 = N1, taus = taus, M = M)
}



# ---------------------------- Algorithm 6 -------------------------------
non_increasing <- function(g, omega = NULL, Cw = 1.0, M = 6,
                           mu_s = NULL, delta_s = NULL,
                           ntilde_s = NULL, eps_s = NULL,
                           eta = 1e-2) {
  n <- length(g)
  if (n <= 1) return(g)
  

  if (is.null(omega)) {
    denom <- ntilde_s^2 * eps_s^2
    temp  <- (mu_s^2) * (M^4) * log(1 / delta_s) * log(M / eta) / denom
    omega <- Cw * sqrt(sum(temp))
  }
  

  u <- pmin(g + omega, 1)   # upper band
  l <- pmax(g - omega, 0)   # lower band
  
  # Right-to-left construction
  f <- numeric(n)
  f[n] <- u[n]
  for (i in (n-1):1) {
    f[i] <- min(u[i], max(l[i], f[i+1]))
  }
  
  f <- pmin(pmax(f, l), u)
  
  f
}



# ---------------------------- Algorithm 3 -------------------------------
# tails and totals from trees
tails_from_tree <- function(tree, j) {
  M <- tree$M
  
  tail_for_group <- function(Nlist) {
    count <- 0
    for (ell in 1:M) {
      K <- 2^ell
      k <- 1:K
      
      left_idx_child   <- (k - 1) * 2^(M - ell) + 1
      left_idx_parent  <- (ceiling(k/2) - 1) * 2^(M - ell + 1) + 1
      
      use <- (left_idx_child >= j) & (left_idx_parent < j)
      if (any(use)) count <- count + sum(Nlist[[ell]][use])
    }
    count
  }
  
  # tails
  tail0 <- tail_for_group(tree$N0)
  tail1 <- tail_for_group(tree$N1)
  
  # totals from the top layer
  Ns0 <- sum(tree$N0[[1]])
  Ns1 <- sum(tree$N1[[1]])
  
  # robust clipping
  clip_pos <- function(x) pmax(x, 0)
  Ns0 <- max(clip_pos(Ns0), 1e-8)
  Ns1 <- max(clip_pos(Ns1), 1e-8)
  tail0 <- min(max(clip_pos(tail0), 0), Ns0)
  tail1 <- min(max(clip_pos(tail1), 0), Ns1)
  
  list(tail0 = tail0, tail1 = tail1, Ns0 = Ns0, Ns1 = Ns1)
}


FDP_Threshold_Search <- function(site_trees, mu_s, Cw = 1.0, omega = NULL, 
                                 eps_s, delta_s, ntilde_s, eta) {
  taus <- site_trees[[1]]$taus
  L <- length(taus) - 1
  DD_hat <- numeric(L)
  
  DD_hat[1] <- 1
  
  for (j in 2:L) {
    site_dd <- vapply(seq_along(site_trees), function(s) {
      temp <- tails_from_tree(site_trees[[s]], j)
      Ns1 <- max(temp$Ns1, 1e-8)
      Ns0 <- max(temp$Ns0, 1e-8)
      (temp$tail1 /Ns1) - (Ns0-temp$tail0) /Ns0
    }, numeric(1))
    DD_hat[j] <- sum(mu_s * site_dd)
  }
  
  f <- non_increasing(g = DD_hat, omega = omega, Cw = Cw, M = site_trees[[1]]$M,
                      mu_s = mu_s, delta_s = delta_s, ntilde_s = ntilde_s,
                      eps_s = eps_s, eta = eta)
  
  list(taus = taus[1:L], DD = DD_hat, DD_smoothed = f)
}



# ---------------------------- Algorithm 2 -------------------------------
select_tau_fdp <- function(taus, DD, alpha, rho = 0, tau_cap = NULL) {
  stopifnot(is.numeric(taus), is.numeric(DD), length(taus) == length(DD))
  ok <- is.finite(taus) & is.finite(DD)
  taus <- taus[ok]; DD <- DD[ok]
  if (!length(taus)) return(list(tau = 0, j = NA_integer_, feasible = FALSE))
  
  if (!is.null(tau_cap) && is.finite(tau_cap) && tau_cap >= 0) {
    taus <- pmax(pmin(taus,  tau_cap), -tau_cap)
  } else {
    tau_cap <- Inf
  }
  
  # locate τ
  j0 <- which.min(abs(taus))
  D0 <- DD[j0]
  
  low  <- if (rho > 0) max(alpha - rho, 0) else alpha
  high <- if (rho > 0) alpha + rho else alpha
  
  # If |DD(0)| <= alpha + rho, accept tau = 0 
  if (abs(D0) <= high) {
    return(list(tau = 0, j = j0, feasible = TRUE))
  }
  
  # Decide side
  if (D0 < -alpha) {
    
    side_idx <- which(taus <=  1e-15)
    target   <- -alpha
    
  } else if (D0 > alpha) {
    
    side_idx <- which(taus >= -1e-15)
    target   <- +alpha
    
  } else {
    
    return(list(tau = 0, j = j0, feasible = TRUE))
    
  }
  if (!length(side_idx)) side_idx <- seq_along(taus)
  
  tg <- taus[side_idx]
  fg <- DD[side_idx]        
  ord <- order(tg); tg <- tg[ord]; fg <- fg[ord]
  
  feas <- which(abs(fg - target) <= (if (rho > 0) rho else 0))
  if (length(feas)) {
    idx <- feas[which.min(abs(tg[feas]))]
    tie <- feas[which(abs(tg[feas]) == abs(tg[idx]))]
    if (length(tie) > 1) idx <- tie[which.min(abs(fg[tie] - target))]
    return(list(tau = tg[idx], j = side_idx[ord][idx], feasible = TRUE))
  }
  
  # No feasible point: pick minimal violation on that side
  viol <- abs(fg - target)
  idx  <- which.min(viol)
  tie  <- which(viol == viol[idx])
  if (length(tie) > 1) idx <- tie[which.min(abs(tg[tie]))]
  tau  <- tg[idx]
  jret <- side_idx[ord][idx]
  
  # Try linear interpolation
  neigh <- intersect(c(idx - 1L, idx + 1L), seq_along(tg))
  if (length(neigh)) {
    k <- neigh[which.min(abs(fg[neigh] - target))]
    if ((fg[idx] - target) * (fg[k] - target) < 0) {
      w <- (target - fg[idx]) / (fg[k] - fg[idx])  # in (0,1)
      tau <- tg[idx] + w * (tg[k] - tg[idx])
      tau <- sign(tau) * min(abs(tau), tau_cap)
      jret <- side_idx[ord][ if (abs(tg[idx]-tau) <= abs(tg[k]-tau)) idx else k ]
    }
  }
  
  list(tau = as.numeric(tau), j = as.integer(jret), feasible = FALSE)
}





FDP_fair_classifier_once <- function(Sites_train, Sites_cal,
                                     alpha, h, eps, delta, mu, nu,
                                     C_K = 1, Cw = 1, rho = NULL, M = 6,
                                     omega = NULL, eta = 1e-2,
                                     pi_method = c("dp","empirical")) {
  pi_method <- match.arg(pi_method)
  S <- length(Sites_train)
  
  per_site <- vector("list", S)
  for (r in seq_len(S)) {
    Xcal_r <- Sites_cal[[r]]$X
    Acal_r <- Sites_cal[[r]]$A
    agg_r <- compute_Z(Sites_train = Sites_train, Xcal = Xcal_r, Acal = Acal_r,
                       h = h, eps = eps, delta = delta, nu = nu,
                       C_K = C_K, d = ncol(Xcal_r), pi_method = pi_method)
    per_site[[r]] <- list(train = Sites_train[[r]], cal = Sites_cal[[r]],
                          Zs = agg_r$Z, pi_agg_on_cal = agg_r$pi_agg, eta_agg_on_cal = agg_r$eta_agg)
  }
  
  site_trees <- lapply(seq_len(S), function(s) {
    FDP_Binary_Tree(Zs = per_site[[s]]$Zs, As = per_site[[s]]$cal$A,
                    eps_s = eps[s], delta_s = delta[s], M = M)
  })
  ntilde_vec <- vapply(Sites_cal, function(x) length(x$Y), numeric(1))
  thr <- FDP_Threshold_Search(site_trees, mu_s = mu, Cw = Cw,
                              omega = omega, eps_s = eps, delta_s = delta,
                              ntilde_s = ntilde_vec, eta = eta)
  taus <- thr$taus
  DD   <- thr$DD_smoothed
  rho0 <- if (is.null(rho)) 0 else rho
  
  pi0_bar <- sum(nu * sapply(per_site, function(z) as.numeric(z$pi_agg_on_cal[['0']])))
  pi1_bar <- sum(nu * sapply(per_site, function(z) as.numeric(z$pi_agg_on_cal[['1']])))
  tau_cap <- 2 * max(min(pi0_bar, pi1_bar), 0)
  sel <- select_tau_fdp(taus = taus, DD = DD, alpha = alpha, rho = rho0, tau_cap = tau_cap)
  tau_hat <- sel$tau
  
  predict <- function(X, A) {
    m <- nrow(X)
    comps_list <- lapply(seq_len(S), function(s) {
      site_private_components(train = Sites_train[[s]], Xeval = X,
                              h = h, eps_s = eps[s], delta_s = delta[s], C_K = C_K, d = ncol(X),
                              kernel_fun = gaussian_kernel, pi_method = pi_method)
    })
    pi_agg <- Reduce(`+`, lapply(seq_len(S), function(s) nu[s] * comps_list[[s]]$pi_tilde))
    p_agg <- list(
      '0' = Reduce(`+`, lapply(seq_len(S), function(s) nu[s] * comps_list[[s]]$p_X_given_A[['0']])),
      '1' = Reduce(`+`, lapply(seq_len(S), function(s) nu[s] * comps_list[[s]]$p_X_given_A[['1']]))
    )
    px1_agg <- list(
      '0' = Reduce(`+`, lapply(seq_len(S), function(s) nu[s] * comps_list[[s]]$p_XY1_given_A[['0']])),
      '1' = Reduce(`+`, lapply(seq_len(S), function(s) nu[s] * comps_list[[s]]$p_XY1_given_A[['1']]))
    )
    eta_agg <- list(
      '0' = pmin(pmax(px1_agg[['0']] / pmax(p_agg[['0']], .Machine$double.eps), 0), 1),
      '1' = pmin(pmax(px1_agg[['1']] / pmax(p_agg[['1']], .Machine$double.eps), 0), 1)
    )
    
    thr0 <- 0.5 - tau_hat / (2 * pi_agg[['0']])
    thr1 <- 0.5 + tau_hat / (2 * pi_agg[['1']])
    
    y <- integer(m)
    i0 <- which(A == 0); if (length(i0)) y[i0] <- as.integer(eta_agg[['0']][i0] >= thr0)
    i1 <- which(A == 1); if (length(i1)) y[i1] <- as.integer(eta_agg[['1']][i1] >= thr1)
    y
  }
  
  list(tau_hat = tau_hat, taus = taus, DD_used = DD, per_site = per_site,
       site_trees = site_trees, predict = predict)
}



# ---------------------------- Algorithm 2 (cross fit) -------------------------------
FDP_fair_classifier_cross_fit <- function(Sites_train, Sites_cal,
                                          alpha, h, eps, delta, mu, nu,
                                          C_K = 1, Cw = 1, rho = NULL, M = 6,
                                          omega = NULL, eta = 1e-2,
                                          pi_method = c("dp","empirical")) {
  pi_method <- match.arg(pi_method)
  
  runA <- FDP_fair_classifier_once(Sites_train, Sites_cal,
                                   alpha, h, eps, delta, mu, nu, C_K, Cw, rho, M, omega, eta,
                                   pi_method = pi_method)
  runB <- FDP_fair_classifier_once(Sites_cal, Sites_train,
                                   alpha, h, eps, delta, mu, nu, C_K, Cw, rho, M, omega, eta,
                                   pi_method = pi_method)
  
  predict <- function(X, A) {
    y1 <- runA$predict(X, A)
    y2 <- runB$predict(X, A)
    n  <- length(y1)
    return(rbinom(n, 1, (y1 + y2) / 2))
  }
  
  list(fold1 = runA, fold2 = runB,
       tau_A = runA$tau_hat, tau_B = runB$tau_hat,
       predict = predict)
}





# ---------------------------- Empirical DD for test dataset -------------------------------
empirical_DD <- function(y_hat, A) {
  A <- as.integer(A)
  p1 <-  mean(y_hat[A == 1] == 1) 
  p0 <-  mean(y_hat[A == 0] == 1) 
  return(list(DD = p1 - p0, abs_DD = abs(p1 - p0), p1 = p1, p0 = p0))
}





# ---------------------------- Cross validation -------------------------------
eta_from_kde <- function(Xtrain, Ytrain, Xeval, h,
                         kernel_fun = gaussian_kernel, tol = 1e-12) {
  X1 <- Xtrain[Ytrain == 1, , drop = FALSE]
  X0 <- Xtrain[Ytrain == 0, , drop = FALSE]
  
  S1 <- if (nrow(X1)) kde_vec(X1, Xeval, h, kernel_fun) else rep(0, nrow(Xeval))
  S0 <- if (nrow(X0)) kde_vec(X0, Xeval, h, kernel_fun) else rep(0, nrow(Xeval))
  
  den <- S1 + S0
  prior <- if (length(Ytrain)) mean(Ytrain) else 0.5
  eta  <- ifelse(!is.finite(den) | den <= tol, prior, S1 / den)
  pmin(pmax(eta, 0), 1)
}




cv_select_h <- function(data, h_grid, n_folds = 3,
                        kernel_fun = gaussian_kernel, tol = 1e-12) {
  S <- length(data$all)
  idx <- sample(1:S,1)
  X <- data$all[[idx]]$X
  Y <- as.integer(data$all[[idx]]$Y)
  n <- nrow(X)
  
  
  fold_ids <- sample(rep(1:n_folds, length.out = n))
  
  mean_err <- numeric(length(h_grid))
  for (j in 1:length(h_grid)) {
    h <- h_grid[j]
    errs <- numeric(n_folds)
    for (k in 1:n_folds) {
      idx_val <- which(fold_ids == k)
      idx_tr  <- which(fold_ids != k)
      eta_vl  <- eta_from_kde(X[idx_tr,,drop=FALSE], Y[idx_tr], X[idx_val,,drop=FALSE],
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










