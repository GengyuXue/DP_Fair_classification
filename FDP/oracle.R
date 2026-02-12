eta_arctan <- function(x, a, w = c(1.5, 0), b = -0.05, gamma = 0.05, k = 4, eps = 1e-3) {
  s <- b + sum(w * (x-0.5)) + gamma * (2*a-1)
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
  Y <- vapply(seq_len(n), function(i) rbinom(1, 1, eta_fun(X[i,], A[i], w,b,gamma,k)), 1)
  
  
  idx_11 <- which(A == 1 & Y == 1)
  idx_10 <- which(A == 1 & Y == 0)
  idx_01 <- which(A == 0 & Y == 1)
  idx_00 <- which(A == 0 & Y == 0)
  
  list(Y = Y, X = X, A = A,
       X_11 = X[idx_11, , drop = FALSE],
       X_10 = X[idx_10, , drop = FALSE],
       X_01 = X[idx_01, , drop = FALSE],
       X_00 = X[idx_00, , drop = FALSE])
}

eta_vec <- function(X, A, w, b, gamma, k, eps) {
  # vectorized version of eta_arctan
  s <- b + (X - 0.5) %*% w + gamma * (2*A - 1)
  p <- 0.5 + atan(k * as.numeric(s)) / pi
  pmin(pmax(p, eps), 1 - eps)
}

# Given eta and A, evaluate disparity DD(tau)
DD_of_tau <- function(tau, eta, A, pi1, pi0) {
  thr1 <- 0.5 + tau/(2*pi1)
  thr0 <- 0.5 - tau/(2*pi0)
  part1 <- mean(eta[A==1] >= thr1)
  part0 <- mean(eta[A==0] >= thr0)
  part1 - part0
}



# Given eta and A, compute 0-1 classification risk for the DP-thresholds
risk_at_tau <- function(tau, eta, A, pi1, pi0) {
  thr1 <- 0.5 + tau/(2*pi1)
  thr0 <- 0.5 - tau/(2*pi0)
  yhat1 <- as.integer(eta[A==1] >= thr1)
  yhat0 <- as.integer(eta[A==0] >= thr0)
  # plug-in risk E[ I(yhat=1)(1-eta) + I(yhat=0)eta ]
  r1 <- mean( yhat1*(1 - eta[A==1]) + (1 - yhat1)*eta[A==1] )
  r0 <- mean( yhat0*(1 - eta[A==0]) + (1 - yhat0)*eta[A==0] )
  pi1 * r1 + pi0 * r0
}

# Find the smallest |tau| s.t. |DD(tau)| <= alpha.
solve_tau <- function(eta, A, pi1, pi0, alpha, B_init = 2) {
  DD0 <- DD_of_tau(0, eta, A, pi1, pi0)
  if (abs(DD0) <= alpha) return(0)
  
  target <- sign(DD0) * alpha
  g <- function(tau) DD_of_tau(tau, eta, A, pi1, pi0) - target
  
  
  B <- B_init * max(pi1, pi0)
  lo <- -B; hi <-  B
  glo <- g(lo); ghi <- g(hi)
  tries <- 0
  while (glo * ghi > 0 && tries < 12) {
    B <- 2*B
    lo <- -B; hi <- B
    glo <- g(lo); ghi <- g(hi)
    tries <- tries + 1
  }
  
  if (glo * ghi < 0) {
    uniroot(g, interval = c(lo, hi), tol = 1e-8)$root
  } else {
    opt <- optimize(function(tau) abs(g(tau)), interval = c(lo, hi))
    opt$minimum
  }
}

oracle_dp_montecarlo <- function(alpha_list,
                                 n_pop = 200000,
                                 d = 2, p_a = 0.3,
                                 alpha1 = 4, beta1 = 2,
                                 alpha0 = 4.5, beta0 = 2,
                                 w = c(1,1), b = 0, gamma = -0.3, k = 12, eps_clip = 1e-3) {
  pop <- generate_data(n_pop, d, p_a,
                       alpha1, beta1, alpha0, beta0,
                       eta_fun = eta_arctan,
                       w = w, b = b, gamma = gamma, k = k, eps = eps_clip)
  A <- as.integer(pop$A); X <- pop$X
  pi1 <- mean(A == 1); pi0 <- 1 - pi1
  eta <- eta_vec(X, A, w, b, gamma, k, eps_clip)
  out <- lapply(alpha_list, function(a) {
    tau_hat <- solve_tau(eta, A, pi1, pi0, a)
    data.frame(alpha = a,
               tau = tau_hat,
               disparity = DD_of_tau(tau_hat, eta, A, pi1, pi0),
               error = risk_at_tau(tau_hat, eta, A, pi1, pi0))
  })
  do.call(rbind, out)
}
