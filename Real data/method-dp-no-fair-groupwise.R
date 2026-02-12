CDP_unfair_classifier_groupwise <- function(test_data, input_data, h, epsilon, delta,kernel_fun = gaussian_kernel,
                                  tol = 1e-8, C_K = 1, pi_mode = c("private","empirical")) {
  pi_mode <- match.arg(pi_mode)
  d <- ncol(input_data$X)
  input_data <- sample_split(input_data, split_ratio = 0.5)
  test_X <- test_data$X
  pis_tr <- compute_pi_tilde(input_data, split = "train",
                             epsilon = epsilon, delta = delta, pi_mode = pi_mode)
  pi_tilde_1 <- pis_tr$pi_tilde_1; pi_tilde_0 <- pis_tr$pi_tilde_0
  n_1_train <- nrow(input_data$X_10_train) + nrow(input_data$X_11_train)
  n_0_train <- nrow(input_data$X_00_train) + nrow(input_data$X_01_train)
  fac <- kernel_eig_factor_trunc(test_X, h, tol = tol, r = 35)
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
  
  eta_output_0 <- compute_for_a(0)
  eta_output_1 <- compute_for_a(1)

  eta_test <- ifelse(test_data$A == 0, eta_output_0$eta_hat, eta_output_1$eta_hat)
  
  y_hat    <- as.integer(eta_test >= 0.5)
  
    return(list(
      pi_tilde_1 = pi_tilde_1, pi_tilde_0 = pi_tilde_0,
      eta_0 = eta_output_0$eta_hat, eta_1 = eta_output_1$eta_hat, y_hat = y_hat
    ))
}



