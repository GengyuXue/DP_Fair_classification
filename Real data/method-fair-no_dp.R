fair_eval_eta_for_a <- function(Xeval, input_data, a, h,
                           d = ncol(Xeval),
                           kernel_fun = gaussian_kernel, C_K = 1, tol = 1e-8,
                           p_floor = 1e-12, use_cal = FALSE, fac = NULL) {
  
  if (!isTRUE(use_cal)) {
    if (a == 0) { Xa1 <- input_data$X_01_train; Xa0 <- input_data$X_00_train }
    else        { Xa1 <- input_data$X_11_train; Xa0 <- input_data$X_10_train }
  } else {
    if (a == 0) { Xa1 <- input_data$X_01_cal;   Xa0 <- input_data$X_00_cal   }
    else        { Xa1 <- input_data$X_11_cal;   Xa0 <- input_data$X_10_cal   }
  }
  n_a <- nrow(Xa1) + nrow(Xa0)
  m   <- nrow(Xeval)
  
  if (n_a == 0 || m == 0) {
    return(list(p_hat = rep(NA_real_, m), eta_hat = rep(NA_real_, m)))
  }
  
  S1 <- kde_vec(Xa1, Xeval, h, kernel_fun)
  S0 <- kde_vec(Xa0, Xeval, h, kernel_fun)
  
  
  p_hat <- (S1 + S0) / n_a
  bad <- !is.finite(p_hat) | (p_hat <= 0)
  if (any(bad)) p_hat[bad] <- p_floor
  
  numer <- S1 / n_a 
  numer[!is.finite(numer)] <- 0
  
  eta_hat <- pmin(pmax(numer / p_hat, 0), 1)
  list(p_hat = p_hat, eta_hat = eta_hat)
}



#' find the cutoff point through binary search
fair_calibrate_tau_bin_search <- function(
    input_data, h, alpha,
    kernel_fun = gaussian_kernel, C_K = 1, tol = 1e-8,
    pi_tilde_1, pi_tilde_0, switch = FALSE){
  
  
  if (switch == FALSE){
    Xcal_0 <- rbind(input_data$X_00_cal, input_data$X_01_cal)  # A=0
    Xcal_1 <- rbind(input_data$X_10_cal, input_data$X_11_cal)  # A=1
  }
  
  else {
    Xcal_0 <- rbind(input_data$X_00_train, input_data$X_01_train)  # A=0
    Xcal_1 <- rbind(input_data$X_10_train, input_data$X_11_train)  # A=1
  }
  
  n0_cal <- nrow(Xcal_0); n1_cal <- nrow(Xcal_1)
  d <- if (n0_cal + n1_cal > 0) ncol(rbind(Xcal_0, Xcal_1)) else 2
  
  eta0_cal <- if (n0_cal) fair_eval_eta_for_a(Xcal_0, input_data, 0, h, d, 
                                         kernel_fun, C_K, tol, use_cal = switch)$eta_hat else numeric(0)
  eta1_cal <- if (n1_cal) fair_eval_eta_for_a(Xcal_1, input_data, 1, h, d, 
                                         kernel_fun, C_K, tol, use_cal = switch)$eta_hat else numeric(0)
  
  if (min(n0_cal, n1_cal) == 0L) {
    warning("One calibration group is empty.")
  }
  

  pi0 <- max(pi_tilde_0, 1e-6)  
  pi1 <- max(pi_tilde_1, 1e-6)
  
  DD_tilde <- function(tau){
    thr1 <- 0.5 + tau / (2 * pi1)
    thr0 <- 0.5 - tau / (2 * pi0)
    part1 <- if (n1_cal) mean(eta1_cal >= thr1) else 0
    part0 <- if (n0_cal) mean(eta0_cal >= thr0) else 0
    part1 - part0 
  }
  
  DD0 <- DD_tilde(0)
  
  if (abs(DD0) <= alpha) {
    return(list(tauhat = 0, DD0 = DD0))
  }
  

  else {
    temp <- max(pi_tilde_1 , pi_tilde_0) + 0.2
    tau_min <- -temp
    tau_max <- temp
    if (DD0>alpha){ 
      tau_min = 0
      tau_old = tau_min
      tau_new = (tau_min+tau_max)/2
      while((tau_max-tau_min)>1e-6){
        Dhat_new = DD_tilde(tau_new)
        tau_old = tau_new
        if(Dhat_new > alpha){
          tau_min = tau_old
        }else{
          tau_max = tau_old
        }
        tau_new = (tau_min+tau_max)/2
      }
      tauhat = tau_old
    } 
    
    else {
      tau_max = 0
      tau_old = tau_max
      tau_new = (tau_min+tau_max)/2
      while((tau_max-tau_min)>1e-6){
        Dhat_new = DD_tilde(tau_new)
        tau_old = tau_new
        if(Dhat_new>(-alpha)){
          tau_min = tau_old
        }else{
          tau_max = tau_old
        }
        tau_new = (tau_min+tau_max)/2
      }
      tauhat = tau_old
    }
  }
  
  
  if(abs(DD_tilde(tauhat)) >= alpha ){ warning("No feasible tau found during bracketing.")}
  
  return(list(tauhat= tauhat, DD = DD_tilde(tauhat)))
  
}



DD_train <- function(test_data, y_hat){
  A <- as.integer(test_data$A)
  
  p1 <- mean(y_hat[A == 1] == 1)
  p0 <- mean(y_hat[A == 0] == 1)
  gap <- p1 - p0
  
  list(
    disparity = gap,        
    abs_disparity = abs(gap) 
  )
}


fair_classifier <- function(test_data, input_data, alpha, h,
                                kernel_fun = gaussian_kernel, tol, C_K = 1, 
                                cross_fit = TRUE){
  d <- ncol(input_data$X)
  input_data <- sample_split(input_data, split_ratio = 0.5)
  test_X <- test_data$X
  
  n_1_train <- nrow(input_data$X_10_train) + nrow(input_data$X_11_train)
  n_0_train <- nrow(input_data$X_00_train) + nrow(input_data$X_01_train)
  n_train <- n_1_train + n_0_train
  
  
  pi_tilde_1 <- n_1_train/n_train 
  pi_tilde_0 <- n_0_train/n_train 
  
  

  # Evaluation of Kernels
  compute_for_a <- function(a) {
    if (a == 0) {
      Xa1 <- input_data$X_01_train; Xa0 <- input_data$X_00_train
      na <- n_0_train
    } else {
      Xa1 <- input_data$X_11_train; Xa0 <- input_data$X_10_train
      na <- n_1_train
    }
    
    S1 <- kde_vec(Xa1, test_X, h, kernel_fun)
    S0 <- kde_vec(Xa0, test_X, h, kernel_fun)
    
    p_hat <- (S1 + S0) / na 
    p_hat <- pmax(p_hat, .Machine$double.eps)
    
    numer <- S1 / na
    eta_hat <- pmin(pmax(numer / p_hat, 0), 1)
    
    list(p_hat = p_hat, eta_hat = eta_hat)
  }
  
  eta_output_0 <- compute_for_a(0)
  eta_output_1 <- compute_for_a(1)
  

  
  tau_data <- fair_calibrate_tau_bin_search(input_data, h, alpha,
                                           kernel_fun = gaussian_kernel, C_K = 1, tol = 1e-8,
                                           pi_tilde_1, pi_tilde_0, switch = FALSE)
  
  
  tau_hat <- tau_data$tauhat
  
  
  thr0 <- 0.5 - tau_hat / (2 * pi_tilde_0)
  thr1 <- 0.5 + tau_hat / (2 * pi_tilde_1)
  
  eta_test <- ifelse(test_data$A == 0, eta_output_0$eta_hat, eta_output_1$eta_hat)
  thr_vec  <- ifelse(test_data$A == 0, thr0, thr1)
  y_hat   <- as.integer(eta_test >= thr_vec)
  DD_hat <- DD_train(test_data, y_hat)
  
  if (cross_fit == FALSE) {
    DD_hat <- DD_train(test_data, y_hat)
    return(list(pi_tilde_1 = pi_tilde_1, pi_tilde_0 =pi_tilde_0,
                eta_0 = eta_output_0$eta_hat, eta_1 = eta_output_1$eta_hat,
                tau_hat = tau_hat, DD_hat = DD_hat$abs_disparity, y_hat = y_hat))
  }
  
  else {
    n_1_cal <- nrow(input_data$X_10_cal) + nrow(input_data$X_11_cal)
    n_0_cal <- nrow(input_data$X_00_cal) + nrow(input_data$X_01_cal)
    n_cal <- n_1_cal + n_0_cal
    
    pi_tilde_1_cal <- n_1_cal/n_cal
    pi_tilde_0_cal <- n_0_cal/n_cal 

    
    # Evaluation of Kernels
    compute_for_a_cal <- function(a) {
      if (a == 0) {
        Xa1 <- input_data$X_01_cal; Xa0 <- input_data$X_00_cal
        na <- n_0_cal
      } else {
        Xa1 <- input_data$X_11_cal; Xa0 <- input_data$X_10_cal
        na <- n_1_cal
      }
      
      S1 <- kde_vec(Xa1, test_X, h, kernel_fun)
      S0 <- kde_vec(Xa0, test_X, h, kernel_fun)
      
      p_hat <- (S1 + S0) / na 
      p_hat <- pmax(p_hat, .Machine$double.eps)
      
      numer <- S1 / na 
      eta_hat <- pmin(pmax(numer / p_hat, 0), 1)
      
      list(p_hat = p_hat, eta_hat = eta_hat)
    }
    
    eta_output_0_cal <- compute_for_a_cal(0)
    eta_output_1_cal <- compute_for_a_cal(1)
    
    
    tau_data_cal <- fair_calibrate_tau_bin_search(input_data, h, alpha,
                                                 kernel_fun = gaussian_kernel, C_K = 1, tol = 1e-8,
                                                 pi_tilde_1_cal, pi_tilde_0_cal, switch = TRUE)
    
    
    tau_hat_cal <- tau_data_cal$tauhat
    DD_hat_cal <- tau_data_cal$DD

    thr0_cal <- 0.5 - tau_hat_cal / (2 * pi_tilde_0_cal)
    thr1_cal <- 0.5 + tau_hat_cal / (2 * pi_tilde_1_cal)
    
    eta_test_cal <- ifelse(test_data$A == 0, eta_output_0_cal$eta_hat, eta_output_1_cal$eta_hat)
    thr_vec_cal  <- ifelse(test_data$A == 0, thr0_cal, thr1_cal)
    y_hat_cal   <- as.integer(eta_test_cal >= thr_vec_cal)
    
    # cross-fit
    y_hat_out <- rbinom(length(y_hat), 1, 0.5*y_hat + 0.5 * y_hat_cal)
    DD_hat_out <- DD_train(test_data, y_hat_out)
    
    return(list(pi_tilde_1_train = pi_tilde_1, pi_tilde_0_train =pi_tilde_0,
                pi_tilde_1_cal = pi_tilde_1_cal, pi_tilde_0_cal =pi_tilde_0_cal,
                eta_0_train = eta_output_0$eta_hat, eta_1_train = eta_output_1$eta_hat,
                eta_0_cal = eta_output_0_cal$eta_hat, eta_1_cal = eta_output_1_cal$eta_hat,
                tau_hat_train = tau_hat, tau_hat_cal=tau_hat_cal,
                DD_hat_train = DD_hat, DD_hat_cal=DD_hat_cal,
                y_hat_train = y_hat, y_hat_cal = y_hat_cal,
                y_hat_out = y_hat_out, DD_hat_out = DD_hat_out$abs_disparity))
    
  }
  
  
  
}








