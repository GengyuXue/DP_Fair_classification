p_kernel_target_only <- function(new_data, train_target,
                                 h , priv = TRUE,
                                 eps, delta,
                                 kernel_f = gaussian_kernel, C_K) {
  n_train <- nrow(train_target$X)
  n_test  <- nrow(new_data$X)
  d <- ncol(new_data$X)
  
  # kernel estimate 
  sum_diff <- numeric(n_test)
  for (j in 1:n_test) {
    kw <- kernel_f(train_target$X, new_data$X[j, ], h)
    sum_diff[j] <- sum((train_target$Y - 0.5) * kw)
  }

  p_estimates <- sum_diff / n_train
  
  if (priv) {
    
    fac <- kernel_eig_factor_trunc(new_data$X, h, tol=1e-8, r = 35)
    
    
    Z <- GP_draw_from_factor(fac)
    
    scl <-  sqrt(6 * C_K * log(2/delta))/(n_train * eps * h^d)
    
    noise <- scl * Z
    
    p_estimates <- p_estimates + noise
  }
  
  preds <- as.integer(p_estimates > 0)
  list(p_estimates = p_estimates, predictions = preds)
}


