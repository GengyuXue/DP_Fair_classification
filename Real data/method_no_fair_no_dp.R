classifier_xy <- function(new_data, train_target,
                                 h, kernel_f = gaussian_kernel) {
  n_train <- nrow(train_target$X)
  n_test  <- nrow(new_data$X)
  
  # kernel estimate 
  sum_diff <- numeric(n_test)
  for (j in 1:n_test) {
    kw <- kernel_f(train_target$X, new_data$X[j, ], h)
    sum_diff[j] <- sum((train_target$Y - 0.5) * kw)
  }
  
  p_estimates <- sum_diff / n_train
  
  preds <- as.integer(p_estimates > 0)
  list(p_estimates = p_estimates, predictions = preds)
}




classifier_xya <- function(train_data, test_data, h,
                          kernel_fun = gaussian_kernel) {

  Xtr <- train_data$X
  Ytr <- train_data$Y
  Atr <- train_data$A
  
  Xte <- test_data$X
  A_te <- test_data$A
  
  yhat <- integer(length(A_te))
  
  for (a in sort(unique(Atr))) {
    # train
    Xa1 <- Xtr[Atr == a & Ytr == 1, , drop = FALSE]
    Xa0 <- Xtr[Atr == a & Ytr == 0, , drop = FALSE]
    
    # test
    idx <- which(A_te == a)
    Xeval <- Xte[idx, , drop = FALSE]
    
    # ---- KDE density estimates ----
    f1 <- kde_vec(Xa1, Xeval, h, kernel_fun)
    f0 <- kde_vec(Xa0, Xeval, h, kernel_fun)
    
    # ---- Compute eta_A, compare to 1/2 ----
    eta <- f1 / (f1 + f0 + 1e-12)
    yhat[idx] <- as.integer(eta >= 0.5)
  }
  return(yhat)
}

