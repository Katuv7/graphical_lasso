coordinate_descent_lasso <- function(W11, s12, rho, max_iter = 10000, tol = 1e-4) {
  p <- length(s12)
  beta <- rep(0, p)
  for (iter in 1:max_iter) {
    beta_old <- beta
    for (j in 1:p) {
      # Soft-thresholding operator
      S <- function(x, lambda) sign(x) * max(abs(x) - lambda, 0)
      
      # Update rule for beta_j
      r_j <- s12[j] - W11[-j, j] %*% beta[-j]
      beta[j] <- S(r_j, rho) / W11[j, j]
    }
    # Check for convergence
    if (sqrt(sum((beta - beta_old)^2)) < tol) {
      break
    }
  }
  return(beta)
}


graphical_lasso <- function(S, rho, threshold = 0.001, max_iter = 10000) {
  p <- ncol(S)
  W <- S + rho * diag(p)
  for (iter in 1:max_iter) {
    W_old <- W
    for (j in 1:p) {
      # Partition W and S
      W11 <- W[-j, -j]
      s12 <- S[-j, j]
      
      # Solve the lasso problem for beta_hat
      beta_hat <- coordinate_descent_lasso(W11, s12, rho)
      
      # Update W
      W[-j, j] <- W11 %*% beta_hat
      W[j, -j] <- W[-j, j]
    }
    
    # Check for convergence
    if (mean(abs(W - W_old)) < threshold) {
      break
    }
  }
  
  # Invert W to get Theta
  Theta <- solve(W)
  colnames(Theta) <- colnames(S)
  rownames(Theta) <- rownames(S)
  return(list(W = W, Theta = Theta))
}