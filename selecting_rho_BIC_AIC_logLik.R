
select_rho <- function(S, rho_range, penalize.diagonal = TRUE, threshold = 0.001, max_iter = 10000) {
  # Initialize a data frame to store the results
  results <- data.frame(rho = numeric(length(rho_range)), 
                        AIC = numeric(length(rho_range)), 
                        BIC = numeric(length(rho_range)), 
                        LogLik = numeric(length(rho_range)))
  
  n <- nrow(S) # Number of observations, assuming S is square (p = ncol(S))
  
  for (i in seq_along(rho_range)) {
    rho <- rho_range[i]
    glasso_result <- graphical_lasso(S, rho, threshold, max_iter)
    Theta <- glasso_result$Theta
    
    # Log-Likelihood calculation
    loglik <- -n/2 * (sum(diag(S %*% Theta)) - log(det(Theta)))
    
    # Number of nonzero parameters in Theta
    if (penalize.diagonal) {
      k <- sum(Theta != 0) # If penalizing diagonal, count all non-zero elements
    } else {
      k <- sum(Theta[upper.tri(Theta)] != 0) + sum(Theta[lower.tri(Theta)] != 0) # Count only off-diagonal non-zeros
    }
    
    # AIC and BIC calculations
    aic <- 2 * k - 2 * loglik
    bic <- log(n) * k - 2 * loglik
    
    # Store results in the data frame
    results[i, ] <- c(rho, aic, bic, loglik)
  }
  
  colnames(results) <- c("Rho", "AIC", "BIC", "LogLik")
  
  return(results)
}



rho_range <- seq(0.0001, 0.5, by = 0.05)


select_rho <- function(S, rho_range, threshold=0.001, max_iter=100) {
  results <- data.frame(rho = numeric(length(rho_range)), 
                        AIC = numeric(length(rho_range)), 
                        BIC = numeric(length(rho_range)), 
                        LogLik = numeric(length(rho_range)))
  
  n <- nrow(S) # Number of observations
  
  for (i in seq_along(rho_range)) {
    rho <- rho_range[i]
    glasso_result <- graphical_lasso(S, rho, threshold, max_iter)
    Theta <- glasso_result$Theta
    
    # Improved Log-Likelihood calculation
    loglik <- -n/2 * (sum(diag(S %*% Theta)) - log(det(Theta)) - ncol(S))
    
    # Count of non-zero off-diagonal elements in Theta for penalized parameters
    k <- sum(Theta[lower.tri(Theta)] != 0)
    
    # AIC and BIC calculations (make sure to use the correct formula and counts)
    aic <- 2 * k - 2 * loglik
    bic <- log(n) * k - 2 * loglik
    
    # Store results
    results[i, ] <- c(rho, aic, bic, loglik)
  }
  
  colnames(results) <- c("Rho", "AIC", "BIC", "LogLik")
  return(results)
}
# Find the best rho value
results_selecting_rho <- select_rho(S, rho_range)

