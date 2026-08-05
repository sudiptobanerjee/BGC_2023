# src/rawSamplerDivideAndConquer.r

#' Exact Divide and Conquer Sampler for Bayesian Linear Regression
#' 
#' Performs sequential Bayesian updating over T blocks of data, updating the 
#' conjugate Normal-Gamma parameters analytically.
#' 
#' @param X_list A list of design matrices (length T).
#' @param y_list A list of response vectors (length T).
#' @param M0_inv Prior precision matrix.
#' @param m0 Prior precision-weighted mean vector.
#' @param a0 Prior shape parameter for Gamma.
#' @param b0 Prior rate parameter for Gamma.
#' @param N Number of posterior samples to draw at the end.
rawSamplerDivideAndConquer <- function(X_list, y_list, M0_inv, m0, a0, b0, N) {
  
  T_blocks <- length(X_list)
  
  # Initialize working parameters with prior values
  M_inv_t <- M0_inv
  m_t <- m0
  a_t <- a0
  b_t <- b0
  
  # Sequential Bayesian updating across blocks
  for (t in 1:T_blocks) {
    X <- X_list[[t]]
    y <- y_list[[t]]
    
    n <- nrow(X)
    p <- ncol(X)
    
    # 1. Compute Cross-Products for current block
    XtX <- crossprod(X)
    XtY <- crossprod(X, y)
    YtY <- crossprod(y)
    
    # 2. Analytic Update for Precision and Mean
    M_inv_new <- M_inv_t + XtX
    m_new <- m_t + XtY
    
    # Update Shape Parameter
    a_new <- a_t + (n / 2)
    
    # Update Rate Parameter (b) avoiding full matrix inversions
    U_M_inv_new <- chol(M_inv_new)
    Mm_new <- backsolve(U_M_inv_new, backsolve(U_M_inv_new, m_new, transpose = TRUE))
    m_new_term <- crossprod(U_M_inv_new %*% Mm_new)
    
    # Prior penalty term from the previous block's parameters
    if (all(m_t == 0)) {
      m_t_term <- 0
    } else {
      U_M_inv_t <- chol(M_inv_t)
      mu_t <- backsolve(U_M_inv_t, backsolve(U_M_inv_t, m_t, transpose = TRUE))
      m_t_term <- crossprod(m_t, mu_t)
    }
    
    b_new <- as.numeric(b_t + 0.5 * (YtY + m_t_term - m_new_term))
    
    # Assign updated parameters for the next iteration t+1
    M_inv_t <- M_inv_new
    m_t <- m_new
    a_t <- a_new
    b_t <- b_new
  }
  
  # --- End of Divide & Conquer. Draw samples from the final posterior ---
  
  U_Minv <- chol(M_inv_t)
  Mm <- backsolve(U_Minv, backsolve(U_Minv, m_t, transpose = TRUE))
  
  # Vectorized sampling for precision/variance
  tau_samples <- rgamma(N, shape = a_t, rate = b_t)
  sigma2_samples <- 1 / tau_samples
  sigma_samples <- sqrt(sigma2_samples)
  
  # Vectorized conditional draw for Beta
  Z <- matrix(rnorm(p * N), nrow = p, ncol = N)
  W <- backsolve(U_Minv, Z)
  scaled_noise <- sweep(W, 2, sigma_samples, `*`)
  
  beta_samples <- t(as.numeric(Mm) + scaled_noise)
  
  if (!is.null(colnames(X_list[[1]]))) {
    colnames(beta_samples) <- colnames(X_list[[1]])
  }
  
  posterior_draws <- cbind(beta_samples, sigma = sigma_samples, sigma2 = sigma2_samples, tau = tau_samples)
  
  return(list(
    beta = beta_samples,      
    sigma2 = sigma2_samples,  
    draws = posterior_draws,
    final_params = list(M_inv = M_inv_t, m = m_t, a = a_t, b = b_t)
  ))
}

