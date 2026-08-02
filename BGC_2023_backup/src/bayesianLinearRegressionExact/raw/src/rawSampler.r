# src/rawSampler.r

rawSampler <- function(y, X, m0, M0_inv, a0, b0, N) {
  
  n <- nrow(X)
  p <- ncol(X)
  
  # 1. Compute Cross-Products[cite: 14]
  XtX <- crossprod(X)
  XtY <- crossprod(X, y)
  
  # 2. Posterior Precision and Cholesky Factorization[cite: 14]
  M_inv <- M0_inv + XtX
  U_Minv <- chol(M_inv) 
  
  # 3. Compute Posterior Mean (Mm) avoiding t() via transpose = TRUE[cite: 14]
  m <- m0 + XtY
  Mm <- backsolve(U_Minv, backsolve(U_Minv, m, transpose = TRUE))
  
  # 4. Exact Marginal Draw for Precision/Variance[cite: 14]
  YtY <- crossprod(y)
  
  # Safely compute the prior term without matrix inversion[cite: 14]
  if (all(m0 == 0)) {
    m0_term <- 0
  } else {
    U_M0inv <- chol(M0_inv)
    mu0 <- backsolve(U_M0inv, backsolve(U_M0inv, m0, transpose = TRUE))
    m0_term <- crossprod(m0, mu0)
  }
  
  m_term <- crossprod(U_Minv %*% Mm)
  
  a_star <- a0 + (n / 2)
  b_star <- as.numeric(b0 + 0.5 * (YtY + m0_term - m_term))
  
  # 5. Vectorized sampling[cite: 14]
  tau_samples <- rgamma(N, shape = a_star, rate = b_star)
  sigma2_samples <- 1 / tau_samples
  sigma_samples <- sqrt(sigma2_samples)
  
  # 6. Vectorized Conditional Draw for Beta[cite: 14]
  Z <- matrix(rnorm(p * N), nrow = p, ncol = N)
  W <- backsolve(U_Minv, Z)
  
  scaled_noise <- sweep(W, 2, sigma_samples, `*`)
  
  # Structural transpose to align to standard N x p format[cite: 14]
  beta_samples <- t(as.numeric(Mm) + scaled_noise)
  
  if (!is.null(colnames(X))) {
    colnames(beta_samples) <- colnames(X)
  }
  
  posterior_draws <- cbind(beta_samples, sigma = sigma_samples, sigma2 = sigma2_samples, tau = tau_samples)
  
  return(list(
    beta = beta_samples,      
    sigma2 = sigma2_samples,  
    draws = posterior_draws   
  ))
}
