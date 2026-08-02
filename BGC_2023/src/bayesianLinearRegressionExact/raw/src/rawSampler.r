# src/rawSampler.r

rawSampler <- function(y, X, m0, M0_inv, a0, b0, N) {
  
  n <- nrow(X)
  p <- ncol(X)
  
  # 1. Compute Cross-Products (Matches XtX and XtY)[cite: 1]
  XtX <- crossprod(X)
  XtY <- crossprod(X, y)
  
  # 2. Posterior Precision and Cholesky Factorization
  M_inv <- M0_inv + XtX
  U_Minv <- chol(M_inv) 
  
  # 3. Compute Posterior Mean (Mm)
  m <- m0 + XtY
  
  # EXACT MATCH to mu_n <- backsolve(U_n, forwardsolve(t(U_n), b))[cite: 1]
  Mm <- backsolve(U_Minv, forwardsolve(t(U_Minv), m))
  
  # 4. Exact Marginal Draw for Precision/Variance
  YtY <- crossprod(y)
  
  # Safely compute the prior term (m0^T M0 m0) without matrix inversion
  if (all(m0 == 0)) {
    m0_term <- 0
  } else {
    # Solve M0_inv * mu0 = m0 to safely extract the implied prior mean
    U_M0inv <- chol(M0_inv)
    mu0 <- backsolve(U_M0inv, forwardsolve(t(U_M0inv), m0))
    # m0^T M0 m0 is strictly equivalent to m0^T mu0
    m0_term <- crossprod(m0, mu0)
  }
  
  # EXACT MATCH to mun_term <- crossprod(U_n %*% mu_n)[cite: 1]
  m_term <- crossprod(U_Minv %*% Mm)
  
  a_star <- a0 + (n / 2)
  b_star <- as.numeric(b0 + 0.5 * (YtY + m0_term - m_term))
  
  # 5. Vectorized sampling[cite: 1]
  tau_samples <- rgamma(N, shape = a_star, rate = b_star)
  sigma2_samples <- 1 / tau_samples
  sigma_samples <- sqrt(sigma2_samples)
  
  # 6. Vectorized Conditional Draw for Beta[cite: 1]
  Z <- matrix(rnorm(p * N), nrow = p, ncol = N)
  W <- backsolve(U_Minv, Z)
  
  scaled_noise <- sweep(W, 2, sigma_samples, `*`)
  beta_samples <- t(as.numeric(Mm) + scaled_noise)
  
  if (!is.null(colnames(X))) {
    colnames(beta_samples) <- colnames(X)
  }
  
  # Bind output to perfectly match the reference sampler format[cite: 1]
  posterior_draws <- cbind(beta_samples, sigma = sigma_samples, sigma2 = sigma2_samples, tau = tau_samples)
  
  return(list(
    beta = beta_samples,      
    sigma2 = sigma2_samples,  
    draws = posterior_draws   
  ))
}

