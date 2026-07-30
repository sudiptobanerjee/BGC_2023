# src/rawSamplerWithQR.r

#' Exact Bayesian Linear Regression Sampler using QR Reparameterization
#'
#' @param y Vector (n x 1) of response values.
#' @param Q Matrix (n x p) with orthonormal columns from QR decomposition of X.
#' @param R Upper-triangular matrix (p x p) from QR decomposition of X.
#' @param m0 Vector (p x 1) prior parameter location vector.
#' @param M0_inv Matrix (p x p) prior precision matrix factor.
#' @param a0 Inverse-Gamma prior shape parameter.
#' @param b0 Inverse-Gamma prior scale parameter.
#' @param N Number of posterior draws requested.
#'
#' @return List containing beta draws (N x p), sigma2 draws (N x 1), and full draws matrix.
rawSamplerWithQR <- function(y, Q, R, m0, M0_inv, a0, b0, N) {
  
  n <- nrow(Q)
  p <- ncol(Q)
  
  # 1. Compute orthogonal cross-product Q^T y
  Qty <- crossprod(Q, y)
  
  # 2. Solve R^T m_theta = m0 using lower-triangular forward substitution
  if (all(m0 == 0)) {
    m_theta <- matrix(0, nrow = p, ncol = 1)
  } else {
    m_theta <- forwardsolve(t(R), m0)
  }
  
  # 3. Compute M_theta_inv = (R^T)^(-1) M0_inv R^(-1) via two triangular solves
  if (all(M0_inv == 0)) {
    M_theta_inv <- matrix(0, nrow = p, ncol = p)
  } else {
    A <- forwardsolve(t(R), M0_inv)
    M_theta_inv <- forwardsolve(t(R), t(A))
  }
  
  # 4. Posterior precision and Cholesky factor: R_theta_n^T R_theta_n = M_theta_n^(-1)
  M_theta_n_inv <- M_theta_inv + diag(p)
  R_theta_n <- chol(M_theta_n_inv)  # Upper-triangular Cholesky factor
  
  # 5. Updated location vector m_theta_n and mean mu_theta_n via triangular solves
  m_theta_n <- m_theta + Qty
  w <- forwardsolve(t(R_theta_n), m_theta_n)
  mu_theta_n <- backsolve(R_theta_n, w)
  
  # 6. Inverse-Gamma Parameter Updates
  YtY <- crossprod(y)
  w_term <- crossprod(w)
  
  # Safely compute prior quadratic term m0^T M0 m0
  if (all(m0 == 0) || all(M0_inv == 0)) {
    m0_term <- 0
  } else {
    U_M0inv <- chol(M0_inv)
    mu0 <- backsolve(U_M0inv, forwardsolve(t(U_M0inv), m0))
    m0_term <- crossprod(m0, mu0)
  }
  
  a_n <- a0 + (n / 2)
  b_n <- as.numeric(b0 + 0.5 * (YtY + m0_term - w_term))
  
  # 7. Vectorized marginal draw for precision/variance
  tau_samples <- rgamma(N, shape = a_n, rate = b_n)
  sigma2_samples <- 1 / tau_samples
  sigma_samples <- sqrt(sigma2_samples)
  
  # 8. Vectorized conditional draw for theta
  Z <- matrix(rnorm(p * N), nrow = p, ncol = N)
  V <- backsolve(R_theta_n, Z)
  scaled_noise_theta <- sweep(V, 2, sigma_samples, `*`)
  theta_samples <- as.numeric(mu_theta_n) + scaled_noise_theta
  
  # 9. Recover original regression parameters beta via backsolve(R, theta)
  beta_draws <- backsolve(R, theta_samples)
  beta_samples <- t(beta_draws)
  
  if (!is.null(colnames(R))) {
    colnames(beta_samples) <- colnames(R)
  } else if (!is.null(colnames(Q))) {
    colnames(beta_samples) <- colnames(Q)
  }
  
  posterior_draws <- cbind(beta_samples, sigma = sigma_samples, sigma2 = sigma2_samples, tau = tau_samples)
  
  return(list(
    beta = beta_samples,
    sigma2 = sigma2_samples,
    draws = posterior_draws
  ))
}

