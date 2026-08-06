# src/flSampling.r

#' Compute local sufficient statistics on the client side
#' @param X The local design matrix
#' @param y The local response vector
compute_local_stats <- function(X, y) {
  n <- nrow(X)
  
  # crossprod(A) computes A^T A, crossprod(A, B) computes A^T B
  delta_M_inv <- crossprod(X)
  delta_m <- crossprod(X, y)
  delta_Q <- crossprod(y)[1, 1] # Extract scalar
  
  return(list(delta_M_inv = delta_M_inv, delta_m = delta_m, delta_Q = delta_Q, n = n))
}

#' Aggregate stats and sample from the global posterior on the server side
#' @param client_stats_list A list of summary statistics from all clients
#' @param M0_inv Prior precision matrix
#' @param m0 Prior natural mean vector
#' @param a0 Prior shape parameter for sigma^2
#' @param b0 Prior scale parameter for sigma^2
#' @param N_samples Number of posterior draws
flSampler <- function(client_stats_list, M0_inv, m0, a0, b0, N_samples) {
  p <- length(m0)
  
  # 1. Server Aggregation without for-loops using Reduce and lapply
  delta_M_inv_all <- lapply(client_stats_list, `[[`, "delta_M_inv")
  M_inv_global <- M0_inv + Reduce(`+`, delta_M_inv_all)
  
  delta_m_all <- lapply(client_stats_list, `[[`, "delta_m")
  m_global <- m0 + Reduce(`+`, delta_m_all)
  
  n_total <- sum(sapply(client_stats_list, `[[`, "n"))
  delta_Q_total <- sum(sapply(client_stats_list, `[[`, "delta_Q"))
  
  # 2. Compute structural matrices via Cholesky (U is upper triangular where U^T U = M_inv_global)
  U <- chol(M_inv_global)
  
  # U_inv is U^{-1}. Since U is upper triangular, U_inv is also upper triangular.
  U_inv <- backsolve(U, diag(p))
  
  # M_global = M_inv_global^{-1} = (U^T U)^{-1} = U^{-1} U^{-T}
  # tcrossprod(A) computes A A^T, so tcrossprod(U_inv) computes U^{-1} U^{-T}
  M_global <- tcrossprod(U_inv)
  
  # 3. Calculate the quadratic forms purely with triangular solvers
  # Handle a non-informative (zero matrix) prior gracefully
  if (sum(abs(M0_inv)) == 0) {
    q_prior <- 0
  } else {
    U0 <- chol(M0_inv)
    # We want m0^T M0 m0 = m0^T (U0^T U0)^{-1} m0 = || U0^{-T} m0 ||^2
    # backsolve with transpose=TRUE solves U0^T x = m0 for x
    x0 <- backsolve(U0, m0, transpose = TRUE)
    q_prior <- crossprod(x0)[1, 1]
  }
  
  # q_post = m_global^T M_global m_global = || U^{-T} m_global ||^2
  x_post <- backsolve(U, m_global, transpose = TRUE)
  q_post <- crossprod(x_post)[1, 1]
  
  # 4. Global Shape and Scale for Inverse-Gamma
  a_T <- a0 + n_total / 2
  b_T <- b0 + 0.5 * (delta_Q_total + q_prior - q_post)
  
  # 5. Global Posterior Sampling (Fully Vectorized)
  sigma2_draws <- 1 / rgamma(N_samples, shape = a_T, rate = b_T)
  
  # beta_mean = M_global * m_global. 
  # Since M_global is symmetric, crossprod(M, m) = M^T m = M m
  beta_mean <- crossprod(M_global, m_global) 
  
  # Draw standard normals in bulk (N_samples rows, p columns)
  Z <- matrix(rnorm(N_samples * p), nrow = N_samples, ncol = p)
  
  # Scale rows of Z by sqrt(sigma2_draws). 
  # Because matrices are column-major in R, vector multiplication natively cascades down rows.
  Z_scaled <- Z * sqrt(sigma2_draws)
  
  # We want beta = beta_mean^T + Z_scaled * U_inv^T. 
  # tcrossprod(A, B) computes A B^T. Thus tcrossprod(Z_scaled, U_inv) gives exactly what we need.
  # rep(..., each = N_samples) ensures beta_mean cascades perfectly into an N x p matrix.
  beta_draws <- tcrossprod(Z_scaled, U_inv) + rep(as.numeric(beta_mean), each = N_samples)
  
  draws <- cbind(beta_draws, sigma2_draws)
  return(list(draws = draws, a_T = a_T, b_T = b_T, M_T = M_global, m_T = m_global))
}
