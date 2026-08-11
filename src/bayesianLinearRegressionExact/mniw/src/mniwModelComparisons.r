# ==============================================================================
# Script: mniwModelComparisons.r
# Purpose: Vectorized Analytical and Sample-Based WAIC (Loop-Free, No solve/t)
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Analytical WAIC Computation (Vectorized)
# ------------------------------------------------------------------------------
compute_analytical_waic <- function(Y, X, Mn, mn, Psi_n, nu_n) {
  q_cols <- ncol(Y)
  I_q <- diag(q_cols)
  
  # 1. Cholesky Factorization of Psi_n (Psi_n = U_Psi^T %*% U_Psi)
  U_Psi <- chol(Psi_n)
  log_det_Psi <- 2 * sum(log(diag(U_Psi)))
  
  # Upper triangular inverse U_Psi^{-1} via backsolve
  U_Psi_inv <- backsolve(U_Psi, I_q)
  
  # 2. Vectorized Predictive Variance Scalars (h_i)
  h <- rowSums((X %*% Mn) * X)
  c_vec <- 1 + h
  
  # 3. Vectorized Residual Matrix E (N x q)
  E <- Y - X %*% mn
  
  # 4. Vectorized Mahalanobis Distances (d_i)
  # e_i^T %*% Psi_n^{-1} %*% e_i = ||e_i^T %*% U_Psi^{-1}||^2
  Z_E <- E %*% U_Psi_inv
  d <- rowSums(Z_E * Z_E)
  
  # --- Vectorized Pointwise lppd ---
  term1 <- lgamma((nu_n + 1) / 2) - lgamma((nu_n - q_cols + 1) / 2)
  term2 <- -(q_cols / 2) * log(pi)
  term3 <- -0.5 * (q_cols * log(c_vec) + log_det_Psi)
  term4 <- -((nu_n + 1) / 2) * log(1 + d / c_vec)
  
  lppd_i <- term1 + term2 + term3 + term4
  lppd_sum <- sum(lppd_i)
  
  # --- Vectorized Pointwise p_waic ---
  var_A <- sum(trigamma((nu_n - 1:q_cols + 1) / 2))
  var_B <- 2 * (h^2) * q_cols + 4 * h * nu_n * d + 2 * nu_n * (d^2)
  cov_A_B <- -2 * d
  
  p_waic_i <- 0.25 * (var_A + var_B + 2 * cov_A_B)
  p_waic_sum <- sum(p_waic_i)
  
  waic <- -2 * lppd_sum + 2 * p_waic_sum
  
  return(list(
    WAIC = waic, 
    lppd = lppd_sum, 
    p_waic = p_waic_sum
  ))
}

# ------------------------------------------------------------------------------
# 2. Sample-Based WAIC Computation (Functional & LOO-Aligned)
# ------------------------------------------------------------------------------
compute_sample_waic <- function(Y, X, B_samples, Sigma_samples) {
  q_cols <- ncol(Y)
  S_samples <- dim(B_samples)[3]
  I_q <- diag(q_cols)
  
  # Vectorized calculation per sample via functional mapping (lapply)
  log_lik_list <- lapply(seq_len(S_samples), function(s) {
    beta_s <- B_samples[,,s]
    sigma_s <- Sigma_samples[,,s]
    
    # Cholesky Factorization of Sigma_s
    U_sigma <- chol(sigma_s)
    log_det_s <- 2 * sum(log(diag(U_sigma)))
    
    # Upper triangular inverse U_sigma^{-1} via backsolve
    U_sigma_inv <- backsolve(U_sigma, I_q)
    
    diff <- Y - X %*% beta_s
    Z_diff <- diff %*% U_sigma_inv
    quad_form <- rowSums(Z_diff * Z_diff)
    
    -0.5 * q_cols * log(2 * pi) - 0.5 * log_det_s - 0.5 * quad_form
  })
  
  # Bind into an S x N log-likelihood matrix
  log_lik_mat <- do.call(rbind, log_lik_list)
  
  # --- LOO-Aligned Aggregation ---
  
  # 1. Pointwise lppd (Vectorized Log-Sum-Exp across columns)
  max_ll <- apply(log_lik_mat, 2, max)
  lppd_i <- max_ll + log(colMeans(exp(sweep(log_lik_mat, 2, max_ll, "-"))))
  lppd_sum <- sum(lppd_i)
  
  # 2. Pointwise p_waic (Column-wise sample variance)
  p_waic_i <- apply(log_lik_mat, 2, stats::var)
  p_waic_sum <- sum(p_waic_i)
  
  waic <- -2 * lppd_sum + 2 * p_waic_sum
  
  return(list(
    WAIC = waic, 
    lppd = lppd_sum, 
    p_waic = p_waic_sum
  ))
}
