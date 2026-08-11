# ==============================================================================
# Script: mniwPosteriorPredictive.r
# Purpose: Exact Posterior Sampling & Predictive Simulations (No solve/t/loops)
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Package-Based MNIW Posterior Sampler
# ------------------------------------------------------------------------------
draw_mniw_samples_pkg <- function(S, mn, Mn_inv, Psi_n, nu_n) {
  # mniw::rmniw accepts Omega = Mn_inv directly
  res <- mniw::rmniw(n = S, Lambda = mn, Omega = Mn_inv, Psi = Psi_n, nu = nu_n)
  return(list(X = res$X, Sigma = res$V))
}

# ------------------------------------------------------------------------------
# 2. Custom Exact MNIW Posterior Sampler
# ------------------------------------------------------------------------------
draw_mniw_samples_custom <- function(S, mn, Mn, Psi_n, nu_n) {
  p_cols <- nrow(mn)
  q_cols <- ncol(mn)
  I_q <- diag(q_cols)
  
  # Precompute Cholesky factors (Upper Triangular)
  U_Mn <- chol(Mn)
  U_Psi <- chol(Psi_n)
  
  # 1. Sample Sigma_draws ~ InvWishart(nu_n, Psi_n)
  Sigma_list <- lapply(seq_len(S), function(s) {
    # Generate upper triangular Bartlett factor matrix W
    W <- matrix(0, nrow = q_cols, ncol = q_cols)
    diag(W) <- sqrt(rchisq(q_cols, df = nu_n - (1:q_cols) + 1))
    
    if (q_cols > 1) {
      lower_indices <- lower.tri(W)
      W[upper.tri(W)] <- rnorm(sum(lower_indices))
    }
    
    # Correct Bartlett transformation: V = W^{-T} U_Psi
    V <- backsolve(W, U_Psi, transpose = TRUE)
    
    # Sigma = V^T V
    crossprod(V)
  })
  Sigma_draws <- array(unlist(Sigma_list), dim = c(q_cols, q_cols, S))
  
  # 2. Sample beta_draws ~ MatrixNormal(mn, Mn, Sigma_s) given drawn Sigma_s
  beta_list <- lapply(seq_len(S), function(s) {
    Sigma_s <- Sigma_draws[,,s]
    U_Sigma_s <- chol(Sigma_s)
    
    # Generate Matrix Normal noise: U_Mn^T %*% Z %*% U_Sigma_s
    Z <- matrix(rnorm(p_cols * q_cols), nrow = p_cols, ncol = q_cols)
    Noise <- crossprod(U_Mn, Z) %*% U_Sigma_s
    
    mn + Noise
  })
  beta_draws <- array(unlist(beta_list), dim = c(p_cols, q_cols, S))
  
  return(list(X = beta_draws, Sigma = Sigma_draws))
}

# ------------------------------------------------------------------------------
# 3. Custom Posterior Predictive Generator
# ------------------------------------------------------------------------------
mniw_posterior_predictive <- function(X_new, beta_draws, Sigma_draws) {
  N_new <- nrow(X_new)
  q_cols <- ncol(Sigma_draws[,,1])
  S_samples <- dim(beta_draws)[3]
  
  # Functional mapping across posterior draws S
  Y_draws_list <- lapply(seq_len(S_samples), function(s) {
    beta_s <- beta_draws[,,s]
    Sigma_s <- Sigma_draws[,,s]
    
    # Upper triangular Cholesky factor U_Sigma (U^T %*% U = Sigma_s)
    U_Sigma_s <- chol(Sigma_s)
    
    # Conditional Mean Matrix (N_new x q)
    Mu_s <- X_new %*% beta_s
    
    # Correlated Normal Noise: Z %*% U_Sigma_s yields covariance U^T U = Sigma_s
    Z <- matrix(rnorm(N_new * q_cols), nrow = N_new, ncol = q_cols)
    Noise_s <- Z %*% U_Sigma_s
    
    Mu_s + Noise_s
  })
  
  # Stack list of (N_new x q) draws into an array of dimension (N_new x q x S)
  array(unlist(Y_draws_list), dim = c(N_new, q_cols, S_samples))
}

# ------------------------------------------------------------------------------
# 4. Package-Based Posterior Predictive Generator (using mniw::rMNorm)
# ------------------------------------------------------------------------------
mniw_posterior_predictive_pkg <- function(X_new, beta_draws, Sigma_draws) {
  N_new <- nrow(X_new)
  q_cols <- ncol(Sigma_draws[,,1])
  S_samples <- dim(beta_draws)[3]
  I_N <- diag(N_new)
  
  # Generate Matrix-Normal predictive draws via mniw::rMNorm
  Y_draws_list <- lapply(seq_len(S_samples), function(s) {
    mniw::rMNorm(
      n = 1,
      Lambda = X_new %*% beta_draws[,,s],
      SigmaR = I_N,
      SigmaC = Sigma_draws[,,s]
    )
  })
  
  array(unlist(Y_draws_list), dim = c(N_new, q_cols, S_samples))
}
