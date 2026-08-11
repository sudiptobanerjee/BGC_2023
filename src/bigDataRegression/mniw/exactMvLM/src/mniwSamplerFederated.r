# ==============================================================================
# File: src/mniwSamplerFederated.r
# Purpose: Exact Federated Learning MNIW Sampler (Server-Client Architecture)
# Rules: Strictly NO solve(), NO t(), NO general inversions, NO for loops.
# ==============================================================================

#' Exact Federated MNIW Sampler (Custom R Implementation)
#'
#' @param client_data_list List of data blocks from T clients, each containing list(X = X_t, Y = Y_t)
#' @param m0 Prior location matrix (p x q)
#' @param M0_inv Prior row precision matrix M_0^{-1} (p x p)
#' @param nu0 Prior Inverse-Wishart degrees of freedom
#' @param Psi0 Prior Inverse-Wishart scale matrix (q x q)
#' @param n_samples Number of joint posterior draws to generate
#' @return List containing exact global posterior parameters and MCMC-free joint draws
mniwSamplerFederated <- function(client_data_list, m0, M0_inv, nu0, Psi0, n_samples = 1000) {
  p <- nrow(m0)
  q <- ncol(m0)
  
  # ----------------------------------------------------------------------------
  # 1. Local Client Computations (Parallelizable across T clients)
  # ----------------------------------------------------------------------------
  compute_local_stats <- function(client_data) {
    X_t <- client_data$X
    Y_t <- client_data$Y
    list(
      delta_M_inv = crossprod(X_t),
      delta_m     = crossprod(X_t, Y_t),
      delta_Q     = crossprod(Y_t),
      n_t         = nrow(X_t)
    )
  }
  
  client_stats <- lapply(client_data_list, compute_local_stats)
  
  # ----------------------------------------------------------------------------
  # 2. Server One-Shot Aggregation
  # ----------------------------------------------------------------------------
  sum_M_inv <- Reduce("+", lapply(client_stats, `[[`, "delta_M_inv"))
  sum_m     <- Reduce("+", lapply(client_stats, `[[`, "delta_m"))
  sum_Q     <- Reduce("+", lapply(client_stats, `[[`, "delta_Q"))
  sum_n     <- sum(sapply(client_stats, `[[`, "n_t"))
  
  # Update precision matrix, location sum, and degrees of freedom
  M_T_inv <- M0_inv + sum_M_inv
  m_T     <- m0 + sum_m
  nu_T    <- nu0 + sum_n
  
  # Compute w_0 and w_T via Cholesky factorizations (avoids solve and t)
  U0 <- chol(M0_inv)
  w0 <- backsolve(U0, m0, transpose = TRUE)
  
  UT <- chol(M_T_inv)
  wT <- backsolve(UT, m_T, transpose = TRUE)
  
  # Global mean matrix mu_T = M_T * m_T
  mu_T <- backsolve(UT, wT)
  
  # One-shot telescopic global scale matrix calculation
  Psi_T <- Psi0 + sum_Q + crossprod(w0) - crossprod(wT)
  
  # ----------------------------------------------------------------------------
  # 3. Server Joint Posterior Sampling (Custom R Composition)
  # ----------------------------------------------------------------------------
  U_Psi_T <- chol(Psi_T)
  
  draw_single <- function(dummy) {
    # Sample Sigma ~ IW(nu_T, Psi_T) via Bartlett Decomposition
    W <- matrix(0, nrow = q, ncol = q)
    diag(W) <- sqrt(rchisq(q, df = nu_T - seq_len(q) + 1))
    
    if (q > 1) {
      W[lower.tri(W)] <- rnorm(q * (q - 1) / 2)
    }
    
    V <- backsolve(W, U_Psi_T, transpose = TRUE)
    Sigma_draw <- crossprod(V)
    
    # Sample Beta ~ MN(mu_T, M_T, Sigma_draw)
    U_Sigma <- chol(Sigma_draw)
    Z       <- matrix(rnorm(p * q), nrow = p, ncol = q)
    
    noise     <- backsolve(UT, Z) %*% U_Sigma
    Beta_draw <- mu_T + noise
    
    list(Beta = Beta_draw, Sigma = Sigma_draw)
  }
  
  draws <- lapply(seq_len(n_samples), draw_single)
  
  Beta_samples  <- simplify2array(lapply(draws, `[[`, "Beta"))
  Sigma_samples <- simplify2array(lapply(draws, `[[`, "Sigma"))
  
  list(
    post_params = list(m = m_T, M_inv = M_T_inv, nu = nu_T, Psi = Psi_T),
    samples     = list(Beta = Beta_samples, Sigma = Sigma_samples)
  )
}

#' Exact Federated MNIW Sampler (C++ Accelerated via mniw::rmniw)
#'
#' @param client_data_list List of data blocks from T clients, each containing list(X = X_t, Y = Y_t)
#' @param m0 Prior location matrix (p x q)
#' @param M0_inv Prior row precision matrix M_0^{-1} (p x p)
#' @param nu0 Prior Inverse-Wishart degrees of freedom
#' @param Psi0 Prior Inverse-Wishart scale matrix (q x q)
#' @param n_samples Number of joint posterior draws to generate
#' @return List containing exact global posterior parameters and MCMC-free joint draws
mniwSamplerFederated_rmniw <- function(client_data_list, m0, M0_inv, nu0, Psi0, n_samples = 1000) {
  p <- nrow(m0)
  q <- ncol(m0)
  
  # ----------------------------------------------------------------------------
  # 1. Local Client Computations
  # ----------------------------------------------------------------------------
  compute_local_stats <- function(client_data) {
    X_t <- client_data$X
    Y_t <- client_data$Y
    list(
      delta_M_inv = crossprod(X_t),
      delta_m     = crossprod(X_t, Y_t),
      delta_Q     = crossprod(Y_t),
      n_t         = nrow(X_t)
    )
  }
  
  client_stats <- lapply(client_data_list, compute_local_stats)
  
  # ----------------------------------------------------------------------------
  # 2. Server One-Shot Aggregation
  # ----------------------------------------------------------------------------
  sum_M_inv <- Reduce("+", lapply(client_stats, `[[`, "delta_M_inv"))
  sum_m     <- Reduce("+", lapply(client_stats, `[[`, "delta_m"))
  sum_Q     <- Reduce("+", lapply(client_stats, `[[`, "delta_Q"))
  sum_n     <- sum(sapply(client_stats, `[[`, "n_t"))
  
  M_T_inv <- M0_inv + sum_M_inv
  m_T     <- m0 + sum_m
  nu_T    <- nu0 + sum_n
  
  U0 <- chol(M0_inv)
  w0 <- backsolve(U0, m0, transpose = TRUE)
  
  UT <- chol(M_T_inv)
  wT <- backsolve(UT, m_T, transpose = TRUE)
  
  mu_T  <- backsolve(UT, wT)
  Psi_T <- Psi0 + sum_Q + crossprod(w0) - crossprod(wT)
  
  # ----------------------------------------------------------------------------
  # 3. Server Joint Posterior Sampling via C++ mniw::rmniw
  # Note: Omega = Row Precision Matrix M_T_inv; draws$V = Inverse-Wishart draws
  # ----------------------------------------------------------------------------
  draws <- mniw::rmniw(n = n_samples, Lambda = mu_T, Omega = M_T_inv, Psi = Psi_T, nu = nu_T)
  
  list(
    post_params = list(m = m_T, M_inv = M_T_inv, nu = nu_T, Psi = Psi_T),
    samples     = list(Beta = draws$X, Sigma = draws$V)
  )
}
