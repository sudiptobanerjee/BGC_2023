# ==============================================================================
# File: src/mniwSamplerDivideAndConquer.r
# Purpose: Exact Matrix-Normal-Inverse-Wishart (MNIW) Streaming Sampler
# Rules: No solve(), no t(), no general inversions, no for loops.
# ==============================================================================

#' Exact Sequential Divide & Conquer MNIW Sampler (Custom R Implementation)
#'
#' @param data_blocks List of data blocks, where each element is list(X = X_t, Y = Y_t)
#' @param m0 Prior location matrix (p x q)
#' @param M0_inv Prior row precision matrix M_0^{-1} (p x p)
#' @param nu0 Prior Inverse-Wishart degrees of freedom
#' @param Psi0 Prior Inverse-Wishart scale matrix (q x q)
#' @param n_samples Number of posterior draws to generate
#' @return List containing exact posterior parameters and MCMC-free joint draws
mniwSamplerDivideAndConquer <- function(data_blocks, m0, M0_inv, nu0, Psi0, n_samples = 1000) {
  p <- nrow(m0)
  q <- ncol(m0)
  
  init_state <- list(
    m = m0,
    M_inv = M0_inv,
    nu = nu0,
    Psi = Psi0
  )
  
  update_block <- function(state, block) {
    X_t <- block$X
    Y_t <- block$Y
    n_t <- nrow(X_t)
    
    U_prev <- chol(state$M_inv)
    w_prev <- backsolve(U_prev, state$m, transpose = TRUE)
    Q_part1 <- crossprod(Y_t) + crossprod(w_prev)
    
    M_inv_new <- state$M_inv + crossprod(X_t)
    m_new     <- state$m + crossprod(X_t, Y_t)
    nu_new    <- state$nu + n_t
    
    U_new <- chol(M_inv_new)
    w_new <- backsolve(U_new, m_new, transpose = TRUE)
    
    Q_t <- Q_part1 - crossprod(w_new)
    Psi_new <- state$Psi + Q_t
    
    list(
      m = m_new,
      M_inv = M_inv_new,
      nu = nu_new,
      Psi = Psi_new
    )
  }
  
  final_state <- Reduce(update_block, data_blocks, init = init_state)
  
  m_N     <- final_state$m
  M_N_inv <- final_state$M_inv
  nu_N    <- final_state$nu
  Psi_N   <- final_state$Psi
  
  U_M_N   <- chol(M_N_inv)
  U_Psi_N <- chol(Psi_N)
  
  w_N  <- backsolve(U_M_N, m_N, transpose = TRUE)
  mu_N <- backsolve(U_M_N, w_N)
  
  draw_single <- function(dummy) {
    W <- matrix(0, nrow = q, ncol = q)
    diag(W) <- sqrt(rchisq(q, df = nu_N - seq_len(q) + 1))
    
    if (q > 1) {
      W[lower.tri(W)] <- rnorm(q * (q - 1) / 2)
    }
    
    V <- backsolve(W, U_Psi_N, transpose = TRUE)
    Sigma_draw <- crossprod(V)
    
    U_Sigma <- chol(Sigma_draw)
    Z       <- matrix(rnorm(p * q), nrow = p, ncol = q)
    
    noise     <- backsolve(U_M_N, Z) %*% U_Sigma
    Beta_draw <- mu_N + noise
    
    list(Beta = Beta_draw, Sigma = Sigma_draw)
  }
  
  draws <- lapply(seq_len(n_samples), draw_single)
  
  Beta_samples  <- simplify2array(lapply(draws, `[[`, "Beta"))
  Sigma_samples <- simplify2array(lapply(draws, `[[`, "Sigma"))
  
  list(
    post_params = list(m = m_N, M_inv = M_N_inv, nu = nu_N, Psi = Psi_N),
    samples     = list(Beta = Beta_samples, Sigma = Sigma_samples)
  )
}

#' Exact Sequential Divide & Conquer MNIW Sampler (C++ Accelerated via mniw::rmniw)
#'
#' @param data_blocks List of data blocks, where each element is list(X = X_t, Y = Y_t)
#' @param m0 Prior location matrix (p x q)
#' @param M0_inv Prior row precision matrix M_0^{-1} (p x p)
#' @param nu0 Prior Inverse-Wishart degrees of freedom
#' @param Psi0 Prior Inverse-Wishart scale matrix (q x q)
#' @param n_samples Number of posterior draws to generate
#' @return List containing exact posterior parameters and MCMC-free joint draws
mniwSamplerDivideAndConquer_rmniw <- function(data_blocks, m0, M0_inv, nu0, Psi0, n_samples = 1000) {
  p <- nrow(m0)
  q <- ncol(m0)
  
  init_state <- list(
    m = m0,
    M_inv = M0_inv,
    nu = nu0,
    Psi = Psi0
  )
  
  update_block <- function(state, block) {
    X_t <- block$X
    Y_t <- block$Y
    n_t <- nrow(X_t)
    
    U_prev <- chol(state$M_inv)
    w_prev <- backsolve(U_prev, state$m, transpose = TRUE)
    Q_part1 <- crossprod(Y_t) + crossprod(w_prev)
    
    M_inv_new <- state$M_inv + crossprod(X_t)
    m_new     <- state$m + crossprod(X_t, Y_t)
    nu_new    <- state$nu + n_t
    
    U_new <- chol(M_inv_new)
    w_new <- backsolve(U_new, m_new, transpose = TRUE)
    
    Q_t <- Q_part1 - crossprod(w_new)
    Psi_new <- state$Psi + Q_t
    
    list(
      m = m_new,
      M_inv = M_inv_new,
      nu = nu_new,
      Psi = Psi_new
    )
  }
  
  final_state <- Reduce(update_block, data_blocks, init = init_state)
  
  m_N     <- final_state$m
  M_N_inv <- final_state$M_inv
  nu_N    <- final_state$nu
  Psi_N   <- final_state$Psi
  
  # Calculate location matrix mu_N = (M_N_inv)^(-1) * m_N via Cholesky decomposition
  U_M_N <- chol(M_N_inv)
  w_N   <- backsolve(U_M_N, m_N, transpose = TRUE)
  mu_N  <- backsolve(U_M_N, w_N)
  
  # Sample using C++ accelerated mniw::rmniw (Omega = precision matrix M_N_inv)
  draws <- mniw::rmniw(n = n_samples, Lambda = mu_N, Omega = M_N_inv, Psi = Psi_N, nu = nu_N)
  
  list(
    post_params = list(m = m_N, M_inv = M_N_inv, nu = nu_N, Psi = Psi_N),
    samples     = list(Beta = draws$X, Sigma = draws$V)
  )
}
