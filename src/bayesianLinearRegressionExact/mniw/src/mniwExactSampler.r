# ==============================================================================
# src/mniwExactSampler.r
# Optimized Matrix Normal Inverse Wishart (MNIW) Sampling and Bayesian Regression
# Strictly loop-free and strictly avoids solve() and t()
# ==============================================================================

#' Optimized Helper: Sample from Inverse Wishart Distribution
sample_iw <- function(nu, S) {
  p <- ncol(S)
  U_S <- chol(S) 
  
  T_mat <- matrix(0, p, p)
  T_mat[upper.tri(T_mat)] <- rnorm(p * (p - 1) / 2)
  diag(T_mat) <- sqrt(rchisq(p, df = nu - 0:(p - 1)))
  
  X <- forwardsolve(T_mat, U_S, upper.tri = TRUE, transpose = TRUE)
  
  return(crossprod(X))
}

#' (i) Composition Sampling Method 1 (Covariance Parameterized)
sample_mniw_method1 <- function(M, V, nu, S) {
  k <- nrow(M)
  p <- ncol(M)
  
  Sigma <- sample_iw(nu, S)
  
  U_Sigma <- chol(Sigma)
  U_V <- chol(V)
  Z <- matrix(rnorm(k * p), nrow = k, ncol = p)
  
  B <- M + crossprod(U_V, Z %*% U_Sigma)
  
  return(list(B = B, Sigma = Sigma))
}

#' (ii) Composition Sampling Method 2 (Precision Parameterized)
sample_mniw_method2 <- function(M, Lambda, nu, S) {
  k <- nrow(M)
  p <- ncol(M)
  
  Sigma <- sample_iw(nu, S)
  
  U_Sigma <- chol(Sigma)
  U_Lambda <- chol(Lambda)
  Z <- matrix(rnorm(k * p), nrow = k, ncol = p)
  
  B <- M + backsolve(U_Lambda, Z %*% U_Sigma)
  
  return(list(B = B, Sigma = Sigma))
}

#' (iii) Bayesian Multivariate Linear Regression
#' @param Y n x p response matrix
#' @param X n x k design matrix
#' @param B0 k x p prior mean matrix
#' @param Lambda0 k x k prior precision matrix
#' @param nu0 prior degrees of freedom
#' @param S0 p x p prior scale matrix
#' @param method Sampling method to use (1 = Covariance, 2 = Precision)
#' @param n_samples Number of MNIW samples to draw
bayes_mvr <- function(Y, X, B0, Lambda0, nu0, S0, method = 2, n_samples = 1) {
  n <- nrow(X)
  k <- ncol(X)
  
  Lambda_n <- crossprod(X) + Lambda0
  U_Lambda_n <- chol(Lambda_n)
  
  XTY <- crossprod(X, Y)
  RHS <- XTY + Lambda0 %*% B0
  
  temp <- forwardsolve(U_Lambda_n, RHS, upper.tri = TRUE, transpose = TRUE)
  B_hat <- backsolve(U_Lambda_n, temp)
  
  Res <- Y - X %*% B_hat
  Diff <- B_hat - B0
  S_n <- S0 + crossprod(Res) + crossprod(Diff, Lambda0 %*% Diff)
  
  nu_n <- nu0 + n
  
  # Functional mapping to eliminate 'for' loops
  if (method == 1) {
    I_k <- diag(k)
    U_inv <- backsolve(U_Lambda_n, I_k)
    V_n <- tcrossprod(U_inv)
    
    samples <- lapply(seq_len(n_samples), function(x) {
      sample_mniw_method1(B_hat, V_n, nu_n, S_n)
    })
  } else if (method == 2) {
    samples <- lapply(seq_len(n_samples), function(x) {
      sample_mniw_method2(B_hat, Lambda_n, nu_n, S_n)
    })
  } else {
    stop("Method must be either 1 or 2.")
  }
  
  # Isolate mapped results optimally 
  B_samples <- lapply(samples, `[[`, "B")
  Sigma_samples <- lapply(samples, `[[`, "Sigma")
  
  return(list(
    post_params = list(B_hat = B_hat, Lambda_n = Lambda_n, nu_n = nu_n, S_n = S_n),
    samples = list(B = B_samples, Sigma = Sigma_samples)
  ))
}
