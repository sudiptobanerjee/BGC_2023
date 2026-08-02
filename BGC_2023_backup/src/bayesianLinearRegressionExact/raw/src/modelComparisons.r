# src/modelComparisons.r

#' Model Assessment and Comparisons (WAIC)
#' Functions for Exact Analytical WAIC, Raw Sample-Based WAIC, and 'loo' Package Wrapper[cite: 12]

#' 1. Compute Pointwise Log-Likelihood Matrix[cite: 12]
#'
#' Computes the (n x N) log-likelihood matrix LL where LL[i, s] = log p(y_i | beta^(s), sigma2^(s)).[cite: 12]
#' Computationally optimized to avoid explicit transpositions and loop constructs.[cite: 12]
#'
#' @param X Matrix (n x p) of covariates.
#' @param y Vector (n x 1) of observed target values.
#' @param beta_draws Matrix (N x p) of posterior draws for regression parameters.
#' @param sigma2_draws Vector/Matrix (N x 1) of posterior draws for error variance sigma^2.
#'
#' @return Matrix (n x N) of pointwise log-likelihood values.
compute_log_likelihood_matrix <- function(X, y, beta_draws, sigma2_draws) {
  if (!is.matrix(X)) X <- as.matrix(X)
  if (!is.matrix(beta_draws)) beta_draws <- as.matrix(beta_draws)
  
  n <- nrow(X)
  sigma2_vec <- as.vector(sigma2_draws)
  N <- length(sigma2_vec)
  
  # Fitted values mu: (n x N) matrix via tcrossprod (X %*% t(beta_draws))[cite: 12]
  mu_matrix <- tcrossprod(X, beta_draws)
  
  # Squared residuals: (n x N)[cite: 12]
  resids_sq <- (y - mu_matrix)^2
  
  # Vectorized scaling across memory-contiguous column-major order[cite: 12]
  sigma2_rep <- rep(sigma2_vec, each = n)
  
  # Pointwise Log-Likelihood matrix (n x N) including the -0.5*log(2*pi) constant[cite: 12]
  LL <- -0.5 * log(2 * pi * sigma2_rep) - (resids_sq / (2 * sigma2_rep))
  
  dim(LL) <- c(n, N)
  return(LL)
}


#' 2. Exact Analytical WAIC for Conjugate Linear Model
#'
#' Computes closed-form lppd_exact, p_WAIC_exact, and WAIC_exact without MCMC error.[cite: 12]
#' Refactored to completely eliminate solve() using Cholesky factorizations.
#'
#' @param X Matrix (n x p) of covariates.
#' @param y Vector (n x 1) of target values.
#' @param m0 Prior mean vector (p x 1).
#' @param M0_inv Prior precision matrix (p x p).
#' @param a0 Prior Inverse-Gamma shape parameter.
#' @param b0 Prior Inverse-Gamma scale parameter.
#'
#' @return List containing lppd, p_waic, waic, and pointwise vectors.
waic_exact <- function(X, y, m0, M0_inv, a0, b0) {
  if (!is.matrix(X)) X <- as.matrix(X)
  if (!is.matrix(y)) y <- as.matrix(y)
  
  n <- nrow(X)
  p <- ncol(X)
  
  # A. Joint Posterior Parameters via Cholesky
  M_inv <- M0_inv + crossprod(X)
  U_Minv <- chol(M_inv)
  
  rhs <- M0_inv %*% m0 + crossprod(X, y)
  
  # Solve M_inv * mu_beta = rhs using Cholesky factor
  mu_beta <- backsolve(U_Minv, backsolve(U_Minv, rhs, transpose = TRUE))
  
  a_star  <- a0 + n / 2
  
  quad_y  <- sum(y^2)
  quad_m0 <- as.numeric(crossprod(m0, M0_inv %*% m0))
  
  # quad_mu: mu_beta^T M_inv mu_beta mapped through Cholesky factor
  U_mu <- U_Minv %*% mu_beta
  quad_mu <- as.numeric(crossprod(U_mu))
  
  b_star  <- b0 + 0.5 * (quad_y + quad_m0 - quad_mu)
  
  # B. Leverage and Residuals
  # Leverage k_i = diag(X M_inv^{-1} X^T) via triangular solve
  W_T <- backsolve(U_Minv, t(X), transpose = TRUE)
  k <- colSums(W_T^2)
  
  h <- 1 + k
  e <- as.vector(y - X %*% mu_beta)
  
  # C. Exact lppd via Marginal Student's t Distribution[cite: 12]
  lppd_i <- lgamma(a_star + 0.5) - lgamma(a_star) - 
            0.5 * log(2 * pi * b_star * h) - 
            (a_star + 0.5) * log(1 + (e^2) / (2 * b_star * h))
  lppd_exact <- sum(lppd_i)
  
  # D. Exact p_WAIC via Law of Total Variance (5-term formula)[cite: 12]
  term1 <- 0.5 * (k^2) + (a_star * k * (e^2)) / b_star
  term2 <- 0.25 * trigamma(a_star) + (a_star * (e^4)) / (4 * (b_star^2)) - (e^2) / (2 * b_star)
  
  p_waic_i <- term1 + term2
  p_waic_exact <- sum(p_waic_i)
  
  # E. Exact WAIC[cite: 12]
  waic_exact_val <- -2 * (lppd_exact - p_waic_exact)
  
  return(list(
    lppd   = lppd_exact,
    p_waic = p_waic_exact,
    waic   = waic_exact_val,
    lppd_i = lppd_i,
    p_waic_i = p_waic_i
  ))
}


#' 3. Sample-Based WAIC via Raw R Code (with Log-Sum-Exp Trick)[cite: 12]
#'
#' Evaluates lppd, p_WAIC, and WAIC empirically from posterior draws using 
#' the Log-Sum-Exp identity to ensure numerical stability against overflow/underflow.[cite: 12]
#'
#' @param LL Matrix (n x N) of pointwise log-likelihoods.
#'
#' @return List containing lppd, p_waic, waic, and pointwise vectors.
waic_sample <- function(LL) {
  n <- nrow(LL)
  N <- ncol(LL)
  
  # A. Numerically Stable lppd via Log-Sum-Exp Trick[cite: 12]
  c_i <- apply(LL, 1, max)
  LL_shifted <- LL - c_i  # LL - c_i is <= 0, preventing exp() overflow/underflow[cite: 12]
  
  lppd_i <- c_i - log(N) + log(rowSums(exp(LL_shifted)))
  lppd_sample <- sum(lppd_i)
  
  # B. Sample Variance Penalty p_WAIC[cite: 12]
  p_waic_i <- apply(LL, 1, var) # Sample variance across N draws for each observation i[cite: 12]
  p_waic_sample <- sum(p_waic_i)
  
  # C. WAIC on Deviance Scale[cite: 12]
  waic_sample_val <- -2 * (lppd_sample - p_waic_sample)
  
  return(list(
    lppd   = lppd_sample,
    p_waic = p_waic_sample,
    waic   = waic_sample_val,
    lppd_i = lppd_i,
    p_waic_i = p_waic_i
  ))
}


#' 4. Wrapper for WAIC via the 'loo' Package[cite: 12]
#'
#' Wrapper function that formats the LL matrix for the external 'loo' package.[cite: 12]
#'
#' @param LL Matrix (n x N) of pointwise log-likelihoods.
#'
#' @return List containing lppd, p_waic, waic, and the raw loo object.
waic_loo_package <- function(LL) {
  if (!requireNamespace("loo", quietly = TRUE)) {
    stop("Package 'loo' is required. Please install it with install.packages('loo').")
  }
  
  # loo::waic expects matrix of size (N x n) [draws x observations][cite: 12]
  loo_obj <- loo::waic(t(LL))
  
  elpd_waic <- loo_obj$estimates["elpd_waic", "Estimate"]
  p_waic    <- loo_obj$estimates["p_waic", "Estimate"]
  waic_val  <- loo_obj$estimates["waic", "Estimate"]
  
  # Reconstruct lppd from elpd_waic = lppd - p_waic[cite: 12]
  lppd_val  <- elpd_waic + p_waic
  
  return(list(
    lppd       = lppd_val,
    p_waic     = p_waic,
    waic       = waic_val,
    loo_object = loo_obj
  ))
}
