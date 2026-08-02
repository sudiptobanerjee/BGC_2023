#' Exact Posterior Predictive Sampling for Linear Regression (Vy = I_m)
#'
#' Computationally optimized routine using BLAS tcrossprod() to eliminate 
#' explicit matrix transpositions t() and memory-heavy allocations.
#'
#' @param X_new Matrix (m x p) of predictor variables for new targets.
#' @param beta_draws Matrix (N x p) of posterior draws for regression parameters.
#' @param sigma2_draws Vector (N x 1) of posterior draws for error variance sigma^2.
#'
#' @return Matrix (m x N) where column i is predictive draw vector y_tilde_(i).
rawPosteriorPredictive <- function(X_new, beta_draws, sigma2_draws) {
  
  # 1. Type verification without unnecessary deep copies
  if (!is.matrix(X_new)) X_new <- as.matrix(X_new)
  if (!is.matrix(beta_draws)) beta_draws <- as.matrix(beta_draws)
  
  sigma <- sqrt(as.vector(sigma2_draws))
  m <- nrow(X_new)
  N <- length(sigma)
  
  if (ncol(X_new) != ncol(beta_draws)) {
    stop("Dimension mismatch: ncol(X_new) must equal ncol(beta_draws).")
  }
  
  # 2. Fast Conditional Mean Computation: Mu = X_new %*% t(beta_draws)
  # tcrossprod(A, B) computes A %*% B^T directly via BLAS without creating t(B)
  Mu <- tcrossprod(X_new, beta_draws)
  
  # 3. Vectorized Noise Generation & Scaling
  # Multiply column j of noise matrix directly by scalar sigma[j] 
  # using contiguous column-major memory recycling
  Y_tilde <- Mu + matrix(rnorm(m * N), nrow = m, ncol = N) * rep(sigma, each = m)
  
  return(Y_tilde)
}
