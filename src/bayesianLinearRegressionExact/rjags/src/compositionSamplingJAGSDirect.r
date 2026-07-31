library(rjags)

run_direct_jags_composition <- function(X, Y = NULL, n_samples = 5000, a0 = 0.01, b0 = 0.01, 
                                        m0 = NULL, M0 = NULL, mu0 = NULL, V0 = NULL, 
                                        y = NULL, M_samples = NULL) {
  # Handle argument aliases
  if (is.null(Y) && !is.null(y)) Y <- y
  if (!is.null(M_samples)) n_samples <- M_samples
  
  N <- nrow(X)
  P <- ncol(X)
  
  # Prior Covariance Factor M0 & Cholesky Decomposition
  if (is.null(M0)) {
    if (!is.null(V0)) M0 <- V0 else M0 <- diag(100, P)
  }
  R0     <- chol(M0)
  M0_inv <- chol2inv(R0) # Prior precision factor via Cholesky inverse
  
  # Precision-weighted prior mean m0 and prior mean vector mu0_vec
  if (is.null(m0)) {
    if (!is.null(mu0)) {
      # Solve M0 %*% m0 = mu0 for m0 using triangular solvers
      z0 <- forwardsolve(R0, mu0, upper.tri = TRUE, transpose = TRUE)
      m0 <- as.vector(backsolve(R0, z0))
    } else {
      m0 <- rep(0, P)
    }
  }
  # mu0_vec = M0 %*% m0 = R0^T %*% (R0 %*% m0)
  mu0_vec <- as.vector(crossprod(R0, R0 %*% m0))
  
  # Posterior Parameter Updates using crossprod and Cholesky
  M_inv <- M0_inv + crossprod(X) # Posterior precision factor
  RN    <- chol(M_inv)           # Cholesky factor of posterior precision
  
  m  <- as.vector(m0 + crossprod(X, Y)) # Precision-weighted posterior mean
  
  # Solve M_inv %*% mu = m for mu using triangular solvers (RN^T RN mu = m)
  zN <- forwardsolve(RN, m, upper.tri = TRUE, transpose = TRUE)
  mu <- as.vector(backsolve(RN, zN))   # Posterior mean vector
  
  a_N <- a0 + N / 2
  b_N <- as.numeric(b0 + 0.5 * (crossprod(Y) + crossprod(mu0_vec, m0) - crossprod(mu, m)))
  
  # JAGS direct sampler
  model_string <- "
  model {
    for (k in 1:N_SAMPLES) {
      tau[k] ~ dgamma(a_N, b_N)
      sigma2[k] <- 1 / tau[k]
      sigma[k] <- sqrt(sigma2[k])
      beta[k, 1:P] ~ dmnorm(mu[1:P], tau[k] * M_inv[1:P, 1:P])
    }
  }
  "
  
  jags_data <- list(
    P = P,
    N_SAMPLES = n_samples,
    mu = mu,
    M_inv = M_inv,
    a_N = a_N,
    b_N = b_N
  )
  
  model <- jags.model(textConnection(model_string), data = jags_data, n.chains = 1, quiet = TRUE)
  samples <- jags.samples(model, variable.names = c("beta", "sigma2", "tau", "sigma"), n.iter = 1)
  
  # Reshape 4D mcarray directly into an (n_samples x P) matrix without t()
  beta_draws <- matrix(samples$beta, nrow = n_samples, ncol = P)
  if (is.null(colnames(X))) {
    colnames(beta_draws) <- paste0("beta[", 1:P, "]")
  } else {
    colnames(beta_draws) <- colnames(X)
  }
  
  sigma2_draws <- as.vector(samples$sigma2)
  tau_draws    <- as.vector(samples$tau)
  sigma_draws  <- as.vector(samples$sigma)
  
  draws <- cbind(beta_draws, sigma2 = sigma2_draws, tau = tau_draws, sigma = sigma_draws)
  
  return(list(draws = draws))
}
