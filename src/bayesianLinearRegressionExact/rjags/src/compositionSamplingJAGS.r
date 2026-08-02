# src/compositionSamplingJAGS.r

library(rjags)

run_replicated_jags_composition <- function(X, Y, 
                                            M_samples = 5000, 
                                            a0 = 0.01, 
                                            b0 = 0.01, 
                                            m0 = NULL, 
                                            M0_inv = NULL,
                                            mu0 = NULL) {
  N <- nrow(X)
  P <- ncol(X)
  
  # -------------------------------------------------------------------
  # Setup Prior Precision Factor M0_inv and Mean Vectors (m0, mu0)[cite: 15]
  # -------------------------------------------------------------------
  if (is.null(M0_inv)) M0_inv <- diag(0.01, P)
  
  if (is.null(m0) && is.null(mu0)) {
    m0  <- rep(0, P)
    mu0 <- rep(0, P)
  } else if (!is.null(mu0) && is.null(m0)) {
    m0 <- as.vector(M0_inv %*% mu0)
  } else if (!is.null(m0) && is.null(mu0)) {
    R0_prec <- chol(M0_inv)
    # Optimized triangular solve using transpose = TRUE
    z0      <- backsolve(R0_prec, m0, transpose = TRUE)
    mu0     <- as.vector(backsolve(R0_prec, z0))
  }
  
  # -------------------------------------------------------------------
  # 1. Compute Analytical Posterior Parameters in R[cite: 15]
  # -------------------------------------------------------------------
  M_inv <- M0_inv + crossprod(X)                  # Posterior precision factor M^-1[cite: 15]
  RN    <- chol(M_inv)                             # Cholesky factor of M^-1[cite: 15]
  
  m_vec <- m0 + crossprod(X, Y)                   # Precision-weighted posterior mean m[cite: 15]
  
  # Optimized triangular solve using transpose = TRUE
  zN    <- backsolve(RN, m_vec, transpose = TRUE)
  mu    <- as.vector(backsolve(RN, zN))           # Posterior mean mu = M * m[cite: 15]
  
  YtY             <- as.numeric(crossprod(Y))
  prior_quad_term <- as.numeric(crossprod(m0, mu0))
  post_quad_term  <- sum(zN^2)
  
  a_N <- a0 + (N / 2)
  b_N <- as.numeric(b0 + 0.5 * (YtY + prior_quad_term - post_quad_term))
  
  # -------------------------------------------------------------------
  # 2. Draw tau directly in R (Exact Marginal Sample)[cite: 15]
  # -------------------------------------------------------------------
  tau_samples    <- rgamma(M_samples, shape = a_N, rate = b_N)
  sigma2_samples <- 1 / tau_samples
  
  # -------------------------------------------------------------------
  # 3. Precompute 3D precision array (Bypasses JAGS DAG build overhead)[cite: 15]
  # -------------------------------------------------------------------
  beta_prec_array <- array(NA, dim = c(M_samples, P, P))
  for (m in 1:M_samples) {
    beta_prec_array[m, , ] <- M0_inv * tau_samples[m] 
  }
  
  # -------------------------------------------------------------------
  # 4. Replicate Data Matrix Y_rep (M_samples x N)[cite: 15]
  # -------------------------------------------------------------------
  Y_rep <- matrix(Y, nrow = M_samples, ncol = N, byrow = TRUE)
  
  # -------------------------------------------------------------------
  # 5. Inline JAGS Model Definition[cite: 15]
  # -------------------------------------------------------------------
  model_string <- "
  model {
    for (m in 1:M_samples) {
      beta[m, 1:P] ~ dmnorm(beta_mean_prior[1:P], beta_prec[m, 1:P, 1:P])
      for (i in 1:N) {
        mu_y[m, i] <- inprod(X[i, 1:P], beta[m, 1:P])
        Y_rep[m, i] ~ dnorm(mu_y[m, i], tau[m])
      }
    }
  }
  "
  
  jags_data <- list(
    N = N,
    P = P,
    M_samples = M_samples,
    X = X,
    Y_rep = Y_rep,
    tau = tau_samples, 
    beta_prec = beta_prec_array,
    beta_mean_prior = as.numeric(mu0)
  )
  
  # -------------------------------------------------------------------
  # 6. Compile and Sample in 1 Iteration without Adaptation[cite: 15]
  # -------------------------------------------------------------------
  jags_model <- jags.model(
    textConnection(model_string), 
    data = jags_data, 
    n.chains = 1,
    n.adapt = 0,
    quiet = TRUE 
  )
  
  jags_out <- jags.samples(
    jags_model, 
    variable.names = c("beta"), 
    n.iter = 1,
    progress.bar = "none"
  )
  
  beta_samples_jags <- jags_out$beta[, , 1, 1]
  
  draws <- cbind(
    beta_samples_jags, 
    sigma2 = sigma2_samples,
    tau    = tau_samples,
    sigma  = sqrt(sigma2_samples)
  )
  
  if (is.null(colnames(X))) {
    colnames(draws)[1:P] <- paste0("beta[", 1:P, "]")
  } else {
    colnames(draws)[1:P] <- colnames(X)
  }
  
  return(list(draws = draws))
}
