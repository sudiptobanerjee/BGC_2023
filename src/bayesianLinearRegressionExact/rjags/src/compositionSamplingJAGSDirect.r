# src/compositionSamplingJAGSDirect.r

library(rjags)

run_direct_jags_composition <- function(X, Y, 
                                        n_samples = 5000, 
                                        a0 = 0.01, 
                                        b0 = 0.01, 
                                        m0 = NULL, 
                                        M0_inv = NULL,
                                        mu0 = NULL) {
  N <- nrow(X)
  P <- ncol(X)
  
  # -------------------------------------------------------------------
  # Setup Prior Precision Factor M0_inv and Mean Vectors (m0, mu0)[cite: 14]
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
  # 1. Compute Analytical Posterior Parameters in R[cite: 14]
  # -------------------------------------------------------------------
  M_inv <- M0_inv + crossprod(X)                  # Posterior precision factor M^-1[cite: 14]
  RN    <- chol(M_inv)                             # Cholesky factor of M^-1[cite: 14]
  
  m_vec <- m0 + crossprod(X, Y)                   # Precision-weighted posterior mean m[cite: 14]
  
  # Optimized triangular solve using transpose = TRUE
  zN    <- backsolve(RN, m_vec, transpose = TRUE)
  mu    <- as.vector(backsolve(RN, zN))           # Posterior mean mu = M * m[cite: 14]
  
  YtY             <- as.numeric(crossprod(Y))
  prior_quad_term <- as.numeric(crossprod(m0, mu0))
  post_quad_term  <- sum(zN^2)
  
  a_N <- a0 + (N / 2)
  b_N <- as.numeric(b0 + 0.5 * (YtY + prior_quad_term - post_quad_term))
  
  # -------------------------------------------------------------------
  # 2. Inline Likelihood-Free JAGS Model Definition[cite: 14]
  # -------------------------------------------------------------------
  model_string <- "
  model {
    # Exact Marginal Draw for Precision
    tau ~ dgamma(a_N, b_N)
    sigma2 <- 1 / tau
    
    # Construct Conditional Precision Matrix for Beta
    for (j in 1:P) {
      for (k in 1:P) {
        Omega[j, k] <- tau * M_inv[j, k]
      }
    }
    
    # Exact Conditional Draw for Beta Coefficients
    beta[1:P] ~ dmnorm(mu[1:P], Omega[1:P, 1:P])
  }
  "
  
  jags_data <- list(
    P = P,
    a_N = a_N,
    b_N = b_N,
    mu = as.numeric(mu),
    M_inv = M_inv
  )
  
  # -------------------------------------------------------------------
  # 3. Compile and Sample Direct Forward Composition[cite: 14]
  # -------------------------------------------------------------------
  jags_model <- jags.model(
    textConnection(model_string),
    data = jags_data,
    n.chains = 1,
    n.adapt = 0,
    quiet = TRUE
  )
  
  coda_samples <- coda.samples(
    model = jags_model,
    variable.names = c("beta", "sigma2", "tau"),
    n.iter = n_samples,
    progress.bar = "none"
  )
  
  draws <- as.matrix(coda_samples)
  
  beta_cols <- grep("beta", colnames(draws))
  if (is.null(colnames(X))) {
    colnames(draws)[beta_cols] <- paste0("beta[", 1:P, "]")
  } else {
    colnames(draws)[beta_cols] <- colnames(X)
  }
  
  draws <- cbind(draws, sigma = sqrt(draws[, "sigma2"]))
  
  # Standardize Column Order: (betas, sigma2, tau, sigma)[cite: 14]
  col_order <- c(colnames(X), "sigma2", "tau", "sigma")
  draws <- draws[, col_order]
  
  return(list(draws = draws))
}
