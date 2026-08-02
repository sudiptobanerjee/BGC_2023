# src/compositionSamplingStanDirect.r

library(rstan)

# Configure Stan options for performance[cite: 18]
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

run_direct_stan_composition <- function(X, Y, 
                                        n_samples = 5000, 
                                        M_samples = NULL,
                                        a0 = 0.01, 
                                        b0 = 0.01, 
                                        m0 = NULL, 
                                        M0_inv = NULL,
                                        mu0 = NULL,
                                        V0_inv = NULL) {
  
  N <- nrow(X)
  P <- ncol(X)
  
  # Handle parameter naming compatibility (m0 vs mu0, M0_inv vs V0_inv, M_samples vs n_samples)[cite: 18]
  if (!is.null(M_samples)) n_samples <- M_samples
  if (is.null(m0)) {
    if (!is.null(mu0)) m0 <- mu0 else m0 <- rep(0, P)
  }
  if (is.null(M0_inv)) {
    if (!is.null(V0_inv)) M0_inv <- V0_inv else M0_inv <- diag(0.01, P)
  }
  
  # 1. Compute Exact Marginal Posterior for tau in R[cite: 18]
  XtX <- crossprod(X)
  XtY <- crossprod(X, Y)
  
  Lambda_n <- M0_inv + XtX
  U_n <- chol(Lambda_n) 
  b_vec <- M0_inv %*% m0 + XtY
  
  # Use triangular solvers for mean vector with transpose = TRUE
  mun <- backsolve(U_n, backsolve(U_n, b_vec, transpose = TRUE))
  
  YtY <- crossprod(Y)
  mu0_term <- crossprod(m0, M0_inv %*% m0)
  mun_term <- crossprod(U_n %*% mun)
  
  an <- a0 + (N / 2)
  bn <- as.numeric(b0 + 0.5 * (YtY + mu0_term - mun_term))
  
  tau_samples <- rgamma(n_samples, shape = an, rate = bn)
  sigma2_samples <- 1 / tau_samples
  
  # 2. Package Data for Stan[cite: 18]
  # Invert precision using the Cholesky factor instead of solve()[cite: 18]
  Sigma_n <- chol2inv(U_n) 
  
  stanData <- list(
    P = P,
    M_samples = n_samples,
    mun = as.numeric(mun),
    Sigma_n = Sigma_n,
    tau = tau_samples
  )

  # 3. Inline Stan Model (Generator only)[cite: 18]
  stan_code <- "
  data {
    int<lower=1> P;
    int<lower=1> M_samples;
    vector[P] mun;
    matrix[P, P] Sigma_n;
    vector[M_samples] tau;
  }

  generated quantities {
    matrix[M_samples, P] beta;
    for (m in 1:M_samples) {
      beta[m] = to_row_vector(multi_normal_rng(mun, Sigma_n / tau[m]));
    }
  }
  "

  cat("Compiling and executing Stan C++ Direct Sampler...\n")
  fit <- stan(
    model_code = stan_code,
    data = stanData,
    algorithm = "Fixed_param",
    iter = 1,
    chains = 1,
    refresh = 0
  )

  # 4. Extract generated beta matrix (M x P)[cite: 18]
  beta_array <- extract(fit, pars = "beta")$beta[1, , ]
  
  joint_posterior <- cbind(
    beta_array, 
    sigma2 = sigma2_samples,
    tau = tau_samples,
    sigma = sqrt(sigma2_samples)
  )
  
  colnames(joint_posterior)[1:P] <- paste0("beta[", 1:P, "]")
  
  return(list(
    draws = joint_posterior
  ))
}
