# src/compositionSamplingNIMBLE.r

library(nimble)

run_replicated_nimble_composition <- function(X, Y, M_samples, a0, b0, m0, M0_inv) {
  N <- nrow(X)
  P <- ncol(X)
  
  # 1. Exact posterior parameters via R (Cholesky/Triangular Solves)
  M_star_inv <- M0_inv + crossprod(X)
  U <- chol(M_star_inv)
  v <- M0_inv %*% m0 + crossprod(X, Y)
  
  # Optimized triangular solve using transpose = TRUE
  y_solve <- backsolve(U, v, transpose = TRUE)
  mu_star <- backsolve(U, y_solve)
  
  a_star <- a0 + N / 2
  b_star <- as.numeric(b0 + 0.5 * (crossprod(Y) + crossprod(m0, M0_inv %*% m0) - crossprod(mu_star, v)))
  
  # Pre-sample tau
  tau_samples <- rgamma(M_samples, shape = a_star, rate = b_star)
  
  # 2. Define a Dummy Model purely for exact simulation (Now including Y_rep)
  direct_code <- nimbleCode({
    prec_beta[1:P, 1:P] <- tau[1] * M_star_inv[1:P, 1:P]
    beta[1:P] ~ dmnorm(mu_star[1:P], prec = prec_beta[1:P, 1:P])
    
    # Simulating a SINGLE replication of size N to avoid memory overload
    for(i in 1:N) {
      mu_y[i] <- inprod(X[i, 1:P], beta[1:P])
      Y_rep[i] ~ dnorm(mu_y[i], tau = tau[1])
    }
  })
  
  # X is now passed as a constant so it is baked into the C++ model
  constants <- list(P = P, N = N, mu_star = as.numeric(mu_star), M_star_inv = M_star_inv, X = X)
  inits <- list(tau = c(tau_samples[1]), beta = rep(0, P), Y_rep = rep(0, N))
  
  # Build and compile model
  model <- nimbleModel(code = direct_code, constants = constants, inits = inits)
  cModel <- compileNimble(model)
  
  # 3. Use nimbleFunction to loop and simulate in C++ (Bypassing MCMC entirely)
  simulate_loop <- nimbleFunction(
    setup = function(model, M, P, N) {},
    run = function(tau_vec = double(1)) {
      # Return a combined matrix to easily pass both beta and Y_rep back to R
      returnType(double(2))
      out <- matrix(0, nrow = M, ncol = P + N)
      
      for(m in 1:M) {
        # 1. Update the scalar tau
        model$tau[1] <<- tau_vec[m]
        
        # 2. Update deterministic precision and simulate exact beta
        model$calculate('prec_beta') 
        model$simulate('beta')       
        
        # 3. Update deterministic mu_y and simulate exact Y_rep
        model$calculate('mu_y')
        model$simulate('Y_rep')
        
        # 4. Store results in the pre-allocated C++ matrix
        out[m, 1:P] <- model$beta[1:P]
        out[m, (P + 1):(P + N)] <- model$Y_rep[1:N]
      }
      return(out)
    }
  )
  
  c_loop <- compileNimble(simulate_loop(model = cModel, M = M_samples, P = P, N = N), project = model)
  
  # 4. Execute the loop
  raw_samples <- c_loop$run(tau_samples)
  
  # 5. Extract and Format Output
  beta_samples <- raw_samples[, 1:P]
  Y_rep_samples <- raw_samples[, (P + 1):(P + N)]
  
  draws <- matrix(NA, nrow = M_samples, ncol = P + 3)
  colnames(draws) <- c(paste0("beta[", 1:P, "]"), "sigma2", "tau", "sigma")
  
  draws[, 1:P] <- beta_samples
  draws[, "tau"] <- tau_samples
  draws[, "sigma2"] <- 1 / tau_samples
  draws[, "sigma"] <- sqrt(draws[, "sigma2"])
  
  return(list(
    draws = draws,
    Y_rep = Y_rep_samples
  ))
}
