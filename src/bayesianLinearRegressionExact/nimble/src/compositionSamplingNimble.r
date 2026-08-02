library(nimble)

run_replicated_nimble_composition <- function(X, Y, M_samples, a0, b0, m0, M0_inv) {
  N <- nrow(X)
  P <- ncol(X)
  
  # 1. Exact posterior parameters via R (Cholesky/Triangular Solves)
  M_star_inv <- M0_inv + crossprod(X)
  U <- chol(M_star_inv)
  v <- M0_inv %*% m0 + crossprod(X, Y)
  y_solve <- forwardsolve(t(U), v)
  mu_star <- backsolve(U, y_solve)
  
  a_star <- a0 + N / 2
  b_star <- as.numeric(b0 + 0.5 * (crossprod(Y) + crossprod(m0, M0_inv %*% m0) - crossprod(mu_star, v)))
  
  # Pre-sample tau
  tau_samples <- rgamma(M_samples, shape = a_star, rate = b_star)
  
  # 2. Define a Dummy Model purely for exact simulation
  direct_code <- nimbleCode({
    # Explicitly index tau[1] so NIMBLE perfectly maps the C++ array dimensions
    prec_beta[1:P, 1:P] <- tau[1] * M_star_inv[1:P, 1:P]
    beta[1:P] ~ dmnorm(mu_star[1:P], prec = prec_beta[1:P, 1:P])
  })
  
  constants <- list(P = P, mu_star = as.numeric(mu_star), M_star_inv = M_star_inv)
  
  # Initialize tau as a vector of length 1
  inits <- list(tau = c(tau_samples[1]), beta = rep(0, P))
  
  # Build and compile model
  model <- nimbleModel(code = direct_code, constants = constants, inits = inits)
  cModel <- compileNimble(model)
  
  # 3. Use nimbleFunction to loop and simulate in C++ (Bypassing MCMC entirely)
  simulate_loop <- nimbleFunction(
    setup = function(model, M, P) {},
    run = function(tau_vec = double(1)) {
      returnType(double(2))
      out <- matrix(0, nrow = M, ncol = P)
      
      for(m in 1:M) {
        # FIX: Explicitly assign the scalar to the first index of the model's array
        model$tau[1] <<- tau_vec[m]
        
        # Update the deterministic precision matrix
        model$calculate('prec_beta') 
        
        # Draw EXACTLY from the conditional distribution
        model$simulate('beta')       
        
        out[m, 1:P] <- model$beta[1:P]
      }
      return(out)
    }
  )
  
  c_loop <- compileNimble(simulate_loop(model = cModel, M = M_samples, P = P), project = model)
  
  # 4. Execute the loop
  beta_samples <- c_loop$run(tau_samples)
  
  # 5. Format Output for your Q-Q plots
  draws <- matrix(NA, nrow = M_samples, ncol = P + 3)
  colnames(draws) <- c(paste0("beta[", 1:P, "]"), "sigma2", "tau", "sigma")
  
  draws[, 1:P] <- beta_samples
  draws[, "tau"] <- tau_samples
  draws[, "sigma2"] <- 1 / tau_samples
  draws[, "sigma"] <- sqrt(draws[, "sigma2"])
  
  return(list(draws = draws))
}
