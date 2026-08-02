library(nimble)

run_direct_nimble_composition <- function(X, Y,
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

  # Handle parameter naming compatibility (m0 vs mu0, M0_inv vs V0_inv, M_samples vs n_samples)
  if (!is.null(M_samples)) n_samples <- M_samples
  if (is.null(m0)) {
    if (!is.null(mu0)) m0 <- mu0 else m0 <- rep(0, P)
  }
  if (is.null(M0_inv)) {
    if (!is.null(V0_inv)) M0_inv <- V0_inv else M0_inv <- diag(0.01, P)
  }

  # 1. Compute Analytical Posterior Parameters in R
  XtX <- crossprod(X)
  XtY <- crossprod(X, Y)

  Lambda_n <- M0_inv + XtX
  U_n <- chol(Lambda_n)

  b_vec <- M0_inv %*% m0 + XtY
  mun <- backsolve(U_n, forwardsolve(t(U_n), b_vec))

  YtY <- crossprod(Y)
  mu0_term <- crossprod(m0, M0_inv %*% m0)
  mun_term <- crossprod(U_n %*% mun)

  an <- a0 + (N / 2)
  bn <- as.numeric(b0 + 0.5 * (YtY + mu0_term - mun_term))

  # 2. Package constants for NIMBLE
  nimbleConstants <- list(
    P = P,
    an = an,
    bn = bn,
    mun = as.numeric(mun),
    Lambda_n = Lambda_n
  )

  # 3. Define inline NIMBLE model code (self-contained, no external txt file needed)
  direct_code <- nimbleCode({
    tau ~ dgamma(an, rate = bn)
    sigma2 <- 1 / tau
    prec_beta[1:P, 1:P] <- tau * Lambda_n[1:P, 1:P]
    beta[1:P] ~ dmnorm(mun[1:P], prec = prec_beta[1:P, 1:P])
  })

  # 4. Build and Compile the NIMBLE Model to C++
  cat("\nBuilding direct NIMBLE model...\n")
  nimble_model <- nimbleModel(code = direct_code, constants = nimbleConstants)

  cat("Compiling model to C++...\n")
  compiled_model <- compileNimble(nimble_model)

  # 5. Configure and Compile Forward Sampler
  mcmc_conf <- configureMCMC(nimble_model, monitors = c("beta", "sigma2", "tau"))
  mcmc <- buildMCMC(mcmc_conf)
  compiled_mcmc <- compileNimble(mcmc, project = nimble_model)

  # 6. Run exact forward sampling (IID draws, zero burn-in)
  cat("\nDrawing", n_samples, "independent composition samples...\n")
  samples <- runMCMC(compiled_mcmc, niter = n_samples, nburnin = 0, nchains = 1, progressBar = FALSE)

  # 7. Format matrix and output
  draws <- as.matrix(samples)

  # Re-order and rename columns to standardize with replicated sampler
  beta_cols <- grep("beta", colnames(draws))
  colnames(draws)[beta_cols] <- paste0("beta[", 1:P, "]")
  
  draws <- draws[, c(paste0("beta[", 1:P, "]"), "sigma2", "tau")]
  draws <- cbind(draws, sigma = sqrt(draws[, "sigma2"]))

  return(list(
    draws = draws
  ))
}
