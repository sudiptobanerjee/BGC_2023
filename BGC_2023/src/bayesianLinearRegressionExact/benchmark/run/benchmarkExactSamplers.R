## ----setup, message=FALSE, warning=FALSE--------------------------------------
# Load required libraries
library(microbenchmark)
library(rjags)
library(nimble)
library(rstan)

# Configure Stan to utilize parallel cores
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# ---------------------------------------------------------
# Source Modular Samplers
# ---------------------------------------------------------
# 1. Raw Base R Samplers
source("../../raw/src/rawSampler.r")
source("../../raw/src/rawSamplerWithQR.r")

# 2. JAGS Samplers
source("../../rjags/src/compositionSamplingJAGSDirect.r")
source("../../rjags/src/compositionSamplingJAGS.r")

# 3. NIMBLE Samplers
source("../../nimble/src/compositionSamplingNimbleDirect.r")
source("../../nimble/src/compositionSamplingNimble.r")

# 4. Stan Samplers
source("../../rstan/src/compositionSamplingStan.r")


## ----data-simulation----------------------------------------------------------
# Set seed for reproducible results
set.seed(42)

# Dimensions
N <- 250  # Number of observations
P <- 10   # Number of predictors

# Simulate Predictors (X) and True Parameters
X <- matrix(rnorm(N * P), nrow = N, ncol = P)
true_beta <- rnorm(P, mean = 0, sd = 2)
true_sigma <- 1.5

# Simulate Response (Y)
Y <- X %*% true_beta + rnorm(N, mean = 0, sd = true_sigma)
Y <- as.numeric(Y)

# Define Shared Prior Parameters
m0 <- rep(0, P)
M0_inv <- diag(0.01, P)
a0 <- 0.01
b0 <- 0.01

# Calculate QR Decomposition for the QR-optimized Raw Sampler
qr_decomp <- qr(X)
Q_mat <- qr.Q(qr_decomp)
R_mat <- qr.R(qr_decomp)

# Specify number of posterior draws requested
n_draws <- 5000


## ----benchmark, warning=FALSE, message=FALSE, cache=TRUE, results='hide'------
# Run the benchmark (times = 5 to account for compilation overhead in C++ methods)
benchmark_results <- microbenchmark(
  
  # Base R Implementations
  Raw_Standard = rawSampler(
    y = Y, X = X, m0 = m0, M0_inv = M0_inv, a0 = a0, b0 = b0, N = n_draws
  ),
  
  Raw_QR = rawSamplerWithQR(
    y = Y, Q = Q_mat, R = R_mat, m0 = m0, M0_inv = M0_inv, a0 = a0, b0 = b0, N = n_draws
  ),
  
  # JAGS Implementations
  JAGS_Direct = run_direct_jags_composition(
    X = X, Y = Y, n_samples = n_draws, a0 = a0, b0 = b0, m0 = m0, M0_inv = M0_inv
  ),
  
  JAGS_Replicated = run_replicated_jags_composition(
    X = X, Y = Y, M_samples = n_draws, a0 = a0, b0 = b0, m0 = m0, M0_inv = M0_inv
  ),
  
  # NIMBLE Implementations
  NIMBLE_Direct = run_direct_nimble_composition(
    X = X, Y = Y, n_samples = n_draws, a0 = a0, b0 = b0, m0 = m0, M0_inv = M0_inv
  ),
  
  NIMBLE_Replicated = run_replicated_nimble_composition(
    X = X, Y = Y, M_samples = n_draws, a0 = a0, b0 = b0, m0 = m0, M0_inv = M0_inv
  ),
  
  # Stan Implementation
  Stan_Direct = run_direct_stan_composition(
    X = X, Y = Y, n_samples = n_draws, a0 = a0, b0 = b0, m0 = m0, M0_inv = M0_inv
  ),
  
  times = 5
)


## ----display-benchmark, echo=FALSE--------------------------------------------
# Print formatted results
print(benchmark_results, unit = "s", signif = 4)


## ----plot-benchmark, echo=FALSE, fig.width=10, fig.height=6-------------------
# Visualize the benchmark distributions
library(ggplot2)
autoplot(benchmark_results) + 
  theme_minimal() +
  labs(
    title = "Execution Time Comparison for Exact Bayesian Sampling",
    subtitle = paste(n_draws, "Posterior Draws |", N, "Observations |", P, "Predictors"),
    x = "Method",
    y = "Time (Seconds)"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12))

