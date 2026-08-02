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


## ----data-loading-------------------------------------------------------------
# Load the empirical dataset
# Assuming a standard tabular format (e.g., space or tab-separated)
dataset <- read.table("../../data/data.txt", header = TRUE)

# Parse Response (Y) and Predictors (X)
# NOTE: This assumes Y is the FIRST column. Change to dataset[, ncol(dataset)] if Y is last.
Y <- as.numeric(dataset[, 1])
X <- as.matrix(dataset[, -1])

# Extract Dimensions
N <- nrow(X)
P <- ncol(X)

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


## ----benchmark, warning=FALSE, message=FALSE, cache=FALSE, results='hide'-----
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

