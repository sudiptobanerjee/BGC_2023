## ----setup, message=FALSE, warning=FALSE--------------------------------------
# Load required libraries
library(rjags)
library(nimble)
library(rstan)
library(ggplot2)

# Restrict Stan to 1 core to prevent OOM kills on barebones hardware
rstan_options(auto_write = TRUE)
options(mc.cores = 1)

# Load Data
bench_dir <- getwd()
data_path <- normalizePath("../../data/data.txt", mustWork = FALSE)
rawData <- read.table(data_path, header = TRUE)

Y <- rawData$Y
X <- cbind(Intercept = 1, scale(as.matrix(rawData[, -1])))
N <- nrow(X)
P <- ncol(X)
n_samples <- 5000

# Shared Priors
mu0 <- rep(0, P)
V0_inv <- diag(1/100, P)
a0 <- 0.01
b0 <- 0.01

# Pre-compute exact marginals needed for Stan and NIMBLE
XtX <- crossprod(X)
XtY <- crossprod(X, Y)
Lambda_n <- V0_inv + XtX
U_n <- chol(Lambda_n) 
b_vec <- V0_inv %*% mu0 + XtY
mun <- backsolve(U_n, backsolve(U_n, b_vec, transpose = TRUE))

YtY <- crossprod(Y)
mu0_term <- crossprod(mu0, V0_inv %*% mu0)
mun_term <- crossprod(U_n %*% mun)
an <- a0 + (N / 2)
bn <- as.numeric(b0 + 0.5 * (YtY + mu0_term - mun_term))

set.seed(123)
tau_samples <- rgamma(n_samples, shape = an, rate = bn)


## ----compilation, message=FALSE, warning=FALSE, results='hide'----------------
# --- STAN PRE-COMPILATION ---
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
stan_compiled <- stan_model(model_code = stan_code)
stan_data <- list(
  P = P, M_samples = n_samples, mun = as.numeric(mun), 
  Sigma_n = chol2inv(U_n), tau = tau_samples
)

# --- NIMBLE PRE-COMPILATION ---
nimbleConstants <- list(P = P, an = an, bn = bn, mun = as.numeric(mun), Lambda_n = Lambda_n)
direct_code <- nimbleCode({
  tau ~ dgamma(an, rate = bn)
  sigma2 <- 1 / tau
  prec_beta[1:P, 1:P] <- tau * Lambda_n[1:P, 1:P]
  beta[1:P] ~ dmnorm(mun[1:P], prec = prec_beta[1:P, 1:P])
})

nimble_model <- nimbleModel(code = direct_code, constants = nimbleConstants)
compiled_model <- compileNimble(nimble_model)
mcmc_conf <- configureMCMC(nimble_model, monitors = c("beta", "sigma2", "tau"))
mcmc <- buildMCMC(mcmc_conf)
compiled_mcmc <- compileNimble(mcmc, project = nimble_model)


## ----benchmark, warning=FALSE, message=FALSE----------------------------------
results_list <- list()

# 1. Base R Standard
source("../../raw/src/rawSampler.r")
time_raw <- system.time({
  res_raw <- rawSampler(y = Y, X = X, m0 = mu0, M0_inv = V0_inv, a0 = a0, b0 = b0, N = n_samples)
})
results_list[["Raw_Standard"]] <- time_raw["elapsed"]

# 2. JAGS Direct Wrapper (Interpreter overhead only)
source("../../rjags/src/compositionSamplingJAGSDirect.r")
time_jags <- system.time({
  res_jags <- run_direct_jags_composition(X = X, Y = Y, n_samples = n_samples, a0 = a0, b0 = b0, mu0 = mu0, M0_inv = V0_inv)
})
results_list[["JAGS_Direct"]] <- time_jags["elapsed"]

# 3. Stan (Execution Only)
time_stan <- system.time({
  fit <- sampling(stan_compiled, data = stan_data, algorithm = "Fixed_param", iter = 1, chains = 1, refresh = 0)
})
results_list[["Stan_Execution"]] <- time_stan["elapsed"]

# 4. NIMBLE (Execution Only)
time_nimble <- system.time({
  compiled_mcmc$run(n_samples)
})
results_list[["NIMBLE_Execution"]] <- time_nimble["elapsed"]


## ----display-benchmark, echo=FALSE, fig.width=10, fig.height=6----------------
# Combine results into a data frame
bench_df <- data.frame(
  Method = names(results_list),
  Time_Seconds = unlist(results_list)
)
rownames(bench_df) <- NULL
bench_df <- bench_df[order(bench_df$Time_Seconds), ]

# Print the numerical table
print(bench_df)

# Plot the results
ggplot(bench_df, aes(x = reorder(Method, Time_Seconds), y = Time_Seconds)) +
  geom_bar(stat = "identity", fill = "steelblue", color = "black", alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "Raw Execution Time for Exact Bayesian Sampling",
    subtitle = paste(n_samples, "Draws | Note: C++ compilation times are completely excluded"),
    x = "Method",
    y = "Time (Seconds)"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12))

