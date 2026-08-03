## ----setup, message=FALSE, warning=FALSE--------------------------------------
# Load required libraries
library(rjags)
library(nimble)
library(rstan)
library(ggplot2)

# Strictly limit Stan to 1 core to prevent Out-Of-Memory (OOM) compilation crashes
rstan_options(auto_write = TRUE)
options(mc.cores = 1)

# Load Data & Process Predictors
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

# Calculate QR Decomposition for the QR-optimized Raw Sampler
qr_decomp <- qr(X)
Q_mat <- qr.Q(qr_decomp)
R_mat <- qr.R(qr_decomp)


## ----benchmark, warning=FALSE, message=FALSE, results='hide'------------------
results_list <- list()

# 1. Base R Standard
source("../../raw/src/rawSampler.r")
time_raw <- system.time({
  res_raw <- rawSampler(y = Y, X = X, m0 = mu0, M0_inv = V0_inv, a0 = a0, b0 = b0, N = n_samples)
})
results_list[["Raw_Standard"]] <- time_raw["elapsed"]
rm(res_raw)
gc()

# 2. Base R QR-Optimized
source("../../raw/src/rawSamplerWithQR.r")
time_qr <- system.time({
  res_qr <- rawSamplerWithQR(y = Y, Q = Q_mat, R = R_mat, m0 = mu0, M0_inv = V0_inv, a0 = a0, b0 = b0, N = n_samples)
})
results_list[["Raw_QR"]] <- time_qr["elapsed"]
rm(res_qr)
gc()

# 3. JAGS Direct Wrapper
source("../../rjags/src/compositionSamplingJAGSDirect.r")
time_jags <- system.time({
  res_jags <- run_direct_jags_composition(X = X, Y = Y, n_samples = n_samples, a0 = a0, b0 = b0, mu0 = mu0, M0_inv = V0_inv)
})
results_list[["JAGS_Direct"]] <- time_jags["elapsed"]
rm(res_jags)
gc()

# 4. NIMBLE Direct Wrapper
source("../../nimble/src/compositionSamplingNimbleDirect.r")
time_nimble <- system.time({
  res_nimble <- run_direct_nimble_composition(X = X, Y = Y, n_samples = n_samples, a0 = a0, b0 = b0, mu0 = mu0, V0_inv = V0_inv)
})
results_list[["NIMBLE_Direct"]] <- time_nimble["elapsed"]
rm(res_nimble)
gc()

# 5. Stan Direct Wrapper
source("../../rstan/src/compositionSamplingStan.r")
time_stan <- system.time({
  res_stan <- run_direct_stan_composition(X = X, Y = Y, n_samples = n_samples, a0 = a0, b0 = b0, mu0 = mu0, V0_inv = V0_inv)
})
results_list[["Stan_Direct"]] <- time_stan["elapsed"]
rm(res_stan)
gc()


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
    title = "End-to-End Execution Time for Exact Bayesian Sampling Wrappers",
    subtitle = paste(n_samples, "Draws | Note: NIMBLE and Stan timings include C++ compilation overhead"),
    x = "Method",
    y = "Time (Seconds)"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12))

