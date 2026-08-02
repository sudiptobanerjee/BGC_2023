## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  echo = TRUE,
  warning = FALSE,
  message = FALSE,
  fig.width = 10,
  fig.height = 8
)


## ----helper_functions---------------------------------------------------------
library(spBayes)
library(rstan)
library(knitr)

# Source translated RStan composition sampler function
source("../src/compositionSamplingStan.r")


## ----data_prep----------------------------------------------------------------
# Load Data from two levels above rstan/run
rawData <- read.table("../../data/data.txt", header = TRUE)

# Prepare response and scaled design matrix
Y <- rawData$Y
X_raw <- as.matrix(rawData[, -1])
X_scaled <- scale(X_raw)
X <- cbind(Intercept = 1, X_scaled)

N <- nrow(X)
P <- ncol(X)
n_samples <- 5000

# Shared Prior Specifications (Prior Precision Factor M0_inv)
M0_inv <- diag(0.01, P) 
m0     <- rep(0, P)     
a0     <- 0.01
b0     <- 0.01

cat("Dataset loaded successfully:", N, "observations and", P, "predictors.\n")


## ----exact_stan_direct--------------------------------------------------------
set.seed(123)

cat("Compiling and executing exact direct sampling in RStan...\n")
runtime_stan <- system.time({
  stan_dir_results <- run_direct_stan_composition(
    X = X, Y = Y, n_samples = n_samples,
    a0 = a0, b0 = b0, m0 = m0, M0_inv = M0_inv
  )
})

print(runtime_stan)
stan_samples <- stan_dir_results$draws

summary_stan <- round(cbind(Mean = colMeans(stan_samples), SD = apply(stan_samples, 2, sd)), 4)
kable(summary_stan, caption = "RStan Direct Composition Posterior Summary")


## ----qqplot_stan_vs_spbayes, fig.height=9, fig.width=12-----------------------
df <- data.frame(Y = Y, X_scaled)

set.seed(123)
spb_model <- bayesLMConjugate(
  formula = Y ~ ., data = df, n.samples = n_samples,
  beta.prior.mean = m0, beta.prior.precision = M0_inv,
  prior.shape = a0, prior.rate = b0
)

spb_raw     <- spb_model$p.beta.tauSq.samples
spb_samples <- cbind(spb_raw[, 1:P], sigma2 = spb_raw[, P + 1], tau = 1/spb_raw[, P + 1], sigma = sqrt(spb_raw[, P + 1]))
colnames(spb_samples) <- colnames(stan_samples)

par(mfrow = c(3, 4), mar = c(4, 4, 2, 1))
param_names <- colnames(stan_samples)

stan_mat <- as.matrix(stan_samples)
spb_mat  <- as.matrix(spb_samples)

for (param in param_names) {
  x_stan_raw <- as.numeric(stan_mat[, param])
  y_spb_raw  <- as.numeric(spb_mat[, param])
  
  # Quantile calculation for consistent plot resolution and alignment
  q_stan <- quantile(x_stan_raw, probs = seq(0, 1, length.out = 1000), na.rm = TRUE)
  q_spb  <- quantile(y_spb_raw, probs = seq(0, 1, length.out = 1000), na.rm = TRUE)
  
  r2 <- cor(q_stan, q_spb)^2

  plot(q_stan, q_spb,
       main = param,
       xlab = "direct composition sampling in RStan", 
       ylab = "bayesLMConjugate sampler in spBayes",
       pch = 20, col = rgb(0.1, 0.2, 0.7, alpha = 0.5))
  
  abline(a = 0, b = 1, col = "red", lwd = 2)
  legend("topleft", legend = bquote(R^2 == .(format(r2, digits = 4))), 
         bty = "n", text.col = "black")
}
par(mfrow = c(1, 1))

