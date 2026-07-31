## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  echo = TRUE,
  warning = FALSE,
  message = FALSE,
  fig.width = 10,
  fig.height = 8
)


## ----data_prep----------------------------------------------------------------
library(spBayes)
library(rjags)
library(knitr)

# Source JAGS wrappers from relative source path
source("../src/compositionSamplingJAGS.r")
source("../src/compositionSamplingJAGSDirect.r")

# Load data from project-level data directory
rawData <- read.table("../../data/data.txt", header = TRUE)

# Prepare response and scaled design matrix
Y <- rawData$Y
X_raw <- as.matrix(rawData[, -1])
X_scaled <- scale(X_raw)
X <- cbind(Intercept = 1, X_scaled)

N <- nrow(X)
P <- ncol(X)
n_samples <- 5000

# Shared prior specification (M0 = prior covariance factor)
M0 <- diag(100, P)
m0 <- rep(0, P) # mu0_vec = M0 %*% m0 = 0
a0 <- 0.01
b0 <- 0.01

cat("Dataset loaded successfully:", N, "observations and", P, "predictors (including Intercept).\n")


## ----exact_jags_replicated----------------------------------------------------
set.seed(123)

cat("Executing parallel universe composition sampling in JAGS...\n")
jags_rep_results <- run_replicated_jags_composition(
  X = X, 
  Y = Y, 
  M_samples = n_samples,
  a0 = a0,
  b0 = b0,
  m0 = m0,
  M0 = M0
)

rep_samples <- jags_rep_results$draws

# Summary statistics (transposition via crossprod against identity matrix)
rep_means     <- colMeans(rep_samples)
rep_sds       <- apply(rep_samples, 2, sd)
rep_quantiles <- crossprod(apply(rep_samples, 2, quantile, probs = c(0.5, 0.025, 0.975)), diag(3))

summary_rep <- round(cbind(Mean = rep_means, SD = rep_sds, rep_quantiles), 4)
kable(summary_rep, caption = "JAGS Replicated Likelihood Posterior Summary")


## ----exact_jags_direct--------------------------------------------------------
set.seed(123)

cat("Executing exact direct sampling in JAGS (Likelihood-Free)...\n")
jags_dir_results <- run_direct_jags_composition(
  X = X,
  Y = Y,
  n_samples = n_samples,
  a0 = a0,
  b0 = b0,
  m0 = m0,
  M0 = M0
)

dir_samples <- jags_dir_results$draws

# Summary statistics (transposition via crossprod against identity matrix)
dir_means     <- colMeans(dir_samples)
dir_sds       <- apply(dir_samples, 2, sd)
dir_quantiles <- crossprod(apply(dir_samples, 2, quantile, probs = c(0.5, 0.025, 0.975)), diag(3))

summary_dir <- round(cbind(Mean = dir_means, SD = dir_sds, dir_quantiles), 4)
kable(summary_dir, caption = "JAGS Direct Posterior Summary (Likelihood-Free)")


## ----compare_rep_vs_spbayes---------------------------------------------------
df <- data.frame(Y = Y, X_scaled)

# 1. Run Replicated JAGS
set.seed(123)
jags_rep <- run_replicated_jags_composition(
  X = X, Y = Y, M_samples = n_samples,
  a0 = a0, b0 = b0, m0 = m0, M0 = M0
)$draws

# 2. Run spBayes (prior precision computed via chol2inv(chol(M0)))
set.seed(123)
spb_model_1 <- bayesLMConjugate(
  formula = Y ~ ., data = df, n.samples = n_samples,
  beta.prior.mean = as.vector(crossprod(chol(M0), chol(M0) %*% m0)), 
  beta.prior.precision = chol2inv(chol(M0)),
  prior.shape = a0, prior.rate = b0
)

# Extract and align spBayes parameter columns
spb_raw1     <- spb_model_1$p.beta.tauSq.samples
spb_var1     <- spb_raw1[, P + 1]
spb_sigma1   <- sqrt(spb_var1)
spb_tau1     <- 1 / spb_var1

spb_samples1 <- cbind(spb_raw1[, 1:P], sigma2 = spb_var1, tau = spb_tau1, sigma = spb_sigma1)
colnames(spb_samples1) <- colnames(jags_rep)

# Display Posterior Summaries
s_jags_rep <- round(cbind(Mean = colMeans(jags_rep), SD = apply(jags_rep, 2, sd),
                           crossprod(apply(jags_rep, 2, quantile, probs = c(0.5, 0.025, 0.975)), diag(3))), 4)

s_spb1     <- round(cbind(Mean = colMeans(spb_samples1), SD = apply(spb_samples1, 2, sd),
                           crossprod(apply(spb_samples1, 2, quantile, probs = c(0.5, 0.025, 0.975)), diag(3))), 4)

kable(s_jags_rep, caption = "Posterior Summary: Replicated JAGS Sampler")
kable(s_spb1, caption = "Posterior Summary: spBayes::bayesLMConjugate")


## ----qqplot_rep_vs_spbayes, fig.height=9, fig.width=12------------------------
# Save PNG copy to disk
png("comparisonsJagsVsSpBayesQ-Qplots.png", width = 1200, height = 900, res = 150)
par(mfrow = c(3, 4), mar = c(4, 4, 2, 1))

param_names <- colnames(jags_rep)

for (param in param_names) {
  x_jags    <- sort(jags_rep[, param])
  y_spbayes <- sort(spb_samples1[, param])
  r2        <- cor(x_jags, y_spbayes)^2

  plot(x_jags, y_spbayes,
       main = param,
       xlab = "composition sampling in rjags", 
       ylab = "bayesLMConjugate in spBayes",
       pch = 20, col = rgb(0.1, 0.2, 0.7, alpha = 0.5))
  
  abline(a = 0, b = 1, col = "red", lwd = 2)
  legend("topleft", legend = bquote(R^2 == .(format(r2, digits = 4))), 
         bty = "n", text.col = "black")
}
par(mfrow = c(1, 1))
dev.off()

# Re-plot inline for the knitted Rmd report
par(mfrow = c(3, 4), mar = c(4, 4, 2, 1))
for (param in param_names) {
  x_jags    <- sort(jags_rep[, param])
  y_spbayes <- sort(spb_samples1[, param])
  r2        <- cor(x_jags, y_spbayes)^2

  plot(x_jags, y_spbayes,
       main = param,
       xlab = "composition sampling in rjags", 
       ylab = "bayesLMConjugate in spBayes",
       pch = 20, col = rgb(0.1, 0.2, 0.7, alpha = 0.5))
  
  abline(a = 0, b = 1, col = "red", lwd = 2)
  legend("topleft", legend = bquote(R^2 == .(format(r2, digits = 4))), 
         bty = "n", text.col = "black")
}
par(mfrow = c(1, 1))


## ----compare_dir_vs_spbayes---------------------------------------------------
# 1. Run Direct JAGS
set.seed(123)
jags_dir <- run_direct_jags_composition(
  X = X, Y = Y, n_samples = n_samples,
  a0 = a0, b0 = b0, m0 = m0, M0 = M0
)$draws

# 2. Run spBayes (prior precision computed via chol2inv(chol(M0)))
set.seed(123)
spb_model_2 <- bayesLMConjugate(
  formula = Y ~ ., data = df, n.samples = n_samples,
  beta.prior.mean = as.vector(crossprod(chol(M0), chol(M0) %*% m0)), 
  beta.prior.precision = chol2inv(chol(M0)),
  prior.shape = a0, prior.rate = b0
)

# Extract and align spBayes parameter columns
spb_raw2     <- spb_model_2$p.beta.tauSq.samples
spb_var2     <- spb_raw2[, P + 1]
spb_sigma2   <- sqrt(spb_var2)
spb_tau2     <- 1 / spb_var2

spb_samples2 <- cbind(spb_raw2[, 1:P], sigma2 = spb_var2, tau = spb_tau2, sigma = spb_sigma2)
colnames(spb_samples2) <- colnames(jags_dir)

# Display Posterior Summaries
s_jags_dir <- round(cbind(Mean = colMeans(jags_dir), SD = apply(jags_dir, 2, sd),
                           crossprod(apply(jags_dir, 2, quantile, probs = c(0.5, 0.025, 0.975)), diag(3))), 4)

s_spb2     <- round(cbind(Mean = colMeans(spb_samples2), SD = apply(spb_samples2, 2, sd),
                           crossprod(apply(spb_samples2, 2, quantile, probs = c(0.5, 0.025, 0.975)), diag(3))), 4)

kable(s_jags_dir, caption = "Posterior Summary: Direct JAGS Sampler")
kable(s_spb2, caption = "Posterior Summary: spBayes::bayesLMConjugate")


## ----qqplot_dir_vs_spbayes, fig.height=9, fig.width=12------------------------
# Save PNG copy to disk
png("comparisonsJagsDirectVsSpBayesQ-Qplots.png", width = 1200, height = 900, res = 150)
par(mfrow = c(3, 4), mar = c(4, 4, 2, 1))

param_names <- colnames(jags_dir)

for (param in param_names) {
  x_jags    <- sort(jags_dir[, param])
  y_spbayes <- sort(spb_samples2[, param])
  r2        <- cor(x_jags, y_spbayes)^2

  plot(x_jags, y_spbayes,
       main = param,
       xlab = "composition sampling in rjags", 
       ylab = "bayesLMConjugate sampler in spBayes",
       pch = 20, col = rgb(0.1, 0.2, 0.7, alpha = 0.5))
  
  abline(a = 0, b = 1, col = "red", lwd = 2)
  legend("topleft", legend = bquote(R^2 == .(format(r2, digits = 4))), 
         bty = "n", text.col = "black")
}
par(mfrow = c(1, 1))
dev.off()

# Re-plot inline for the knitted Rmd report
par(mfrow = c(3, 4), mar = c(4, 4, 2, 1))
for (param in param_names) {
  x_jags    <- sort(jags_dir[, param])
  y_spbayes <- sort(spb_samples2[, param])
  r2        <- cor(x_jags, y_spbayes)^2

  plot(x_jags, y_spbayes,
       main = param,
       xlab = "composition sampling in rjags", 
       ylab = "bayesLMConjugate sampler in spBayes",
       pch = 20, col = rgb(0.1, 0.2, 0.7, alpha = 0.5))
  
  abline(a = 0, b = 1, col = "red", lwd = 2)
  legend("topleft", legend = bquote(R^2 == .(format(r2, digits = 4))), 
         bty = "n", text.col = "black")
}
par(mfrow = c(1, 1))

