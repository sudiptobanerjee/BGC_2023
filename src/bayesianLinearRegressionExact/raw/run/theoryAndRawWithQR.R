## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)


## ----setup-libs, message=FALSE, warning=FALSE---------------------------------
# Source raw sampler implementations from ../src/
source("../src/rawSampler.r")
source("../src/rawSamplerWithQR.r")
source("../src/rawPosteriorPredictive.r")

library(knitr)


## ----simulate-data------------------------------------------------------------
set.seed(2026)
n <- 250
p <- 3
N_samples <- 10000

# Predictors with extreme scale differences
X_raw <- cbind(
  Intercept = 1,
  X1 = runif(n, 0, 1),
  X2 = runif(n, 1e4, 1e5)
)

true_beta <- c(12.5, -4.2, 0.00085)
true_sigma2 <- 4.0

# Generate response y
y <- as.numeric(X_raw %*% true_beta + rnorm(n, mean = 0, sd = sqrt(true_sigma2)))


## ----fit-models---------------------------------------------------------------
# Set flat prior parameters
m0 <- rep(0, p)
M0_inv <- matrix(0, p, p)
a0 <- -p / 2
b0 <- 0

# -----------------------------------------------------------------------------
# Method 1: Standardized Sampler (rawSampler.r)
# -----------------------------------------------------------------------------
X_center <- colMeans(X_raw[, -1])
X_scale  <- apply(X_raw[, -1], 2, sd)
X_scaled <- cbind(Intercept = 1, scale(X_raw[, -1]))

res_scaled <- rawSampler(y, X_scaled, m0, M0_inv, a0, b0, N_samples)

# Centered intercept draws (alpha ~ y_bar)
alpha_scaled_draws <- res_scaled$beta[, 1]

# Unscale Beta draws back to raw data units (Vectorized)
beta_scaled_draws <- res_scaled$beta
beta_unscaled_draws <- beta_scaled_draws
beta_unscaled_draws[, 2:3] <- sweep(beta_scaled_draws[, 2:3], 2, X_scale, "/")
beta_unscaled_draws[, 1]   <- beta_scaled_draws[, 1] - 
  rowSums(sweep(beta_scaled_draws[, 2:3], 2, X_center / X_scale, "*"))

# -----------------------------------------------------------------------------
# Method 2: QR Sampler on Centered Predictors (rawSamplerWithQR.r)
# -----------------------------------------------------------------------------
X_centered <- cbind(Intercept = 1, scale(X_raw[, -1], center = TRUE, scale = FALSE))
qr_X <- qr(X_centered)
Q <- qr.Q(qr_X)
R <- qr.R(qr_X)

res_qr <- rawSamplerWithQR(y, Q, R, m0, M0_inv, a0, b0, N_samples)

# Centered intercept draws (alpha ~ y_bar)
alpha_qr_draws <- res_qr$beta[, 1]

# Convert centered QR Beta draws to raw unscaled predictors (Vectorized)
beta_qr_draws <- res_qr$beta
beta_unscaled_qr <- beta_qr_draws
beta_unscaled_qr[, 1] <- beta_qr_draws[, 1] - rowSums(sweep(beta_qr_draws[, 2:3], 2, X_center, "*"))


## ----posterior-table, echo=FALSE----------------------------------------------
calc_stats <- function(draws_matrix, sigma2_vec) {
  means <- colMeans(draws_matrix)
  sds   <- apply(draws_matrix, 2, sd)
  q025  <- apply(draws_matrix, 2, quantile, probs = 0.025)
  q975  <- apply(draws_matrix, 2, quantile, probs = 0.975)
  
  data.frame(
    Mean = sprintf("%.6f", means),
    SD   = sprintf("%.6f", sds),
    `2.5%` = sprintf("%.6f", q025),
    `97.5%` = sprintf("%.6f", q975),
    check.names = FALSE
  )
}

draws_scaled_full <- cbind(Centered_Alpha = alpha_scaled_draws, beta_unscaled_draws, sigma2 = res_scaled$sigma2)
draws_qr_full     <- cbind(Centered_Alpha = alpha_qr_draws, beta_unscaled_qr, sigma2 = res_qr$sigma2)

stats_scaled <- calc_stats(draws_scaled_full, res_scaled$sigma2)
stats_qr     <- calc_stats(draws_qr_full, res_qr$sigma2)

param_names <- c("alpha (Centered Intercept ~ y_bar)", "beta_0 (Raw Intercept at X=0)", "beta_1 (X1)", "beta_2 (X2)", "sigma^2")

comp_table <- data.frame(
  Parameter = rep(param_names, each = 2),
  Method    = rep(c("rawSampler (Standardized X)", "rawSamplerWithQR (Centered X)"), times = length(param_names)),
  Mean      = c(rbind(stats_scaled$Mean, stats_qr$Mean)),
  SD        = c(rbind(stats_scaled$SD, stats_qr$SD)),
  `2.5%`    = c(rbind(stats_scaled$`2.5%`, stats_qr$`2.5%`)),
  `97.5%`   = c(rbind(stats_scaled$`97.5%`, stats_qr$`97.5%`)),
  check.names = FALSE
)

kable(comp_table, caption = "Posterior Estimation Comparison: rawSampler vs. rawSamplerWithQR", align = "lccccc")


## ----posterior-qqplots, fig.width=9, fig.height=8, fig.align='center'---------
par(mfrow = c(2, 2), mar = c(4.5, 4.5, 2.5, 1))

qqplot(alpha_scaled_draws, alpha_qr_draws,
       xlab = "rawSampler Centered Intercept", ylab = "rawSamplerWithQR Centered Intercept",
       main = expression("Q-Q Plot for Centered Intercept (" * alpha * " ~ " * bar(y) * ")"))
abline(0, 1, col = "firebrick", lwd = 2)

qqplot(beta_unscaled_draws[, 1], beta_unscaled_qr[, 1],
       xlab = "rawSampler Raw Intercept", ylab = "rawSamplerWithQR Raw Intercept",
       main = expression("Q-Q Plot for Raw Intercept (" * beta[0] * ")"))
abline(0, 1, col = "firebrick", lwd = 2)

qqplot(beta_unscaled_draws[, 2], beta_unscaled_qr[, 2],
       xlab = "rawSampler Draw Quantiles", ylab = "rawSamplerWithQR Draw Quantiles",
       main = expression("Q-Q Plot for " * beta[1]))
abline(0, 1, col = "firebrick", lwd = 2)

qqplot(res_scaled$sigma2, res_qr$sigma2,
       xlab = "rawSampler Draw Quantiles", ylab = "rawSamplerWithQR Draw Quantiles",
       main = expression("Q-Q Plot for " * sigma^2))
abline(0, 1, col = "firebrick", lwd = 2)


## ----predictive-sampling------------------------------------------------------
m_test <- 100
X_new_raw <- cbind(
  Intercept = 1,
  X1 = runif(m_test, 0, 1),
  X2 = runif(m_test, 1e4, 1e5)
)

# 1. Predictions using rawSamplerWithQR Draws
Y_pred_qr <- rawPosteriorPredictive(X_new_raw, beta_unscaled_qr, res_qr$sigma2)

# 2. Predictions using rawSampler Draws (unscaled parameters)
Y_pred_scaled <- rawPosteriorPredictive(X_new_raw, beta_unscaled_draws, res_scaled$sigma2)


## ----predictive-qqplots, fig.width=9, fig.height=8, fig.align='center'--------
par(mfrow = c(2, 2), mar = c(4.5, 4.5, 2.5, 1))

test_indices <- c(1, 25, 50, 100)
invisible(lapply(test_indices, function(idx) {
  qqplot(Y_pred_scaled[idx, ], Y_pred_qr[idx, ],
         xlab = "rawSampler Predictive Quantiles",
         ylab = "rawSamplerWithQR Predictive Quantiles",
         main = paste("Posterior Predictive Q-Q Plot: Obs #", idx))
  abline(0, 1, col = "dodgerblue3", lwd = 2)
}))


## ----empirical-analysis-------------------------------------------------------
# Read full empirical dataset
data_path <- if (file.exists("../../data/data.txt")) "../../data/data.txt" else "data.txt"
full_data <- read.table(data_path, header = TRUE)

y_emp <- full_data[[1]]
X_emp_covariates <- as.matrix(full_data[, -1])

n_emp <- length(y_emp)
p_emp <- ncol(X_emp_covariates) + 1

X_center_emp <- colMeans(X_emp_covariates)
X_scale_emp  <- apply(X_emp_covariates, 2, sd)

# Flat Jeffreys prior setup
m0_emp <- rep(0, p_emp)
M0_inv_emp <- matrix(0, p_emp, p_emp)
a0_emp <- -p_emp / 2
b0_emp <- 0

# -----------------------------------------------------------------------------
# Method 1: Standardized Sampler on Empirical Data (rawSampler.r)
# -----------------------------------------------------------------------------
X_scaled_emp <- cbind(Intercept = 1, scale(X_emp_covariates))

res_scaled_emp <- rawSampler(y_emp, X_scaled_emp, m0_emp, M0_inv_emp, a0_emp, b0_emp, N_samples)

# Centered intercept draws (alpha ~ y_bar ~ 11.00)
alpha_scaled_emp <- res_scaled_emp$beta[, 1]

# Unscale Beta draws back to natural physical units (Vectorized)
beta_scaled_draws_emp <- res_scaled_emp$beta
beta_unscaled_draws_emp <- beta_scaled_draws_emp
beta_unscaled_draws_emp[, -1] <- sweep(beta_scaled_draws_emp[, -1], 2, X_scale_emp, "/")
beta_unscaled_draws_emp[, 1]  <- beta_scaled_draws_emp[, 1] - 
  rowSums(sweep(beta_scaled_draws_emp[, -1], 2, X_center_emp / X_scale_emp, "*"))

# -----------------------------------------------------------------------------
# Method 2: QR Sampler on Centered Empirical Predictors (rawSamplerWithQR.r)
# Mean-centering predictors anchors centered intercept alpha at y_bar (~11.00)
# -----------------------------------------------------------------------------
X_centered_emp <- cbind(Intercept = 1, scale(X_emp_covariates, center = TRUE, scale = FALSE))
qr_X_emp <- qr(X_centered_emp)
Q_emp <- qr.Q(qr_X_emp)
R_emp <- qr.R(qr_X_emp)

res_qr_emp <- rawSamplerWithQR(y_emp, Q_emp, R_emp, m0_emp, M0_inv_emp, a0_emp, b0_emp, N_samples)

# Centered intercept draws (alpha ~ y_bar ~ 11.00)
alpha_qr_emp <- res_qr_emp$beta[, 1]

# Convert centered QR Beta draws to raw unscaled predictors (Vectorized)
beta_qr_draws_emp <- res_qr_emp$beta
beta_unscaled_qr_emp <- beta_qr_draws_emp
beta_unscaled_qr_emp[, 1] <- beta_qr_draws_emp[, 1] - 
  rowSums(sweep(beta_qr_draws_emp[, -1], 2, X_center_emp, "*"))


## ----empirical-table, echo=FALSE----------------------------------------------
draws_scaled_emp_full <- cbind(Centered_Alpha = alpha_scaled_emp, beta_unscaled_draws_emp, sigma2 = res_scaled_emp$sigma2)
draws_qr_emp_full     <- cbind(Centered_Alpha = alpha_qr_emp, beta_unscaled_qr_emp, sigma2 = res_qr_emp$sigma2)

stats_scaled_emp <- calc_stats(draws_scaled_emp_full, res_scaled_emp$sigma2)
stats_qr_emp     <- calc_stats(draws_qr_emp_full, res_qr_emp$sigma2)

param_names_emp <- c("alpha (Centered Intercept ~ y_bar)", "beta_0 (Raw Intercept at X=0)", colnames(X_emp_covariates), "sigma^2")

comp_table_emp <- data.frame(
  Parameter = rep(param_names_emp, each = 2),
  Method    = rep(c("rawSampler (Standardized X)", "rawSamplerWithQR (Centered X)"), times = length(param_names_emp)),
  Mean      = c(rbind(stats_scaled_emp$Mean, stats_qr_emp$Mean)),
  SD        = c(rbind(stats_scaled_emp$SD, stats_qr_emp$SD)),
  `2.5%`    = c(rbind(stats_scaled_emp$`2.5%`, stats_qr_emp$`2.5%`)),
  `97.5%`   = c(rbind(stats_scaled_emp$`97.5%`, stats_qr_emp$`97.5%`)),
  check.names = FALSE
)

kable(comp_table_emp, caption = "Empirical Posterior Estimation Comparison: rawSampler vs. rawSamplerWithQR (data.txt)", align = "lccccc")


## ----empirical-qqplots, fig.width=10, fig.height=9, fig.align='center'--------
par(mfrow = c(3, 4), mar = c(4, 4, 2.5, 1))

# Centered Intercept Alpha
qqplot(alpha_scaled_emp, alpha_qr_emp,
       xlab = "rawSampler Draw Quantiles", 
       ylab = "rawSamplerWithQR Draw Quantiles",
       main = expression("Q-Q Plot: " * alpha * " (~ " * bar(y) * ")"))
abline(0, 1, col = "firebrick", lwd = 2)

# Raw Intercept Beta_0 and Slopes
invisible(lapply(1:p_emp, function(i) {
  p_name <- if (i == 1) "beta_0 (Raw Intercept)" else colnames(X_emp_covariates)[i - 1]
  qqplot(beta_unscaled_draws_emp[, i], beta_unscaled_qr_emp[, i],
         xlab = "rawSampler Draw Quantiles", 
         ylab = "rawSamplerWithQR Draw Quantiles",
         main = paste("Q-Q Plot:", p_name))
  abline(0, 1, col = "firebrick", lwd = 2)
}))

# Sigma^2
qqplot(res_scaled_emp$sigma2, res_qr_emp$sigma2,
       xlab = "rawSampler Draw Quantiles", 
       ylab = "rawSamplerWithQR Draw Quantiles",
       main = expression("Q-Q Plot for " * sigma^2))
abline(0, 1, col = "firebrick", lwd = 2)


## ----empirical-predictive-sampling--------------------------------------------
# Construct raw unscaled design matrix for empirical observations
X_emp_raw <- cbind(Intercept = 1, X_emp_covariates)

# 1. Generate posterior predictive distributions using QR Sampler draws
Y_pred_qr_emp <- rawPosteriorPredictive(X_emp_raw, beta_unscaled_qr_emp, res_qr_emp$sigma2)

# 2. Generate posterior predictive distributions using Standard Sampler draws
Y_pred_scaled_emp <- rawPosteriorPredictive(X_emp_raw, beta_unscaled_draws_emp, res_scaled_emp$sigma2)

# Compute posterior predictive point estimates and 95% interval bounds
pred_means_qr <- rowMeans(Y_pred_qr_emp)
pred_q025_qr  <- apply(Y_pred_qr_emp, 1, quantile, probs = 0.025)
pred_q975_qr  <- apply(Y_pred_qr_emp, 1, quantile, probs = 0.975)

# Calculate empirical 95% predictive coverage and Root Mean Squared Error (RMSE)
emp_coverage <- mean(y_emp >= pred_q025_qr & y_emp <= pred_q975_qr) * 100
emp_rmse     <- sqrt(mean((y_emp - pred_means_qr)^2))


## ----empirical-predictive-plots, fig.width=10, fig.height=8, fig.align='center'----
par(mfrow = c(2, 2), mar = c(4.5, 4.5, 2.5, 1))

# Q-Q Plots across selected empirical test observations
emp_test_indices <- c(1, 500, 1000)
invisible(lapply(emp_test_indices, function(idx) {
  qqplot(Y_pred_scaled_emp[idx, ], Y_pred_qr_emp[idx, ],
         xlab = "rawSampler Predictive Quantiles",
         ylab = "rawSamplerWithQR Predictive Quantiles",
         main = paste("Empirical Predictive Q-Q: Obs #", idx))
  abline(0, 1, col = "dodgerblue3", lwd = 2)
}))

# Posterior Predictive Check: Observed Response vs Posterior Predictive Means
plot(y_emp, pred_means_qr,
     xlab = "Observed Response (y_emp)", 
     ylab = "Posterior Predictive Mean",
     main = expression("PPC: Observed vs. Predictive Mean (" * R^2 * " Check)"),
     pch = 16, col = rgb(0.1, 0.4, 0.8, 0.3))
abline(0, 1, col = "firebrick", lwd = 2)


## ----empirical-predictive-summary, echo=FALSE---------------------------------
pred_metrics <- data.frame(
  Metric = c("Posterior Predictive RMSE", "95% Credible Interval Coverage", "Max Predictive Discrepancy (QR vs Std)"),
  Value  = c(sprintf("%.4f", emp_rmse), sprintf("%.2f%%", emp_coverage), sprintf("%.2e", max(abs(rowMeans(Y_pred_qr_emp) - rowMeans(Y_pred_scaled_emp))))),
  check.names = FALSE
)

kable(pred_metrics, caption = "Empirical Posterior Predictive Metrics (data.txt)", align = "lc")


## ----benchmark-timing---------------------------------------------------------
# Benchmark 1: Simulated Data Runtime (n = 250, p = 3)
time_sim_standard <- system.time({
  res_sim_std <- rawSampler(y, X_scaled, m0, M0_inv, a0, b0, N_samples)
})["elapsed"]

time_sim_qr <- system.time({
  res_sim_qr <- rawSamplerWithQR(y, Q, R, m0, M0_inv, a0, b0, N_samples)
})["elapsed"]

# Benchmark 2: Empirical Data Runtime (n = 1369, p = 9)
time_emp_standard <- system.time({
  res_emp_std <- rawSampler(y_emp, X_scaled_emp, m0_emp, M0_inv_emp, a0_emp, b0_emp, N_samples)
})["elapsed"]

time_emp_qr <- system.time({
  res_emp_qr <- rawSamplerWithQR(y_emp, Q_emp, R_emp, m0_emp, M0_inv_emp, a0_emp, b0_emp, N_samples)
})["elapsed"]

# Compute relative speedup factors
speedup_sim <- time_sim_standard / time_sim_qr
speedup_emp <- time_emp_standard / time_emp_qr

# Construct summary benchmark table
benchmark_df <- data.frame(
  Dataset = c("Simulated Data (n = 250, p = 3)", "Empirical Data (n = 1369, p = 9)"),
  `Standard Sampler (s)` = sprintf("%.4f", c(time_sim_standard, time_emp_standard)),
  `QR Sampler (s)` = sprintf("%.4f", c(time_sim_qr, time_emp_qr)),
  `Speedup Factor` = sprintf("%.2fx", c(speedup_sim, speedup_emp)),
  check.names = FALSE
)

kable(benchmark_df, caption = "Execution Time Benchmarking for 10,000 Draws", align = "lccc")

