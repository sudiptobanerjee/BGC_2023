## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(knitr)
library(mniw)

# Source the custom MNIW sampler for parameter inference
source("../src/mniwExactSampler.r")


## ----initialize_workspace, message=FALSE, warning=FALSE-----------------------
# 1. Data Ingestion
Y <- as.matrix(read.table("../data/Y.txt", header = FALSE))
X_raw <- as.matrix(read.table("../data/X.txt", header = FALSE))

n <- nrow(X_raw)
X <- cbind(1, X_raw)
k <- ncol(X)
p <- ncol(Y)

# Extract 20 points from the training set for prediction targets
set.seed(8675309)
m_in_sample <- 20
sample_indices <- sample(1:nrow(X_raw), m_in_sample)
X_raw_tilde <- X_raw[sample_indices, , drop = FALSE]
X_tilde <- cbind(1, X_raw_tilde)

# 2. Prior Specification
B0 <- matrix(0, nrow = k, ncol = p)
Lambda0 <- diag(1, k)
nu0 <- p + 2
S0 <- diag(1, p)
n_samples <- 1000

# 3. Custom Posterior Sampling (Method 2)
set.seed(42)
custom_results <- bayes_mvr(Y = Y, X = X, B0 = B0, Lambda0 = Lambda0, 
                            nu0 = nu0, S0 = S0, method = 2, n_samples = n_samples)

beta_samples <- custom_results$samples$B
if(is.list(beta_samples)) beta_samples <- simplify2array(beta_samples)

sigma_samples <- custom_results$samples$Sigma
if(is.list(sigma_samples)) sigma_samples <- simplify2array(sigma_samples)

# 4. Native `mniw` Posterior Sampling
Mn_inv <- Lambda0 + crossprod(X)
U_inv <- chol(Mn_inv) 
mn <- Lambda0 %*% B0 + crossprod(X, Y)

y_tmp <- backsolve(U_inv, mn, transpose = TRUE)
post_mean <- backsolve(U_inv, y_tmp)

nu_n <- nu0 + nrow(X)
Psi_n <- S0 + crossprod(Y) + crossprod(B0, Lambda0 %*% B0) - crossprod(mn, post_mean)

set.seed(42)
mniw_draws <- rmniw(n = n_samples, Lambda = post_mean, Omega = Mn_inv, Psi = Psi_n, nu = nu_n)
beta_mniw_array <- mniw_draws$X
sigma_mniw_array <- mniw_draws$V


## ----custom_predictive_sampler------------------------------------------------
N_draws <- dim(beta_samples)[3]
q_cols <- dim(beta_samples)[2]

set.seed(123)
# Vectorize the sampling over all posterior draws
Y_tilde_custom_list <- lapply(1:N_draws, function(i) {
  B_i <- beta_samples[, , i]
  Sigma_i <- sigma_samples[, , i]

  M_i <- X_tilde %*% B_i
  Z <- matrix(rnorm(m_in_sample * q_cols), nrow = m_in_sample, ncol = q_cols)
  U_Sigma <- chol(Sigma_i)

  return(M_i + Z %*% U_Sigma)
})

# Collapse the list into a 3D array
Y_tilde_custom <- simplify2array(Y_tilde_custom_list)


## ----predictive_diagnostics_exact, fig.width=11, fig.height=5.5---------------
# Extract predictions for the primary outcome 
Y_tilde <- Y_tilde_custom[, 1, ] 
Y_true  <- Y[sample_indices, 1]

# Compute posterior predictive metrics
pred_mean <- rowMeans(Y_tilde)
PI_lower  <- apply(Y_tilde, 1, quantile, probs = 0.025)
PI_upper  <- apply(Y_tilde, 1, quantile, probs = 0.975)

# Calculate Empirical Coverage Probability
covered <- (Y_true >= PI_lower) & (Y_true <= PI_upper)
coverage_prob <- mean(covered)

cat(sprintf("**Empirical In-Sample 95%% Coverage:** %.2f%% (%d / %d observations covered)\n", 
            coverage_prob * 100, sum(covered), length(Y_true)))

# Define color schemes for covered vs uncovered observations
pt_cols  <- ifelse(covered, "black", "firebrick")
bar_cols <- ifelse(covered, rgb(0.3, 0.6, 0.85, 0.6), rgb(1, 0.4, 0.4, 0.7))

par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))

# Plot 1: Ordered 95% PIs vs True In-Sample Y
sort_idx     <- order(Y_true)
Y_true_sort  <- Y_true[sort_idx]
pred_sort    <- pred_mean[sort_idx]
lower_sort   <- PI_lower[sort_idx]
upper_sort   <- PI_upper[sort_idx]
pt_cols_sort <- pt_cols[sort_idx]
bar_cols_sort<- bar_cols[sort_idx]

y_lims <- range(c(lower_sort, upper_sort, Y_true_sort))

plot(1:length(Y_true), Y_true_sort, type = "n", ylim = y_lims,
     xlab = "In-Sample Observation Index (Sorted by True Y)", 
     ylab = "Target Variable Y",
     main = sprintf("In-Sample 95%% PIs (Coverage: %.1f%%)", coverage_prob * 100))

segments(x0 = 1:length(Y_true), y0 = lower_sort, 
         x1 = 1:length(Y_true), y1 = upper_sort, 
         col = bar_cols_sort, lwd = 1.5)
points(1:length(Y_true), pred_sort, pch = 18, col = "mediumblue", cex = 0.9)
points(1:length(Y_true), Y_true_sort, pch = 16, col = pt_cols_sort, cex = 0.9)

legend("topleft", 
       legend = c("Actual Y (Covered)", "Actual Y (Uncovered)", "Predictive Mean", "95% Predictive Interval"),
       col = c("black", "firebrick", "mediumblue", rgb(0.3, 0.6, 0.85)),
       pch = c(16, 16, 18, NA), lty = c(NA, NA, NA, 1), lwd = c(NA, NA, NA, 2),
       bty = "o", bg = "white", cex = 0.75)

# Plot 2: Observed vs. Predicted Scatter
plot(Y_true, pred_mean, type = "n", ylim = y_lims, xlim = y_lims,
     xlab = "Observed In-Sample Y", 
     ylab = "Predicted Y (Posterior Mean)",
     main = "Observed vs. Predicted In-Sample Y")

abline(0, 1, col = "red", lty = 2, lwd = 2)

segments(x0 = Y_true, y0 = PI_lower, 
         x1 = Y_true, y1 = PI_upper, 
         col = bar_cols, lwd = 1.2)
points(Y_true, pred_mean, pch = 16, col = pt_cols, cex = 0.9)

legend("topleft", 
       legend = c("1:1 Reference Line", "Covered Observation", "Uncovered Observation"),
       col = c("red", "black", "firebrick"),
       lty = c(2, NA, NA), pch = c(NA, 16, 16), lwd = c(2, NA, NA),
       bty = "o", bg = "white", cex = 0.75)

par(mfrow = c(1, 1))


## ----data_split---------------------------------------------------------------
# Load raw design and response data
X_raw <- as.matrix(read.table("../data/X.txt"))
Y     <- as.matrix(read.table("../data/Y.txt"))

set.seed(42)

# Define 20% hold-out proportion
n_total <- nrow(X_raw)
m <- floor(0.20 * n_total)

# Partition indices
test_idx  <- sample(1:n_total, size = m)
train_idx <- setdiff(1:n_total, test_idx)

# Training Set (Used exclusively for model estimation)
X_train <- cbind(1, X_raw[train_idx, , drop = FALSE])
Y_train <- Y[train_idx, , drop = FALSE]

# Held-out Test Set (Used exclusively for validation)
X_test  <- cbind(1, X_raw[test_idx, , drop = FALSE]) # X_tilde
Y_test  <- Y[test_idx, , drop = FALSE]                # Actual held-out targets


## ----fit_training_model-------------------------------------------------------
library(mniw)

# 1. Training dimensions
n_tr   <- nrow(X_train)
p_cols <- ncol(X_train)
q_cols <- ncol(Y_train)

# 2. BLAS-optimized OLS Statistics via Cholesky decomposition
R_X <- chol(crossprod(X_train))
XtY <- crossprod(X_train, Y_train)

# Compute B_hat = (X^T X)^{-1} X^T Y using forward/backward substitution
B_hat <- backsolve(R_X, backsolve(R_X, XtY, transpose = TRUE))

# Compute Residual Sum of Squares (E^T E) using crossprod
E_residuals <- Y_train - X_train %*% B_hat
S_residuals <- crossprod(E_residuals)

# 3. Conjugate Posterior Hyperparameters
Psi_post <- S_residuals
nu_post  <- n_tr - p_cols
N_draws  <- 2000

# 4. Draw Inverse-Wishart covariance matrices using mniw::riwish
set.seed(123)
sigma_samples <- riwish(N_draws, Psi = Psi_post, nu = nu_post)

# 5. Draw Beta parameters using functional vectorization (lapply)
beta_list <- lapply(1:N_draws, function(i) {
  Sigma_i <- sigma_samples[, , i]
  U_Sigma <- chol(Sigma_i)
  Z       <- matrix(rnorm(p_cols * q_cols), nrow = p_cols, ncol = q_cols)

  return(B_hat + backsolve(R_X, Z) %*% U_Sigma)
})

# Collapse list into a (p x q x N_draws) array
beta_samples <- simplify2array(beta_list)


## ----custom_out_of_sample-----------------------------------------------------
set.seed(123)

Y_tilde_custom_list <- lapply(1:N_draws, function(i) {
  B_i     <- beta_samples[, , i]
  Sigma_i <- sigma_samples[, , i]
  
  M_i     <- X_test %*% B_i
  Z       <- matrix(rnorm(m * q_cols), nrow = m, ncol = q_cols)
  U_Sigma <- chol(Sigma_i)
  
  return(M_i + Z %*% U_Sigma)
})

Y_tilde_custom <- simplify2array(Y_tilde_custom_list)


## ----mniw_out_of_sample-------------------------------------------------------
set.seed(123)

Y_tilde_mniw_list <- lapply(1:N_draws, function(i) {
  B_i     <- beta_samples[, , i]
  Sigma_i <- sigma_samples[, , i]
  
  M_i <- X_test %*% B_i
  return(drop(rMNorm(n = 1, Lambda = M_i, SigmaR = diag(m), SigmaC = Sigma_i)))
})

Y_tilde_mniw <- simplify2array(Y_tilde_mniw_list)


## ----validation_diagnostics, fig.width=11, fig.height=5.5---------------------
# Extract predictions for primary outcome Y[,1]
Y_tilde <- Y_tilde_custom[, 1, ] # Dimension: (m x N_draws)
Y_true  <- Y_test[, 1]

# Posterior predictive summary metrics
pred_mean <- rowMeans(Y_tilde)
PI_lower  <- apply(Y_tilde, 1, quantile, probs = 0.025)
PI_upper  <- apply(Y_tilde, 1, quantile, probs = 0.975)

# Calculate empirical out-of-sample coverage
covered <- (Y_true >= PI_lower) & (Y_true <= PI_upper)
coverage_prob <- mean(covered)

cat(sprintf("**Empirical Out-of-Sample 95%% Coverage:** %.2f%% (%d / %d test points covered)\n", 
            coverage_prob * 100, sum(covered), m))

# Styling setup
pt_cols  <- ifelse(covered, "black", "firebrick")
bar_cols <- ifelse(covered, rgb(0.3, 0.6, 0.85, 0.6), rgb(1, 0.4, 0.4, 0.7))

par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))

# --- Plot 1: Ordered 95% PIs vs True Hold-out Y ---
sort_idx    <- order(Y_true)
Y_true_sort <- Y_true[sort_idx]
pred_sort   <- pred_mean[sort_idx]
lower_sort  <- PI_lower[sort_idx]
upper_sort  <- PI_upper[sort_idx]
pt_sort     <- pt_cols[sort_idx]
bar_sort    <- bar_cols[sort_idx]

y_lims <- range(c(lower_sort, upper_sort, Y_true_sort))

plot(1:m, Y_true_sort, type = "n", ylim = y_lims,
     xlab = "Hold-out Observation Index (Sorted by True Y)", ylab = "Target Variable Y",
     main = sprintf("Out-of-Sample 95%% PIs (Coverage: %.1f%%)", coverage_prob * 100))

segments(x0 = 1:m, y0 = lower_sort, x1 = 1:m, y1 = upper_sort, col = bar_sort, lwd = 1.5)
points(1:m, pred_sort, pch = 18, col = "mediumblue", cex = 0.8)
points(1:m, Y_true_sort, pch = 16, col = pt_sort, cex = 0.8)

legend("topleft", 
       legend = c("Actual Y (Covered)", "Actual Y (Uncovered)", "Predictive Mean", "95% Predictive Interval"),
       col = c("black", "firebrick", "mediumblue", rgb(0.3, 0.6, 0.85)),
       pch = c(16, 16, 18, NA), lty = c(NA, NA, NA, 1), lwd = c(NA, NA, NA, 2), bty = "o", bg = "white", cex = 0.75)

# --- Plot 2: Observed vs. Predicted Scatter ---
plot(Y_true, pred_mean, type = "n", ylim = y_lims, xlim = y_lims,
     xlab = "Observed Hold-out Y", ylab = "Predicted Y (Posterior Mean)",
     main = "Observed vs. Predicted Hold-out Y")

abline(0, 1, col = "red", lty = 2, lwd = 2)
segments(x0 = Y_true, y0 = PI_lower, x1 = Y_true, y1 = PI_upper, col = bar_cols, lwd = 1.2)
points(Y_true, pred_mean, pch = 16, col = pt_cols, cex = 0.8)

legend("topleft", 
       legend = c("1:1 Reference Line", "Covered Observation", "Uncovered Observation"),
       col = c("red", "black", "firebrick"), lty = c(2, NA, NA), pch = c(NA, 16, 16), lwd = c(2, NA, NA),
       bty = "o", bg = "white", cex = 0.75)

par(mfrow = c(1, 1))


## ----modular_custom_predictive, fig.width=11, fig.height=5.5------------------
# Source updated helper functions
source("../src/mniwPosteriorPredictive.r")

# Generate posterior predictive draws via mniw_posterior_predictive
set.seed(123)
Y_tilde_mod_custom <- mniw_posterior_predictive(X_test, beta_samples, sigma_samples)

# Extract metrics for primary target variable Y[,1]
Y_pred_draws <- Y_tilde_mod_custom[, 1, ] # (m x N_draws)
Y_true       <- Y_test[, 1]

pred_mean <- rowMeans(Y_pred_draws)
PI_lower  <- apply(Y_pred_draws, 1, quantile, probs = 0.025)
PI_upper  <- apply(Y_pred_draws, 1, quantile, probs = 0.975)

covered       <- (Y_true >= PI_lower) & (Y_true <= PI_upper)
coverage_prob <- mean(covered)

# --- Summary Table of Test Set Predictions ---
summary_df <- data.frame(
  Test_Obs  = 1:m,
  True_Y    = round(Y_true, 3),
  Pred_Mean = round(pred_mean, 3),
  Lower_95  = round(PI_lower, 3),
  Upper_95  = round(PI_upper, 3),
  Covered   = covered
)

knitr::kable(head(summary_df, 10), caption = "First 10 Test Set Predictions (mniw_posterior_predictive)")

cat(sprintf("\n**Empirical Out-of-Sample 95%% Coverage (`mniw_posterior_predictive`):** %.2f%% (%d / %d test points covered)\n", 
            coverage_prob * 100, sum(covered), m))

# --- Diagnostic Plotting ---
pt_cols  <- ifelse(covered, "black", "firebrick")
bar_cols <- ifelse(covered, rgb(0.3, 0.6, 0.85, 0.6), rgb(1, 0.4, 0.4, 0.7))

par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))

# Plot 1: Ordered 95% PIs
sort_idx    <- order(Y_true)
Y_true_sort <- Y_true[sort_idx]
pred_sort   <- pred_mean[sort_idx]
lower_sort  <- PI_lower[sort_idx]
upper_sort  <- PI_upper[sort_idx]
pt_sort     <- pt_cols[sort_idx]
bar_sort    <- bar_cols[sort_idx]

y_lims <- range(c(lower_sort, upper_sort, Y_true_sort))

plot(1:m, Y_true_sort, type = "n", ylim = y_lims,
     xlab = "Hold-out Observation Index (Sorted by True Y)", ylab = "Target Variable Y",
     main = sprintf("mniw_posterior_predictive 95%% PIs (Coverage: %.1f%%)", coverage_prob * 100))

segments(x0 = 1:m, y0 = lower_sort, x1 = 1:m, y1 = upper_sort, col = bar_sort, lwd = 1.5)
points(1:m, pred_sort, pch = 18, col = "mediumblue", cex = 0.8)
points(1:m, Y_true_sort, pch = 16, col = pt_sort, cex = 0.8)

legend("topleft", 
       legend = c("Actual Y (Covered)", "Actual Y (Uncovered)", "Predictive Mean", "95% Predictive Interval"),
       col = c("black", "firebrick", "mediumblue", rgb(0.3, 0.6, 0.85)),
       pch = c(16, 16, 18, NA), lty = c(NA, NA, NA, 1), lwd = c(NA, NA, NA, 2), bty = "o", bg = "white", cex = 0.75)

# Plot 2: Observed vs. Predicted Scatter
plot(Y_true, pred_mean, type = "n", ylim = y_lims, xlim = y_lims,
     xlab = "Observed Hold-out Y", ylab = "Predicted Y (Posterior Mean)",
     main = "Observed vs. Predicted Hold-out Y")

abline(0, 1, col = "red", lty = 2, lwd = 2)
segments(x0 = Y_true, y0 = PI_lower, x1 = Y_true, y1 = PI_upper, col = bar_cols, lwd = 1.2)
points(Y_true, pred_mean, pch = 16, col = pt_cols, cex = 0.8)

legend("topleft", 
       legend = c("1:1 Reference Line", "Covered Observation", "Uncovered Observation"),
       col = c("red", "black", "firebrick"), lty = c(2, NA, NA), pch = c(NA, 16, 16), lwd = c(2, NA, NA),
       bty = "o", bg = "white", cex = 0.75)

par(mfrow = c(1, 1))


## ----modular_pkg_predictive, fig.width=11, fig.height=5.5---------------------
# Generate posterior predictive draws via mniw_posterior_predictive_pkg
set.seed(123)
Y_tilde_mod_pkg <- mniw_posterior_predictive_pkg(X_test, beta_samples, sigma_samples)

# Extract metrics for primary target variable Y[,1]
Y_pred_draws_pkg <- Y_tilde_mod_pkg[, 1, ] # (m x N_draws)
Y_true           <- Y_test[, 1]

pred_mean_pkg <- rowMeans(Y_pred_draws_pkg)
PI_lower_pkg  <- apply(Y_pred_draws_pkg, 1, quantile, probs = 0.025)
PI_upper_pkg  <- apply(Y_pred_draws_pkg, 1, quantile, probs = 0.975)

covered_pkg       <- (Y_true >= PI_lower_pkg) & (Y_true <= PI_upper_pkg)
coverage_prob_pkg <- mean(covered_pkg)

# --- Summary Table of Test Set Predictions ---
summary_pkg_df <- data.frame(
  Test_Obs  = 1:m,
  True_Y    = round(Y_true, 3),
  Pred_Mean = round(pred_mean_pkg, 3),
  Lower_95  = round(PI_lower_pkg, 3),
  Upper_95  = round(PI_upper_pkg, 3),
  Covered   = covered_pkg
)

knitr::kable(head(summary_pkg_df, 10), caption = "First 10 Test Set Predictions (mniw_posterior_predictive_pkg)")

cat(sprintf("\n**Empirical Out-of-Sample 95%% Coverage (`mniw_posterior_predictive_pkg`):** %.2f%% (%d / %d test points covered)\n", 
            coverage_prob_pkg * 100, sum(covered_pkg), m))

# --- Diagnostic Plotting ---
pt_cols  <- ifelse(covered_pkg, "black", "firebrick")
bar_cols <- ifelse(covered_pkg, rgb(0.3, 0.6, 0.85, 0.6), rgb(1, 0.4, 0.4, 0.7))

par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))

# Plot 1: Ordered 95% PIs
sort_idx    <- order(Y_true)
Y_true_sort <- Y_true[sort_idx]
pred_sort   <- pred_mean_pkg[sort_idx]
lower_sort  <- PI_lower_pkg[sort_idx]
upper_sort  <- PI_upper_pkg[sort_idx]
pt_sort     <- pt_cols[sort_idx]
bar_sort    <- bar_cols[sort_idx]

y_lims <- range(c(lower_sort, upper_sort, Y_true_sort))

plot(1:m, Y_true_sort, type = "n", ylim = y_lims,
     xlab = "Hold-out Observation Index (Sorted by True Y)", ylab = "Target Variable Y",
     main = sprintf("mniw_posterior_predictive_pkg 95%% PIs (Coverage: %.1f%%)", coverage_prob_pkg * 100))

segments(x0 = 1:m, y0 = lower_sort, x1 = 1:m, y1 = upper_sort, col = bar_sort, lwd = 1.5)
points(1:m, pred_sort, pch = 18, col = "mediumblue", cex = 0.8)
points(1:m, Y_true_sort, pch = 16, col = pt_sort, cex = 0.8)

legend("topleft", 
       legend = c("Actual Y (Covered)", "Actual Y (Uncovered)", "Predictive Mean", "95% Predictive Interval"),
       col = c("black", "firebrick", "mediumblue", rgb(0.3, 0.6, 0.85)),
       pch = c(16, 16, 18, NA), lty = c(NA, NA, NA, 1), lwd = c(NA, NA, NA, 2), bty = "o", bg = "white", cex = 0.75)

# Plot 2: Observed vs. Predicted Scatter
plot(Y_true, pred_mean_pkg, type = "n", ylim = y_lims, xlim = y_lims,
     xlab = "Observed Hold-out Y", ylab = "Predicted Y (Posterior Mean)",
     main = "Observed vs. Predicted Hold-out Y")

abline(0, 1, col = "red", lty = 2, lwd = 2)
segments(x0 = Y_true, y0 = PI_lower_pkg, x1 = Y_true, y1 = PI_upper_pkg, col = bar_cols, lwd = 1.2)
points(Y_true, pred_mean_pkg, pch = 16, col = pt_cols, cex = 0.8)

legend("topleft", 
       legend = c("1:1 Reference Line", "Covered Observation", "Uncovered Observation"),
       col = c("red", "black", "firebrick"), lty = c(2, NA, NA), pch = c(NA, 16, 16), lwd = c(2, NA, NA),
       bty = "o", bg = "white", cex = 0.75)

par(mfrow = c(1, 1))

