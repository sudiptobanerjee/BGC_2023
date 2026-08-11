## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(knitr)

# Source the optimized MNIW sampler functions
source("../src/mniwExactSampler.r")


## ----load_data----------------------------------------------------------------
# Read the continuous numerical data into standard matrices
Y <- as.matrix(read.table("../data/Y.txt", header = FALSE))
X_raw <- as.matrix(read.table("../data/X.txt", header = FALSE))

# Prepend an intercept to the design matrix
n <- nrow(X_raw)
X <- cbind(1, X_raw)

# Extract dimensions for prior specification
k <- ncol(X)
p <- ncol(Y)


## ----specify_priors-----------------------------------------------------------
# B0: Prior mean for regression coefficients (k x p)
B0 <- matrix(0, nrow = k, ncol = p)

# Lambda0: Prior precision matrix for rows (k x k)
Lambda0 <- diag(1, k)

# nu0: Prior degrees of freedom (must be > p - 1)
nu0 <- p + 2

# S0: Prior scale matrix (p x p)
S0 <- diag(1, p)


## ----execute_model------------------------------------------------------------
# Execute the Bayesian Multivariate Linear Regression
# We utilize Method 2 (Precision Parameterized) to efficiently sample
set.seed(42)
n_samples <- 1000

model_results <- bayes_mvr(
  Y = Y, 
  X = X, 
  B0 = B0, 
  Lambda0 = Lambda0, 
  nu0 = nu0, 
  S0 = S0, 
  method = 2, 
  n_samples = n_samples
)

# Extract samples and coerce to 3D arrays if necessary
beta_samples <- model_results$samples$B
if(is.list(beta_samples)) beta_samples <- simplify2array(beta_samples)

sigma_samples <- model_results$samples$Sigma
if(is.list(sigma_samples)) sigma_samples <- simplify2array(sigma_samples)


## ----summarize_posterior, echo=FALSE, results='asis'--------------------------
# Helper function to compute summaries for a 3D array of posterior samples
generate_summary_table <- function(samples_array, param_prefix) {
  n_rows <- dim(samples_array)[1]
  n_cols <- dim(samples_array)[2]
  
  # Create a grid of matrix indices to avoid looping
  grid <- expand.grid(Row = 1:n_rows, Col = 1:n_cols)
  
  # Process summaries efficiently over the index grid
  summary_list <- lapply(1:nrow(grid), function(idx) {
    i <- grid$Row[idx]
    j <- grid$Col[idx]
    samps <- samples_array[i, j, ]
    
    data.frame(
      Parameter = paste0(param_prefix, "[", i, ",", j, "]"),
      Mean = round(mean(samps), 4),
      SD = round(sd(samps), 4),
      `2.5%` = round(quantile(samps, 0.025), 4),
      `50%` = round(median(samps), 4),
      `97.5%` = round(quantile(samps, 0.975), 4),
      check.names = FALSE,
      row.names = NULL
    )
  })
  
  df_final <- do.call(rbind, summary_list)
  # Explicitly clear row names to prevent artifact columns in markdown
  rownames(df_final) <- NULL
  return(df_final)
}

# Generate dataframes
df_beta <- generate_summary_table(beta_samples, "Beta")
df_sigma <- generate_summary_table(sigma_samples, "Sigma")

# Print beautifully formatted tables using knitr::kable
cat("## Posterior Summary for $\\beta$\n\n")
print(kable(df_beta, format = "markdown", align = "c", row.names = FALSE))

cat("\n## Posterior Summary for $\\Sigma$\n\n")
print(kable(df_sigma, format = "markdown", align = "c", row.names = FALSE))


## ----bayesm_benchmark, message=FALSE, warning=FALSE, results='asis'-----------
# Load the bayesm package
library(bayesm)

# 1. Execute the model
# rmultireg takes arguments directly and draws 1 sample per call, 
# Using lapply instead of replicate/sapply to avoid hidden loops
set.seed(42)
bayesm_draws <- lapply(1:n_samples, function(x) {
  rmultireg(Y = Y, X = X, Bbar = B0, A = Lambda0, nu = nu0, V = S0)
})

# 2. Extract and format the samples into 3D arrays to match our earlier structure
beta_bayesm_list <- lapply(bayesm_draws, function(x) x$B)
sigma_bayesm_list <- lapply(bayesm_draws, function(x) x$Sigma)

beta_bayesm_array <- simplify2array(beta_bayesm_list)
sigma_bayesm_array <- simplify2array(sigma_bayesm_list)

# 3. Generate summaries using the helper function defined in Section 6
df_beta_bayesm <- generate_summary_table(beta_bayesm_array, "Beta (bayesm)")
df_sigma_bayesm <- generate_summary_table(sigma_bayesm_array, "Sigma (bayesm)")

# 4. Print the tables for comparison
cat("### `bayesm` Posterior Summary for $\\beta$\n\n")
print(kable(df_beta_bayesm, format = "markdown", align = "c", row.names = FALSE))

cat("\n### `bayesm` Posterior Summary for $\\Sigma$\n\n")
print(kable(df_sigma_bayesm, format = "markdown", align = "c", row.names = FALSE))


## ----qqplots_bayesm, fig.width=10, fig.height=8, echo=FALSE, message=FALSE, warning=FALSE----
# Determine dimensions for plotting
k_dim <- dim(beta_samples)[1]
p_dim <- dim(beta_samples)[2]

# --- QQ-Plots for Beta ---
# Set up plotting area dynamically based on matrix dimensions
par(mfrow = c(k_dim, p_dim), mar = c(4, 4, 2, 1))
grid_beta <- expand.grid(Row = 1:k_dim, Col = 1:p_dim)

invisible(lapply(1:nrow(grid_beta), function(idx) {
  i <- grid_beta$Row[idx]
  j <- grid_beta$Col[idx]
  qqplot(beta_samples[i, j, ], beta_bayesm_array[i, j, ],
         main = paste0("QQ-Plot: Beta[", i, ",", j, "]"),
         xlab = "Custom Sampler Quantiles", 
         ylab = "bayesm Quantiles",
         pch = 16, col = rgb(0.2, 0.4, 0.6, 0.5))
  abline(0, 1, col = "red", lwd = 2)
}))

# --- QQ-Plots for Sigma ---
# Set up plotting area for the p x p covariance matrix
par(mfrow = c(p_dim, p_dim), mar = c(4, 4, 2, 1))
grid_sigma <- expand.grid(Row = 1:p_dim, Col = 1:p_dim)

invisible(lapply(1:nrow(grid_sigma), function(idx) {
  i <- grid_sigma$Row[idx]
  j <- grid_sigma$Col[idx]
  qqplot(sigma_samples[i, j, ], sigma_bayesm_array[i, j, ],
         main = paste0("QQ-Plot: Sigma[", i, ",", j, "]"),
         xlab = "Custom Sampler Quantiles", 
         ylab = "bayesm Quantiles",
         pch = 16, col = rgb(0.2, 0.4, 0.6, 0.5))
  abline(0, 1, col = "red", lwd = 2)
}))

# Reset graphics parameters
par(mfrow = c(1, 1))


## ----mniw_benchmark, message=FALSE, warning=FALSE, results='asis'-------------
# Load the mniw package
library(mniw)

# 1. Compute posterior parameters analytically
# Utilizing highly efficient, numerically stable linear algebra 
# (avoiding solve() and t() via Cholesky decomposition and cross products)

# Evaluate precision matrix: Mn_inv = Lambda0 + X'X
Mn_inv <- Lambda0 + crossprod(X)

# Cholesky decomposition of the precision matrix
U_inv <- chol(Mn_inv) 

# Posterior linear anchor: mn = Lambda0 %*% B0 + X'Y
mn <- Lambda0 %*% B0 + crossprod(X, Y)

# Solve for posterior mean: post_mean = (Mn_inv)^{-1} mn
# We solve the system t(U_inv) %*% U_inv %*% post_mean = mn using stable triangular solvers
y_tmp <- backsolve(U_inv, mn, transpose = TRUE)
post_mean <- backsolve(U_inv, y_tmp)

# Solve for posterior scale matrix: Mn = (Mn_inv)^{-1} = U_inv^{-1} * (U_inv')^{-1}
#U_inv_inv <- backsolve(U_inv, diag(nrow(U_inv)))
#Mn <- tcrossprod(U_inv_inv)

# Posterior degrees of freedom
nu_n <- nu0 + nrow(X)

# Posterior scale matrix for covariance: Psi_n
# Expanded efficiently using crossprod to avoid standard transpositions
Psi_n <- S0 + crossprod(Y) + crossprod(B0, Lambda0 %*% B0) - crossprod(mn, post_mean)

# 2. Execute the mniw sampler
set.seed(42)
# rmniw draws from the MNIW distribution and returns a list containing:
# X (the Matrix-Normal draws) and V (the Inverse-Wishart draws)
# NOTE: In mniw, the row-covariance matrix is passed as Sigma, not Omega
mniw_draws <- rmniw(n = n_samples, Lambda = post_mean, Omega = Mn_inv, Psi = Psi_n, nu = nu_n)

# 3. Extract and format the samples
# mniw conveniently outputs 3D arrays natively (p x q x n)
beta_mniw_array <- mniw_draws$X
sigma_mniw_array <- mniw_draws$V  # FIXED: mniw returns covariance as V, not Sigma

# 4. Generate summaries using the helper function defined in Section 6
df_beta_mniw <- generate_summary_table(beta_mniw_array, "Beta (mniw)")
df_sigma_mniw <- generate_summary_table(sigma_mniw_array, "Sigma (mniw)")

# 5. Print the tables for comparison
cat("### `mniw` Posterior Summary for $\\beta$\n\n")
print(kable(df_beta_mniw, format = "markdown", align = "c", row.names = FALSE))

cat("\n### `mniw` Posterior Summary for $\\Sigma$\n\n")
print(kable(df_sigma_mniw, format = "markdown", align = "c", row.names = FALSE))


## ----qqplots_mniw, fig.width=10, fig.height=8, echo=FALSE, message=FALSE, warning=FALSE----
# --- QQ-Plots for Beta ---
# Set up plotting area dynamically based on matrix dimensions
par(mfrow = c(k_dim, p_dim), mar = c(4, 4, 2, 1))

invisible(lapply(1:nrow(grid_beta), function(idx) {
  i <- grid_beta$Row[idx]
  j <- grid_beta$Col[idx]
  qqplot(beta_samples[i, j, ], beta_mniw_array[i, j, ],
         main = paste0("QQ-Plot: Beta[", i, ",", j, "]"),
         xlab = "Custom Sampler Quantiles", 
         ylab = "mniw Quantiles",
         pch = 16, col = rgb(0.2, 0.6, 0.4, 0.5))
  abline(0, 1, col = "red", lwd = 2)
}))

# --- QQ-Plots for Sigma ---
# Set up plotting area for the p x p covariance matrix
par(mfrow = c(p_dim, p_dim), mar = c(4, 4, 2, 1))

invisible(lapply(1:nrow(grid_sigma), function(idx) {
  i <- grid_sigma$Row[idx]
  j <- grid_sigma$Col[idx]
  qqplot(sigma_samples[i, j, ], sigma_mniw_array[i, j, ],
         main = paste0("QQ-Plot: Sigma[", i, ",", j, "]"),
         xlab = "Custom Sampler Quantiles", 
         ylab = "mniw Quantiles",
         pch = 16, col = rgb(0.2, 0.6, 0.4, 0.5))
  abline(0, 1, col = "red", lwd = 2)
}))

# Reset graphics parameters
par(mfrow = c(1, 1))


## ----benchmark_execution, message=FALSE, warning=FALSE------------------------
# Load the microbenchmark package for precise timing
if (!requireNamespace("microbenchmark", quietly = TRUE)) {
  install.packages("microbenchmark")
}
library(microbenchmark)

# Define a smaller sample size for the benchmark to run quickly
bench_samples <- 500

# Run the benchmark
set.seed(42)
timing_results <- microbenchmark(
  "Custom (Method 1)" = {
    bayes_mvr(Y = Y, X = X, B0 = B0, Lambda0 = Lambda0, 
              nu0 = nu0, S0 = S0, method = 1, n_samples = bench_samples)
  },
  
  "Custom (Method 2)" = {
    bayes_mvr(Y = Y, X = X, B0 = B0, Lambda0 = Lambda0, 
              nu0 = nu0, S0 = S0, method = 2, n_samples = bench_samples)
  },
  
  "bayesm::rmultireg" = {
    lapply(1:bench_samples, function(x) {
      rmultireg(Y = Y, X = X, Bbar = B0, A = Lambda0, nu = nu0, V = S0)
    })
  },
  
  "mniw::rmniw" = {
    # Include the analytical parameter preparation time to be fair
    Mn_inv <- Lambda0 + crossprod(X)
    U_inv <- chol(Mn_inv) 
    mn <- Lambda0 %*% B0 + crossprod(X, Y)
    
    y_tmp <- backsolve(U_inv, mn, transpose = TRUE)
    post_mean <- backsolve(U_inv, y_tmp)
    
 #   U_inv_inv <- backsolve(U_inv, diag(nrow(U_inv)))
 #   Mn <- tcrossprod(U_inv_inv)
    nu_n <- nu0 + nrow(X)
    Psi_n <- S0 + crossprod(Y) + crossprod(B0, Lambda0 %*% B0) - crossprod(mn, post_mean)
    
    rmniw(n = bench_samples, Lambda = post_mean, Omega = Mn_inv, Psi = Psi_n, nu = nu_n)
  },
  times = 50 # Run each method 50 times for a stable average
)


## ----benchmark_plot, echo=FALSE, fig.width=9, fig.height=6--------------------
# Print the textual summary table
print(timing_results)

# Generate a traditional base R boxplot of the computational times
boxplot(
  timing_results, 
  unit = "ms", 
  main = "Execution Time Comparison",
  ylab = "Time (milliseconds)",
  col = c("lightblue", "lightgreen", "lightpink", "wheat"),
  las = 1,       # Keep y-axis labels horizontal
  cex.axis = 0.9 # Slightly scale down axis text to fit labels
)


## ----prep_predictive_data-----------------------------------------------------
# Set a seed for reproducible row sampling
set.seed(8675309)

# Define the number of new predictions to generate
m <- 20

# Randomly sample row indices from the raw predictor matrix
sample_indices <- sample(1:nrow(X_raw), m)

# Extract the corresponding rows to form the new raw predictors
X_raw_tilde <- X_raw[sample_indices, , drop = FALSE]

# Prepend the intercept to match the trained model structure
X_tilde <- cbind(1, X_raw_tilde)


## ----custom_predictive_sampler------------------------------------------------
# Ensure the number of predictive samples matches our parameter draws
N_draws <- dim(beta_samples)[3]
q_cols <- dim(beta_samples)[2]

set.seed(123)
# We use lapply to vectorize the sampling over all posterior draws,
# completely avoiding explicit for-loops.
Y_tilde_custom_list <- lapply(1:N_draws, function(i) {
  # Extract the i-th posterior parameter realizations
  B_i <- beta_samples[, , i]
  Sigma_i <- sigma_samples[, , i]

  # Compute the predicted mean location matrix
  M_i <- X_tilde %*% B_i

  # Draw the standard normal anchor matrix Z
  Z <- matrix(rnorm(m * q_cols), nrow = m, ncol = q_cols)

  # Calculate the upper Cholesky factor of the column covariance
  U_Sigma <- chol(Sigma_i)

  # Map to the Matrix-Normal likelihood and return
  # Z %*% U_Sigma mathematically equals Z %*% t(L_Sigma)
  Y_draw <- M_i + Z %*% U_Sigma
  return(Y_draw)
})

# Collapse the list of matrices into a 3D array (m x q x N_draws)
Y_tilde_custom <- simplify2array(Y_tilde_custom_list)


## ----mniw_predictive_sampler, message=FALSE, warning=FALSE--------------------
set.seed(123)

# Vectorize the predictive draws using the mniw parameter arrays
Y_tilde_mniw_list <- lapply(1:N_draws, function(i) {
  # Extract the mniw posterior parameter realizations
  B_i <- beta_mniw_array[, , i]
  Sigma_i <- sigma_mniw_array[, , i]

  # Compute the predicted mean location matrix
  M_i <- X_tilde %*% B_i

  # Draw directly from the Matrix-Normal distribution using mniw
  # rMNorm returns an array, so we use drop() to cast it to a matrix
  Y_draw <- rMNorm(n = 1, Lambda = M_i, SigmaR = diag(m), SigmaC = Sigma_i)
  return(drop(Y_draw))
})

# Collapse into a 3D array
Y_tilde_mniw <- simplify2array(Y_tilde_mniw_list)


## ----verify_predictive, echo=FALSE, fig.width=6, fig.height=5-----------------
# Extract the predictive distributions for the first row, first column response
pred_custom_11 <- Y_tilde_custom[1, 1, ]
pred_mniw_11 <- Y_tilde_mniw[1, 1, ]

# Print quick comparative summary statistics
cat("**Predictive Summary for Y_tilde[1,1] (Custom Sampler):**\n")
print(summary(pred_custom_11))

cat("\n**Predictive Summary for Y_tilde[1,1] (`mniw` Sampler):**\n")
print(summary(pred_mniw_11))

# Generate a QQ-Plot to visually verify distributional alignment
qqplot(
  pred_custom_11,
  pred_mniw_11,
  main = expression(paste("Posterior Predictive QQ-Plot for ", tilde(Y)[1:1])),
  xlab = "Custom Sampler Quantiles",
  ylab = "mniw Quantiles",
  pch = 16,
  col = rgb(0.6, 0.2, 0.4, 0.5)
)
abline(0, 1, col = "red", lwd = 2)


## ----predictive_diagnostics_exact, fig.width=11, fig.height=5.5---------------
# Extract (m x N_draws) matrix for primary outcome and hold-out observations
Y_tilde <- Y_tilde_custom[, 1, ] # Rows = observations, Columns = draws
Y_true  <- Y[sample_indices, 1]

# 1. Compute posterior predictive metrics
pred_mean <- rowMeans(Y_tilde)
PI_lower  <- apply(Y_tilde, 1, quantile, probs = 0.025)
PI_upper  <- apply(Y_tilde, 1, quantile, probs = 0.975)

# 2. Calculate Empirical Coverage Probability
covered <- (Y_true >= PI_lower) & (Y_true <= PI_upper)
coverage_prob <- mean(covered)

# Print explicit summary metric to markdown output
cat(sprintf("**Empirical 95%% Coverage Probability:** %.2f%% (%d / %d observations covered)\n", 
            coverage_prob * 100, sum(covered), length(Y_true)))

# Define color schemes for covered vs uncovered observations
pt_cols  <- ifelse(covered, "black", "firebrick")
bar_cols <- ifelse(covered, rgb(0.3, 0.6, 0.85, 0.6), rgb(1, 0.4, 0.4, 0.7))

# Set up side-by-side diagnostic layout
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))

# ==============================================================================
# Plot 1: 95% Posterior Predictive Intervals vs Hold-out Y (Sorted by True Y)
# ==============================================================================
sort_idx     <- order(Y_true)
Y_true_sort  <- Y_true[sort_idx]
pred_sort    <- pred_mean[sort_idx]
lower_sort   <- PI_lower[sort_idx]
upper_sort   <- PI_upper[sort_idx]
pt_cols_sort <- pt_cols[sort_idx]
bar_cols_sort<- bar_cols[sort_idx]

y_lims <- range(c(lower_sort, upper_sort, Y_true_sort))

plot(1:length(Y_true), Y_true_sort, type = "n", ylim = y_lims,
     xlab = "Hold-out Observation Index (Sorted by True Y)", 
     ylab = "Target Variable Y",
     main = sprintf("95%% Predictive Intervals (Coverage: %.1f%%)", coverage_prob * 100))

# Draw predictive intervals, means, and true values
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

# ==============================================================================
# Plot 2: Observed vs. Predicted Y with 95% PIs
# ==============================================================================
plot(Y_true, pred_mean, type = "n", ylim = y_lims, xlim = y_lims,
     xlab = "Observed Hold-out Y", 
     ylab = "Predicted Y (Posterior Predictive Mean)",
     main = "Observed vs. Predicted Y with 95% PIs")

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

# Reset graphics layout
par(mfrow = c(1, 1))

