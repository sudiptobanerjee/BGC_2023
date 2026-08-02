# -----------------------------------------------------------------------------
# 1. Source Optimized Samplers
# -----------------------------------------------------------------------------
source("../src/rawSampler.r")
source("../src/rawPosteriorPredictive.r")

# -----------------------------------------------------------------------------
# 2. Load Data and Partition into Training (80%) and Hold-out/Test (20%)
# -----------------------------------------------------------------------------
dataset <- read.table("../../data/data.txt", header = TRUE)

set.seed(42)
n_total   <- nrow(dataset)
train_idx <- sample.int(n_total, size = floor(0.8 * n_total))

train_data <- dataset[train_idx, ]
test_data  <- dataset[-train_idx, ]

# -----------------------------------------------------------------------------
# 3. Train Covariate Matrix & Scaling Attributes
# -----------------------------------------------------------------------------
y_train     <- train_data[, 1]
X_train_raw <- as.matrix(train_data[, -1])

# Center and scale training predictors, saving transformation parameters
X_train_scaled <- scale(X_train_raw)
train_center   <- attr(X_train_scaled, "scaled:center")
train_scale    <- attr(X_train_scaled, "scaled:scale")

X_train   <- cbind(Intercept = 1, X_train_scaled)
p         <- ncol(X_train)
N_samples <- 5000

# -----------------------------------------------------------------------------
# 4. Fit Model Posterior via rawSampler
# -----------------------------------------------------------------------------
posterior <- rawSampler(
  y = y_train, 
  X = X_train, 
  m0 = matrix(0, p, 1), 
  M0_inv = matrix(0, p, p), 
  a0 = -p/2, 
  b0 = 0, 
  N = N_samples
)

beta_draws   <- posterior$draws[, 1:p]
sigma2_draws <- posterior$draws[, "sigma2"]

# -----------------------------------------------------------------------------
# 5. Prepare Hold-out Test Matrix (Re-using Training Scale Parameters)
# -----------------------------------------------------------------------------
X_test_raw    <- as.matrix(test_data[, -1])
X_test_scaled <- scale(X_test_raw, center = train_center, scale = train_scale)
X_test        <- cbind(Intercept = 1, X_test_scaled)
y_test_actual <- test_data[, 1]

# -----------------------------------------------------------------------------
# 6. Generate Posterior Predictive Draws (High Performance via tcrossprod)
# -----------------------------------------------------------------------------
Y_tilde_draws <- rawPosteriorPredictive(
  X_new = X_test, 
  beta_draws = beta_draws, 
  sigma2_draws = sigma2_draws
)

# -----------------------------------------------------------------------------
# 7. Summarise Predictive Draws & Calculate Performance Metrics
# -----------------------------------------------------------------------------

# Point predictions (Posterior Predictive Mean)
y_pred_mean <- rowMeans(Y_tilde_draws)

# 95% Posterior Predictive Intervals (2.5% and 97.5% quantiles)
y_pred_pi   <- t(apply(Y_tilde_draws, 1, quantile, probs = c(0.025, 0.975)))

# Coverage Check: Determine if actual hold-out Y falls within the 95% PI
in_interval   <- (y_test_actual >= y_pred_pi[, 1]) & (y_test_actual <= y_pred_pi[, 2])
coverage_rate <- mean(in_interval) * 100

# Error metrics
rmspe <- sqrt(mean((y_test_actual - y_pred_mean)^2))
mae   <- mean(abs(y_test_actual - y_pred_mean))

cat("====================================================\n")
cat("       POSTERIOR PREDICTIVE PERFORMANCE             \n")
cat("====================================================\n")
cat(sprintf("Hold-out Sample Size (m)    : %d\n", length(y_test_actual)))
cat(sprintf("RMSPE                      : %.4f\n", rmspe))
cat(sprintf("MAE                        : %.4f\n", mae))
cat(sprintf("Empirical 95%% Coverage Rate: %.2f%%\n", coverage_rate))
cat("====================================================\n")

# -----------------------------------------------------------------------------
# 8. Diagnostic Plots: 95% Predictive Intervals vs Hold-out Data
# -----------------------------------------------------------------------------

par(mfrow = c(1, 2), mar = c(5, 4.5, 4, 2) + 0.1)

# --- PLOT 1: Sorted Hold-out Observations vs 95% Predictive Intervals ---
sort_idx        <- order(y_test_actual)
y_test_sorted   <- y_test_actual[sort_idx]
y_mean_sorted   <- y_pred_mean[sort_idx]
pi_lower_sorted <- y_pred_pi[sort_idx, 1]
pi_upper_sorted <- y_pred_pi[sort_idx, 2]
in_cov_sorted   <- in_interval[sort_idx]

m_test <- length(y_test_actual)

plot(1:m_test, y_test_sorted, type = "n",
     ylim = range(c(pi_lower_sorted, pi_upper_sorted, y_test_sorted)),
     xlab = "Hold-out Observation Index (Sorted by True Y)",
     ylab = "Target Variable Y",
     main = "95% Posterior Predictive Intervals vs Hold-out Y",
     cex.lab = 1.1, cex.main = 1.2)

# Draw 95% predictive interval segments (Blue = Covered, Red = Uncovered)
segments(1:m_test, pi_lower_sorted, 1:m_test, pi_upper_sorted, 
         col = ifelse(in_cov_sorted, rgb(0.2, 0.5, 0.8, 0.6), rgb(0.9, 0.2, 0.2, 0.8)), 
         lwd = 1.8)

# Add posterior predictive mean
points(1:m_test, y_mean_sorted, pch = 18, col = "navy", cex = 1.1)

# Add true actual hold-out y values
points(1:m_test, y_test_sorted, 
       pch = 16, 
       col = ifelse(in_cov_sorted, "black", "firebrick3"), 
       cex = 1.1)

legend("topleft", 
       legend = c("Actual Y (Covered)", "Actual Y (Uncovered)", 
                  "Predictive Mean", "95% Predictive Interval"),
       col = c("black", "firebrick3", "navy", rgb(0.2, 0.5, 0.8, 0.8)),
       pch = c(16, 16, 18, NA), 
       lty = c(NA, NA, NA, 1), 
       lwd = c(NA, NA, NA, 2), 
       bty = "n", cex = 0.9)

# --- PLOT 2: Observed vs. Predicted Scatter Plot ---
plot(y_test_actual, y_pred_mean, 
     xlim = range(c(y_test_actual, y_pred_pi)),
     ylim = range(c(y_test_actual, y_pred_pi)),
     xlab = "Observed Hold-out Y", 
     ylab = "Predicted Y (Posterior Predictive Mean)",
     main = "Observed vs. Predicted Y with 95% PIs",
     pch = 16, 
     col = ifelse(in_interval, rgb(0.1, 0.1, 0.1, 0.7), "firebrick3"), 
     cex = 1.2, cex.lab = 1.1, cex.main = 1.2)

# Draw vertical error bars representing the 95% predictive interval
segments(y_test_actual, y_pred_pi[, 1], y_test_actual, y_pred_pi[, 2], 
         col = ifelse(in_interval, rgb(0.2, 0.5, 0.8, 0.35), rgb(0.9, 0.2, 0.2, 0.5)), 
         lwd = 1.2)

# 45-degree reference line (Ideal perfect prediction)
abline(0, 1, col = "red", lty = 2, lwd = 2)

legend("topleft", 
       legend = c("1:1 Reference Line", "Covered Observation", "Uncovered Observation"),
       col = c("red", "black", "firebrick3"), 
       lty = c(2, NA, NA), 
       pch = c(NA, 16, 16), 
       bty = "n", cex = 0.9)

# Reset graphics layout
par(mfrow = c(1, 1))
