params <-
list(sampler_option = 2L)

## ----analytical_waic, warning=FALSE, message=FALSE----------------------------
# Source the optimized WAIC function
source("../src/mniwModelComparisons.r")

# 1. Ingest the Data
Y <- as.matrix(read.table("../data/Y.txt", header = FALSE))
X_raw <- as.matrix(read.table("../data/X.txt", header = FALSE))

# Prepend intercept to design matrix
X <- cbind(1, X_raw)

# Extract dimensions
n_obs <- nrow(Y)
q_cols <- ncol(Y)
p_cols <- ncol(X)

# 2. Specify Conjugate Priors
M0_inv <- diag(1, p_cols)
m0 <- matrix(0, nrow = p_cols, ncol = q_cols)
nu0 <- q_cols + 2
Psi0 <- diag(1, q_cols)

# 3. Compute Posterior Parameters (Optimized Cholesky Solvers)
Mn_inv <- M0_inv + crossprod(X)
U_Mn_inv <- chol(Mn_inv)

mn_unscaled <- M0_inv %*% m0 + crossprod(X, Y)

tmp_mn <- backsolve(U_Mn_inv, mn_unscaled, transpose = TRUE)
mn <- backsolve(U_Mn_inv, tmp_mn)

nu_n <- nu0 + n_obs

Psi_n <- Psi0 + crossprod(Y) + crossprod(m0, M0_inv %*% m0) - crossprod(tmp_mn)

U_Mn_inv_inv <- backsolve(U_Mn_inv, diag(p_cols))
Mn <- tcrossprod(U_Mn_inv_inv)

# 4. Execute the analytical WAIC calculation
waic_analytical <- compute_analytical_waic(Y, X, Mn, mn, Psi_n, nu_n)


## ----sample_based_waic, warning=FALSE, message=FALSE--------------------------
# Load required package and source external sampling scripts
library(mniw)
source("../src/mniwPosteriorPredictive.r")

# Set seed for reproducibility
set.seed(42)
S_samples <- 5000

# Select sampling method based on parameter option
if (params$sampler_option == 1) {
  cat("Evaluating using Option #1 (mniw package)\n")
  samples <- draw_mniw_samples_pkg(S = S_samples, mn = mn, Mn_inv = Mn_inv, Psi_n = Psi_n, nu_n = nu_n)
} else if (params$sampler_option == 2) {
  cat("Evaluating using Option #2 (Custom Exact Sampler)\n")
  samples <- draw_mniw_samples_custom(S = S_samples, mn = mn, Mn = Mn, Psi_n = Psi_n, nu_n = nu_n)
} else {
  stop("Invalid sampler_option. Please choose either 1 or 2.")
}

# Compute WAIC using custom sampling function
waic_sample <- compute_sample_waic(Y, X, samples$X, samples$Sigma)


## ----loo_waic, warning=FALSE, message=FALSE-----------------------------------
library(loo)

# 1. Build the S x N pointwise log-likelihood matrix across posterior draws
log_lik_mat <- do.call(rbind, lapply(seq_len(S_samples), function(s) {
  beta_s <- samples$X[,,s]
  sigma_s <- samples$Sigma[,,s]
  
  inv_sigma_s <- solve(sigma_s)
  log_det_s <- as.numeric(determinant(sigma_s, logarithm = TRUE)$modulus)
  
  diff <- Y - X %*% beta_s
  quad_form <- rowSums((diff %*% inv_sigma_s) * diff)
  
  -0.5 * q_cols * log(2 * pi) - 0.5 * log_det_s - 0.5 * quad_form
}))

# 2. Evaluate WAIC using loo::waic()
loo_obj <- loo::waic(log_lik_mat)

# Extract elpd, p_waic, and calculate corresponding lppd
waic_loo <- list(
  WAIC = loo_obj$estimates["waic", "Estimate"],
  lppd = loo_obj$estimates["elpd_waic", "Estimate"] + loo_obj$estimates["p_waic", "Estimate"],
  p_waic = loo_obj$estimates["p_waic", "Estimate"]
)


## ----comparison_summary, warning=FALSE, message=FALSE-------------------------
library(knitr)

# -----------------------------------------------------------------------------
# Print Comparative Performance Summary
# -----------------------------------------------------------------------------
comparison_df <- data.frame(
  Metric = c("lppd (Log Predictive Density)", "p_WAIC (Effective Parameters)", "WAIC (Deviance Scale)"),
  `Exact_Analytical` = c(waic_analytical$lppd, waic_analytical$p_waic, waic_analytical$WAIC),
  `Custom_Sampled`   = c(waic_sample$lppd, waic_sample$p_waic, waic_sample$WAIC),
  `LOO_Package`      = c(waic_loo$lppd, waic_loo$p_waic, waic_loo$WAIC),
  `Abs_Diff_(Exact_vs_Sample)` = c(
    abs(waic_analytical$lppd - waic_sample$lppd),
    abs(waic_analytical$p_waic - waic_sample$p_waic),
    abs(waic_analytical$WAIC - waic_sample$WAIC)
  ),
  check.names = FALSE
)

# Display styled summary table
kable(
  comparison_df, 
  digits = 4,
  col.names = c("Metric", "Exact Analytical", "Custom Sampled", "LOO Package", "Abs Diff (Exact vs Sample)"),
  caption = "Comparative WAIC Metrics across Evaluation Methods"
)


## ----waic-sequential-eval, message=FALSE, warning=FALSE, fig.width=10, fig.height=7----
library(ggplot2)
library(tidyr)
library(knitr)

# Source the analytical WAIC computation script
source("../src/mniwModelComparisons.r")

# Extract matrix dimensions
n_obs  <- nrow(Y)
q_cols <- ncol(Y)
p_cols <- ncol(X)

# Assign standard prior matrices if not present in active environment
M0_inv <- if (exists("M0_inv")) M0_inv else diag(1, p_cols)
m0     <- if (exists("m0")) m0 else matrix(0, nrow = p_cols, ncol = q_cols)
nu0    <- if (exists("nu0")) nu0 else q_cols + 2
Psi0   <- if (exists("Psi0")) Psi0 else diag(1, q_cols)

# Sequentially evaluate models via lapply without solve(), t(), or for loops
waic_seq_results <- do.call(rbind, lapply(1:p_cols, function(k) {
  # Sub-matrices for column subset k
  X_k      <- X[, 1:k, drop = FALSE]
  M0_inv_k <- M0_inv[1:k, 1:k, drop = FALSE]
  m0_k     <- m0[1:k, , drop = FALSE]
  
  # Posterior parameters derived using Cholesky solvers (no solve, no t)
  Mn_inv_k   <- M0_inv_k + crossprod(X_k)
  U_Mn_inv_k <- chol(Mn_inv_k)
  
  mn_unscaled_k <- M0_inv_k %*% m0_k + crossprod(X_k, Y)
  tmp_mn_k      <- backsolve(U_Mn_inv_k, mn_unscaled_k, transpose = TRUE)
  mn_k          <- backsolve(U_Mn_inv_k, tmp_mn_k)
  
  nu_n_k  <- nu0 + n_obs
  Psi_n_k <- Psi0 + crossprod(Y) + crossprod(m0_k, M0_inv_k %*% m0_k) - crossprod(tmp_mn_k)
  
  U_Mn_inv_inv_k <- backsolve(U_Mn_inv_k, diag(k))
  Mn_k           <- tcrossprod(U_Mn_inv_inv_k)
  
  # Evaluate analytical WAIC directly using compute_analytical_waic
  waic_k <- compute_analytical_waic(
    Y     = Y, 
    X     = X_k, 
    Mn    = Mn_k, 
    mn    = mn_k, 
    Psi_n = Psi_n_k, 
    nu_n  = nu_n_k
  )
  
  # Model sequence label: M0 (Intercept), M1 (+X1), M2 (+X2), ...
  model_label <- if (k == 1) "M0 (Intercept)" else paste0("M", k - 1, " (+X", k - 1, ")")
  
  data.frame(
    k          = k,
    Model      = model_label,
    lppd       = waic_k$lppd,
    p_waic     = waic_k$p_waic,
    waic       = waic_k$WAIC,
    stringsAsFactors = FALSE
  )
}))

# Display summary comparison table
kable(
  waic_seq_results[, c("Model", "lppd", "p_waic", "waic")], 
  digits = 3,
  col.names = c("Model Sequence", "lppd", "p_WAIC", "WAIC"),
  caption = "MNIW Model Criteria vs. Number of Included Covariates"
)

# Reshape data and format factors for stacked panel display
waic_long <- pivot_longer(
  waic_seq_results, 
  cols = c("waic", "lppd", "p_waic"), 
  names_to = "Metric", 
  values_to = "Value"
)

# Set factor ordering for x-axis models
waic_long$Model <- factor(waic_long$Model, levels = waic_seq_results$Model)

# Set factor ordering and labels for stacked panels
waic_long$Metric_Label <- factor(
  waic_long$Metric,
  levels = c("waic", "lppd", "p_waic"),
  labels = c("WAIC (Lower is Better)", "lppd (Higher Fit is Better)", "p_WAIC (Penalty)")
)

# Match exact color palette from NIG profile
metric_colors <- c(
  "WAIC (Lower is Better)"     = "#c85200",
  "lppd (Higher Fit is Better)" = "#2b5c8f",
  "p_WAIC (Penalty)"           = "#7570b3"
)

# Generate stacked profile plot
ggplot(waic_long, aes(x = Model, y = Value, color = Metric_Label, group = Metric_Label)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 3) +
  facet_wrap(~ Metric_Label, ncol = 1, scales = "free_y") +
  scale_color_manual(values = metric_colors) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "grey30", size = 11, margin = margin(b = 10)),
    strip.text = element_text(face = "bold", size = 11),
    axis.text.x = element_text(angle = 30, hjust = 1, color = "black"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey92")
  ) +
  labs(
    title = "Sequential Model Comparison Profile",
    subtitle = "Tracking Model Fit (lppd), Complexity Penalty (p_WAIC), and Overall Criterion (WAIC)",
    x = "Model Sequence",
    y = "Model Value"
  )

