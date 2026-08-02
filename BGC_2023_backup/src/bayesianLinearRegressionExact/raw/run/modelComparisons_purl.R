# -----------------------------------------------------------------------------
# 1. Source Libraries & Routines
# -----------------------------------------------------------------------------
source("../src/rawSampler.r")
source("../src/modelComparisons.r")

# -----------------------------------------------------------------------------
# 2. Load Data and Setup Model Matrices
# -----------------------------------------------------------------------------
dataset <- read.table("../../data/data.txt", header = TRUE)

y_train <- dataset[, 1]
X_train_raw <- as.matrix(dataset[, -1])

X_train <- cbind(Intercept = 1, scale(X_train_raw))
n <- nrow(X_train)
p <- ncol(X_train)

# Prior specifications
m0     <- matrix(0, p, 1)
M0_inv <- matrix(0, p, p)  # Non-informative prior precision
a0     <- -p / 2
b0     <- 0
N_draws <- 10000

# -----------------------------------------------------------------------------
# 3. Fit Model & Generate Posterior Draws
# -----------------------------------------------------------------------------
set.seed(42)
posterior <- rawSampler(
  y = y_train, 
  X = X_train, 
  m0 = m0, 
  M0_inv = M0_inv, 
  a0 = a0, 
  b0 = b0, 
  N = N_draws
)

beta_draws   <- posterior$draws[, 1:p]
sigma2_draws <- posterior$draws[, "sigma2"]

# -----------------------------------------------------------------------------
# 4. Method (i): Exact Analytical WAIC
# -----------------------------------------------------------------------------
res_exact <- waic_exact(
  X = X_train, 
  y = y_train, 
  m0 = m0, 
  M0_inv = M0_inv, 
  a0 = a0, 
  b0 = b0
)

# -----------------------------------------------------------------------------
# 5. Method (ii): Raw Sample-Based WAIC (with Log-Sum-Exp Trick)
# -----------------------------------------------------------------------------
LL_matrix <- compute_log_likelihood_matrix(
  X = X_train, 
  y = y_train, 
  beta_draws = beta_draws, 
  sigma2_draws = sigma2_draws
)

res_sample <- waic_sample(LL_matrix)

# -----------------------------------------------------------------------------
# 6. Method (iii): 'loo' Package WAIC
# -----------------------------------------------------------------------------
res_loo <- waic_loo_package(LL_matrix)

# -----------------------------------------------------------------------------
# 7. Print Comparative Performance Summary
# -----------------------------------------------------------------------------
comparison_df <- data.frame(
  Metric = c("lppd (Log Predictive Density)", "p_WAIC (Effective Parameters)", "WAIC (Deviance Scale)"),
  `Exact_Analytical` = c(res_exact$lppd, res_exact$p_waic, res_exact$waic),
  `Raw_Sample_MCMC`  = c(res_sample$lppd, res_sample$p_waic, res_sample$waic),
  `LOO_Package`       = c(res_loo$lppd, res_loo$p_waic, res_loo$waic),
  `Abs_Diff_(Exact_vs_Sample)` = c(
    abs(res_exact$lppd - res_sample$lppd),
    abs(res_exact$p_waic - res_sample$p_waic),
    abs(res_exact$waic - res_sample$waic)
  )
)

knitr::kable(
  comparison_df, 
  digits = 4, 
  col.names = c("Metric", "Exact Analytical", "Raw Sample (Log-Sum-Exp)", "loo Package", "Abs Diff (Exact vs Raw)"),
  caption = "Comparison of WAIC Calculations across Analytical, Sample-Based, and 'loo' Methods"
)

library(ggplot2)
library(dplyr)
library(tidyr)

# -----------------------------------------------------------------------------
# 1. Fit Sequential Models M_0 to M_p
# -----------------------------------------------------------------------------
num_features <- ncol(X_train_raw)
model_results <- vector("list", num_features + 1)

for (k in 0:num_features) {
  if (k == 0) {
    X_k <- X_train[, 1, drop = FALSE] # Intercept only
    model_label <- "M0 (Intercept)"
  } else {
    X_k <- X_train[, 1:(k + 1), drop = FALSE] # Intercept + X1 + ... + Xk
    model_label <- paste0("M", k, " (+X", k, ")")
  }
  
  p_k <- ncol(X_k)
  
  # Hyperparameters for subset model k
  m0_k     <- matrix(0, p_k, 1)
  M0_inv_k <- matrix(0, p_k, p_k)
  a0_k     <- -p_k / 2
  b0_k     <- 0
  
  # Exact WAIC calculation
  res_k <- waic_exact(
    X = X_k, 
    y = y_train, 
    m0 = m0_k, 
    M0_inv = M0_inv_k, 
    a0 = a0_k, 
    b0 = b0_k
  )
  
  model_results[[k + 1]] <- data.frame(
    Model_Index = k,
    Model       = model_label,
    Num_Params  = p_k + 1,  # p_k regression coefficients + 1 variance
    lppd        = res_k$lppd,
    p_waic      = res_k$p_waic,
    WAIC        = res_k$waic
  )
}

seq_df <- bind_rows(model_results)

# -----------------------------------------------------------------------------
# 2. Display Model Assessment Table
# -----------------------------------------------------------------------------
best_model_idx <- which.min(seq_df$WAIC)

knitr::kable(
  seq_df,
  digits = 3,
  col.names = c("Index", "Model Structure", "True Params (p+1)", "lppd (Fit)", "p_WAIC (Penalty)", "WAIC (Deviance)"),
  caption = "Sequential Bayesian Model Assessment (Closed-Form Exact WAIC)"
)

# -----------------------------------------------------------------------------
# 3. Plot A: Faceted Metrics Overview
# -----------------------------------------------------------------------------
df_long <- seq_df %>%
  pivot_longer(
    cols = c("lppd", "p_waic", "WAIC"),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Metric_Label = factor(
      Metric,
      levels = c("WAIC", "lppd", "p_waic"),
      labels = c("WAIC (Lower is Better)", "lppd (Higher Fit is Better)", "p_WAIC (Penalty)")
    ),
    Model = factor(Model, levels = seq_df$Model)
  )

p1 <- ggplot(df_long, aes(x = Model, y = Value, group = Metric_Label, color = Metric_Label)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 3) +
  facet_wrap(~ Metric_Label, scales = "free_y", ncol = 1) +
  scale_color_manual(values = c("WAIC (Lower is Better)" = "#d95f02", 
                               "lppd (Higher Fit is Better)" = "#2b5c8f", 
                               "p_WAIC (Penalty)" = "#7570b3")) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 11),
    axis.text.x = element_text(angle = 35, hjust = 1, color = "black"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 14)
  ) +
  labs(
    title = "Sequential Model Comparison Profile",
    subtitle = "Tracking Model Fit (lppd), Complexity Penalty (p_WAIC), and Overall Criterion (WAIC)",
    x = "Model Sequence",
    y = "Metric Value"
  )

print(p1)

# -----------------------------------------------------------------------------
# 4. Plot B: Focused WAIC Trajectory Plot with Minimum Highlight
# -----------------------------------------------------------------------------
p2 <- ggplot(seq_df, aes(x = factor(Model, levels = Model), y = WAIC, group = 1)) +
  geom_line(color = "#2b5c8f", linewidth = 1.2) +
  geom_point(color = "#2b5c8f", size = 3.5) +
  # Highlight optimal model point
  geom_point(
    data = seq_df[best_model_idx, ],
    aes(x = factor(Model, levels = Model), y = WAIC),
    color = "#d95f02", size = 5
  ) +
  annotate(
    "text",
    x = best_model_idx,
    y = seq_df$WAIC[best_model_idx] + (max(seq_df$WAIC) - min(seq_df$WAIC)) * 0.06,
    label = paste0("Optimal Model: ", seq_df$Model[best_model_idx]),
    color = "#d95f02", fontface = "bold", size = 3.8
  ) +
  # Add vertical expansion & prevent text clipping
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.18))) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, color = "black"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    plot.margin = margin(t = 15, r = 20, b = 10, l = 10)
  ) +
  labs(
    title = "Out-of-Sample Predictive Trajectory (WAIC)",
    subtitle = "Orange highlight indicates the model minimizing out-of-sample deviance",
    x = "Model Sub-Structure",
    y = "WAIC Score"
  )

print(p2)
