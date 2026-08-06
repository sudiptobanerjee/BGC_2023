## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)
library(dplyr)
library(knitr)
library(ggplot2)
library(tidyr)

# Source the custom divide and conquer script
source("../src/rawSamplerDivideAndConquer.r")


## ----partition_data-----------------------------------------------------------
# 1. Load Data
data <- read.table("../data/data.txt", header = TRUE)

# Prepare the Model Matrix and Response Vector
X_full <- as.matrix(cbind(Intercept = 1, data[, c("X1", "X2", "X3", "X4", "X5", "X6", "Lat", "Longi")]))
y_full <- as.numeric(data$Y)

# 2. Partition into Blocks using vectorized operations
T_blocks <- 10
block_size <- 130
total_n <- nrow(X_full)

# Create a grouping vector for splitting
group_indices <- c(
  rep(1:(T_blocks - 1), each = block_size), 
  rep(T_blocks, total_n - (T_blocks - 1) * block_size)
)

# Split data into lists efficiently
X_list <- lapply(split(as.data.frame(X_full), group_indices), as.matrix)
y_list <- unname(split(y_full, group_indices))

cat("Successfully created", length(X_list), "data blocks.\n")
cat("Rows in block 1-9:", nrow(X_list[[1]]), "| Rows in block 10:", nrow(X_list[[10]]), "\n")


## ----seq_update---------------------------------------------------------------
# Setup Priors (Weakly Informative)
p <- ncol(X_full)
M0_inv <- matrix(0, nrow = p, ncol = p)
m0 <- rep(0, p)
a0 <- 0.01
b0 <- 0.01
N_samples <- 5000

# Run Divide and Conquer Sampler (T = 10)
set.seed(123)
res_T10 <- rawSamplerDivideAndConquer(
  X_list = X_list, 
  y_list = y_list, 
  M0_inv = M0_inv, m0 = m0, a0 = a0, b0 = b0, 
  N = N_samples
)

# Summarize T=10
draws_T10 <- as.data.frame(res_T10$draws)
summary_T10 <- data.frame(
  Parameter = colnames(draws_T10),
  Mean_T10 = colMeans(draws_T10),
  SD_T10 = apply(draws_T10, 2, sd)
)
kable(summary_T10, digits = 4, caption = "Posterior Summaries (Divide and Conquer, T=10)")


## ----exact_comparison---------------------------------------------------------
# Run Single Block Sampler (T = 1)
set.seed(123) # Use same seed to demonstrate distributional equivalence
res_T1 <- rawSamplerDivideAndConquer(
  X_list = list(X_full),
  y_list = list(y_full),
  M0_inv = M0_inv, m0 = m0, a0 = a0, b0 = b0,
  N = N_samples
)

draws_T1 <- as.data.frame(res_T1$draws)
summary_T1 <- data.frame(
  Parameter = colnames(draws_T1),
  Mean_T1 = colMeans(draws_T1),
  SD_T1 = apply(draws_T1, 2, sd)
)

# Merge summaries
comparison_df <- merge(summary_T1, summary_T10, by = "Parameter")
kable(comparison_df, digits = 4, caption = "Posterior Comparison: T=1 vs T=10")


## ----qq_plots, fig.width=12, fig.height=10------------------------------------
params_to_plot <- c("Intercept", "X1", "X2", "X3", "X4", "X5", "X6", "Lat", "Longi", "sigma2")

# Vectorized alignment and R-squared calculation
plot_data_list <- lapply(params_to_plot, function(param) {
  q_T1 <- sort(draws_T1[[param]])
  q_T10 <- sort(draws_T10[[param]])
  r_squared <- cor(q_T1, q_T10)^2

  data.frame(
    Facet_Label = paste0(param, " (R-squared = ", round(r_squared, 6), ")"),
    Exact_T1 = q_T1,
    Sequential_T10 = q_T10
  )
})

df_plot_all <- do.call(rbind, plot_data_list)

# Convert Facet_Label to a factor to enforce the specific plotting order
df_plot_all$Facet_Label <- factor(df_plot_all$Facet_Label, levels = unique(df_plot_all$Facet_Label))

ggplot(df_plot_all, aes(x = Exact_T1, y = Sequential_T10)) +
  geom_point(alpha = 0.3, color = "steelblue") +
  geom_abline(slope = 1, intercept = 0, color = "darkred", linetype = "dashed", linewidth = 1) +
  facet_wrap(~ Facet_Label, scales = "free", ncol = 3) +
  labs(
    title = "Q-Q Plots: Exact Posterior (T=1) vs Sequential (T=10)",
    x = "Exact Posterior (T = 1)",
    y = "Sequential Posterior (T = 10)"
  ) +
  theme_bw()


## ----scaled_analysis----------------------------------------------------------
# Scale covariates (exclude intercept to prevent variance collapse)
covariates <- data[, c("X1", "X2", "X3", "X4", "X5", "X6", "Lat", "Longi")]
X_scaled <- as.matrix(cbind(Intercept = 1, scale(covariates)))

# Split the scaled data efficiently
X_list_scaled <- lapply(split(as.data.frame(X_scaled), group_indices), as.matrix)

# Run Divide and Conquer Sampler on Scaled Data (T = 10)
set.seed(123)
res_scaled_T10 <- rawSamplerDivideAndConquer(
  X_list = X_list_scaled, 
  y_list = y_list, 
  M0_inv = M0_inv, m0 = m0, a0 = a0, b0 = b0, 
  N = N_samples
)

# Summarize Scaled T=10
draws_scaled_T10 <- as.data.frame(res_scaled_T10$draws)
summary_scaled_T10 <- data.frame(
  Parameter = colnames(draws_scaled_T10),
  Mean_Scaled_T10 = colMeans(draws_scaled_T10),
  SD_Scaled_T10 = apply(draws_scaled_T10, 2, sd)
)

# Run Single Block Sampler on Scaled Data (T = 1)
set.seed(123)
res_scaled_T1 <- rawSamplerDivideAndConquer(
  X_list = list(X_scaled), 
  y_list = list(y_full), 
  M0_inv = M0_inv, m0 = m0, a0 = a0, b0 = b0, 
  N = N_samples
)

# Summarize Scaled T=1
draws_scaled_T1 <- as.data.frame(res_scaled_T1$draws)
summary_scaled_T1 <- data.frame(
  Parameter = colnames(draws_scaled_T1),
  Mean_Scaled_T1 = colMeans(draws_scaled_T1),
  SD_Scaled_T1 = apply(draws_scaled_T1, 2, sd)
)

# Merge summaries for the scaled comparison
comparison_scaled_df <- merge(summary_scaled_T1, summary_scaled_T10, by = "Parameter")
kable(comparison_scaled_df, digits = 4, caption = "Posterior Comparison (Scaled): T=1 vs T=10")


## ----qq_plots_scaled, fig.width=12, fig.height=10-----------------------------
# Vectorized alignment and R-squared calculation for scaled data
plot_data_list_scaled <- lapply(params_to_plot, function(param) {
  q_T1_scaled <- sort(draws_scaled_T1[[param]])
  q_T10_scaled <- sort(draws_scaled_T10[[param]])
  r_squared_scaled <- cor(q_T1_scaled, q_T10_scaled)^2
  
  data.frame(
    Facet_Label = paste0(param, " (R-squared = ", round(r_squared_scaled, 6), ")"),
    Exact_T1 = q_T1_scaled,
    Sequential_T10 = q_T10_scaled
  )
})

df_plot_all_scaled <- do.call(rbind, plot_data_list_scaled)

# Convert Facet_Label to a factor to enforce the specific plotting order
df_plot_all_scaled$Facet_Label <- factor(df_plot_all_scaled$Facet_Label, levels = unique(df_plot_all_scaled$Facet_Label))

ggplot(df_plot_all_scaled, aes(x = Exact_T1, y = Sequential_T10)) +
  geom_point(alpha = 0.3, color = "seagreen") +
  geom_abline(slope = 1, intercept = 0, color = "darkred", linetype = "dashed", linewidth = 1) +
  facet_wrap(~ Facet_Label, scales = "free", ncol = 3) +
  labs(
    title = "Q-Q Plots (Scaled): Exact Posterior (T=1) vs Sequential (T=10)",
    x = "Exact Posterior (Scaled, T = 1)",
    y = "Sequential Posterior (Scaled, T = 10)"
  ) +
  theme_bw()

