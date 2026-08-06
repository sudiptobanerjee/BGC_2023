## ----setup, include=TRUE------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)
library(knitr)
library(ggplot2)
library(dplyr)

# Source the federated learning sampler
source("../src/flSampling.r")


## ----setup_priors-------------------------------------------------------------
# Define universal priors and parameters
p <- 9 # Intercept + 8 predictors
M0_inv <- matrix(0, nrow = p, ncol = p)
m0 <- rep(0, p)
a0 <- 0.01
b0 <- 0.01
N_samples <- 5000
params_names <- c("Intercept", "X1", "X2", "X3", "X4", "X5", "X6", "Lat", "Longi", "sigma2")


## ----unscaled_analysis--------------------------------------------------------
# 1. Load Data
data <- read.table("../data/data.txt", header = TRUE)
X_full <- as.matrix(cbind(Intercept = 1, data[, c("X1", "X2", "X3", "X4", "X5", "X6", "Lat", "Longi")]))
y_full <- as.numeric(data$Y)

# 2. Partition into 10 clients
T_clients <- 10
block_size <- 130
total_n <- nrow(X_full)
group_indices <- c(rep(1:(T_clients - 1), each = block_size), 
                   rep(T_clients, total_n - (T_clients - 1) * block_size))

X_list <- lapply(split(as.data.frame(X_full), group_indices), as.matrix)
y_list <- unname(split(y_full, group_indices))

# 3. CLIENT SIDE: Compute local stats using mapply
client_stats <- mapply(compute_local_stats, X_list, y_list, SIMPLIFY = FALSE)

# 4. SERVER SIDE: Aggregate Federated Learning (T = 10)
set.seed(123)
fl_results <- flSampler(client_stats_list = client_stats, 
                        M0_inv = M0_inv, m0 = m0, a0 = a0, b0 = b0, N_samples = N_samples)
draws_fl <- as.data.frame(fl_results$draws)
colnames(draws_fl) <- params_names

summary_fl <- data.frame(
  Parameter = colnames(draws_fl),
  Mean_FL = colMeans(draws_fl),
  SD_FL = apply(draws_fl, 2, sd)
)

# 5. EXACT CENTRALIZED: Compute local stats on the entire dataset at once (T = 1)
centralized_stats <- list(compute_local_stats(X_full, y_full))
set.seed(123) # Match seed for exact draw comparison
central_results <- flSampler(client_stats_list = centralized_stats, 
                             M0_inv = M0_inv, m0 = m0, a0 = a0, b0 = b0, N_samples = N_samples)
draws_central <- as.data.frame(central_results$draws)
colnames(draws_central) <- params_names

summary_central <- data.frame(
  Parameter = colnames(draws_central),
  Mean_Central = colMeans(draws_central),
  SD_Central = apply(draws_central, 2, sd)
)

# 6. Compare Summaries
comparison_df <- merge(summary_central, summary_fl, by = "Parameter")
kable(comparison_df, digits = 4, caption = "Unscaled Posterior Comparison: Centralized vs Federated")


## ----qq_plots_unscaled, fig.width=12, fig.height=10---------------------------
# Vectorized alignment and R-squared calculation
plot_data_list <- lapply(params_names, function(param) {
  q_central <- sort(draws_central[[param]])
  q_fl <- sort(draws_fl[[param]])
  r_squared <- cor(q_central, q_fl)^2
  
  data.frame(
    Facet_Label = paste0(param, " (R-squared = ", round(r_squared, 6), ")"),
    Centralized = q_central,
    Federated = q_fl
  )
})

df_plot_all <- do.call(rbind, plot_data_list)
df_plot_all$Facet_Label <- factor(df_plot_all$Facet_Label, levels = unique(df_plot_all$Facet_Label))

ggplot(df_plot_all, aes(x = Centralized, y = Federated)) +
  geom_point(alpha = 0.3, color = "steelblue") +
  geom_abline(slope = 1, intercept = 0, color = "darkred", linetype = "dashed", linewidth = 1) +
  facet_wrap(~ Facet_Label, scales = "free", ncol = 3) +
  labs(
    title = "Q-Q Plots: Centralized vs Federated (Unscaled)",
    x = "Centralized Posterior (T = 1)",
    y = "Federated Posterior (T = 10)"
  ) +
  theme_bw()


## ----scaled_analysis----------------------------------------------------------
# 1. Center and Scale covariates (exclude intercept)
covariates <- data[, c("X1", "X2", "X3", "X4", "X5", "X6", "Lat", "Longi")]
X_scaled_full <- as.matrix(cbind(Intercept = 1, scale(covariates)))

# 2. Partition scaled data into 10 clients
X_list_scaled <- lapply(split(as.data.frame(X_scaled_full), group_indices), as.matrix)

# 3. CLIENT SIDE: Compute local stats on scaled data
client_stats_scaled <- mapply(compute_local_stats, X_list_scaled, y_list, SIMPLIFY = FALSE)

# 4. SERVER SIDE: Aggregate Federated Learning (T = 10)
set.seed(123)
fl_results_scaled <- flSampler(client_stats_list = client_stats_scaled, 
                               M0_inv = M0_inv, m0 = m0, a0 = a0, b0 = b0, N_samples = N_samples)
draws_fl_scaled <- as.data.frame(fl_results_scaled$draws)
colnames(draws_fl_scaled) <- params_names

summary_fl_scaled <- data.frame(
  Parameter = colnames(draws_fl_scaled),
  Mean_FL_Scaled = colMeans(draws_fl_scaled),
  SD_FL_Scaled = apply(draws_fl_scaled, 2, sd)
)

# 5. EXACT CENTRALIZED: Compute local stats on entire scaled dataset at once (T = 1)
centralized_stats_scaled <- list(compute_local_stats(X_scaled_full, y_full))
set.seed(123)
central_results_scaled <- flSampler(client_stats_list = centralized_stats_scaled, 
                                    M0_inv = M0_inv, m0 = m0, a0 = a0, b0 = b0, N_samples = N_samples)
draws_central_scaled <- as.data.frame(central_results_scaled$draws)
colnames(draws_central_scaled) <- params_names

summary_central_scaled <- data.frame(
  Parameter = colnames(draws_central_scaled),
  Mean_Central_Scaled = colMeans(draws_central_scaled),
  SD_Central_Scaled = apply(draws_central_scaled, 2, sd)
)

# 6. Compare Summaries
comparison_scaled_df <- merge(summary_central_scaled, summary_fl_scaled, by = "Parameter")
kable(comparison_scaled_df, digits = 4, caption = "Scaled Posterior Comparison: Centralized vs Federated")


## ----qq_plots_scaled, fig.width=12, fig.height=10-----------------------------
# Vectorized alignment and R-squared calculation for scaled data
plot_data_list_scaled <- lapply(params_names, function(param) {
  q_central_scaled <- sort(draws_central_scaled[[param]])
  q_fl_scaled <- sort(draws_fl_scaled[[param]])
  r_squared_scaled <- cor(q_central_scaled, q_fl_scaled)^2
  
  data.frame(
    Facet_Label = paste0(param, " (R-squared = ", round(r_squared_scaled, 6), ")"),
    Centralized = q_central_scaled,
    Federated = q_fl_scaled
  )
})

df_plot_all_scaled <- do.call(rbind, plot_data_list_scaled)
df_plot_all_scaled$Facet_Label <- factor(df_plot_all_scaled$Facet_Label, levels = unique(df_plot_all_scaled$Facet_Label))

ggplot(df_plot_all_scaled, aes(x = Centralized, y = Federated)) +
  geom_point(alpha = 0.3, color = "seagreen") +
  geom_abline(slope = 1, intercept = 0, color = "darkred", linetype = "dashed", linewidth = 1) +
  facet_wrap(~ Facet_Label, scales = "free", ncol = 3) +
  labs(
    title = "Q-Q Plots (Scaled): Centralized vs Federated",
    x = "Centralized Posterior (Scaled, T = 1)",
    y = "Federated Posterior (Scaled, T = 10)"
  ) +
  theme_bw()

