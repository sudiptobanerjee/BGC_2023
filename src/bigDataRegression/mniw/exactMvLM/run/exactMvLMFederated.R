## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)
library(knitr)
library(mniw)
source("../src/mniwSamplerFederated.r")


## ----execute_federated_custom-------------------------------------------------
# Ingest dataset files
X_raw  <- as.matrix(read.table("../data/X.txt"))
Y_raw  <- as.matrix(read.table("../data/Y.txt"))
X_full <- cbind(1, X_raw)
Y_full <- Y_raw

p_dim <- ncol(X_full)
q_dim <- ncol(Y_full)

# Global prior parameters
m0     <- matrix(0, nrow = p_dim, ncol = q_dim)
M0_inv <- diag(1, p_dim)
nu0    <- q_dim + 2
Psi0   <- diag(1, q_dim)

# Partition data across T = 5 clients
T_clients   <- 5
client_idxs <- split(seq_len(nrow(X_full)), cut(seq_len(nrow(X_full)), T_clients, labels = FALSE))
client_data_list <- lapply(client_idxs, function(idx) {
  list(X = X_full[idx, ], Y = Y_full[idx, ])
})

# 1. Federated Learning Sampler (Custom R)
set.seed(42)
res_fed_custom <- mniwSamplerFederated(client_data_list, m0, M0_inv, nu0, Psi0, n_samples = 2000)

# 2. Single-Server Full Batch Sampler
full_server_data <- list(list(X = X_full, Y = Y_full))
set.seed(42)
res_batch_custom <- mniwSamplerFederated(full_server_data, m0, M0_inv, nu0, Psi0, n_samples = 2000)


## ----compare_parameters_custom------------------------------------------------
diff_m     <- max(abs(res_fed_custom$post_params$m - res_batch_custom$post_params$m))
diff_M_inv <- max(abs(res_fed_custom$post_params$M_inv - res_batch_custom$post_params$M_inv))
diff_Psi   <- max(abs(res_fed_custom$post_params$Psi - res_batch_custom$post_params$Psi))
diff_nu    <- abs(res_fed_custom$post_params$nu - res_batch_custom$post_params$nu)

param_comp <- data.frame(
  Parameter_Matrix = c("Location Vector (m)", "Precision Matrix (M^{-1})", "Scale Matrix (Psi)", "Degrees of Freedom (nu)"),
  Max_Absolute_Difference = c(diff_m, diff_M_inv, diff_Psi, diff_nu)
)

kable(param_comp, digits = 16, caption = "Numerical Equivalence: Federated vs. Single Server (Custom Sampler)")


## ----posterior_summary_custom, echo=FALSE, results='asis'---------------------
calc_summaries <- function(arr, name) {
  grid <- expand.grid(Row = 1:dim(arr)[1], Col = 1:dim(arr)[2])
  res <- lapply(seq_len(nrow(grid)), function(k) {
    i <- grid$Row[k]; j <- grid$Col[k]
    vals <- arr[i, j, ]
    data.frame(
      Param = paste0(name, "[", i, ",", j, "]"),
      Mean = mean(vals), SD = sd(vals),
      `2.5%` = quantile(vals, 0.025), `97.5%` = quantile(vals, 0.975),
      check.names = FALSE
    )
  })
  do.call(rbind, res)
}

df_beta_fed  <- calc_summaries(res_fed_custom$samples$Beta, "Beta")
df_sigma_fed <- calc_summaries(res_fed_custom$samples$Sigma, "Sigma")

cat("### Regression Coefficients (Beta) Summary (Custom Federated Sampler)\n\n")
print(kable(df_beta_fed, digits = 4, row.names = FALSE))

cat("\n### Covariance Matrix (Sigma) Summary (Custom Federated Sampler)\n\n")
print(kable(df_sigma_fed, digits = 4, row.names = FALSE))


## ----qqplot_beta_custom, fig.width=9, fig.height=7, echo=FALSE----------------
par(mfrow = c(p_dim, q_dim), mar = c(3, 3, 2, 1))
grid_beta <- expand.grid(Row = 1:p_dim, Col = 1:q_dim)

invisible(lapply(seq_len(nrow(grid_beta)), function(k) {
  i <- grid_beta$Row[k]; j <- grid_beta$Col[k]
  qqplot(res_fed_custom$samples$Beta[i, j, ], res_batch_custom$samples$Beta[i, j, ],
         main = paste0("Beta[", i, ",", j, "]"),
         xlab = "Federated Quantiles", ylab = "Single Server Quantiles",
         pch = 16, col = rgb(0.2, 0.4, 0.8, 0.4))
  abline(0, 1, col = "firebrick", lwd = 2)
}))
par(mfrow = c(1, 1))


## ----qqplot_sigma_custom, fig.width=8, fig.height=8, echo=FALSE---------------
par(mfrow = c(q_dim, q_dim), mar = c(3, 3, 2, 1))
grid_sigma <- expand.grid(Row = 1:q_dim, Col = 1:q_dim)

invisible(lapply(seq_len(nrow(grid_sigma)), function(k) {
  i <- grid_sigma$Row[k]; j <- grid_sigma$Col[k]
  qqplot(res_fed_custom$samples$Sigma[i, j, ], res_batch_custom$samples$Sigma[i, j, ],
         main = paste0("Sigma[", i, ",", j, "]"),
         xlab = "Federated Quantiles", ylab = "Single Server Quantiles",
         pch = 16, col = rgb(0.8, 0.4, 0.2, 0.4))
  abline(0, 1, col = "firebrick", lwd = 2)
}))
par(mfrow = c(1, 1))


## ----execute_federated_rmniw--------------------------------------------------
# 1. Federated Learning Sampler (rmniw)
set.seed(42)
res_fed_rmniw <- mniwSamplerFederated_rmniw(client_data_list, m0, M0_inv, nu0, Psi0, n_samples = 2000)

# 2. Single-Server Full Batch Sampler (rmniw)
set.seed(42)
res_batch_rmniw <- mniwSamplerFederated_rmniw(full_server_data, m0, M0_inv, nu0, Psi0, n_samples = 2000)


## ----compare_parameters_rmniw-------------------------------------------------
diff_m_rmniw     <- max(abs(res_fed_rmniw$post_params$m - res_batch_rmniw$post_params$m))
diff_M_inv_rmniw <- max(abs(res_fed_rmniw$post_params$M_inv - res_batch_rmniw$post_params$M_inv))
diff_Psi_rmniw   <- max(abs(res_fed_rmniw$post_params$Psi - res_batch_rmniw$post_params$Psi))
diff_nu_rmniw    <- abs(res_fed_rmniw$post_params$nu - res_batch_rmniw$post_params$nu)

param_comp_rmniw <- data.frame(
  Parameter_Matrix = c("Location Vector (m)", "Precision Matrix (M^{-1})", "Scale Matrix (Psi)", "Degrees of Freedom (nu)"),
  Max_Absolute_Difference = c(diff_m_rmniw, diff_M_inv_rmniw, diff_Psi_rmniw, diff_nu_rmniw)
)

kable(param_comp_rmniw, digits = 16, caption = "Numerical Equivalence: Federated vs. Single Server (rmniw Sampler)")


## ----posterior_summary_rmniw, echo=FALSE, results='asis'----------------------
df_beta_fed_rmniw  <- calc_summaries(res_fed_rmniw$samples$Beta, "Beta")
df_sigma_fed_rmniw <- calc_summaries(res_fed_rmniw$samples$Sigma, "Sigma")

cat("### Regression Coefficients (Beta) Summary (rmniw Federated Sampler)\n\n")
print(kable(df_beta_fed_rmniw, digits = 4, row.names = FALSE))

cat("\n### Covariance Matrix (Sigma) Summary (rmniw Federated Sampler)\n\n")
print(kable(df_sigma_fed_rmniw, digits = 4, row.names = FALSE))


## ----qqplot_beta_rmniw, fig.width=9, fig.height=7, echo=FALSE-----------------
par(mfrow = c(p_dim, q_dim), mar = c(3, 3, 2, 1))

invisible(lapply(seq_len(nrow(grid_beta)), function(k) {
  i <- grid_beta$Row[k]; j <- grid_beta$Col[k]
  qqplot(res_fed_rmniw$samples$Beta[i, j, ], res_batch_rmniw$samples$Beta[i, j, ],
         main = paste0("Beta[", i, ",", j, "]"),
         xlab = "Federated rmniw Quantiles", ylab = "Single Server rmniw Quantiles",
         pch = 16, col = rgb(0.2, 0.4, 0.8, 0.4))
  abline(0, 1, col = "firebrick", lwd = 2)
}))
par(mfrow = c(1, 1))


## ----qqplot_sigma_rmniw, fig.width=8, fig.height=8, echo=FALSE----------------
par(mfrow = c(q_dim, q_dim), mar = c(3, 3, 2, 1))

invisible(lapply(seq_len(nrow(grid_sigma)), function(k) {
  i <- grid_sigma$Row[k]; j <- grid_sigma$Col[k]
  qqplot(res_fed_rmniw$samples$Sigma[i, j, ], res_batch_rmniw$samples$Sigma[i, j, ],
         main = paste0("Sigma[", i, ",", j, "]"),
         xlab = "Federated rmniw Quantiles", ylab = "Single Server rmniw Quantiles",
         pch = 16, col = rgb(0.8, 0.4, 0.2, 0.4))
  abline(0, 1, col = "firebrick", lwd = 2)
}))
par(mfrow = c(1, 1))


## ----timing_comparison--------------------------------------------------------
n_bench <- 10000

t_custom <- system.time({
  res_bench_custom <- mniwSamplerFederated(client_data_list, m0, M0_inv, nu0, Psi0, n_samples = n_bench)
})["elapsed"]

t_rmniw <- system.time({
  res_bench_rmniw <- mniwSamplerFederated_rmniw(client_data_list, m0, M0_inv, nu0, Psi0, n_samples = n_bench)
})["elapsed"]

timing_df <- data.frame(
  Sampler_Implementation = c("Custom R Sampler (mniwSamplerFederated)", 
                             "C++ Accelerated Sampler (mniwSamplerFederated_rmniw)"),
  Elapsed_Time_Sec = c(as.numeric(t_custom), as.numeric(t_rmniw)),
  Speedup_Factor   = c(1.0, as.numeric(t_custom / t_rmniw))
)

kable(timing_df, digits = 4, caption = paste0("Federated Sampling Execution Time Comparison (", n_bench, " Samples)"))

