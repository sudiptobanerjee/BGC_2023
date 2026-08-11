## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)
library(knitr)
library(mniw)
library(microbenchmark)

source("../src/mniwSamplerDivideAndConquer.r")
source("../src/mniwSamplerFederated.r")


## ----flops_table, echo=FALSE--------------------------------------------------
flops_df <- data.frame(
  Phase = c("Local/Block Crossproducts", "Intermediate Overhead (per block)", "Server Aggregation", "Final Cholesky & Backsolve", "Total Asymptotic FLOPs"),
  Sequential_Divide_and_Conquer = c("O(N(p + q)^2)", "O(p^3 + p^2 q + p q^2)", "N/A (Streaming)", "O(p^3 + p^2 q + p q^2)", "O(N(p + q)^2 + T(p^3 + p^2 q + p q^2))"),
  Federated_Learning = c("O((N/T)(p + q)^2) [Parallel]", "None (Zero local matrix solves)", "O(T(p^2 + pq + q^2))", "O(p^3 + p^2 q + p q^2) [One-Shot]", "O((N/T)(p + q)^2 + T(p^2 + pq + q^2) + p^3 + p^2 q + p q^2)")
)

kable(flops_df, caption = "Theoretical FLOPs Comparison (Sequential D&C vs. Federated Learning)")


## ----storage_table, echo=FALSE------------------------------------------------
storage_df <- data.frame(
  Metric = c("Memory per Local Node", "Transmitted Network Payload per Client", "Data Privacy Preserved", "Server Memory Footprint"),
  Sequential_Divide_and_Conquer = c("O(n_t (p + q))", "O(n_t (p + q)) [Sequential Stream]", "No (Raw data streamed)", "O(p^2 + pq + q^2)"),
  Federated_Learning = c("O(n_t (p + q))", "p^2 + pq + q^2 + 1 scalars [Fixed]", "Yes (Zero raw data transfer)", "O(T(p^2 + pq + q^2))")
)

kable(storage_df, caption = "Storage and Communication Bandwidth Complexity")


## ----execute_real_benchmark---------------------------------------------------
# Ingest raw dataset files
X_raw  <- as.matrix(read.table("../data/X.txt"))
Y_raw  <- as.matrix(read.table("../data/Y.txt"))
X_full <- cbind(1, X_raw)
Y_full <- Y_raw

p_dim <- ncol(X_full)
q_dim <- ncol(Y_full)

# Global prior setup
m0     <- matrix(0, nrow = p_dim, ncol = q_dim)
M0_inv <- diag(1, p_dim)
nu0    <- q_dim + 2
Psi0   <- diag(1, q_dim)

# Partition data into T = 5 blocks/clients
T_blocks    <- 5
block_sizes <- split(seq_len(nrow(X_full)), cut(seq_len(nrow(X_full)), T_blocks, labels = FALSE))
data_blocks <- lapply(block_sizes, function(idx) {
  list(X = X_full[idx, ], Y = Y_full[idx, ])
})

n_draws <- 10000

# 1. Custom Sampler Benchmark (Sequential D&C vs. Federated)
set.seed(42)
t_dc_custom <- system.time({
  res_dc_custom <- mniwSamplerDivideAndConquer(data_blocks, m0, M0_inv, nu0, Psi0, n_samples = n_draws)
})["elapsed"]

set.seed(42)
t_fed_custom <- system.time({
  res_fed_custom <- mniwSamplerFederated(data_blocks, m0, M0_inv, nu0, Psi0, n_samples = n_draws)
})["elapsed"]

# 2. C++ Accelerated Sampler Benchmark (mniw::rmniw)
set.seed(42)
t_dc_rmniw <- system.time({
  res_dc_rmniw <- mniwSamplerDivideAndConquer_rmniw(data_blocks, m0, M0_inv, nu0, Psi0, n_samples = n_draws)
})["elapsed"]

set.seed(42)
t_fed_rmniw <- system.time({
  res_fed_rmniw <- mniwSamplerFederated_rmniw(data_blocks, m0, M0_inv, nu0, Psi0, n_samples = n_draws)
})["elapsed"]


## ----real_data_results_table, echo=FALSE--------------------------------------
bench_real_df <- data.frame(
  Sampler_Type = c("Custom R Composition Sampler", "Custom R Composition Sampler", 
                   "C++ Accelerated Sampler (mniw::rmniw)", "C++ Accelerated Sampler (mniw::rmniw)"),
  Method = c("Sequential Divide & Conquer", "Federated Learning", 
             "Sequential Divide & Conquer", "Federated Learning"),
  Elapsed_Time_Sec = c(as.numeric(t_dc_custom), as.numeric(t_fed_custom), 
                       as.numeric(t_dc_rmniw), as.numeric(t_fed_rmniw)),
  Speedup_vs_DC = c(1.0, as.numeric(t_dc_custom / t_fed_custom), 
                    1.0, as.numeric(t_dc_rmniw / t_fed_rmniw))
)

kable(bench_real_df, digits = 4, caption = paste0("Empirical Real-Data Benchmark Results (", n_draws, " Posterior Draws)"))


## ----simulation_benchmark-----------------------------------------------------
run_simulation_benchmark <- function(N_total, T_count, p = 8, q = 4, n_draws = 5000) {
  # Generate synthetic dataset
  X_sim <- cbind(1, matrix(rnorm(N_total * (p - 1)), nrow = N_total, ncol = p - 1))
  Beta_true <- matrix(rnorm(p * q), nrow = p, ncol = q)
  Y_sim <- X_sim %*% Beta_true + matrix(rnorm(N_total * q), nrow = N_total, ncol = q)

  m0_sim     <- matrix(0, nrow = p, ncol = q)
  M0_inv_sim <- diag(1, p)
  nu0_sim    <- q + 2
  Psi0_sim   <- diag(1, q)

  # Partition into T blocks (Fixed: cut() parenthesis)
  idxs <- split(seq_len(N_total), cut(seq_len(N_total), T_count, labels = FALSE))
  sim_blocks <- lapply(idxs, function(idx) list(X = X_sim[idx, ], Y = Y_sim[idx, ]))

  # Timing: Custom Samplers
  t_dc_c  <- system.time(mniwSamplerDivideAndConquer(sim_blocks, m0_sim, M0_inv_sim, nu0_sim, Psi0_sim, n_samples = n_draws))["elapsed"]
  t_fed_c <- system.time(mniwSamplerFederated(sim_blocks, m0_sim, M0_inv_sim, nu0_sim, Psi0_sim, n_samples = n_draws))["elapsed"]

  # Timing: C++ rmniw Samplers
  t_dc_r  <- system.time(mniwSamplerDivideAndConquer_rmniw(sim_blocks, m0_sim, M0_inv_sim, nu0_sim, Psi0_sim, n_samples = n_draws))["elapsed"]
  t_fed_r <- system.time(mniwSamplerFederated_rmniw(sim_blocks, m0_sim, M0_inv_sim, nu0_sim, Psi0_sim, n_samples = n_draws))["elapsed"]

  data.frame(
    N = N_total,
    T = T_count,
    Custom_DC_Sec  = as.numeric(t_dc_c),
    Custom_FED_Sec = as.numeric(t_fed_c),
    Custom_Speedup = as.numeric(t_dc_c / t_fed_c),
    rmniw_DC_Sec   = as.numeric(t_dc_r),
    rmniw_FED_Sec  = as.numeric(t_fed_r),
    rmniw_Speedup  = as.numeric(t_dc_r / t_fed_r)
  )
}

# Run benchmark grid over varying N and T
set.seed(123)
grid_params <- list(
  list(N = 10000, T = 5),
  list(N = 50000, T = 10),
  list(N = 100000, T = 20)
)

sim_results <- do.call(rbind, lapply(grid_params, function(g) {
  run_simulation_benchmark(N_total = g$N, T_count = g$T, n_draws = 5000)
}))


## ----simulation_results_table, echo=FALSE-------------------------------------
kable(sim_results, digits = 4, caption = "Simulation Experiment Scaling Performance (5,000 Posterior Draws)")

