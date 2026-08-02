# 1. Source the external sampling function and load required packages
source("../src/rawSampler.r")
library(spBayes) # Required for bayesLMConjugate

# 2. Load the dataset
dataset <- read.table("../../data/data.txt", header = TRUE)
y <- dataset[, 1]

# Extract raw predictors and scale them
X_raw <- as.matrix(dataset[, -1])
X_scaled <- scale(X_raw)

# Explicitly add a column of 1s for the intercept to the scaled predictors
X <- cbind(Intercept = 1, X_scaled)

# Ensure column names are clean and explicitly labeled
if (!is.null(colnames(dataset))) {
  colnames(X)[-1] <- colnames(dataset)[-1]
} else {
  colnames(X)[-1] <- paste0("X", 1:ncol(X_scaled))
}

n <- nrow(X)
p <- ncol(X)

# 3. Set prior parameters (Flat, uninformative prior)
a0 <- -p/2
b0 <- 0.0
M0_inv <- matrix(0, nrow = p, ncol = p)
m0 <- matrix(0, nrow = p, ncol = 1)

# 4. Perform Data Analysis Using rawSampler
set.seed(42)
N_samples <- 5000

posterior_samples <- rawSampler(y = y, X = X, m0 = m0, M0_inv = M0_inv, a0 = a0, b0 = b0, N = N_samples)
raw_draws <- posterior_samples$draws

# 5. Perform Data Analysis Using bayesLMConjugate
# Create a data frame for the formula interface
df <- data.frame(y = y, X_scaled)

spB_model <- bayesLMConjugate(
  y ~ ., 
  data = df, 
  n.samples = N_samples, 
  beta.prior.mean = rep(0, p),
  beta.prior.precision = M0_inv,
  prior.shape = a0,
  prior.rate = b0
)

# Extract and align spBayes samples using the correct list tag
spB_mcmc <- as.matrix(spB_model$p.beta.tauSq.samples)

# Robustly extract the variance parameter (it is always the last column)
spB_sigma2 <- spB_mcmc[, ncol(spB_mcmc)] 
spB_sigma <- sqrt(spB_sigma2)
spB_tau <- 1 / spB_sigma2

# Bind spBayes output into a matching matrix structure
spB_draws <- cbind(spB_mcmc[, 1:p], sigma = spB_sigma, sigma2 = spB_sigma2, tau = spB_tau)
colnames(spB_draws)[1:p] <- colnames(X)

# 6. Generate Q-Q Plots comparing the two samplers
num_params <- p + 2 
cols <- 3
rows <- ceiling(num_params / cols)

# Set graphical parameters
par(mfrow = c(rows, cols), pty = "s", mar = c(4, 4, 3, 1), oma = c(0, 0, 4, 0))
plot_params <- c(colnames(X), "sigma", "tau")

for (param in plot_params) {
  raw_vals <- raw_draws[, param]
  spB_vals <- spB_draws[, param]
  
  qq_data <- qqplot(raw_vals, spB_vals, plot.it = FALSE)
  r_squared <- cor(qq_data$x, qq_data$y)^2
  
  plot(qq_data$x, qq_data$y, 
       main = param, 
       xlab = "composition sampling in raw r code", 
       ylab = "bayesLMConjugate sampler in spBayes",
       col = rgb(0.2, 0.4, 0.8, 0.5), 
       pch = 16,
       cex = 1.2,          
       cex.main = 1.5,     
       cex.lab = 1.2,      
       cex.axis = 1.1)     
  
  abline(0, 1, col = "red", lwd = 2)
  
  legend("topleft", bty = "n", legend = bquote(R^2 == .(round(r_squared, 4))), cex = 1.2)
}

# Add a main title to the overall plot window
mtext("Comparison of Posterior Draws: rawSampler vs bayesLMConjugate", outer = TRUE, font = 2, cex = 1.8)

# Reset plotting parameters
par(mfrow = c(1, 1), pty = "m", mar = c(5, 4, 4, 2) + 0.1, oma = c(0, 0, 0, 0))
