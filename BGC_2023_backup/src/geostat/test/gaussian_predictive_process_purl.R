knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)

# Check if spBayes is installed, install if not
if (!requireNamespace("spBayes", quietly = TRUE)) {
  install.packages("spBayes")
}
# Check if ggplot2 is installed, install if not
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
# Check if fields is installed, install if not
if (!requireNamespace("fields", quietly = TRUE)) {
  install.packages("fields")
}

# Load libraries
library(spBayes)
library(ggplot2)
library(fields)

set.seed(123)
n <- 500  # Number of observed data points
coords <- cbind(runif(n, 0, 1), runif(n, 0, 1))
x <- runif(n, 0, 1)
beta <- c(1, 2)
phi <- 6  # Spatial decay parameter
sigma.sq <- 1
tau.sq <- 0.1

# Generate spatial process w
dist.mat <- as.matrix(dist(coords))
R <- exp(-phi * dist.mat)
w <- t(chol(R)) %*% rnorm(n)

# Generate response variable y
y <- beta[1] + beta[2] * x + w * sqrt(sigma.sq) + rnorm(n, 0, sqrt(tau.sq))

# Combine into a dataframe
data.df <- data.frame(y = y, x = x, coords = coords)

# Choose a smaller number of knots
n.knots <- 50
knot.coords <- fields::cover.design(coords, n.knots)$design

# Plot data locations and knots
plot(coords, main = "Data locations and knots", xlab = "Easting", ylab = "Northing")
points(knot.coords, pch = 19, col = "red")
legend("topright", c("Data locations", "Knots"), pch = c(1, 19), col = c("black", "red"))

# Define model formula
formula <- y ~ x

# Define starting values and priors
starting <- list(beta = beta, sigma.sq = sigma.sq, tau.sq = tau.sq, phi = phi)
tuning <- list(phi = 0.5, sigma.sq=0.01, tau.sq=0.01) 
priors <- list(beta.Norm = list(rep(0, 2), diag(1000, 2)),
               phi.Unif = c(3/1, 3/0.1), 
               sigma.sq.IG = c(2, 1),
               tau.sq.IG = c(2, 1))

# Fit the model
n.samples <- 5000
model <- spBayes::spLM(formula,
                 coords = coords,
                 knots = knot.coords,
                 data = data.df,
                 starting = starting,
                 tuning = tuning,
                 priors = priors,
                 cov.model = "exponential",
                 n.samples = n.samples)

# Burn-in and thin the MCMC samples
burn.in <- 0.5 * n.samples
model.thin <- spBayes::spRecover(model, start = burn.in, thin = 10, verbose = FALSE)


# Create new prediction locations
n.pred <- 100
pred.coords <- cbind(runif(n.pred, 0, 1), runif(n.pred, 0, 1))
pred.x <- runif(n.pred, 0, 1)
pred.covars <- cbind(1, pred.x)

# Predict at new locations
pred.samples <- spBayes::spPredict(
  model.thin,
  pred.coords = pred.coords,
  pred.covars = pred.covars,
  verbose = FALSE
)

# Summarize the predictive samples
pred.summary <- apply(pred.samples$p.y.predictive.samples, 1, quantile, probs = c(0.025, 0.5, 0.975))

# Create a data frame for plotting predictions
pred.df <- data.frame(
  x = pred.coords[, 1],
  y = pred.coords[, 2],
  median = pred.summary[2, ],
  lower = pred.summary[1, ],
  upper = pred.summary[3, ]
)

#Carry out spatial interpolation using the "akima" package
library(akima)
akima_interp <- interp(x = pred.df$x, y = pred.df$y, z = pred.df$median, nx = 100, ny = 100)
akima_df <- data.frame(
  x = rep(akima_interp$x, times = length(akima_interp$y)),
  y = rep(akima_interp$y, each = length(akima_interp$x)),
  median = as.vector(akima_interp$z)
)

#Carry out spatial interpolation using the "mba" package
library(MBA)

# Create a 3-column data frame or matrix for mba.surf
# The columns must be in the order x, y, and the value to be interpolated (z)
mba_data <- pred.df[, c("x", "y", "median")]
# Perform the interpolation
# `mba.surf` takes the data and the desired resolution (e.g., 100x100).
# `extend = TRUE` ensures the entire grid is covered, not just the convex hull of the points.
mba_surface <- mba.surf(mba_data, 100, 100, extend = TRUE)$xyz.est

# Create a heatmap of the median predicted values using the akima interpolator
library(ggplot2)
ggplot(akima_df, aes(x = x, y = y)) +
  geom_tile(aes(fill = median)) +
  scale_fill_viridis_c(option = "plasma") +
  labs(title = "Median Predicted Surface Using Akima",
       x = "Easting",
       y = "Northing",
       fill = "Predicted Value") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Create a heatmap of the median predicted values using the mba interpolator
x.res <- 100
y.res <- 100
image.plot(mba_surface, xaxs = "r", yaxs = "r", xlab = "Easting (m)", ylab = "Northing (m)", main="Median Predicted Surface Using MBA")
