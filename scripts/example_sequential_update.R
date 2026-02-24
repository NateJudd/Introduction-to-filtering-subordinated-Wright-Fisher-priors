# ============================================================
# Example: Sequential Bayesian Updating under SWF Prior
# ============================================================

here::here()
source("R/swf_core.R")

library(ggplot2)
library(ggridges)
library(reshape2)

# ---- Initial parameters ----
theta0 <- c(alpha = 1.5, beta = 1.5)

# Binomial experiment settings
n_trials <- 5
y_obs <- c(3, 4, 2, 2, 5)
time_steps <- rep(0.001, length(y_obs))

# ---- Grid for densities ----
x_grid <- seq(0.01, 0.99, by = 0.01)

# ---- Prior density ----
prior_density <- dbeta(x_grid, theta0[1], theta0[2])

# ---- First conjugate update (no prediction step) ----
theta_current <- c(
  theta0[1] + y_obs[1],
  theta0[2] + n_trials - y_obs[1]
)

posterior_list <- list(prior_density)

posterior_list[[2]] <- dbeta(x_grid, theta_current[1], theta_current[2])

# ---- Sequential SWF updates ----
for (k in 2:length(y_obs)) {
  
  out <- Phi(
    y = y_obs[k],
    nt = n_trials,
    t = time_steps[k],
    theta = theta_current
  )
  
  posterior_list[[k + 1]] <- out$density
  theta_current <- out$theta
}

# ---- Combine for plotting ----
posterior_matrix <- do.call(rbind, posterior_list)
rownames(posterior_matrix) <- paste0("pi_", seq_along(posterior_list) - 1)

posterior_melt <- melt(posterior_matrix)

# create plot object
p <- ggplot(posterior_melt, aes(x = value, y = Var1)) +
  geom_density_ridges2() +
  theme_minimal()

# create output directory if it doesn't already exist
dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)

# save figure to outpute directory
ggsave("outputs/figures/sequential_update.png", plot = p, width = 8, height = 6)

