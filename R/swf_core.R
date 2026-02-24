# ============================================================
# swf_core.R
# Sequential Wright–Fisher (SWF) Bayesian updating utilities
# ============================================================

# ------------------------------------------------------------
# Dependency management
# ------------------------------------------------------------

load_dependencies <- function(pkgs = c("ggplot2", "ggridges", "reshape2")) {
  
  installed <- rownames(installed.packages())
  
  for (pkg in pkgs) {
    if (!pkg %in% installed) {
      message("Installing missing package: ", pkg)
      install.packages(pkg, dependencies = TRUE)
    }
    suppressPackageStartupMessages(
      library(pkg, character.only = TRUE)
    )
  }
}

# ------------------------------------------------------------
# Core rate functions
# ------------------------------------------------------------

swf_lambda <- function(n, pars) {
  n * (2 * (n - 1) + sum(pars))
}

swf_integrand <- function(z, u) {
  zdens <- dnorm(z, 0, 1)
  1 - exp(u * zdens)
}

swf_psi <- function(n, pars, rho = 1, lim = 4) {
  rate <- swf_lambda(n, pars)
  expectation <- integrate(
    f = swf_integrand,
    lower = -lim,
    upper = lim,
    u = rate
  )
  rho * expectation$value
}

# ------------------------------------------------------------
# Weight computation
# ------------------------------------------------------------

swf_weights <- function(n_succ, n_total, t, pars) {
  
  i <- n_succ
  j <- n_total
  
  nr <- i + 1
  nc <- j - i + 1
  
  weights <- matrix(0, nrow = nr, ncol = nc)
  
  for (k in 0:i) {
    for (l in 0:(j - i)) {
      
      if (k == 0 && l == 0) {
        weights[i - k + 1, j - i - l + 1] <-
          exp(swf_psi(j, pars) * t)
        next
      }
      
      hyp <- choose(i, k) * choose(j - i, l) / choose(j, k + l)
      
      rates_prod <- 1
      for (m in 0:(k + l - 1)) {
        rates_prod <- rates_prod *
          swf_lambda(j - m, pars)
      }
      
      sum_comp <- 0
      for (m1 in 0:(k + l)) {
        
        denom <- 1
        for (m2 in 0:(k + l)) {
          if (m2 != m1) {
            denom <- denom *
              abs(
                swf_lambda(j - m1, pars) -
                  swf_lambda(j - m2, pars)
              )
          }
        }
        
        sum_comp <- sum_comp +
          exp(swf_psi(j - m1, pars) * t) *
          (-1)^m1 / denom
      }
      
      weights[i - k + 1, j - i - l + 1] <-
        hyp * rates_prod * (-1)^(k + l) * sum_comp
    }
  }
  
  weights
}

# ------------------------------------------------------------
# Sequential update step
# ------------------------------------------------------------

Phi <- function(y, nt, t, theta,
                grid = seq(0.01, 0.99, by = 0.01)) {
  
  weights <- swf_weights(
    n_succ = y,
    n_total = nt,
    t = t,
    pars = theta
  )
  
  theta_new <- c(
    theta[1] + y,
    theta[2] + nt - y
  )
  
  density <- rep(0, length(grid))
  
  nr <- nrow(weights) - 1
  
  for (i in 0:nr) {
    for (j in 0:(nt - nr)) {
      density <- density +
        weights[i + 1, 1] *
        dbeta(grid,
              theta_new[1] + i,
              theta_new[2] + j)
    }
  }
  
  list(
    density = density,
    theta = theta_new
  )
}
