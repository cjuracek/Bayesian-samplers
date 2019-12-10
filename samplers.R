library(MASS)

# 2 helper functions for MVN sampler
compute_mu_n <- function(lambda_0, n, sigma, mu_0, y_bar) {
  lambda_0_inv <- solve(lambda_0)
  sigma_inv <- solve(sigma)
  
  first_factor <- (lambda_0_inv + n * sigma_inv) %>% solve()
  second_factor <- lambda_0_inv %*% mu_0 + n * sigma_inv %*% y_bar
  
  return(first_factor %*% second_factor)
}

compute_S_n <- function(Y, theta, S_0) {
  S_theta <- matrix(0, ncol(Y), ncol(Y))
  for(i in 1:nrow(Y)) {
    S_theta <- S_theta + (Y[i, ] - theta) %*% t(Y[i, ] - theta)
  }
  return(S_0 + S_theta)
}

# Gibbs sampler for MCMC approximation of mean, covariance of multivariate normal
# Based off of algorithm in Hoff, page 112
#   nu_0, S_0: hyperparameters for covariance
#   Lambda_0, mu_0: hyperparameters for mean
# Returns named list of thetas (matrix w/each row a predictor mean) and
#   sigmas (nested list of p x p covariance matrices)
gibbs_mvn <- function(Y, num_iter = 10000,
                      nu_0 = 4, S_0 = cov(Y),
                      Lambda_0 = cov(Y), mu_0 = matrix(apply(Y, 2, mean), ncol(Y), 1)) {
  
  p <- ncol(Y)
  n <- nrow(Y)
  
  # https://stackoverflow.com/questions/19340079/how-to-declare-list-object-with-m-elements
  thetas <- matrix(NA, num_iter, p)
  sigmas <- rep(list(matrix(NA, p, p)), num_iter)
  y_bar <- matrix(apply(Y, 2, mean), p, 1)
  
  # Initial theta value for chain
  thetas[1,] <- mu_0
  for(i in seq_len(num_iter - 1)) {
    
    # Update Sigma via full conditional
    S_n <- compute_S_n(Y, thetas[i,], S_0)
    sigmas[[i + 1]] <- solve(rWishart(1, nu_0 + n, solve(S_n))[,,1])
    
    # Update theta via full conditional
    mu_n <- compute_mu_n(Lambda_0, n, sigmas[[i + 1]], mu_0, y_bar)
    Lambda_n <- solve(solve(Lambda_0) + n * solve(sigmas[[i + 1]]))
    thetas[i + 1,] <- mvrnorm(1, mu_n, Lambda_n)
    
  }
  
  return(list(theta = thetas, sigma = sigmas))
}

##################################################################################################################

# Sample a value of mu from its full conditional (normal) in the hierarchical normal model
sample_mu_full <- function(m, theta_bar, tau_sq, mu_0, gamma_0_sq) {
  mu_mean <- ((m * theta_bar) / tau_sq + mu_0 / gamma_0_sq) / ((m / tau_sq) + (1 / gamma_0_sq))
  mu_var <- 1 / ((m / tau_sq) + (1 / gamma_0_sq))
  return(rnorm(1, mu_mean, sqrt(mu_var)))
}

# Sample a value of mu from its full conditional (gamma) in the hierarchical normal model
# Remember, the precision follows a gamma, so take reciprocal at the end.
sample_tau_sq_full <- function(eta_0, m, tau_0_sq, theta, mu) {
  a <- (eta_0 + m) / 2
  b <- (eta_0 * tau_0_sq + sum((theta - mu) ^ 2)) / 2
  return(1 / rgamma(1, a, b))
}

# Sample a value of sigma squared from its full conditional in the hierarchical normal model.
# Remember, the precision follows a gamma, so take reciprocal at the end
sample_sigma_sq_full <- function(data, nu_0, num_obs, sigma_0_sq, theta) {
  total_ssr <- 0
  for(i in seq_along(data)) {
    total_ssr <- total_ssr + sum((data[[i]] - theta[i]) ^ 2)
  }
  a <- (nu_0 + num_obs) / 2
  b <- (nu_0 * sigma_0_sq + total_ssr) / 2
  return(1 / rgamma(1, a, b))
}

# Sample theta's from their respective full conditionals in the hierarchical normal model.
sample_theta_full <- function(data, sizes, avgs, 
                              sigma_sq, mu, tau_sq) {
  theta_mean <- ((sizes * avgs) / sigma_sq + mu / tau_sq) / (sizes / sigma_sq + 1 / tau_sq)
  theta_var <- 1 / (sizes / sigma_sq + 1 / tau_sq)
  return(rnorm(length(data), theta_mean, sqrt(theta_var)))
}

# Uses Gibbs sampling to produce MCMC approximations 
# to full joint posterior distribution of hierarchical normal model
# Reference: Hoff, Chapter 8
# Data: List of groups
# Takes an initial state:
#   - theta_1 - theta_m: group means
#   - mu: mean of sampling distribution of theta's
#   - tau_sq: between school variance
#   - sigma_sq: within school variance
# and prior parameters. Uses Gibbs sampling to produce MCMC approximations 
# to full joint posterior distribution.
# Returns: named list of posterior approximations for theta's,
#   sigma_sq, mu, and tau
gibbs_hierarchical <- function(data, mu_0, gamma_0_sq, tau_0_sq,
                               eta_0, sigma_0_sq, nu_0,
                               num_iter = 5000) {
  
  # Compute convenient invariant information
  num_obs <- length(unlist(data))
  group_sizes <- sapply(data, length)
  group_avgs <- sapply(data, mean)
  
  m <- length(data)
  theta <- matrix(0, num_iter, m); sigma_sq <- numeric(num_iter)
  mu <- numeric(num_iter); tau_sq <- numeric(num_iter)
  
  # Define 1st value in chain
  theta[1,] <- group_avgs
  sigma_sq[1] <- sigma_0_sq
  mu[1] <- mu_0
  tau_sq[1] <- tau_0_sq
  
  for(s in seq_len(num_iter - 1)) {
    
    # Update mu from full conditional
    mu[s + 1] <- sample_mu_full(m, mean(theta[s,]), 
                            tau_sq[s], mu_0, gamma_0_sq)
    
    # Update tau_sq from full conditional
    tau_sq[s + 1] <- sample_tau_sq_full(eta_0, m, tau_0_sq,
                                    theta[s,], mu[s + 1])
    
    # Update sigma_sq...
    sigma_sq[s + 1] <- sample_sigma_sq_full(data, nu_0, num_obs,
                                        sigma_0_sq, theta[s,])
    
    # Update thetas...
    theta[s + 1,] <- sample_theta_full(data = data, sizes = group_sizes, avgs = group_avgs,
                                    sigma_sq[s + 1], mu[s + 1], tau_sq[s + 1])
  }
  
  return(list(theta = theta, sigma_sq = sigma_sq, mu = mu, tau_sq = tau_sq))
}

################################################################################

### log-density of the multivariate normal distribution (Fan Bu)
ldmvnorm <- function(X, mu, Sigma, Sigma_inv = solve(Sigma), det_Sigma = det(Sigma)) {
  Y <- t(t(X) - mu)
  
  return(sum(diag(-0.5 * t(Y) %*% Y %*% Sigma_inv)) - 
         0.5 * (prod(dim(X)) * log(2 * pi) +
         dim(X)[1] * log(det_Sigma))) 
}

# Uses a Metropolis-Hastings algorithm to produce approximations
# of the joint posterior (theta, sigma, and beta's).
# Reference: Hoff, page 203
# Parameters:
#   group: Which variable should be used for groups?
#   response: What is the response variable?
#   mu_0, Lambda_0: Hyperparameters for theta
#   eta_0, S_0: Hyperparameters for Sigma
#   tau: Proposal variance for betas
mh_mixed_logistic <- function(data, group, response, mu_0, Lambda_0, 
                              eta_0, S_0, 
                              tau = 1.5, num_iter = 10000) {
  
  # Assuming 1 col for response, 1 col for group, rest for prediction
  groups <- data[group] %>% unique()
  m <- length(groups)
  p <- ncol(data) - 2
  
  thetas <- matrix(nrow = num_iter, ncol = p + 1) 
  Sigmas <- rep(list(matrix(NA, p + 1, p + 1)), num_iter)
  betas <- matrix(nrow = m, ncol = p + 1)
  
  # Useful invariant information
  Lambda_0_inv <- solve(Lambda_0)
  
  for(s in seq_len(num_iter - 1)) {
    
    # Update theta via its full conditonal
    Lambda_m <- solve(Lambda_0_inv + m * solve(Sigmas[[s]]))
    mu_m <- Lambda_m %*% (Lambda_0_inv %*% mu_0 + m * solve(Sigmas[[s]]) %*% matrix(colMeans(betas)))
    thetas[s + 1] <- mvrnorm(1, mu_m, Lambda_m)
    
    # Update Sigma via its full conditional
    S_theta <- apply(betas, 1, function(beta) (beta - thetas[s + 1,]) %*% t(beta - thetas[s + 1,])) %>% sum()
    Sigma_inv <- rWishart(1, eta_0 + m, solve(S_0 + S_theta))[,,1]
    Sigmas[[s + 1]] <- solve(Sigma_inv)
    
    det_Sigma <- det(Sigmas[[s + 1]])
    # Update Beta's
    for(j in seq_len(m)) {
      y_j <- data %>% filter(group == groups[j]) %>% select(response)
      x_j <- data %>% filter(group == group[j]) %>% select("percentms")
      
      # Proposal
      beta_star <- mvrnorm(1, betas[j,], tau * Sigmas[[s + 1]])
      
      # Calculate acceptance rate (log ratio!)
      # dbinom is actually the sum of several Bernouli variables
      log_r <- sum(dbinom(y_j, 1, plogis(beta_star[1] + beta_star[2] * x_j), log = T) -
                   dbinom(y_j, 1, plogis(betas[j, 1] + betas[j, 2] * x_j), log = T)) +
               ldmvnorm(t(beta_star), thetas[s + 1], Sigmas[[s + 1]], Sigma_inv = Sigma_inv, det_Sigma = det_Sigma) - 
               ldmvnorm(t(betas[j,]), thetas[s + 1], Sigmas[[s + 1]], Sigma_inv = Sigma_inv, det_Sigma = det_Sigma)
      
      # Update if appropriate
      betas[s + 1,] <- ifelse(log(runif(1)) < log_r, beta_star, betas[s,])
    }
  }
  
  return(list(theta = thetas, sigma = Sigmas, beta = betas))
}
