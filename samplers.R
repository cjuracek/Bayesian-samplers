library(MASS)

# 3 helper functions for MVN sampler
compute_mu_n <- function(lambda_0, n, sigma, mu_0, y_bar) {
  lambda_0_inv <- solve(lambda_0)
  sigma_inv <- solve(sigma)
  
  first_factor <- (lambda_0_inv + n * sigma_inv) %>% solve()
  second_factor <- lambda_0_inv %*% mu_0 + n * sigma_inv %*% y_bar
  
  return(first_factor %*% second_factor)
}

compute_lambda_n <- function(lambda_0, n, sigma) {
  lambda_n <- (solve(lambda_0) + n * solve(sigma)) %>% solve()
  return(lambda_n)
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
  
  # Initial 
  thetas[1,] <- mu_0
  for(i in seq_len(num_iter - 1)) {
    
    # Update Sigma via full conditional
    S_n <- compute_S_n(Y, thetas[i,], S_0)
    sigmas[[i + 1]] <- solve(rWishart(1, nu_0 + n, solve(S_n))[,,1])
    
    # Update theta via full conditional
    mu_n <- compute_mu_n(Lambda_0, n, sigmas[[i + 1]], mu_0, y_bar)
    Lambda_n <- compute_lambda_n(Lambda_0, n, sigmas[[i + 1]])
    thetas[i + 1,] <- mvrnorm(1, mu_n, Lambda_n)
    
  }
  
  return(list(theta = thetas, sigma = sigmas))
}