gibbs_mvn <- function(Y, p = ncol(Y), n = nrow(Y), num_iter = 10000,
                       nu_0 = 4, S_0 = cov(Y), Lambda_0 = cov(Y),
                       mu_0 = matrix(apply(Y, 2, mean), p, 1)) {
  
  # https://stackoverflow.com/questions/19340079/how-to-declare-list-object-with-m-elements
  thetas <- rep(list(matrix(NA, p, 1)), num_iter)
  sigmas <- rep(list(matrix(NA, p, p)), num_iter)
  y_bar <- matrix(apply(Y, 2, mean), p, 1)
  
  for(i in 1:num_iter) {
    
    # Full conditional for Theta
    if(i == 1) {  # Generate intial Theta_0
      thetas[[i]] <- mu_0
    }
    else{
      mu_n <- compute_mu_n(Lambda_0, n, sigmas[[i - 1]], mu_0, y_bar)
      Lambda_n <- compute_lambda_n(Lambda_0, n, sigmas[[i - 1]])
      thetas[[i]] <- mvrnorm(1, mu_n, Lambda_n)
    }
    
    # Full conditional for Sigma
    S_n <- compute_S_n(Y, thetas[[i]], S_0)
    sigmas[[i]] <- solve(rWishart(1, nu_0 + n, solve(S_n))[,,1])
    
  }
  
  #return(thetas)
  return(c(thetas, sigmas))
}