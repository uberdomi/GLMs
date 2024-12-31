library(tidyverse)

family_functions <- function(family) {
  functions <- list(b = function(x){return(0)}, # canonical function
                    b_p = function(x){return(0)}, # first derivative
                    b_pp = function(x){return(0)}, # second derivative
                    
                    g = function(x){return(0)}, # link function
                    g_inv = function(x){return(0)}, # inverse link function
                    
                    h = function(x){return(0)}) # inverse mean function
  
  if(family == "gaussian") {
    
    functions$b <- function(x){return(x^2/2)}
    functions$b_p <- function(x){return(x)}
    functions$b_pp <- function(x){return(rep(1, length(x)))}
    
    functions$g <- function(x){return(x)}
    functions$g_inv <- function(x){return(x)}
    
    functions$h <- function(x){return(x)}
    
  } else if(family == "binomial") {
    
    functions$b <- function(x){return(log(1+exp(x)))}
    functions$b_p <- function(x){return(exp(x)/(1+exp(x)))}
    functions$b_pp <- function(x){return(exp(x)/(1+exp(x))^2)}
    
    functions$g <- function(x){return(exp(x)/(1+exp(x)))}
    functions$g_inv <- function(x){return(log(x/(1-x)))}
    
    functions$h <- function(x){return(log(x/(1-x)))}
    
  } else if(family == "poisson") {
    
    functions$b <- function(x){return(exp(x))}
    functions$b_p <- function(x){return(exp(x))}
    functions$b_pp <- function(x){return(exp(x))}
    
    functions$g <- function(x){return(exp(x))}
    functions$g_inv <- function(x){return(log(x))}
    
    functions$h <- function(x){return(log(x))}
    
  } else if(family == "gamma") {
    
    functions$b <- function(x){return(-log(x))}
    functions$b_p <- function(x){return(-1/x)}
    functions$b_pp <- function(x){return(1/x^2)}
    
    functions$g <- function(x){return(log(x))}
    functions$g_inv <- function(x){return(exp(x))}
    
    functions$h <- function(x){return(-1/x)}
    
  } else {
    stop("Invalid family")
  }

  return(functions)
}

# Implementing the Fisher scoring iterative method:
fisher_scoring <- function(X, y, family, beta, max_iter = 100, tol = 1e-6) {
  # Initialize the number of observations
  n <- nrow(X)
  if(length(y) != n) {
    stop("The length of the response vector is not equal to the number of observations")
  }
  
  # Initialize the number of predictors
  p <- ncol(X)
  if(length(beta) != p) {
    stop("The length of the beta vector is not equal to the number of predictors")
  }
  
  # Initialize the weights
  w <- rep(1, n)
  
  # Get the family functions
  functions <- family_functions(family)
  
  b <- functions$b
  b_p <- functions$b_p
  b_pp <- functions$b_pp
  g <- functions$g
  g_inv <- functions$g_inv
  h <- functions$h
  
  # Initialize the iteration counter
  iter <- 0
  # Initialize the difference between the current and previous beta
  diff_beta <- Inf
  # Initialize the beta vector
  beta_new <- beta
  
  # Start the iteration
  while(iter < max_iter & diff_beta > tol) {
    iter <- iter + 1
    print(str_glue("Iteration: {iter} out of {max_iter}"))

    # Store the current beta
    beta_old <- beta_new
    # Compute the linear predictor
    eta <- X %*% beta_old
    # Compute the mean
    mu <- g_inv(eta)
    # Compute the canonical parameter
    theta <- h(mu)
    
    # Compute the weights
    w <- b_pp(theta)
    # Compute the diagonal weight matrix
    W <- diag(w)
    
    # print(str_glue("Dimensions: n = {n}, p = {p}, mu = {length(mu)}, theta = {length(theta)}, w = {length(w)}, W1 = {nrow(W)}, W2 = {ncol(W)}"))
    

    # Compute the gradient
    gradient <- t(X) %*% (y - mu)
    # Compute the Hessian
    hessian <- t(X) %*% W %*% X
    
    # Compute the update
    update <- solve(hessian) %*% gradient
    print("Hessian")
    print(hessian)
    
    print("Gradient")
    print(gradient)
    
    print("Update")
    print(update)
    # Update the beta vector
    beta_new <- beta_old + update
    
    # Compute the difference between the current and previous beta
    diff_beta <- sum((update)^2)
    print(str_glue("Norm squared of the update: {diff_beta}"))
  }
  
  str_glue("The algorithm converged after {iter} iterations") %>% print()
  return(beta_new)
}

set.seed(123)
# Test the functions
beta <- c(1, 2, -3)

p <- length(beta)
n <- 50

X <- matrix(rnorm(n*p), ncol = p)

y <- X %*% beta + rnorm(n, sd = 1)

beta_est <- fisher_scoring(X, y, "gaussian", rep(0, p))


functions <- family_functions("gaussian")

b <- functions$b
b_p <- functions$b_p
b_pp <- functions$b_pp
g <- functions$g
g_inv <- functions$g_inv
h <- functions$h

y_est <- g_inv(X %*% beta_est)

# Compare the estimated and true 'y' values:
p <- ggplot() +
  geom_point(aes(x = y, y = y_est)) +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  labs(title = "Estimated vs True 'y' values") +
  theme_minimal()

print(p)

