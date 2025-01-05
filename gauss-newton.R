library(tidyverse)
exam_data <- read_csv("data/higher_education_data.csv")

#' Helper function to turn vector of parameters into list
#' @param Vector of parameters
#' @return List of parameters
param_extractor <- function(param) {
  mew <- param[1]
  if(mew < 0) {
    stop("mew must be greater than 0")
  }
  
  beta <- param[2]
  if(beta < 0) {
    stop("beta must be greater than 0")
  }
  
  gamma <- param[3]
  
  delta <- param[4]
  if(delta < 0) {
    stop("delta must be greater than 0")
  }
  
  return(
    list(
      mew = mew,
      beta = beta,
      gamma = gamma,
      delta = delta
    )
  )
}

#' Helper function to convert parameters to vector
#' @param List of parameters
#' @return Vector of parameters
param_as_vector <- function(param) {
  return(
    c(param$mew, param$beta, param$gamma, param$delta)
  )
}

#' Function to produce predictions
#' @param Model parameters
#' @param Data
#' @param Whether conversion to parameter list is necessary
#' @return Prediction
f <- function(param, data, extract_param = FALSE) {
  # Extract parameters if necessary
  if(extract_param) {
    param <- param_extractor(param)
  }
  
  # Plug in parameters and data into model
  return(
    param$mew + param$beta * exp(-(param$delta + param$gamma * data$s) * data$t)
  )
}

#' Objective function: Sum of squared errors
#' @param List of parameters
#' @param Data
#' @return Sum of squared errors
g <- function(param, data) {
  # Calculate predictions
  y_hat <- f(param, data)
  
  # Return sum of squared errors
  return(
    sum((data$y - y_hat) ^ 2)
  )
}

#' Jacobian of the objective function with respect to the parameters
#' @param List of parameters
#' @param Data
#' @return Jacobian matrix
jacobian_matrix_generator <- function(param, data) {
  # Residual is y_i - f
  # = y_i - mew - beta * exp(-(delta + gamma * s_i) * t_i)
  # Return partial derivatives of above with respect to each parameter
  
  # Wrt mew is -1
  d_mu <- rep(-1, nrow(data))
  
  # Remaining partial have a common term
  exp_term <- exp(-(param$delta + param$gamma * data$s) * data$t)
  
  # Wrt beta is -exp_term
  d_beta <- -exp_term
  
  # Wrt gamma is beta * s * t * exp_term
  d_gamma <- param$beta * data$s * data$t * exp_term
  
  # Wrt delta is beta * t * exp_term
  d_delta <- param$beta * data$t * exp_term
  
  # Combine partials into matrix
  return(
    cbind(d_mu, d_beta, d_gamma, d_delta)
  )
}

#' Gauss-Newton minimiser
#' @param Initial parameters
#' @param Data
#' @param Tolerance
#' @param Minimum number of iterations
#' @param Maximum number of iterations
#' @return Optimised parameters
gauss_newton_minimiser <- function(param, data, tolerance = 1e-6, min_iteration = 2, max_iteration = 1000) {
  iteration <- 1
  # starting values of parameters
  param <- param_extractor(param)
  all_errors <- numeric(max_iteration)
  
  # If maximum number of iterations exceeded, break
  while(iteration - 1 < max_iteration) {
    # Get jacobian matrix for all parameters
    jacobian <- jacobian_matrix_generator(param, data)
    
    # Compute residuals
    residuals <- data$y - f(param, data, extract_param = FALSE)
    
    # Update parameters
    param_vector <- param_as_vector(param)
    param_vector <- param_vector +
      solve(
        t(jacobian) %*% jacobian,
        t(jacobian) %*% residuals
      )
    param <- param_extractor(param_vector)
    
    # Objective function value for these parameters
    sum_squared_errors <- g(param, data)
    all_errors[iteration] <- sum_squared_errors
    
    print(jacobian)
    print(residuals)
    print(iteration)
    print(param_vector)
    print(all_errors[iteration])
    print("\n")
    # Stopping condition: if improvement in objective function is less than tolerance, break
    if (iteration > min_iteration && abs(all_errors[iteration] - all_errors[iteration - 1]) < tolerance) {      break
      break
    }
    
    iteration <- iteration + 1
  }
  
  return(param)
}

# What are reasonable starting values for the parameters, assuming beta = 1.1?

# Assuming beta = 1.1
starting_beta <- 1.1

# mew: baseline improvement given a beta of zero (intercept)
# starting value: mean of y at t = 1, where beta is likely to be zero
# = 2.070135

starting_mew <- mean(exam_data$y[exam_data$t == 1])

# gamma: effect of higher mathematics
# starting value: difference in rates of decline between no higher and higher mathematics
# = 0.02530709
starting_gamma <- abs(mean(exam_data$y[exam_data$s == 0]) - mean(exam_data$y[exam_data$s == 1]))

# delta: baseline rate of decline, assuming no higher education
# starting value: rate of decline between t = 1 and t = 10

decline <- abs(mean(exam_data$y[exam_data$t == 10]) - mean(exam_data$y[exam_data$t == 1]))

# decline = 1.1 * exp(-delta * (10 - 1))
# decline = 1.1 * exp(-delta * 9)
# decline / 1.1 = exp(-delta * 9)
# log(decline / 1.1, base = e) = -delta * 9
# -log(decline / 1.1, base = e) / 9 = delta

starting_delta <- -log(decline / 1.1) / 9

starting_param <- c(starting_mew, starting_beta, starting_gamma, starting_delta)

# Inputting starting parameters into Gauss-Newton minimiser
gauss_newton_minimiser(starting_param, exam_data, max_iteration = 2)

