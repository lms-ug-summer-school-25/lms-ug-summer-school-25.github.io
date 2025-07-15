# # # # # # # # # # # # # # # # # # # # 
# # # # # #   Part 1.3   # # # # # # #
# # # # # # # # # # # # # # # # # # # 


# This is a comment. R ignores anything after a '#' on a line.

# Assigning values to variables
x <- 10
y = 5 # ' = ' also works for assignment, but '<-' is preferred in R
print(x + y)

# Vectors (ordered collections of values of the same type)
my_vector <- c(1, 2, 3, 4, 5) # 'c()' combines values into a vector
print(my_vector)
print(my_vector * 2)
print(sum(my_vector))

# Functions
my_function <- function(a, b) {
  result <- a * b
  return(result)
}
print(my_function(x, y))

# Data frames (like tables or spreadsheets)
my_data <- data.frame(
  Name = c("Alice", "Bob", "Charlie"),
  Age = c(25, 30, 22),
  City = c("NY", "LA", "Chicago")
)
print(my_data)
print(my_data$Age) # Accessing a column by name
print(my_data[1, ]) # Accessing the first row



# # # # # # # # # # # # # # # # # # # # 
# # # # # #   Part 1.4   # # # # # # #
# # # # # # # # # # # # # # # # # # # 

# Install packages (only need to do this once per R installation)
install.packages("ggplot2") # For nice plots
install.packages("coda")    # For MCMC diagnostics

# Load packages (need to do this every time you start a new R session and want to use them)
library(ggplot2)
library(coda)


# # # # # # # # # # # # # # # # # # # # 
# # # # # #   Part 1.5   # # # # # # #
# # # # # # # # # # # # # # # # # # # 

# Simple plot using R's base plotting system
# 'plot() ' is a versatile function, 'type = "l"' means draw lines.
# 'main ', 'xlab ', 'ylab ' are for setting the title and axis labels.
plot(my_vector, type = "l", main = "My Simple Line Plot",
     xlab = "Index", ylab = "Value")

# # # # # # # # # # # # # # # # # # # # 
# # # # # #   Part 2.2   # # # # # # #
# # # # # # # # # # # # # # # # # # # 

# Define the Stochastic SIS model simulation function
# This function takes:
# N: Total population size
# beta: Transmission rate
# gamma: Recovery rate
# I0: Initial number of infected individuals
# T: Total number of time steps to simulate

simulate_SIS <- function(N, beta, gamma, I0, T) {
  # Initialize vectors to store the number of Susceptible (S) 
  # and Infected (I) individuals
  S <- numeric(T)
  I <- numeric(T)
  
  # Initialize vectors to store the *number of new infections* and 
  # *new recoveries* that occur at each time step. These will be  
  # crucial for calculating the likelihood function.
  new_inf_history <- numeric(T)
  new_rec_history <- numeric(T)
  
  # Set initial conditions for S and I at time t=0
  S[1] <- N - I0
  I[1] <- I0
  new_inf_history[1] <- 0 # No new events at time 0
  new_rec_history[1] <- 0 # No new events at time 0
  
  # Loop through each time step from t=1 to T-1
  for (t in 1:(T-1)) {
    # Calculate the probability of infection for a susceptible 
    # individual based on the current number of infected individuals
    p_inf <- 1 - exp(-beta*I[t]/N)
    
    # Calculate the probability of recovery for an infected individual
    p_rec <- 1 - exp(-gamma)
    
    # Simulate number of events in this time step using binomial 
    # distribution
    # rbinom(n, size, prob) generates n random values from a binomial 
    # distribution size is the number of trials, prob is the 
    #  probability of  success on each trial
    new_inf <- rbinom(1, S[t], p_inf)
    new_rec <- rbinom(1, I[t], p_rec)
    
    # Store the simulated number of events for this time step
    new_inf_history[t+1] <- new_inf
    new_rec_history[t+1] <- new_rec
    
    # Update the number of Susceptible and Infected individuals for the 
    # next time step (t+1)
    # S decreases by new infections and increases by new recoveries
    S[t + 1] <- S[t] - new_inf + new_rec
    # I increases by new infections and decreases by new recoveries
    I[t + 1] <- I[t] + new_inf - new_rec
  }
  # Return a data frame
  return(data.frame(time = 0:(T-1), S = S, I = I,
                    new_inf_obs = new_inf_history, new_rec_obs = new_rec_history))
}


# # # # # # # # # # # # # # # # # # # # 
# # # # # #   Part 2.3   # # # # # # #
# # # # # # # # # # # # # # # # # # # 

# Set a seed for reproducibility of stochastic results.
# This ensures that if you run the code multiple times with the 
# same seed, you will get the exact same sequence of random 
# numbers, making your results reproducible.
set.seed(3)

# Simulate data using the 'simulate_SIS' function.
# N: Total population size
# beta: True transmission rate
# gamma: True recovery rate
# I0: Initial number of infected individuals
# T: Total number of time steps (e.g., days) for the simulation
SISdata <- simulate_SIS(N = 1000, beta = 0.3, gamma = 0.1, I0 = 10, 
                        T = 150)

# View the first few rows of the simulation output.
head(SISdata)

# Plotting the simulation results
plot(SISdata$time, SISdata$I, type = "l", col = "blue", 
     ylab = "Number of individuals", xlab = "Time", ylim = c(0, 1000))
lines(SISdata$time, SISdata$S, col = "green")

# Add legend
legend("topright", legend = c("Infectious", "Susceptible"),
       lty = c(1, 1), col = c("blue", "green"))


# # # # # # # # # # # # # # # # # # # # 
# # # # # #   Part 3.1   # # # # # # #
# # # # # # # # # # # # # # # # # # # 

# True parameters (these are what we want to estimate)
true_beta <- 0.3
true_gamma <- 0.1

# Initial conditions and time for data generation
N_data <- 1000
I0_data <- 5
S0_data <- N_data - I0_data
T_data <- 150

# Simulate the true SIS stochastic model to generate observed data
set.seed(123) # For reproducibility
observed_data <-  simulate_SIS(N = N_data, beta = true_beta,  
                               gamma = true_gamma, I0 = I0_data,  T = T_data)


# # # # # # # # # # # # # # # # # # # # 
# # # # # #   Part 3.2   # # # # # # #
# # # # # # # # # # # # # # # # # # # 

# Define the log-likelihood Function for the SIS model
# This function takes:
# params: A vector of parameters (e.g., c(beta, gamma))
# data: The observed data.

loglik_SIS <- function(params, data) {
  beta = params[1]
  gamma = params[2]
  
  # Ensure parameters are valid (e.g., positive).
  # If any parameter is non-positive, return -Inf for the 
  # log-likelihood.
  # This tells the optimizer that these parameters are invalid.
  if (beta <= 0 || gamma <= 0) {
    return(-Inf) 
  }
  
  loglik <- 0 # Initialize the total log-likelihood to zero
  # Calculate the total population size from the initial states
  N <- data$S[1] + data$I[1] 
  
  # Loop through each time step in the observed data
  for (t in 1:(nrow(data) - 1)) {
    
    # Get the number of Susceptible and Infected individuals 
    # at time t and t+1
    St <- data$S[t]
    It <- data$I[t]
    
    # Get the observed number of new infections and new recoveries 
    # that occurred between time t and time t+1. These come directly 
    # from our observed data.
    new_inf <-  data$new_inf_obs[t + 1] # Events observed at t+1, 
    # from state at t
    new_rec <- data$new_rec_obs[t + 1] # Events observed at t+1, from 
    # state at t
    
    # Calculate probabilities
    p_inf <- 1 - exp(-beta*It/N)
    p_rec <- 1 - exp(-gamma)
    
    # Avoid probabilities of exactly 0 or 1, which can cause log(0) or 
    # log(1-1) issues in hte dbinom function.
    # We clip them to be slightly away from the boundaries.
    p_inf <- min(max(p_inf, 1e-10), 1 - 1e-10)
    p_rec <- min(max(p_rec, 1e-10), 1 - 1e-10)
    
    # Add log-likelihood contribution for this step
    # dbinom(x, size, prob, log = TRUE) gives the log of the 
    # binomial probability mass function.
    loglik <- loglik +
      dbinom(new_inf, St, p_inf, log = TRUE) +
      dbinom(new_rec, It, p_rec, log = TRUE)
  }
  return(loglik) # Return the total log-likelihood
} 


# # # # # # # # # # # # # # # # # # # # 
# # # # # #   Part 3.3   # # # # # # #
# # # # # # # # # # # # # # # # # # # 

?optim

# Set Initial guesses for parameters
# These are the starting points for the optimizer's search.
initial_params  <- c(beta = 0.2, gamma = 0.05)

# Run the optimization using 'optim()'.
# 'par': The initial guesses for the parameters.
# 'fn': The function to be optimized (our log-likelihood function).
# 'data': Our simulated 'SISdata' that contains the observed S, I, 
# new_inf_obs, new_rec_obs.
# 'control = list(fnscale = -1)': This is crucial! It tells 'optim' 
# to maximize 'fn' instead of minimizing it.
mle_result <- optim(par = initial_params,
                    fn = loglik_SIS, data = observed_data,
                    control = list(fnscale = -1)) 
print("MLE Results (initial run):")
print(mle_result)
# Look for 'convergence = 0' which indicates successful convergence.
# 'par' will give you the estimated parameter values.
# 'value' is the maximum log-likelihood found.

# Extract estimated parameters
estimated_beta_mle <- mle_result$par[1]
estimated_gamma_mle <- mle_result$par[2]

cat("\nEstimated Beta (MLE):", estimated_beta_mle, "\n")
cat("Estimated Gamma (MLE):", estimated_gamma_mle, "\n")
cat("True Beta:", true_beta, "\n")
cat("True Gamma:", true_gamma, "\n")


# # # # # # # # # # # # # # # # # # # # 
# # # # # #   Part 3.4   # # # # # # #
# # # # # # # # # # # # # # # # # # # 

# Define the number of times to repeat the MLE process
num_repetitions <- 25

# Create empty vectors to store the estimated beta and gamma from 
#  each repetition
estimated_betas <- numeric(num_repetitions)
estimated_gammas <- numeric(num_repetitions)

# Set the fixed parameters for data generation for this exercise
fixed_N <- 1000
fixed_I0 <- 10
fixed_T <- 150
fixed_true_beta <- 0.3
fixed_true_gamma <- 0.1

# Loop to repeat the simulation and MLE estimation
for (i in 1:num_repetitions) {
  # Set a unique seed for each repetition to get different stochastic 
  # data realizations
  set.seed(100 + i)
  
  # Simulate new SIS data for this repetition
  current_SISdata <- simulate_SIS(N = fixed_N, beta = fixed_true_beta, 
                                  gamma = fixed_true_gamma,  I0 = fixed_I0, T = fixed_T)
  
  # Run MLE for the current simulated data
  # Use the same initial guesses and bounds as before
  current_mle_result <- optim(par = initial_params, 
                              fn = loglik_SIS, data = current_SISdata,
                              control = list(fnscale = -1))
  
  # Store the estimated parameters
  estimated_betas[i] <- current_mle_result$par[1]
  estimated_gammas[i] <- current_mle_result$par[2]
  
}

# Now, visualize the distribution of your 25 estimates using boxplots.
# Boxplots are great for showing the median, quartiles, and outliers 
# of a distribution.
# Set up a 1x2 plotting layout
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1)) 

# Boxplot for estimated Beta values
boxplot(estimated_betas, main = "Distribution of Estimated Beta",
        ylab = "Estimated Beta", col = "lightblue", border = "black")
# Add a line for the true beta value
abline(h = fixed_true_beta, col = "red", lwd = 2, lty = 2) 

# Boxplot for estimated Gamma values
boxplot(estimated_gammas, main = "Distribution of Estimated Gamma",
        ylab = "Estimated Gamma", col = "lightgreen", border = "black")
# Add a line for the true gamma value
abline(h = fixed_true_gamma, col = "red", lwd = 2, lty = 2) 

# Reset plotting layout
par(mfrow = c(1, 1))

# # # # # # # # # # # # # # # # # # # # 
# # # # # #   Part 4.2   # # # # # # #
# # # # # # # # # # # # # # # # # # # 


# Log-Prior Function (assuming exponential priors)
log_prior <- function(beta, gamma, lambda_beta = 1, lambda_gamma = 1) {
  dexp(beta, rate = lambda_beta, log = TRUE) +
    dexp(gamma, rate = lambda_gamma, log = TRUE)
}


# Metropolis-Hastings Algorithm
MH_sampler_SIS <- function(data, n_iter, beta_init, gamma_init,
                           lambda_beta, lambda_gamma, proposal_sd) {
  
  # Initialize current parameters
  beta <- beta_init
  gamma <- gamma_init
  samples <- matrix(NA, n_iter, 2)
  
  # Counter for accepted proposals
  accepted_count <- 0
  
  for (i in 1:n_iter) {
    beta_prop <- rnorm(1, beta, proposal_sd)
    gamma_prop <- rnorm(1, gamma, proposal_sd)
    
    if (beta_prop > 0 && gamma_prop > 0) {
      loglik_curr <- loglik_SIS(c(beta, gamma), data)
      loglik_prop <- loglik_SIS(c(beta_prop, gamma_prop), data)
      
      logprior_curr <- log_prior(beta, gamma, lambda_beta, 
                                 lambda_gamma)
      logprior_prop <- log_prior(beta_prop, gamma_prop, lambda_beta, 
                                 lambda_gamma)
      
      log_accept_ratio <- (loglik_prop + logprior_prop) - 
        (loglik_curr + logprior_curr)
      
      if (log(runif(1)) < log_accept_ratio) { 
        # accept
        beta <- beta_prop
        gamma <- gamma_prop
        accepted_count <- accepted_count + 1
      }
    }
    
    # Store the current (accepted or re-used) parameters in the chain
    samples[i, ] <- c(beta, gamma)
  }
  
  cat("Acceptance Rate:", accepted_count / n_iter, "\n")
  
  colnames(samples) <- c("beta", "gamma")
  return(as.data.frame(samples))
}


# # # # # # # # # # # # # # # # # # # # 
# # # # # #   Part 4.3   # # # # # # #
# # # # # # # # # # # # # # # # # # # 

# True parameters (these are what we want to estimate)
true_beta <- 0.3
true_gamma <- 0.1

# Initial conditions and time for data generation
N_data <- 1000
I0_data <- 5
S0_data <- N_data - I0_data
T_data <- 150

# Simulate the true SIS stochastic model to generate observed data
set.seed(123) # For reproducibility
observed_data <-  simulate_SIS(N = N_data, beta = true_beta,  
                               gamma = true_gamma, I0 = I0_data,  T = T_data)


# Run the MH sampler
mh_chain <- MH_sampler_SIS(data = observed_data, 
                           n_iter = 10000, beta_init = 0.2,  gamma_init = 0.2, 
                           lambda_beta = 1, lambda_gamma = 1,  proposal_sd = 0.01)



# # # # # # # # # # # # # # # # # # # # 
# # # # # #   Part 4.4   # # # # # # #
# # # # # # # # # # # # # # # # # # # 

# Discard burn-in period (e.g., first 10% of iterations)
burn_in <- n_iter * 0.1
mh_chain_post_burnin <- mh_chain[-(1:burn_in), ]

#install.packages("coda")
library("coda")

# Convert to 'mcmc' object for coda package functions
mcmc_object <- as.mcmc(mh_chain_post_burnin)

# 1. Trace Plots: Show the values of parameters over iterations. 
# Should look like "fuzzy caterpillars".
plot(mcmc_object)

# 2. Summary Statistics: Mean, median, credible intervals
# (similar to confidence intervals).
summary(mcmc_object)

# 3. Autocorrelation Plots: Show correlation between 
# samples at different lags. Should drop quickly.
autocorr.plot(mcmc_object)

# 4. Effective Sample Size (ESS): How many independent samples you 
# effectively have. Higher is better.
effectiveSize(mcmc_object)

# Compare MH estimates (mean of posterior) with MLE and true values
cat("\nMH Posterior Means:\n")
print(colMeans(mh_chain_post_burnin))
cat("\nMLE Estimates:\n")
print(mle_result$par)
cat("\nTrue Parameters:\n")
print(c(beta = true_beta, gamma = true_gamma))

