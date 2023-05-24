# Discrete Ricker model implementation in R

# Function to simulate the discrete Ricker model
ricker <- function(K, N0, num_generations, b0_m = 0.01, b1_m = 0.1, b0_g = 0.1, b1_g = 0.1, b0_r = 0.2, distP=0.05) {
  N <- numeric(num_generations)  # Create a vector to store population sizes
  
  N[1] <- N0  # Set initial population size
  out_df <- data.frame(
    t=integer(), N=integer(), x1=numeric(), r=numeric(), m=numeric(), 
    g=numeric(), R=numeric(), disturbance=numeric()
    )
  for (t in 2:num_generations) {
    x1 = runif(1,0,1)
    m = 1-(b0_m + b1_m*x1)
    g = b0_g + b1_g*x1
    r = rpois(1, max(exp(b0_r * (1 - N[t - 1] / K))-1, 0))
    R = log(g+m)
    disturbance = rbinom(1,1,1-distP)
    N[t] <-  N[t - 1] * exp(R * (1 - N[t - 1] / K))*disturbance + r
    
    out_df <- rbind(
      out_df,
      data.frame(t=t, N=N[t], x1=x1, r=r, m=m, g=g, R=R, disturbance=disturbance)
    )
  }
  return(out_df)
}

# Set model parameters
K <- 1000    # Carrying capacity
N0 <- 10    # Initial population size
num_generations <- 500 # Number of generations to simulate

# Simulate the Ricker model
sim_df <- ricker(K, N0, num_generations, distP = 0.005,
                           b0_m = 0.1, b0_r = 0.1, b0_g = 0.2, 
                           b1_m = 0, b1_g = 0)
head(sim_df)
par(mfrow =c(2,1))
plot(N~t, data = sim_df, type = "l")
plot(R~t, data = sim_df, type = "l", ylim = c(0,2))
lines(m~t, data = sim_df, col = "red")
lines(g~t, data = sim_df, col = "blue")
lines(r~t, data = sim_df, col = "green")
legend("topright",legend = c("mortality","growth","recruitment"), col = c("red","blue","green"), lty = 1)
# Print the population sizes at each generation
plot(NULL, NULL, xlim = c(0, 20), ylim = c(0, K+20))
for (t in 1:num_generations) {
  cat("Generation", t, ":", population_sizes[t], "\n")
  points(t,population_sizes[t])
}



as.vector(matrix(1:10, ncol = 5))
