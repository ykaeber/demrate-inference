# Beverton-Holt Model

inv.plogis <- binomial()$linkfun

# Function to simulate Beverton-Holt model
beverton_holt <- function(
    N0, timesteps,
    b0_e = 0.0, b1_e = 0.0,
    b0_k = 100, b1_k = 0.3, b1_g = 0.0, 
    b0_r = 1, b1_recr = 0.1, b2_recr = -0.001, 
    distP = 0.01
    ) {
  # Create vectors to store N sizes and time steps
  N <- numeric(timesteps)
  
  # Set initial N size
  N[1] <- N0
  recr = 1
  
  out_df <- data.frame(
    t = integer(), N=numeric(), x1=numeric(), K=numeric(), g=numeric(), 
    recr = numeric(), recr_mean = numeric(), disturbance = integer(), actualR = numeric())
  
  for (i in 2:timesteps) {
    x1 = b0_e + b1_e*i+runif(1,0,1)
    K = b0_k + b1_k * x1
    g = plogis(b1_g * x1)*2
  
    #recr_mean = b0_r*exp(b1_recr*x1 + b2_recr*N[i] - 1)

    K_r = b0_r + b1_recr*x1
    g_recr = plogis(b2_recr * N[i]+0.51)*2
    recr_mean = (g_recr * recr) / (1 + (recr / K_r))
    recr = recr_mean#rpois(1, recr_mean)    
    
    # Simulate N dynamics
    disturbance = rbinom(1,1,1-distP)
    
    N[i] <- (g * N[i-1]) / (1 + (N[i-1] / K))*disturbance + recr
    
    out_df <- rbind(
      out_df,
      data.frame(t=i, N=N[i], x1, K, g, recr_mean, recr,disturbance, actualR = N[i-1]/N[i])
    )
    
  }
  

  # Return the N sizes and time steps
  return(out_df)
}

# Set the parameters
r <- 1.4   # Intrinsic growth rate
K <- 100   # Carrying capacity
N0 <- 10   # Initial N size
timesteps <- 50   # Number of time steps to simulate

# Run the Beverton-Holt model simulation
sim_df <- beverton_holt(N0 = 10, timesteps = 200, 
  b0_e = 0.0, b1_e = 0.0,
  b0_k = 5.1, b1_k = 0.3, b1_g = 1.4, b0_r = 3, b1_recr = 2.01, b2_recr = 0.001, 
  distP = 0.01
)

par(mfrow =c(2,1),mar=c(4,4,0.5,0.5))
plot(N~t, data = sim_df, type = "l")
abline(h = 100)
abline(v = sim_df[sim_df$disturbance == 0,]$t, col = "orange", lty = 2)
legend("bottomright",legend = c("disturbance"), 
       col = c("orange"), lty = 2, lwd = c(2), cex = 0.7)
# abline(v = c(seq(0,nrow(sim_df), 100)), col = "grey50", lty = 2)
plot(actualR~t, data = sim_df, type = "l", ylim = c(-0.3,3), lwd = 2)
# abline(v = c(seq(0,nrow(sim_df), 100)), col = "grey50", lty = 2)
abline(v = sim_df[sim_df$disturbance == 0,]$t, col = "orange", lty = 2)
lines(g~t, data = sim_df, col = "blue")
lines(recr_mean~t, data = sim_df, col = "green")
lines(K~t, data = sim_df, col = "red")
legend("topright",legend = c("actual rate of change","growth","recr","K"), 
       col = c("black","blue","green","red"), lty = 1, lwd = c(2,1,1,1), cex = 0.7)
legend("bottomright",legend = c("disturbance"), 
       col = c("orange"), lty = 2, lwd = c(2), cex = 0.7)

