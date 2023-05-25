# Beverton-Holt Model

inv.plogis <- binomial()$linkfun

# Function to simulate Beverton-Holt model
beverton_holt <- function(
    N0, timesteps,
    b0_e = 0.0, b1_e = 0.0,
    b0_k = 100, b1_k = 0.3, b1_g = 0.0, 
    b0_r = 1, b1_recr = 0.1, b2_recr = -0.001, opt_x1 = 1,
    distP = 0.01, Nrep = 1
    ) {
  
  out_df <- data.frame(
    rep = integer(), t = integer(), N=numeric(), x1=numeric(), K=numeric(), g=numeric(), 
    recr = numeric(), recr_mean = numeric(), disturbance = integer(), 
    actualR = numeric())
  
  for(rep_i in 1:Nrep){
    
    # Create vectors to store N sizes and time steps
    N <- numeric(timesteps)
    
    # Set initial N size
    N[1] <- N0
    recr = 1

    for (t in 2:timesteps) {
      #if(t < cc_range[1]) 
      
      x1 = plogis(b0_e + b1_e*(t - opt_x1)^2+runif(1,0,1))
      K = b0_k + b1_k * x1
      g = plogis(b1_g * x1)*2
    
      #recr_mean = b0_r*exp(b1_recr*x1 + b2_recr*N[t] - 1)
      
      # Simulate N dynamics
      disturbance = rbinom(1,1,1-distP)
      
      N[t] <- (g * N[t-1]) / (1 + (N[t-1] / K))*disturbance
      
      recr = max(recr, 1)
      K_r = 10
      g_recr = plogis(b0_r + b1_recr*x1+b2_recr*N[t])*2
      recr_mean = (g_recr * recr) / (1 + (recr / K_r))
      recr = rpois(1, recr_mean)  
      
      N[t] = N[t] + recr
      
      out_df <- rbind(
        out_df,
        data.frame(rep=rep_i, t=t, N=N[t], x1, K, g, recr_mean, recr,disturbance, actualR = N[t-1]/N[t])
      )
      
    }
  }

  # Return the N sizes and time steps
  return(out_df)
}

# Run the Beverton-Holt model simulation
sim_df <- beverton_holt(N0 = 10, timesteps = 500, 
  b0_e = 1, b1_e = 0.0,
  b0_k = 100, b1_k = 4.9, b1_g = 1.4, b0_r = 30, b1_recr = 8.1, b2_recr = -3.1, 
  distP = 0.02, Nrep = 1, opt_x1 = 250
)

#sim_df <- aggregate(. ~ t, data = sim_df, mean, na.rm = TRUE)

par(mfrow =c(3,1),mar=c(4,4,0.5,0.5))
plot(N~t, data = sim_df, type = "l")
abline(h = 100)
abline(v = sim_df[sim_df$disturbance == 0,]$t, col = "orange", lty = 2)
legend("bottomright",legend = c("disturbance"), 
       col = c("orange"), lty = 2, lwd = c(2), cex = 0.7)
# abline(v = c(seq(0,nrow(sim_df), 100)), col = "grey50", lty = 2)
plot(actualR~t, data = sim_df, type = "l", ylim = c(-0.3,30), lwd = 2)
# abline(v = c(seq(0,nrow(sim_df), 100)), col = "grey50", lty = 2)
abline(v = sim_df[sim_df$disturbance == 0,]$t, col = "orange", lty = 2)
lines(g~t, data = sim_df, col = "blue")
lines(recr_mean~t, data = sim_df, col = "green")
lines(K~t, data = sim_df, col = "red")

legend("topright",legend = c("actual rate of change","growth","recr","K"), 
       col = c("black","blue","green","red"), lty = 1, lwd = c(2,1,1,1), cex = 0.7)
legend("bottomright",legend = c("disturbance"), 
       col = c("orange"), lty = 2, lwd = c(2), cex = 0.7)

plot(x1~t, data = sim_df, type = "l", lwd = 2, ylim = c(0,1))

