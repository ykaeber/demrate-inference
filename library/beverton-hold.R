# Beverton-Holt Model

inv.plogis <- binomial()$linkfun

# Function to simulate Beverton-Holt model
beverton_holt <- function(
    N0, timesteps, spinup = 0, sample_interval = 1,
    b0_e = 0.0, b1_e = 0.0,
    b0_k = 100, b1_k = 0.3, 
    b1_g = 0.0, b2_g = 0,
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
      x1 = b0_e + rnorm(1,0,.02)
      x1[x1>1] = 1
      x1[x1<0] = 0
      #x1 = plogis(b0_e + b1_e*(t - opt_x1)^2+runif(1,0,1))
      K = b0_k + b1_k * x1
      g = plogis(b1_g * x1 + b2_g*N[t-1])*2
    
      #recr_mean = b0_r*exp(b1_recr*x1 + b2_recr*N[t] - 1)
      
      # Simulate N dynamics
      disturbance = rbinom(1,1,1-distP)
      
      N[t] <- (g * N[t-1]) / (1 + (N[t-1] / K))*disturbance
      
      recr = max(recr, 1)
      K_r = 10
      g_recr = plogis(b0_r + b1_recr*x1+b2_recr*N[t])*2
      if(g_recr < 0.0001) g_recr = 0
      recr_mean = (g_recr * recr) / (1 + (recr / K_r))
      
      recr = rpois(1, recr_mean)  
      
      N[t] = N[t] + recr
      
      out_df <- rbind(
        out_df,
        data.frame(rep=rep_i, t=t, N=N[t], x1, K, g, recr_mean, recr,disturbance, actualR = N[t-1]/N[t])
      )
      
    }
  }

  sample_vec <- seq(spinup, timesteps-1, sample_interval)
  # Return the N sizes and time steps
  return(out_df[sample_vec,])
}

# Run the Beverton-Holt model simulation
sim_df <- beverton_holt(N0 = 10, timesteps = 500, 
  b0_e = 1, b1_e = 0.0,
  b0_k = 100, b1_k = 4.9, b1_g = 1.4, b0_r = 30, b1_recr = 8.1, b2_recr = -3.1, 
  distP = 0.02, Nrep = 1, opt_x1 = 250
)

#sim_df <- aggregate(. ~ t, data = sim_df, mean, na.rm = TRUE)

plot_bh <- function(sim_df){
  par(mfrow =c(3,1),mar=c(2,4,1.2,1))
  plot(N~t, data = sim_df, type = "l")
  title("a)", adj = 0, line = 0.4, cex = 10)
  abline(h = 100)
  abline(v = sim_df[sim_df$disturbance == 0,]$t, col = "orange", lty = 2)
  legend("bottomright",legend = c("disturbance"), 
         col = c("orange"), lty = 2, lwd = c(2), cex = 0.7)
  # abline(v = c(seq(0,nrow(sim_df), 100)), col = "grey50", lty = 2)
  plot(actualR~t, data = sim_df, type = "l", ylim = c(-0.3,3), lwd = 2, ylab = "rate")
  title("b)", adj = 0, line = 0.4, cex = 10)
  # abline(v = c(seq(0,nrow(sim_df), 100)), col = "grey50", lty = 2)
  abline(v = sim_df[sim_df$disturbance == 0,]$t, col = "orange", lty = 2)
  lines(g~t, data = sim_df, col = "blue")
  lines(recr_mean~t, data = sim_df, col = "green")
  #lines(K~t, data = sim_df, col = "red")
  
  legend("topright",legend = c("actual rate of change","growth","r mean"), 
         col = c("black","blue","green"), lty = 1, lwd = c(2,1,1,1), cex = 0.7)
  legend("bottomright",legend = c("disturbance"), 
         col = c("orange"), lty = 2, lwd = c(2), cex = 0.7)
  par(mar=c(4,4,1.2,1))
  plot(recr~t, data = sim_df, type = "l", lwd = 2, col = "green", ylab = "rpois(r)")
  title("c)", adj = 0, line = 0.4, cex = 10)
}
data.frame(
  par = c("distP", "b0_e","b1_recr", "b2_recr", "b1_g","b2_g"),
  default = c(0.01, 0.5, 20, -0.2, 3, -0.1),
  min = c(0, 0, 0, -5, 0, -5),
  max = c(0.2, 1, 50, 0, 10,0)
)

set.seed(1)
sim_df <- beverton_holt(
  N0 = 0,
  timesteps = 250,
  spinup = 10L,
  b0_e = 0.8,
  b1_e = 0.0,
  b0_k = 100,
  b1_k = 0,
  b1_g = 5,
  b2_g = -0.1,
  b0_r = 30,
  b1_recr = 0.1,
  b2_recr = -0.3,
  distP = 0.02,
  Nrep = 1,
  opt_x1 = 250)
# sim_df <- beverton_holt(
#   N0 = 0, timesteps = 250, spinup = 0,
#   b0_e = 0.5, b1_e = 0,
#   b0_k = 100, b1_k = 0, 
#   b1_g = 10, b2_g = -0.1,
#   b0_r = inv.plogis(0.0001), b1_recr = 20, b2_recr = -0.2,
#   distP = 0.02, Nrep = 1, opt_x1 = 250)
plot_bh(sim_df)
tiff("figures/bh-simple.tiff", units = "in", width = 5, height = 5, res = 300, compression = "lzw")
plot_bh(sim_df)
dev.off()

inv.plogis(0)
head(1)
