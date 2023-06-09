
ricker <- function(
    K, N0, timesteps, distP=0.05, env=0, env_trend = 0, 
    spinup = 0, sample_interval = 1,
    b0_m = 0.01, b1_m = 0, b2_m = 0,  
    b0_g = 0, b1_g = 0, b2_g = 0,
    b0_r = 0.2,  b1_r = 0, b2_r = 0  
    ) {
  N <- numeric(timesteps)  # Create a vector to store population sizes
  
  N[1] <- N0  # Set initial population size
  
  out_df <- data.frame(
    t=integer(), N=integer(), env=numeric(), r=numeric(), actualr=integer(), m=numeric(), 
    g=numeric(), R=numeric(), disturbance=numeric()
  )
  
  for (t in 2:timesteps) {
    env = rnorm(1, env + env_trend*t, 0.1)
    g = b0_g + b1_g*env + b2_g * N[t-1]
    m = b0_m * plogis(b1_m*env + b2_m * N[t-1])
    r = b0_r + b1_r*env + b2_r * N[t-1]
    #r = b0_r * exp(b1_r*env + b2_r * N[t-1]/K)
    
    G = exp(g * (1 - N[t - 1] / K))
    M = m * N[t-1]
    R = rpois(1, exp(r * (1 - N[t - 1] / K)))
    disturbance = rbinom(1,1,1-distP)
    
    N[t] = (N[t-1] * G - M)*disturbance + R
    N[t] = max(N[t], 0)
    
    out_df <- rbind(
      out_df,
      data.frame(t=t, N=N[t], env=env, env_trend=env_trend*t, disturbance=disturbance,
                 r=r, R=R, g=g, G=G, m=m, M=M,
                 actualR = N[t-1]/N[t])
    )
  }
  
  sample_vec <- seq(spinup, timesteps-1, sample_interval)
  return(out_df[sample_vec,])
}


# Set model parameters
K <- 1000    # Carrying capacity
N0 <- 10    # Initial population size
timesteps <- 500 # Number of generations to simulate



sim_df <- ricker(K = 100, N0 = 10, timesteps = 10000, distP = 0.005, env = 10,
                 b0_m = 0.01, b1_m = -0.001, b2_m = 0.01,
                 b0_g = 0.02, b1_g = 0.02, b2_g = -0.01,
                 b0_r = 0.001, b1_r = 0.02, b2_r = -0.1)

str(sim_df)
head(sim_df,50)
sim_df$actualR <- c(NA, sim_df$N[-1]/sim_df$N[-length(sim_df$N)])


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
lines(M~t, data = sim_df, col = "red")
lines(G~t, data = sim_df, col = "blue")
lines(exp(r)~t, data = sim_df, col = "green")
legend("topright",legend = c("actual rate of change","mortality","growth","recruitment"), 
       col = c("black","red","blue","green"), lty = 1, lwd = c(2,1,1,1), cex = 0.7)
legend("bottomright",legend = c("disturbance"), 
       col = c("orange"), lty = 2, lwd = c(2), cex = 0.7)

sim_df[6000:6100,]
