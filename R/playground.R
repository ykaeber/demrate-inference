
ricker <- function(
    K, N0, num_generations, distP=0.05,
    b0_m = 0.01, b1_m = 0, b2_m = 0,
    b0_g = 0.1,  b1_g = 0, b2_g = 0,
    b0_r = 0.2,  b1_r = 0, b2_r = 0  
    ) {
  N <- numeric(num_generations)  # Create a vector to store population sizes
  
  N[1] <- N0  # Set initial population size
  
  out_df <- data.frame(
    t=integer(), N=integer(), x1=numeric(), r=numeric(), actualr=integer(), m=numeric(), 
    g=numeric(), R=numeric(), disturbance=numeric()
  )
  for (t in 2:num_generations) {
    x1 = runif(1,0,1)
    m = b0_m + b1_m*x1 + b2_m*N[t]
    g = b0_g + b1_g*x1 + b2_g*N[t]
    r = b0_r + b1_r*x1 + b2_r*N[t]
    actualr = rpois(1, max(exp(r * (1 - N[t - 1] / K))-1, 0))
    R = g
    disturbance = rbinom(1,1,1-distP)
    N[t - 1] = N[t - 1]*(1-m)
    #print(paste0("x=",x1," r=",r," m=",m," g=",g," R=",R, " disturbance=",disturbance))
    N[t] <-  N[t - 1] *  exp(R * (1 - N[t - 1] / K))*disturbance + actualr
    
    out_df <- rbind(
      out_df,
      data.frame(t=t, N=N[t], x1=x1, r=r, actualr=actualr, m=m, g=g, R=R, disturbance=disturbance)
    )
  }
  return(out_df)
}


# Set model parameters
K <- 1000    # Carrying capacity
N0 <- 10    # Initial population size
num_generations <- 500 # Number of generations to simulate



sim_df <- ricker(K = 100, N0 = 10L, num_generations = 1000L, distP = 0.01,
                 b0_m = 0.00, b1_m = 0, b2_m = 0,
                 b0_g = .1,  b1_g = 0, b2_g = -.3,
                 b0_r = 0.01,  b1_r = 0, b2_r = 0)

head(sim_df)
sim_df$actualR <- c(NA, sim_df$N[-1]/sim_df$N[-length(sim_df$N)])


par(mfrow =c(2,1),mar=c(4,4,0.5,0.5))
plot(N~t, data = sim_df, type = "l", ylim = c(0, 120))
abline(h = 100)
abline(v = sim_df[sim_df$disturbance == 0,]$t, col = "orange", lty = 2)
legend("bottomright",legend = c("disturbance"), 
       col = c("orange"), lty = 2, lwd = c(2), cex = 0.7)
# abline(v = c(seq(0,nrow(sim_df), 100)), col = "grey50", lty = 2)
plot(actualR~t, data = sim_df, type = "l", ylim = c(-0.3,3), lwd = 2)
# abline(v = c(seq(0,nrow(sim_df), 100)), col = "grey50", lty = 2)
abline(v = sim_df[sim_df$disturbance == 0,]$t, col = "orange", lty = 2)
lines(m~t, data = sim_df, col = "red")
lines(g~t, data = sim_df, col = "blue")
lines(r~t, data = sim_df, col = "green")
legend("topright",legend = c("actual rate of change","mortality","growth","recruitment"), 
       col = c("black","red","blue","green"), lty = 1, lwd = c(2,1,1,1), cex = 0.7)
legend("bottomright",legend = c("disturbance"), 
       col = c("orange"), lty = 2, lwd = c(2), cex = 0.7)


