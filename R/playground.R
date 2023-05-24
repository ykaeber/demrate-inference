ricker <- function(K, N0, num_generations, b0_m = 0.01, b1_m = 0.1, b0_g = 0.1, b1_g = 0.1, b0_r = 0.2, distP=0.05) {
  N <- numeric(num_generations)  # Create a vector to store population sizes
  
  N[1] <- N0  # Set initial population size
  
  for (t in 2:num_generations) {
    x1 = runif(1,0,1)
    x2 = runif(1, 0, 1)
    m = 1-(b0_m + b1_m*x1)
    g = b0_g + b1_g*x2
    r = rpois(1, max(exp(b0_r * (1 - N[t - 1] / K))-1, 0))
    R = g
    disturbance = rbinom(1,1,1-distP)
    N[t - 1] = N[t - 1]*m
    #print(paste0("x=",x1," r=",r," m=",m," g=",g," R=",R, " disturbance=",disturbance))
    N[t] <-  N[t - 1] *  exp(R * (1 - N[t - 1] / K))*disturbance + r
    
  }
  return(N)
}

ricker(K = 100, N0 = 10L, num_generations = 1000L, distP = 0.01,
       b0_m = 0.0, b0_r = 0.1, b0_g = 0.0, 
       b1_m = 0.2, b1_g = 0.4) |> plot(type = "l")



tmp = res$`Growth+Mortality`
tmp$data
plot(tmp$results$obs[,1], tmp$results$pred[,1], main = "Growth")
plot(tmp$results$obs[,2], tmp$results$pred[,2], main = "Mortality")
cor(tmp$results$obs[,1], tmp$results$pred[,1])
cor(tmp$results$obs[,2], tmp$results$pred[,2])


tmp = res$`Recruitment+Disturbance`
tmp$data
plot(tmp$results$obs[,1], tmp$results$pred[,1], main = "Recruitment")
plot(tmp$results$obs[,2], tmp$results$pred[,2], main = "Disturbance")


tmp = res$`Growth+Mortality+Recruitment+Disturbance`
tmp$data
plot(tmp$results$obs[,1], tmp$results$pred[,1], main = "Growth")
plot(tmp$results$obs[,2], tmp$results$pred[,2], main = "Mortality")
plot(tmp$results$obs[,1], tmp$results$pred[,1], main = "Recruitment")
plot(tmp$results$obs[,2], tmp$results$pred[,2], main = "Disturbance")

