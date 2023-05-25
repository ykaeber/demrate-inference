library(sjSDM)
library(reticulate)
torch = sjSDM:::pkg.env$fa$torch

RNN_model = 
  reticulate::PyClass("model", 
                      inherit = torch$nn$Module,
                      defs = list(
                        `__init__` = function(self, n_response = 1L, n_times = 200L) {
                          super()$`__init__`()
                          self$RNN = torch$nn$GRU(input_size = c(2L),  hidden_size = 32L, num_layers = 2L, batch_first = TRUE)
                          self$linear2 = torch$nn$Linear(32L*as.integer(n_times), 20L)
                          self$linear1 = torch$nn$Linear(20L, n_response)
                          NULL
                        },
                        forward = function(self, x) {
                          x = self$RNN(x)
                          x = x[[1]]$reshape(c(x[[1]]$shape[0], -1L))
                          x = self$linear2(x)
                          x = torch$torch$nn$functional$selu(x)
                          x = self$linear1(x)
                          return(x)
                        }
                      ))

train_function = function(data, split = 0.8, epochs = 1000L, n_response = 1L, device = 0L) {
  device = paste0('cuda:', device)
  # Scaling
  X = (data[,-(1:n_response)])
  N_t = X[,1:200]
  X_t = X[,201:400]
  N_t = (N_t - mean(as.vector(N_t)))/sd(as.vector(N_t))
  X = cbind(N_t, X_t)
  Y = data[, 1:n_response, drop=FALSE]
  split_border = round(nrow(data)*split)
  XX = array(X[1:split_border,], dim = c(split_border, 200L,2L))
  model = RNN_model(n_response = as.integer(n_response) )$to(device)
  optimizer = torch$optim$Adamax(model$parameters(), lr = 0.001)
  lambda1 = torch$tensor(0.0008)$to(device)
  lambda2 = torch$tensor(0.001)$to(device)
  DT = 
    torch$utils$data$TensorDataset(torch$tensor(XX, dtype = torch$float32)$to(device), 
                                   torch$tensor(Y[1:split_border,,drop=FALSE], dtype = torch$float32)$to(device))
  DL = torch$utils$data$DataLoader(DT, batch_size = 50L, shuffle = TRUE, pin_memory = FALSE)
  
  for(e in 1:epochs) {
    DATA = reticulate::iterate(DL)
    for(b in 1:length(DATA)){
      x = DATA[[b]][[0]]
      y = DATA[[b]][[1]]
      pred = model(x)
      loss = torch$nn$functional$mse_loss(pred, y, reduce = FALSE)$mean(0L)$sum()
      K = reticulate::iterate( model$parameters() )
      loss = loss$add(torch$cat(sapply(K, function(k) torch$norm( k ,p = 1L)$reshape(c(1L, 1L))))$sum()$mul(lambda1))
      loss = loss$add(torch$cat(sapply(K, function(k) torch$norm( k ,p = 2L)$reshape(c(1L, 1L))))$sum()$mul(lambda2))
      loss$backward()
      optimizer$step()
      optimizer$zero_grad()
    }
    if(e %% 10 ==0) cat(paste0("Epoch: ",e," Loss: ", loss$item(), "\n"))
  }
  XX = array(X[(split_border+1L):nrow(X),], dim = c(length((split_border+1L):nrow(X)), 200L,2L))
  XT = torch$tensor(XX, 
                    dtype = torch$float32)$to(device)
  Pred = model(XT)$data$cpu()$numpy()
  torch$cuda$empty_cache()
  return(list(obs = Y[(split_border+1L):nrow(X),,drop=FALSE], pred = Pred))
}


# Discrete Ricker model implementation in R

# Function to simulate the discrete Ricker model
ricker <- function(K, N0, num_generations, b0_m = 0.01, b1_m = 0.1, b0_g = 0.1, b1_g = 0.1, b0_r = 0.2, distP=0.05) {
  N <- numeric(num_generations)  # Create a vector to store population sizes
  
  N[1] <- N0  # Set initial population size
  
  for (t in 2:num_generations) {
    x1 = runif(1,0,1)
    m = 1-(b0_m + b1_m*x1)
    g = b0_g + b1_g*x1
    r = rpois(1, max(exp(b0_r * (1 - N[t - 1] / K))-1, 0))
    R = g*m
    disturbance = rbinom(1,1,1-distP)
    #print(paste0("x=",x1," r=",r," m=",m," g=",g," R=",R, " disturbance=",disturbance))
    N[t] <-  N[t - 1] * abs(rnorm(1, exp(R * (1 - N[t - 1] / K)), sd = 0.05))*disturbance + r
    
  }
  return(N)
}

inv.plogis <- binomial()$linkfun

# Function to simulate Beverton-Holt model
beverton_holt_R <- function(
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
  
  sample_vec <- seq(spinup+1, timesteps, sample_interval)
  # Return the N sizes and time steps
  return(out_df[sample_vec,])
}
