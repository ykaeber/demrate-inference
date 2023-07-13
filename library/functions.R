library(sjSDM)
library(reticulate)
torch = sjSDM:::pkg.env$fa$torch

RNN_model = 
  reticulate::PyClass("model", 
                      inherit = torch$nn$Module,
                      defs = list(
                        `__init__` = function(self, n_response = 1L, input_size = 2L, n_times = 200L, device) {
                          super()$`__init__`()
                          self$RNN = torch$nn$GRU(input_size = input_size,  hidden_size = 32L, num_layers = 2L, batch_first = TRUE)$to(device)
                          self$linear2 = torch$nn$Linear(32L*as.integer(n_times), 20L)$to(device)
                          self$linear1 = torch$nn$Linear(20L, n_response)$to(device)
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

train_function = function(data, split = 0.8, epochs = 1000L, n_response = 1L, device = 0L, ntimes = 500L) {
  device = paste0('cuda:', device)
  # Scaling
  X = (data[,-(1:n_response)])
  N_t = X[,1:ntimes]
  X_t = X[,(ntimes+1):(ntimes*2)]
  N_t = (N_t - mean(as.vector(N_t)))/sd(as.vector(N_t))
  X = cbind(N_t, X_t)
  Y = data[, 1:n_response, drop=FALSE]
  split_border = round(nrow(data)*split)
  XX = array(X[1:split_border,], dim = c(split_border, ntimes,2L))
  torch$cuda$set_device(device)
  model = RNN_model(n_response = as.integer(n_response), n_times = ntimes, device = device )$cuda()
  optimizer = torch$optim$Adamax(model$parameters(), lr = 0.001)
  lambda1 = torch$tensor(0.0008)$to(device)
  lambda2 = torch$tensor(0.001)$to(device)
  DT = 
    torch$utils$data$TensorDataset(torch$tensor(XX, dtype = torch$float32, device = device), 
                                   torch$tensor(Y[1:split_border,,drop=FALSE], dtype = torch$float32, device = device))
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
  XX = array(X[(split_border+1L):nrow(X),], dim = c(length((split_border+1L):nrow(X)), ntimes,2L))
  XT = torch$tensor(XX, 
                    dtype = torch$float32)$to(device)
  Pred = model(XT)$data$cpu()$numpy()
  torch$cuda$empty_cache()
  return(list(obs = Y[(split_border+1L):nrow(X),,drop=FALSE], pred = Pred))
}


train_function_demFor = function(N, samples, split = 0.8, epochs = 1000L, n_response = 1L, device = 0L, ntimes = 500L) {
  device = paste0('cuda:', device)
  # Scaling
  X = N
  Y = samples
  split_border = round(nrow(N)*split)
  XX = X[1:split_border,,]
  torch$cuda$set_device(device)
  model = RNN_model(n_response = as.integer(n_response),input_size = 7L, n_times = ntimes, device = device )$cuda()
  optimizer = torch$optim$Adamax(model$parameters(), lr = 0.001)
  lambda1 = torch$tensor(0.0008)$to(device)
  lambda2 = torch$tensor(0.001)$to(device)
  DT = 
    torch$utils$data$TensorDataset(torch$tensor(XX, dtype = torch$float32, device = device), 
                                   torch$tensor(as.matrix(Y[1:split_border,,drop=FALSE]), dtype = torch$float32, device = device))
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
  XX = X[(split_border+1L):nrow(X),,]
  XT = torch$tensor(XX, 
                    dtype = torch$float32)$to(device)
  Pred = model(XT)$data$cpu()$numpy()
  torch$cuda$empty_cache()
  return(list(obs = Y[(split_border+1L):nrow(X),,drop=FALSE], pred = Pred))
}


selected_species_pars <- 
  data.frame(
    stringsAsFactors = FALSE,
    spID = c(0L,1L,2L,3L,4L,5L,6L,7L,
             8L,9L,10L,11L,12L,13L,14L,15L,16L,17L,18L,19L,
             20L,21L,22L,23L,24L,25L,26L,27L,28L,29L),
    species = c("Abies alba","Larix decidua",
                "Picea abies","Pinus cembra","Pinus montana",
                "Pinus sylvestris","Taxus baccata","Acer campestre",
                "Acer platanoides","Acer pseudoplatanus","Alnus glutinosa",
                "Alnus incana","Alnus viridis","Betula pendula",
                "Carpinus betulus","Castanea sativa","Corylus avellana",
                "Fagus sylvatica","Fraxinus excelsior","Populus nigra",
                "Populus tremula","Quercus petraea","Quercus pubescens",
                "Quercus robur","Salix alba","Sorbus aria","Sorbus aucuparia",
                "Tilia cordata","Tilia platyphyllos","Ulmus glabra"),
    kHMax = c(60L,54L,63L,26L,23L,48L,
              22L,25L,35L,40L,40L,25L,6L,30L,35L,35L,15L,52L,
              42L,40L,42L,50L,25L,52L,35L,23L,27L,40L,40L,
              43L),
    kDMax = c(200L,250L,200L,200L,50L,
              150L,350L,150L,200L,250L,150L,100L,25L,100L,150L,
              350L,50L,250L,250L,250L,150L,350L,150L,350L,250L,
              100L,100L,350L,350L,300L),
    kG = c(296L,400L,342L,198L,239L,
           393L,175L,210L,360L,338L,380L,218L,476L,448L,360L,
           375L,245L,307L,363L,394L,390L,378L,226L,376L,
           403L,230L,205L,365L,365L,361L),
    kDDMin = c(641L,323L,385L,323L,436L,
               610L,1011L,1062L,1042L,898L,898L,610L,272L,610L,
               898L,1237L,898L,723L,980L,662L,610L,785L,1011L,
               1042L,1062L,898L,498L,1339L,1339L,1062L),
    kDrTol = c(0.23,0.25,0.15,
               0.3,0.37,0.37,0.23,0.33,0.25,0.25,0.08,0.08,
               0.16,0.16,0.25,0.33,0.33,0.25,0.16,0.08,
               0.25,0.33,0.33,0.25,0.08,0.33,0.33,0.33,0.25,
               0.25),
    kLy = c(0.03,0.5,0.05,0.075,0.5,
            0.4,0.03,0.2,0.075,0.05,0.2,0.2,0.3,0.5,0.075,0.1,
            0.2,0.03,0.075,0.3,0.3,0.2,0.3,0.2,0.3,0.3,
            0.1,0.075,0.075,0.075),
    kLa = c(1L,9L,5L,5L,9L,8L,1L,6L,
            4L,3L,6L,7L,7L,9L,4L,5L,6L,1L,5L,7L,7L,7L,
            8L,7L,7L,8L,6L,5L,5L,5L),
    kA1 = c(0.23,0.1,0.23,0.23,0.23,
            0.17,0.23,0.1,0.06,0.06,0.1,0.1,0.1,0.08,0.06,
            0.06,0.06,0.06,0.1,0.1,0.1,0.06,0.06,0.06,0.08,0.1,
            0.08,0.06,0.06,0.06),
    kA2 = c(1.56,1.43,1.56,1.56,1.56,
            1.4,1.56,1.43,1.7,1.7,1.43,1.43,1.43,1.43,1.7,
            1.7,1.7,1.7,1.43,1.43,1.43,1.7,1.7,1.7,1.43,1.43,
            1.43,1.7,1.7,1.7),
    kC1 = c(0.45,0.35,0.45,0.45,0.45,
            0.45,0.45,0.35,0.35,0.35,0.35,0.35,0.35,0.35,0.35,
            0.35,0.35,0.35,0.35,0.35,0.35,0.35,0.35,0.35,
            0.35,0.35,0.35,0.35,0.35,0.35),
    kC2 = c(6L,12L,6L,6L,6L,6L,6L,
            12L,12L,12L,12L,12L,12L,12L,12L,12L,12L,12L,12L,
            12L,12L,12L,12L,12L,12L,12L,12L,12L,12L,12L)
  )

selected_species_pars <- data.table(selected_species_pars)
selected_species_pars[, env := scale(kDDMin), ]
selected_species_pars[, env := env+abs(min(env)), ]
selected_species_pars[, env := env/max(env), ]



transform_output = function(outMat1, parsModel) {
  outDT <- data.table(outMat1)
  names(outDT) <- c(parsModel$outVars,"t","p")
  out1 <- outDT
  out1 <- merge(out1, selected_species_pars[,c("spID","species")], by = "spID")
  p_dat <- out1[, .(
    lai = mean(lai),
    ba = mean(ba),
    dbh = mean(dbh),
    nTrs = mean(nTrs)
  ), by = .(t, species)]
  p_dat <- melt(p_dat, id.vars = c("t", "species"))
  p_dat = p_dat[p_dat$variable == "ba",]
  out = 
    p_dat %>% 
    pivot_wider(names_from = "species", values_from = "value") %>% 
    select(-t, -variable) %>% as.matrix()
  return(out)  
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
