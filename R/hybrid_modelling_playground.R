library(sjSDM)
torch = sjSDM:::pkg.env$torch
dt = torch$float32
as_ft = function(x) torch$tensor(x, dtype = dt)
as_ftg = function(x) torch$tensor(x, dtype = dt, requires_grad = TRUE)


beverton_holt_torch <- function(
    N0, timesteps, spinup = 0, sample_interval = 1,
    b0_e = as_ftg(0.0), 
    b1_e =as_ftg( 0.0),
    b0_k =as_ftg( 100),
    b1_k =as_ftg( 0.3),
    b1_g =as_ftg( 0.0),
    b2_g =as_ftg( 0),
    b0_r =as_ftg( 1),
    b1_recr =as_ftg( 0.1),
    b2_recr = as_ftg(0.001),
    opt_x1 =as_ftg( 1),
    distP =as_ftg( 0.01),
    noise_env = 0.002,
    Nrep = 1
) {
    out_df = list()
    N = list()
    Nrep = as.integer(Nrep)
    N[[1]] = as_ft(matrix(N0, Nrep, 1L))
    recr = as_ft(matrix(1.0, Nrep, 1L))
    one = as_ft(1.0)
    plogis = function(x) one$div(one+x$neg()$exp())
    
    for (t in 2:timesteps) {
      x1 = b0_e$add(torch$distributions$Normal(0., noise_env)$sample(list(Nrep, 1L)))
      x1 = torch$clamp(x1, 0.0, 1.0)
      K = b0_k$add(b1_k$multiply( x1))
      g = plogis(b1_g$multiply(x1)$add(b2_g$multiply(N[[t-1]])))$multiply(as_ft(2))
      (disturbance = torch$distributions$Binomial(1L, probs = one$sub(distP))$sample(list(Nrep, 1L)))
      N[[t]] <- g$multiply(N[[t-1]])$div( one$add( N[[t-1]]$div(K)  ))$multiply(disturbance) 
      recr = torch$clamp(recr, 0., 1.)
      K_r = as_ft(10.)
      g_recr = plogis( b0_r$add( b1_recr$multiply(x1))$add( b2_recr$multiply(N[[t]])) )$multiply(2.0)
      recr_mean = g_recr$multiply(recr)$div( one$add(recr$div(K_r)))
      recr = torch$exp(recr_mean$add(torch$distributions$Normal( 0., 0.02)$sample(list(Nrep, 1L)))) #rpois(1, recr_mean) 
      N[[t]] = N[[t]]$add( recr)
      out_df[[t]] <- list(rep=rep_i, t=t, N=N[[t]], x1, K, g, recr_mean, recr,disturbance, actualR = N[[t-1]]$div(N[[t]]))
      
    }
  return(out_df)
}

beverton_holt_DNN <- function(
    N0, timesteps, spinup = 0, sample_interval = 1,
    b0_e = as_ftg(0.0), 
    b1_e =as_ftg( 0.0),
    b0_k =as_ftg( 100),
    b1_k =as_ftg( 0.3),
    b1_g =as_ftg( 0.0),
    b2_g = NULL,
    b0_r =as_ftg( 1),
    b1_recr =as_ftg( 0.1),
    b2_recr = as_ftg(0.001),
    opt_x1 =as_ftg( 1),
    distP =as_ftg( 0.01),
    Nrep = 1
) {
  out_df = list()
  N = list()
  Nrep = as.integer(Nrep)
  N[[1]] = as_ft(matrix(N0, Nrep, 1L))
  recr = as_ft(matrix(1.0, Nrep, 1L))
  one = as_ft(1.0)
  plogis = function(x) one$div(one+x$neg()$exp())
  
  for (t in 2:timesteps) {
    x1 = b0_e$add(torch$distributions$Normal(0., 0.002)$sample(list(Nrep, 1L)))
    x1 = torch$clamp(x1, 0.0, 1.0)
    K = b0_k$add(b1_k$multiply( x1))
    g = plogis(b1_g$multiply(x1)$add(b2_g(N[[t-1]])))$multiply(as_ft(2))
    (disturbance = torch$distributions$Binomial(1L, probs = one$sub(distP))$sample(list(Nrep, 1L)))
    N[[t]] <- g$multiply(N[[t-1]])$div( one$add( N[[t-1]]$div(K)  ))$multiply(disturbance) 
    recr = torch$clamp(recr, 0., 1.)
    K_r = as_ft(10.)
    g_recr = plogis( b0_r$add( b1_recr$multiply(x1))$add( b2_recr$multiply(N[[t]])) )$multiply(2.0)
    recr_mean = g_recr$multiply(recr)$div( one$add(recr$div(K_r)))
    recr = torch$exp(recr_mean$add(torch$distributions$Normal( 0., 0.02)$sample(list(Nrep, 1L)))) #rpois(1, recr_mean) 
    N[[t]] = N[[t]]$add( recr)
    out_df[[t]] <- list(rep=rep_i, t=t, N=N[[t]], x1, K, g, recr_mean, recr,disturbance, actualR = N[[t-1]]$div(N[[t]]))
    
  }
  return(out_df)
}

sim_df <- beverton_holt_torch(N0 = 10, 
                              timesteps = 200, 
                              b0_e =as_ftg( 0.3),
                              b1_e =as_ftg( 0.0),
                              b0_k =as_ftg( 100),
                              b1_k =as_ftg( 4.9),
                              b1_g =as_ftg( 0.3),
                              b0_r =as_ftg( 30),
                              b1_recr =as_ftg(20),
                              b2_recr =as_ftg(-0.2),
                              b2_g =as_ftg(-0.2), # -0.1
                              distP =as_ftg( 0.05),
                              Nrep = 1, 
                              opt_x1 = 250
)
plot(sapply(2:length(sim_df), function(x) reticulate::py_to_r(sim_df[[x]]$N$data$cpu()$numpy())), type = "l")

N0 = 10
timesteps = 200
b0_e =as_ftg( 1)
b1_e =as_ftg( 0.0)
b0_k =as_ftg( 100)
b1_k =as_ftg( 4.9)
b1_g =as_ftg( 1.4)
b0_r =as_ftg( 30)
b1_recr =as_ftg( 8.1)
b2_recr =as_ftg(3.1)
b2_g =as_ftg( 0)
distP =as_ftg( 0.02)
Nrep = 1
opt_x1 = 25



Y = as_ft(sapply(2:length(sim_df), function(x) reticulate::py_to_r(sim_df[[x]]$N$data$cpu()$numpy())))


par =as_ftg( -1)
optim = torch$optim$Adamax(params = list(par), lr = 0.1 )

for(i in 1:100) {
  optim$zero_grad()
  sim_df <- beverton_holt_torch(N0 = 10, 
                                timesteps = 200, 
                                b0_e =as_ftg( 0.3),
                                b1_e =as_ftg( 0.0),
                                b0_k =as_ftg( 100),
                                b1_k =as_ftg( 4.9),
                                b1_g =as_ftg( 0.3),
                                b0_r =as_ftg( 30),
                                b1_recr =as_ftg(20),
                                b2_recr =as_ftg(-0.2),
                                b2_g =par, # -0.1
                                distP =as_ftg( 0.05),
                                Nrep = 200, 
                                opt_x1 = 250
  )
  pred = torch$cat(sapply(2:length(sim_df), function(x) sim_df[[x]]$N), 1L)
  loss = (pred - Y$reshape(list(1L, -1L)))$pow(2.0)$mean()
  #loss = torch$nn$functional$mse_loss(Y, pred)
  #K = torch$autograd$grad(loss, b1_recr)
  loss$backward()
  optim$step()
  print(par)
}

model = torch$nn$Sequential(
  torch$nn$Linear(1L, 20L),
  torch$nn$ReLU(),
  torch$nn$Linear(20L, 1L)
)
optim = torch$optim$Adamax(params = model$parameters(), lr = 0.003 )

for(i in 1:100) {
  optim$zero_grad()
  sim_df <- beverton_holt_DNN(N0 = 10, 
                                timesteps = 200, 
                                b0_e =as_ftg( 0.3),
                                b1_e =as_ftg( 0.0),
                                b0_k =as_ftg( 100),
                                b1_k =as_ftg( 4.9),
                                b1_g =as_ftg( 0.3),
                                b0_r =as_ftg( 30),
                                b1_recr =as_ftg(20),
                                b2_recr =as_ftg(-0.2),
                                b2_g =function(N_density) model$forward(N_density), # -0.1
                                distP =as_ftg( 0.05),
                                Nrep = 200, 
                                opt_x1 = 250
  )
  pred = torch$cat(sapply(2:length(sim_df), function(x) sim_df[[x]]$N), 1L)
  loss = (pred - Y$reshape(list(1L, -1L)))$pow(2.0)$mean()
  #loss = torch$nn$functional$mse_loss(Y, pred)
  #K = torch$autograd$grad(loss, b1_recr)
  loss$backward()
  optim$step()
  print(loss)
}

pred = model$forward(as_ft(matrix(seq(1, 100, length.out = 100), 100, 1)))
plot(seq(1, 100, length.out = 100), reticulate::py_to_r(pred$data$cpu()$numpy()))
points(seq(1, 100, length.out = 100), seq(1, 100, length.out = 100)*(-0.2), col = "red")

# distP = ifelse(tmp$distP,     runif(1, 0, 0.2), 0.01)
# b0_e = ifelse(tmp$b0_e,       runif(1, 0, 1.0), 0.5)
# b1_recr = ifelse(tmp$b1_recr, runif(1, 0, 50), 20)
# b2_recr = ifelse(tmp$b2_recr, runif(1, -5, 0.0), -0.2)
# b1_g = ifelse(tmp$b1_g,       runif(1, 0, 1.0), 0.3)
# b2_g = ifelse(tmp$b2_g,       runif(1, -5, 0), -0.1)
# 

# rr = 
# sapply(seq(-1, -0.05, length.out = 25), function(XX) {
#   sim_df <- beverton_holt_torch(N0 = 10, 
#                                 timesteps = 200, 
#                                 b0_e =as_ftg( 0.3),
#                                 b1_e =as_ftg( 0.0),
#                                 b0_k =as_ftg( 100),
#                                 b1_k =as_ftg( 4.9),
#                                 b1_g =as_ftg( 0.3),
#                                 b0_r =as_ftg( 30),
#                                 b1_recr =as_ftg(20),
#                                 b2_recr =as_ftg(-0.2),
#                                 b2_g =as_ftg( XX), # -0.1
#                                 distP =as_ftg( 0.05),
#                                 Nrep = 1000, 
#                                 opt_x1 = 250
#   )
#   pred = torch$cat(sapply(2:length(sim_df), function(x) sim_df[[x]]$N), 1L)
#   loss = (pred - Y$reshape(list(1L, -1L)))$pow(2.0)$mean()
#   return(reticulate::py_to_r(loss$data$cpu()$numpy()))
# })
# plot(seq(-1, -0.05, length.out = 25), rr)

beverton_holt_DNN_growth <- function(
    N0, timesteps, spinup = 0, sample_interval = 1,
    b0_e = as_ftg(0.0), 
    b1_e =as_ftg( 0.0),
    b0_k =as_ftg( 100),
    b1_k =as_ftg( 0.3),
    b0_r =as_ftg( 1),
    b1_recr =as_ftg( 0.1),
    b2_recr = as_ftg(0.001),
    opt_x1 =as_ftg( 1),
    distP =as_ftg( 0.01),
    Nrep = 1,
    growth = NULL
) {
  out_df = list()
  N = list()
  Nrep = as.integer(Nrep)
  N[[1]] = as_ft(matrix(N0, Nrep, 1L))
  recr = as_ft(matrix(1.0, Nrep, 1L))
  one = as_ft(1.0)
  plogis = function(x) one$div(one+x$neg()$exp())
  
  for (t in 2:timesteps) {
    x1 = b0_e$add(torch$distributions$Normal(0., 0.002)$sample(list(Nrep, 1L)))
    x1 = torch$clamp(x1, 0.0, 1.0)
    K = b0_k$add(b1_k$multiply( x1))
    #g = plogis(b1_g$multiply(x1)$add(b2_g(N[[t-1]])))$multiply(as_ft(2))
    g = growth(x1, N[[t-1]] )$multiply(as_ft(2))
    (disturbance = torch$distributions$Binomial(1L, probs = one$sub(distP))$sample(list(Nrep, 1L)))
    N[[t]] <- g$multiply(N[[t-1]])$div( one$add( N[[t-1]]$div(K)  ))$multiply(disturbance) 
    recr = torch$clamp(recr, 0., 1.)
    K_r = as_ft(10.)
    g_recr = plogis( b0_r$add( b1_recr$multiply(x1))$add( b2_recr$multiply(N[[t]])) )$multiply(2.0)
    recr_mean = g_recr$multiply(recr)$div( one$add(recr$div(K_r)))
    recr = torch$exp(recr_mean$add(torch$distributions$Normal( 0., 0.02)$sample(list(Nrep, 1L)))) #rpois(1, recr_mean) 
    N[[t]] = N[[t]]$add( recr)
    out_df[[t]] <- list(rep=rep_i, t=t, N=N[[t]], x1, K, g, recr_mean, recr,disturbance, actualR = N[[t-1]]$div(N[[t]]))
    
  }
  return(out_df)
}

sim_df <- beverton_holt_torch(N0 = 10, 
                              timesteps = 200, 
                              b0_e =as_ftg( 0.3),
                              b1_e =as_ftg( 0.0),
                              b0_k =as_ftg( 100),
                              b1_k =as_ftg( 4.9),
                              b1_g =as_ftg( 0.4),
                              b0_r =as_ftg( 30),
                              b1_recr =as_ftg(20),
                              b2_recr =as_ftg(-0.2),
                              b2_g =as_ftg(-0.1), # -0.1
                              distP =as_ftg( 0.05),
                              Nrep = 1, 
                              noise_env = 0.3,
                              opt_x1 = 250
)
plot(sapply(2:length(sim_df), function(x) reticulate::py_to_r(sim_df[[x]]$N$data$cpu()$numpy())), type = "l")


Y = as_ft(sapply(2:length(sim_df), function(x) reticulate::py_to_r(sim_df[[x]]$N$data$cpu()$numpy())))

model = torch$nn$Sequential(
  torch$nn$Linear(2L, 20L),
  torch$nn$ReLU(),
  torch$nn$Linear(20L, 1L),
  torch$nn$Sigmoid()
)
optim = torch$optim$Adamax(params = model$parameters(), lr = 0.01 )

for(i in 1:100) {
  optim$zero_grad()
  sim_df <- beverton_holt_DNN_growth(N0 = 10, 
                              timesteps = 200, 
                              b0_e =as_ftg( 0.3),
                              b1_e =as_ftg( 0.0),
                              b0_k =as_ftg( 100),
                              b1_k =as_ftg( 4.9),
                              #b1_g =as_ftg( 0.3),
                              b0_r =as_ftg( 30),
                              b1_recr =as_ftg(20),
                              b2_recr =as_ftg(-0.2),
                              #b2_g =function(N_density) model$forward(N_density), # -0.1
                              growth = function(x1, N_density) model$forward(torch$cat(list(x1, N_density), 1L)),
                              distP =as_ftg( 0.05),
                              Nrep = 200, 
                              opt_x1 = 250
  )
  pred = torch$cat(sapply(2:length(sim_df), function(x) sim_df[[x]]$N), 1L)
  loss = (pred - Y$reshape(list(1L, -1L)))$pow(2.0)$mean()
  #loss = torch$nn$functional$mse_loss(Y, pred)
  #K = torch$autograd$grad(loss, b1_recr)
  loss$backward()
  optim$step()
  print(loss)
}
pred = model$forward(as_ft(cbind(0.5,matrix(seq(1, 100, length.out = 100), 100, 1))))
plot(seq(1, 100, length.out = 100), binomial()$linkfun(reticulate::py_to_r(pred$data$cpu()$numpy())))
points(seq(1, 100, length.out = 100), seq(1, 100, length.out = 100)*(-0.2), col = "red")


pred = model$forward(as_ft(cbind(matrix(0.5+rnorm(100, 0, 0.002), 100, 1), 60) ) )
plot(seq(1, 100, length.out = 100), binomial()$linkfun(reticulate::py_to_r(pred$data$cpu()$numpy())))
points(seq(1, 100, length.out = 100), seq(1, 100, length.out = 100)*(-0.2), col = "red")


marginalEffectsGeneric(model, 
                       data = cbind(matrix(0.5+rnorm(100, 0, 0.3), 100, 1), matrix(sample(5:20, 100, replace = TRUE), 100, 1)),
                       predict_func = function(model, newdata) binomial()$linkfun(reticulate::py_to_r(model$forward(as_ft(newdata))$data$cpu()$numpy()))
                        )



library(BayesianTools)
ll = function(par) {
  sim_df <- beverton_holt_torch(N0 = 10, 
                                timesteps = 200, 
                                b0_e =as_ftg( 0.3),
                                b1_e =as_ftg( 0.0),
                                b0_k =as_ftg( 100),
                                b1_k =as_ftg( 4.9),
                                b1_g =as_ft(par[1]),
                                b0_r =as_ftg( 30),
                                b1_recr =as_ftg(20),
                                b2_recr =as_ftg(-0.2),
                                b2_g =as_ft(par[2]), # -0.1
                                distP =as_ftg( 0.05),
                                Nrep = 100, 
                                noise_env = 0.3,
                                opt_x1 = 250
  )
  pred = torch$cat(sapply(2:length(sim_df), function(x) sim_df[[x]]$N), 1L)
  loss = (pred - Y$reshape(list(1L, -1L)))$pow(2.0)$mean()
  return(-reticulate::py_to_r(loss$data$cpu()$numpy()))
}
pr = createUniformPrior(c(-2, 0), c(0, 2))
bs = createBayesianSetup(ll, pr)
run = runMCMC(bs, settings = list(iterations = 1000L))
plot(run)
