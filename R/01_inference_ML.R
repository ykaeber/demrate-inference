library(sjSDM)
library(reticulate)
library(parallel)
library(Rcpp)
sourceCpp("library/beverton-hold.cpp")
source("library/functions.R")


cl = makeCluster(12L)
nodes = unlist(parallel::clusterEvalQ(cl, paste(Sys.info()[['nodename']], Sys.getpid(), sep='-')))
clusterEvalQ(cl, {library(sjSDM);source("library/functions.R")})
clusterEvalQ(cl, {library(Rcpp);sourceCpp("library/beverton-hold.cpp")})

parameter = 
  data.frame(
    par = c("distP", "b0_e","b1_recr", "b2_recr", "b1_g","b2_g"),
    default = c(0.01, 0.5, 20, -0.2, 3, -0.1),
    min = c(0, 0, 0, -5, 0, -5),
    max = c(0.2, 1, 50, 0, 10,0)
  )


data_simulate = expand.grid(factor(c(0,1)), 
                            factor(c(0, 1)), 
                            factor(c(0, 1)), 
                            factor(c(0, 1)), 
                            factor(c(0, 1)), 
                            factor(c(0,1)))[-1,]
data_simulate = sapply(as.data.frame(data_simulate), function(r) as.integer(r) - 1L)
colnames(data_simulate) = c("distP",
                            "b0_e",
                            "b1_recr",
                            "b2_recr",
                            "b1_g",
                            "b2_g")
data_simulate = as.data.frame(data_simulate)

parallel::clusterExport(cl, varlist = ls(envir = .GlobalEnv))

results = parSapply(cl, 1:nrow(data_simulate), function(KK) {
  sourceCpp("R/beverton-hold.cpp")
  
  myself = paste(Sys.info()[['nodename']], Sys.getpid(), sep='-')
  dist = cbind(nodes,0:3)
  dev = as.integer(as.numeric(dist[which(dist[,1] %in% myself, arr.ind = TRUE), 2]))
  
  source("R/functions.R")
  
  tmp = (data_simulate[KK, ])
  
  data =
    lapply(1:2000, function(i) {
      
      distP = ifelse(tmp$distP,     runif(1, 0, 0.2), 0.01)
      b0_e = ifelse(tmp$b0_e,       runif(1, 0, 1.0), 0.5)
      b1_recr = ifelse(tmp$b1_recr, runif(1, 0, 50), 20)
      b2_recr = ifelse(tmp$b2_recr, runif(1, -5, 0.0), -0.2)
      b1_g = ifelse(tmp$b1_g,       runif(1, 0, 1.0), 3)
      b2_g = ifelse(tmp$b2_g,       runif(1, -5, 0), -0.1)
      
      data_tmp =  beverton_holt(N0 = 10,
                                timesteps = 510,
                                spinup = 10L,
                                b0_e = b0_e,
                                b1_e = 0.0,
                                b0_k = 100,
                                b1_k = 4.9,
                                b1_g = b1_g,
                                b2_g = b2_g,
                                b0_r = 30,
                                b1_recr = b1_recr,
                                b2_recr = b2_recr,
                                distP = distP,
                                Nrep = 1,
                                opt_x1 = 250)
      N = data_tmp$N
      E = data_tmp$x1
      N[is.na(N)] = 0
      return(c(distP,
               b0_e,
               b1_recr,
               b2_recr,
               b1_g,
               b2_g,
               N,
               E
      ))
    })
  data = (abind::abind(data, along = 0L))
  samples = data[, c(which(unlist(tmp) > 0.5, arr.ind = TRUE)),drop=FALSE]
  samples = data.frame(samples)
  colnames(samples) = names(which(unlist(tmp) == 1, arr.ind = TRUE))
  
  if(sum(unlist(tmp)) != 6) data = data[, -c(which(unlist(tmp) != 1, arr.ind = TRUE))]
  results = train_function(data, n_response = length(which(unlist(tmp) == 1, arr.ind = TRUE)), epochs = 1000L, device = as.integer(dev), split =0.7, ntimes = 500L)
  result = list(list(results = results, data = samples))
  names(result) = paste0(names(which(unlist(tmp) == 1, arr.ind = TRUE)), collapse = "+")
  return(result)
})
saveRDS(results, file = "results/results_2000.RDS")

torch$cuda$empty_cache()

results = parSapply(cl, 1:nrow(data_simulate), function(KK) {
  sourceCpp("R/beverton-hold.cpp")
  
  myself = paste(Sys.info()[['nodename']], Sys.getpid(), sep='-')
  dist = cbind(nodes,0:3)
  dev = as.integer(as.numeric(dist[which(dist[,1] %in% myself, arr.ind = TRUE), 2]))
  
  source("R/functions.R")
  
  tmp = (data_simulate[KK, ])
  
  data = 
    lapply(1:10000, function(i) {
      
      distP = ifelse(tmp$distP,     runif(1, 0, 0.2), 0.01)
      b0_e = ifelse(tmp$b0_e,       runif(1, 0, 1.0), 0.5)
      b1_recr = ifelse(tmp$b1_recr, runif(1, 0, 50), 20)
      b2_recr = ifelse(tmp$b2_recr, runif(1, -5, 0.0), -0.2)
      b1_g = ifelse(tmp$b1_g,       runif(1, 0, 1.0), 3)
      b2_g = ifelse(tmp$b2_g,       runif(1, -5, 0), -0.1)
      
      data_tmp =  beverton_holt(N0 = 10, 
                                timesteps = 510, 
                                spinup = 10L, 
                                b0_e = b0_e, 
                                b1_e = 0.0,
                                b0_k = 100, 
                                b1_k = 4.9, 
                                b1_g = b1_g, 
                                b2_g = b2_g,
                                b0_r = 30, 
                                b1_recr = b1_recr, 
                                b2_recr = b2_recr, 
                                distP = distP, 
                                Nrep = 1, 
                                opt_x1 = 250)
      N = data_tmp$N
      E = data_tmp$x1
      N[is.na(N)] = 0
      return(c(distP,
               b0_e,
               b1_recr,
               b2_recr,
               b1_g,
               b2_g,
               N,
               E
      ))
    })
  data = (abind::abind(data, along = 0L))
  samples = data[, c(which(unlist(tmp) > 0.5, arr.ind = TRUE)),drop=FALSE]
  samples = data.frame(samples)
  colnames(samples) = names(which(unlist(tmp) == 1, arr.ind = TRUE))
  
  if(sum(unlist(tmp)) != 6) data = data[, -c(which(unlist(tmp) != 1, arr.ind = TRUE))]
  results = train_function(data, n_response = length(which(unlist(tmp) == 1, arr.ind = TRUE)), epochs = 1000L, device = as.integer(dev), split =0.7, ntimes = 500L)
  result = list(list(results = results, data = samples))
  names(result) = paste0(names(which(unlist(tmp) == 1, arr.ind = TRUE)), collapse = "+")
  return(result)
})
saveRDS(results, file = "results/results_10000.RDS")
