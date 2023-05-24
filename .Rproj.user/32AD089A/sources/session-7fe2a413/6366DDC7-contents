library(sjSDM)
library(reticulate)
library(parallel)
source("R/functions.R")


cl = makeCluster(8L)
nodes = unlist(parallel::clusterEvalQ(cl, paste(Sys.info()[['nodename']], Sys.getpid(), sep='-')))

data_simulate = data.frame(
  Growth = c(1L, 0L, 0L, 0L, 1L, 0L, 0L, rep(1L, 5L)),
  Mortality = c(0L, 1L, rep(0L, 3L), 1L, 0L, 1L, 1L, 0L, 1L, 1L),
  Recruitment = c(0L, 0L, 1L, rep(0L, 3L), 1L, 0L, 0L, rep(1L, 3L)),
  Disturbance = c(rep(0L, 3L), rep(1L, 4L), 0L, 1L, 0L, 0L, 1L)
)

clusterEvalQ(cl, {library(sjSDM);source("R/functions.R")})
parallel::clusterExport(cl, varlist = ls(envir = .GlobalEnv))

results = parSapply(cl, 1:12, function(KK) {

  myself = paste(Sys.info()[['nodename']], Sys.getpid(), sep='-')
  dist = cbind(nodes,0:3)
  dev = as.integer(as.numeric(dist[which(dist[,1] %in% myself, arr.ind = TRUE), 2]))
  
  source("R/functions.R")
  
  tmp = (data_simulate[KK, ])
  
  data = 
    lapply(1:1000, function(i) {
      
      b0_g = ifelse(tmp$Growth,      runif(1, 0.01, 1.0), 0.3) # Growth
      b0_m = ifelse(tmp$Mortality,   runif(1, 0.01, 1.0), 0.8) # Mort
      b0_r = ifelse(tmp$Recruitment, runif(1, 0.05, 1.1), 0.1) # Recruit
      disP = ifelse(tmp$Disturbance, runif(1, 0.001, 0.2), 0.01) # Disturbance

      return(c(b0_g,
               b0_m,
               b0_r,
               disP, 
               ricker(K = 100, N0 = 10L, num_generations = 500L, distP = disP,
                      b0_m = b0_m, b0_r = b0_r, b0_g = b0_g, 
                      b1_m = 0, b1_g = 0)))
    })
  data = (abind::abind(data, along = 0L))
  samples = data[, c(which(unlist(tmp) > 0.5, arr.ind = TRUE)),drop=FALSE]
  samples = data.frame(samples)
  colnames(samples) = names(which(unlist(tmp) == 1, arr.ind = TRUE))
  
  if(sum(unlist(tmp)) != 4) data = data[, -c(which(unlist(tmp) != 1, arr.ind = TRUE))]
  results = train_function(data, n_response = length(which(unlist(tmp) == 1, arr.ind = TRUE)), epochs = 1000L, device = as.integer(dev), split =0.7)
  result = list(list(results = results, data = samples))
  names(result) = paste0(names(which(unlist(tmp) == 1, arr.ind = TRUE)), collapse = "+")
  return(result)
})
saveRDS(results, file = "results/results_1000.RDS")

torch$cuda$empty_cache()

results = parSapply(cl, 1:12, function(KK) {
  
  myself = paste(Sys.info()[['nodename']], Sys.getpid(), sep='-')
  dist = cbind(nodes,0:3)
  dev = as.integer(as.numeric(dist[which(dist[,1] %in% myself, arr.ind = TRUE), 2]))
  
  source("R/functions.R")
  
  tmp = (data_simulate[KK, ])
  
  data = 
    lapply(1:8000, function(i) {
      
      b0_g = ifelse(tmp$Growth,      runif(1, 0.01, 1.0), 0.3) # Growth
      b0_m = ifelse(tmp$Mortality,   runif(1, 0.01, 1.0), 0.8) # Mort
      b0_r = ifelse(tmp$Recruitment, runif(1, 0.05, 1.1), 0.1) # Recruit
      disP = ifelse(tmp$Disturbance, runif(1, 0.001, 0.2), 0.01) # Disturbance
      
      return(c(b0_g,
               b0_m,
               b0_r,
               disP, 
               ricker(K = 100, N0 = 10L, num_generations = 500L, distP = disP,
                      b0_m = b0_m, b0_r = b0_r, b0_g = b0_g, 
                      b1_m = 0, b1_g = 0)))
    })
  data = (abind::abind(data, along = 0L))
  samples = data[, c(which(unlist(tmp) > 0.5, arr.ind = TRUE)),drop=FALSE]
  samples = data.frame(samples)
  colnames(samples) = names(which(unlist(tmp) == 1, arr.ind = TRUE))
  
  if(sum(unlist(tmp)) != 4) data = data[, -c(which(unlist(tmp) != 1, arr.ind = TRUE))]
  results = train_function(data, n_response = length(which(unlist(tmp) == 1, arr.ind = TRUE)), epochs = 1000L, device = as.integer(dev), split =0.7)
  result = list(list(results = results, data = samples))
  names(result) = paste0(names(which(unlist(tmp) == 1, arr.ind = TRUE)), collapse = "+")
  return(result)
})
saveRDS(results, file = "results/results_8000.RDS")
