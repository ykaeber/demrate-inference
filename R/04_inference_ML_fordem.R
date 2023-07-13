library(sjSDM)
library(reticulate)
library(parallel)
library(Rcpp)
library(data.table)
library(tidyverse)
sourceCpp("library/fordem2.cpp")
source("library/functions.R")


cl = makeCluster(50L)
nodes = unlist(parallel::clusterEvalQ(cl, paste(Sys.info()[['nodename']], Sys.getpid(), sep='-')))
clusterEvalQ(cl, {library(sjSDM);library(data.table);library(tidyverse);source("library/functions.R")})
clusterEvalQ(cl, {library(Rcpp)})

# distP = 0.01, # range 0 to 0.2
# distInt = 0.01, # (fraction of patches affected) range 0 to 1
# regEnvEff = 1, # range 0 to 10
# regShadeEff = 1, # range 0 to 10
# mortEnvEff = 0.2, # range 0 to 10
# mortShadeEff = 1, # range 0 to 10
# env = 0.4, # range 0 to 1

parameter =
  data.frame(
    par = c("distP", "distInt","regEnvEff", "regShadeEff", "mortEnvEff","mortShadeEff", "env"),
    default = c(0.01, 0.01, 1, 1, 0.2, 1, 0.4),
    min = 0,
    max = c(0.2, 1, 10, 10, 10,10,1)
  )


data_simulate = expand.grid(factor(c(0,1)), 
                            factor(c(0, 1)), 
                            factor(c(0, 1)), 
                            factor(c(0, 1)), 
                            factor(c(0, 1)), 
                            factor(c(0,1)),
                            factor(c(0,1)))[-1,]
data_simulate = sapply(as.data.frame(data_simulate), function(r) as.integer(r) - 1L)
colnames(data_simulate) = c("distP", "distInt","regEnvEff", "regShadeEff", "mortEnvEff","mortShadeEff", "env")
data_simulate = as.data.frame(data_simulate)

# parallel::clusterExport(cl, varlist = ls(envir = .GlobalEnv))

results = list()
for(KK in 1:nrow(data_simulate)) {
  print(KK)
  # myself = paste(Sys.info()[['nodename']], Sys.getpid(), sep='-')
  # dist = cbind(nodes,0:3)
  # dev = as.integer(as.numeric(dist[which(dist[,1] %in% myself, arr.ind = TRUE), 2]))
  # 
  tmp = (data_simulate[KK, ])
  parallel::clusterExport(cl, varlist = list("tmp"), envir = environment())
  
  data =
    parLapply(cl, 1:2000, function(i) {
      source("library/functions.R")
      sourceCpp("library/fordem2.cpp")
      distP =          ifelse(tmp$distP,       runif(1, 0, 0.2), 0.01)
      distInt =        ifelse(tmp$distInt,     runif(1, 0, 1.0), 0.01)
      regEnvEff =      ifelse(tmp$regEnvEff,   runif(1, 0, 10), 1.0)
      regShadeEff =    ifelse(tmp$regShadeEff, runif(1, 0, 10), 1.0)
      mortEnvEff =     ifelse(tmp$mortEnvEff,  runif(1, 0, 10), 0.2) # 0.2
      mortShadeEff =   ifelse(tmp$mortShadeEff,runif(1, 0, 10), 1.0) 
      env =            ifelse(tmp$env,         runif(1, 0, 1.0), 0.4)
      print(distP)
      parsModel <- list(
        timesteps = 500,
        sampleSteps = 1,
        outVars = c("spID", "lai", "nTrs", "dbh", "ba"),
        actualSpecies = c(0,2,13,1,17,5, 27),
        baseReg = 50/0.1, # range: 1/0.1 to 5000/0.1
        baseRegP = 0.1, # (optional for inference) range 0.01 to 1
        regEnvEff = regEnvEff, # range 0 to 10
        regShadeEff = regShadeEff, # range 0 to 10
        mortEnvEff = mortEnvEff, # range 0 to 10
        mortShadeEff = mortShadeEff, # range 0 to 10
        bgMort = 2.3,
        distP = distP, # range 0 to 0.2
        distInt = distInt, # (fraction of patches affected) range 0 to 1
        env = env, # range 0 to 1
        nicheWidth = 4,
        patchesN = 100,
        areaHectar = 0.1,
        heightClassesN = 50,
        initPop = NULL,
        speciesPars = selected_species_pars
      )
      outMat1 = runModel(pars = parsModel, speciesPars = selected_species_pars)
      data_tmp =  transform_output(outMat1, parsModel = parsModel)
      return(list(pars = c(distP,
               distInt,
               regEnvEff,
               regShadeEff,
               mortEnvEff,
               mortShadeEff,
               env
      ), N = data_tmp ))
    })
  N = (abind::abind(lapply(data, function(d) d$N), along = 0L))
  N = (N-max(N))/max(N)
  pars = do.call(rbind, lapply(data, function(d) d$pars))
  samples = pars[, c(which(unlist(tmp) > 0.5, arr.ind = TRUE)),drop=FALSE]
  samples = data.frame(samples)
  colnames(samples) = names(which(unlist(tmp) == 1, arr.ind = TRUE))
  
  #if(sum(unlist(tmp)) != 6) samples = samples[, -c(which(unlist(tmp) != 1, arr.ind = TRUE)),drop=FALSE]
  result = train_function_demFor(N, samples, n_response = ncol(samples), epochs = 1000L, device = as.integer(1L), split =0.7, ntimes = 500L)
  
  results[[KK]]  = list(results = result, data = samples)
  names(results)[KK] = paste0(names(which(unlist(tmp) == 1, arr.ind = TRUE)), collapse = "+")
  saveRDS(results, file = "results/results_2000_forDem.RDS")
}
saveRDS(results, file = "results/results_2000_forDem_tot.RDS")
parallel::stopCluster(cl)
# 
# cl = makeCluster(12L)
# nodes = unlist(parallel::clusterEvalQ(cl, paste(Sys.info()[['nodename']], Sys.getpid(), sep='-')))
# clusterEvalQ(cl, {library(sjSDM);source("library/functions.R")})
# clusterEvalQ(cl, {library(Rcpp);sourceCpp("library/beverton-hold.cpp")})
# 
# parameter = 
#   data.frame(
#     par = c("distP", "b0_e","b1_recr", "b2_recr", "b1_g","b2_g"),
#     default = c(0.01, 0.5, 20, -0.2, 3, -0.1),
#     min = c(0, 0, 0, -5, 0, -5),
#     max = c(0.2, 1, 50, 0, 10,0)
#   )
# 
# 
# data_simulate = expand.grid(factor(c(0,1)), 
#                             factor(c(0, 1)), 
#                             factor(c(0, 1)), 
#                             factor(c(0, 1)), 
#                             factor(c(0, 1)), 
#                             factor(c(0,1)))[-1,]
# data_simulate = sapply(as.data.frame(data_simulate), function(r) as.integer(r) - 1L)
# colnames(data_simulate) = c("distP",
#                             "b0_e",
#                             "b1_recr",
#                             "b2_recr",
#                             "b1_g",
#                             "b2_g")
# data_simulate = as.data.frame(data_simulate)
# 
# parallel::clusterExport(cl, varlist = ls(envir = .GlobalEnv))
# 
# gc()
# torch$cuda$empty_cache()
# 
# results = parSapply(cl, 1:nrow(data_simulate), function(KK) {
#   sourceCpp("library/beverton-hold.cpp")
#   
#   myself = paste(Sys.info()[['nodename']], Sys.getpid(), sep='-')
#   dist = cbind(nodes,0:3)
#   dev = as.integer(as.numeric(dist[which(dist[,1] %in% myself, arr.ind = TRUE), 2]))
#   
#   source("library/functions.R")
#   
#   tmp = (data_simulate[KK, ])
#   
#   data = 
#     lapply(1:10000, function(i) {
#       
#       distP = ifelse(tmp$distP,     runif(1, 0, 0.2), 0.01)
#       b0_e = ifelse(tmp$b0_e,       runif(1, 0, 1.0), 0.5)
#       b1_recr = ifelse(tmp$b1_recr, runif(1, 0, 50), 20)
#       b2_recr = ifelse(tmp$b2_recr, runif(1, -5, 0.0), -0.2)
#       b1_g = ifelse(tmp$b1_g,       runif(1, 0, 1.0), 0.3)
#       b2_g = ifelse(tmp$b2_g,       runif(1, -5, 0), -0.1)
#       
#       data_tmp =  beverton_holt(N0 = 10, 
#                                 timesteps = 510, 
#                                 spinup = 10L, 
#                                 b0_e = b0_e, 
#                                 b1_e = 0.0,
#                                 b0_k = 100, 
#                                 b1_k = 4.9, 
#                                 b1_g = b1_g, 
#                                 b2_g = b2_g,
#                                 b0_r = 30, 
#                                 b1_recr = b1_recr, 
#                                 b2_recr = b2_recr, 
#                                 distP = distP, 
#                                 Nrep = 1, 
#                                 opt_x1 = 250)
#       N = data_tmp$N
#       E = data_tmp$x1
#       N[is.na(N)] = 0
#       return(c(distP,
#                b0_e,
#                b1_recr,
#                b2_recr,
#                b1_g,
#                b2_g,
#                N,
#                E
#       ))
#     })
#   data = (abind::abind(data, along = 0L))
#   samples = data[, c(which(unlist(tmp) > 0.5, arr.ind = TRUE)),drop=FALSE]
#   samples = data.frame(samples)
#   colnames(samples) = names(which(unlist(tmp) == 1, arr.ind = TRUE))
#   
#   if(sum(unlist(tmp)) != 6) data = data[, -c(which(unlist(tmp) != 1, arr.ind = TRUE))]
#   results = train_function(data, n_response = length(which(unlist(tmp) == 1, arr.ind = TRUE)), epochs = 1000L, device = as.integer(dev), split =0.7, ntimes = 500L)
#   result = list(list(results = results, data = samples))
#   names(result) = paste0(names(which(unlist(tmp) == 1, arr.ind = TRUE)), collapse = "+")
#   return(result)
# })
# saveRDS(results, file = "results/results_10000.RDS")
