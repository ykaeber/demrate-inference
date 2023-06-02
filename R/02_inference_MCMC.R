library(parallel)
library(Rcpp)
library(BayesianTools)
sourceCpp("library/beverton-hold.cpp")

data =
  lapply(1:1000, function(i) {
    
    distP = runif(1, 0, 0.2)
    b0_e = runif(1, 0, 1.0) 
    b1_recr = runif(1, 0, 50)
    b2_recr = runif(1, -5, 0.0)
    b1_g = runif(1, 0, 1.0)
    b2_g = runif(1, -5, 0)
    
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
             N
    ))
  })
data = (abind::abind(data, along = 0L))


cl = makeCluster(10L)
parallel::clusterExport(cl, varlist = ls(envir = .GlobalEnv))
clusterEvalQ(cl, {library(BayesianTools);Rcpp::sourceCpp("library/beverton-hold.cpp")})


results = parLapply(cl, 1:nrow(data), function(KK) {
  
  obs = data[KK, -(1:6)]
  true_pars = data[KK, 1:6]
  
  logLik = function(par) {
    ll = 
      replicate(50, {
        data_tmp =  beverton_holt(N0 = 10, 
                                  timesteps = 510, 
                                  spinup = 10L, 
                                  b0_e = par[2], 
                                  b1_e = 0.0,
                                  b0_k = 100, 
                                  b1_k = 4.9, 
                                  b1_g = par[5], 
                                  b2_g = par[6],
                                  b0_r = 30, 
                                  b1_recr = par[3], 
                                  b2_recr = par[4], 
                                  distP = par[1], 
                                  Nrep = 1, 
                                  opt_x1 = 250)$N
        sum((data_tmp-obs)**2)/500
      })
    return(-mean(ll))
  }
  
  bs = createBayesianSetup(logLik, prior = createUniformPrior(c(0,0,0,-5,0,-5), c(0.2, 1.0, 50, 0.0, 1.0, 0)))
  res = runMCMC(bs, settings = list(iterations = 10000))
  inferred = MAP(res, start = 2000)$par
  
  return(cbind(inferred, true_pars))
})
saveRDS(results, file = "results/MCMC_1000.RDS")
