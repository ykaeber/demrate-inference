library(sjSDM)
library(reticulate)
library(parallel)
library(Rcpp)
library(BayesianTools)
sourceCpp("R/beverton-hold.cpp")
source("R/functions.R")

res_b2_recr = 
  sapply(seq(-5, -0.01, length.out = 20), function(P) {
    
    obs =  beverton_holt(N0 = 10, 
                         timesteps = 210, 
                         spinup = 10L, 
                         b0_e = 0.5, 
                         b1_e = 0.0,
                         b0_k = 100, 
                         b1_k = 4.9, 
                         b1_g = 3, 
                         b2_g = -0.1,
                         b0_r = 30, 
                         b1_recr = 20, 
                         b2_recr = P, 
                         distP = 0.01, 
                         Nrep = 1, 
                         opt_x1 = 250)
    
    logLik = function(par) {
      ll = 
        replicate(10, {
          data_tmp =  beverton_holt(N0 = 10, 
                                    timesteps = 510, 
                                    spinup = 10L, 
                                    b0_e = 0.5, 
                                    b1_e = 0.0,
                                    b0_k = 100, 
                                    b1_k = 4.9, 
                                    b1_g = 3, 
                                    b2_g = -0.1,
                                    b0_r = 30, 
                                    b1_recr = 20, 
                                    b2_recr = par, 
                                    distP = 0.01, 
                                    Nrep = 1, 
                                    opt_x1 = 250)$N
          sum((data_tmp-obs$N)**2)/500
        })
      return(-mean(ll))
    }
    
    bs = createBayesianSetup(logLik, prior = createUniformPrior(-7, 0))
    res = runMCMC(bs, settings = list(iterations = 3000))
    return(MAP(res)$par)
  })
plot(y = res_b2_recr, seq(-5, -0.01, length.out = 20))



obs =  beverton_holt(N0 = 10, timesteps = 510, spinup = 10L, b0_e = 0.5, b1_e = 0.0,
                     b0_k = 100, 
                     b1_k = 4.9, 
                     b1_g = 3, 
                     b2_g = -0.1,
                     b0_r = 30, 
                     b1_recr = 20, 
                     b2_recr = -0.2, 
                     distP = 0.01, 
                     Nrep = 1, 
                     opt_x1 = 250)

# distP = ifelse(tmp$distP,     runif(1, 0, 0.2), 0.01)
# b0_e = ifelse(tmp$b0_e,       runif(1, 0, 1.0), 0.5)
# b1_recr = ifelse(tmp$b1_recr, runif(1, 0, 50), 20)
# b2_recr = ifelse(tmp$b2_recr, runif(1, -5, 0.0), -0.2)
# b1_g = ifelse(tmp$b1_g,       runif(1, 0, 1.0), 3)
# b2_g = ifelse(tmp$b2_g,       runif(1, -5, 0), -0.1)

logLik = function(par) {
  ll = 
    replicate(10, {
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
      sum((data_tmp-obs$N)**2)/500
    })
  return(-mean(ll))
}

bs = createBayesianSetup(logLik, prior = createUniformPrior(c(0,0,0,-5,0,-5), c(0.2, 1.0, 50, 0.0, 1.0, 0)))
res = runMCMC(bs, settings = list(iterations = 50000))
plot(res, start = 8000)
MAP(res, start = 5000)$par
