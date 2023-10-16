

envMfunctions = list(
  f1 = function(par, env) {
    return(apply(1*(sapply(1:ncol(env), function(i) env[,i] > par[,i])),1, sum))
  },
  f2 = function(par, env) {
    if(is.matrix(par)) {
      maxDens = dnorm(par, mean = par, sd = 0.1)
      pDens = lapply(1:nrow(par), function(j) {sapply(1:ncol(env), function(i) dnorm(env[,i,drop=FALSE], mean = par[j,i], sd = 0.1))})
      pDens = abind::abind(pDens, along = 0L)
      return(apply(pDens/maxDens, 1, sum))
    } else {
      maxDens = dnorm(par, mean = par, sd = 0.1)
      pDens =  dnorm(env, mean = par, sd = 0.1)
      return(pDens/maxDens)
    }
  }
)

stateF = list(
  height = function(cohortMat, par) {
    return(exp((cohortMat[,2]*par)*0.03)) # speciesPars[cohortMat[,4],4]
  },
  BA = function(cohortMat) {
    return(pi*(cohortMat[,2]/100/2)**2)
  }
)

compF = function(cohortMat, minLight = 50, height) {
  BA = stateF$BA(cohortMat)/0.1
  height = sapply(height, function(i ) sum(height[height > i]))
  BA_height = height
  AL = 1-BA_height/minLight
  AL[AL < 0] = 0
  AL = AL/max(AL)
  return(AL)
}

regF = function(cohortMat, timestep, parReg) {
  AL = compF(cohortMat, height = 0)
  regP = 1*(AL >abs( parReg[,1] ))
  environment = envMfunctions$f2( par = parReg[,2:3], env = envM[timestep,,drop=FALSE])
  regeneration = rpois(nrow(parReg), exp(regP + environment))
  return(regeneration)
}

mortF = function(cohortMat, timestep, parMort) {
  # shade
  AL = compF(cohortMat, height = stateF$height(cohortMat, parGlobal[cohortMat[,4],1]))
  Shade = envMfunctions$f2(par = parMort[cohortMat[,4],1], env = AL)
  # Umwelt
  environment = envMfunctions$f2(env = envM[timestep,,drop=FALSE], par = parMort[cohortMat[,4],2:3])
  # size
  gPSize = 0.1*(cohortMat[,2]/parMort[cohortMat[,4],4])^2.3
  mortP = sapply(gPSize+Shade+environment, function(x) min(c(x, 1)))
  
  return(rbinom(nrow(cohortMat), cohortMat[,3] , mortP))
}

# last parameter % of maximal Growth
growthF = function(cohortMat, timestep, parGrowth) {
  # Shade
  AL = compF(cohortMat, height = stateF$height(cohortMat, parGlobal[cohortMat[,4],1]))
  shade = envMfunctions$f2(env =  AL, par = parGrowth[cohortMat[,4], 1])
  # Environment
  environment = envMfunctions$f2(env = envM[timestep,,drop=FALSE], par = parGrowth[cohortMat[,4],2:3])
  growth = (1-(1-environment+shade)^4) * parGrowth[cohortMat[,4], 4]
  growth = sapply(growth, function(x) max(c(x, 0)))
  return(growth)
}




runModel = function(envM, 
                    envMfunctions = NULL,
                    cohortMat = NULL, 
                    compF = NULL,
                    stateF = NULL,
                    growthF = NULL,
                    mortF = NULL,
                    regF = NULL,
                    parGlobal = NULL,
                    parReg = NULL,
                    parMort = NULL,
                    parGrowth = NULL,
                    timesteps = 100) {
  
  results = vector("list", timesteps+1)
  results[[1]] = cohortMat
  for(tstep in 1:timesteps) {
    g = growthF(cohortMat, tstep, parGrowth)
    cohortMat[,2] = cohortMat[,2] + g
    m = mortF(cohortMat, tstep, parMort)
    cohortMat[,3] = cohortMat[,3] - m
    cohortMat = cohortMat[cohortMat[,3] > 0,]
    r = regF(cohortMat, tstep, parReg)
    if(length(r) > 0) {
      newCohortMat = data.frame(ID=NA, 
                                dbh = 1, 
                                nTree = r[r>0], 
                                Species = which(r>0, arr.ind = TRUE))
      cohortMat = rbind(cohortMat, newCohortMat)
      
    }
    results[[tstep+1]] = cohortMat
  }
  return(results)
}
cohortMat = cbind(ID = 1:10, 
                  dbh = runif(10, 30, 100), 
                  nTree = rpois(10, 10), 
                  Species = sample.int(5, 5, replace = TRUE))

## envM
# matrix, row = Zeit, columns = Predictors
envM = matrix(rnorm(500*2), 500, 2)

# Parameter
parGlobal = matrix(runif(5), 5, 1) 
parReg = matrix(runif(5*3,0, 3), 5, 3)
parMort = matrix(runif(5*4),5,4)
parMort[,4] = runif(5, 250, 350)
parGrowth = matrix(runif(5*4, 0, 7), 5, 4)


res = runModel(envM = envM, 
               envMfunctions = envMfunctions, 
               cohortMat = cohortMat,
               compF = compF, 
               stateF = stateF,
               growthF = growthF,
               mortF = mortF,
               regF = regF,
               parGlobal = parGlobal,
               parReg = parReg,
               parMort = parMort,
               parGrowth = parGrowth, timesteps = 500L)
library(tidyverse)
RR = lapply(res, function(i) i %>% as.data.frame %>%  mutate(Species = as.character(Species) )%>%  group_by(Species) %>% summarise(dbh = sum(dbh)))
plot(unlist(sapply(RR, function(r) r[r$Species == "2", 2])), type = "l")


