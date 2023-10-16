library(data.table)
source("library/ForDem-species.R")

################################################################################
## define model default variables
################################################################################

# foliage function from ForClim
gFolA_f <- function(kA1, kA2, kC1, kC2, dbh, areaHectar = 0.01) {
  gFolW = kA1 * exp( kA2*log(dbh)) * kC1;
  gFolA = kC2 * gFolW / kC1
  gFolA = gFolA/areaHectar/10000 # 0.01 is the patch size in ha
  return(gFolA/2)
}

# function for calculating the lai from the current cohorts
lai_f <- function(cohorts, pars){
  speciesPars = pars$speciesPars
  lai = 0
  for(iCohort in cohorts){
    iSpID = iCohort[["spID"]]
    kA1 = speciesPars[spID == iSpID,]$kA1
    kA2 = speciesPars[spID == iSpID,]$kA2
    kC1 = speciesPars[spID == iSpID,]$kC1
    kC2 = speciesPars[spID == iSpID,]$kC2
    dbh = iCohort[["dbh"]]
    nTrs = iCohort[["nTrs"]]
    lai = lai + gFolA_f(kA1, kA2, kC1, kC2, dbh)*nTrs
  }
  return(lai)
}

envF = function(envMean, envSD=0.1, env){
  maxDens = dnorm(envMean, mean = envMean, sd = envSD)
  pDens = dnorm(env, mean = envMean, sd = envSD)
  pDens/maxDens
}

AL_F = function(LAI) exp(-0.25*LAI)

shadeF = function(LAI, shadeMean, shadeSD = 0.1){
  AL = AL_F(LAI)
  minP = pnorm(0, mean = shadeMean, sd = shadeSD, lower.tail = F)
  shadeCond = pnorm(AL, mean = shadeMean, sd = shadeSD, lower.tail = F)
  shadeCond = (minP-shadeCond)/minP
  return(shadeCond)
}
#x = seq(0,50,0.1)
#plot(x, sapply(x, function(x1) shadeF(x1, shadeMean = 0.01, shadeSD = 0.3)), ylim = c(0,1))

#weather = function() runif(1,0,1)

# weather = function(speciesPars) speciesPars[spID == 17]$kDDMin/1100

regeneration_f = function(cohorts, LAI, env, pars){
  speciesPars = pars$speciesPars
  regP = c()
  baseReg = pars$baseReg
  for(i in 1:length(pars[["actualSpecies"]])){
    iSpID= pars[["actualSpecies"]][i]
    shadeMean = speciesPars[spID == iSpID]$kLy
    shadeCond = shadeF(LAI, shadeMean = shadeMean)
    
    envMean = speciesPars[spID == iSpID]$kDDMin/1100
    envCond = envF(envMean = envMean, env = env)
    
    if(shadeCond*envCond != 0){
      regP[i] = min(c(shadeCond+envCond,1))
      
    }else(
      regP[i] = 0
    )
  }
  
  baseReg = 10
  if(sum(regP) == 0) regTrsMeans = regP else regTrsMeans = regP/sum(regP)*baseReg
  
  if(!is.null(cohorts)) newID = cohorts[length(cohorts)][[1]][["cohortID"]] else newID = 0
  for(i in 1:length(pars[["actualSpecies"]])){
    iSpID= pars[["actualSpecies"]][i]
    if(regTrsMeans[i] > 0){
      newTrs = rpois(1, regTrsMeans)
      if(newTrs > 0){
        newID = newID + 1
        cohorts[[newID]] = list(
          cohortID = newID,
          spID = iSpID,
          nTrs = newTrs,
          dbh = 1
        )
      }
    }
  }
  return(cohorts)
}

parabulaF <- function(x, b0, b1){
  (b0*x + b1*(x)^2)/(b0^2 + b1*(b0)^2)
}


growth_f = function(cohorts, LAI, env, pars){
  aliveCohorts = which(sapply(cohorts, function(x) x$nTrs>0))
  for (i in aliveCohorts) {
    speciesPars = pars$speciesPars
    D = cohorts[[i]]$dbh
    iSpID = cohorts[[i]]$spID
    kLa = speciesPars[spID == iSpID]$kLa
    kHMax = speciesPars[spID == iSpID]$kHMax*100
    kG = speciesPars[spID == iSpID]$kG
    envMean = speciesPars[spID == iSpID]$kDDMin/1100
    
    #kHMax = kHMax*envF(envMean = envMean, env = env)
    
    kB1 = 137
    #kSMin = 1.3 * kLa + 39.5            #Rasche et al. 2012
    #kE1 = 14 * kLa + 13                 #Rasche et al. 2012
    # alternative formulation in ForClim Variant from Nica Huber
    kE1 = 14 * (kLa / 3 + 3) + 13       #NH 2018       
    kSMin = 1.3 * (kLa / 3 + 3) + 39.5  #NH 2018
    
    kSIn = kSMin + 0.75 * kE1           #kS-value for height initialisation
    
    H = kB1 + (kHMax - kB1) * (1 - exp(-kSIn * D / (kHMax - kB1)))
    
    AL = AL_F(LAI)
    #Calculation of gS and related variables [LR]
    gS = kSMin + kE1 * (1.0-AL);
    gFun = max(c(0,  gS*(1-(H-kB1)/(kHMax-kB1))))
    
    gGRF = envF(envMean = envMean, env = env)
    #New diameter growth rate [LR]
    gRateD = max(c(0, gGRF*kG*D*(1-H/kHMax)/(2*H+gFun*D)))
    
    cohorts[[i]]$dbh = cohorts[[i]]$dbh + gRateD
  }
  return(cohorts)
}

mortality_f = function(cohorts, LAI, env, tDist, pars){
  aliveCohorts = which(sapply(cohorts, function(x) x$nTrs>0))
  for (i in aliveCohorts) {
    speciesPars = pars$speciesPars
    kAlpha = pars$bgMort
    D = cohorts[[i]]$dbh
    iSpID= cohorts[[i]]$spID
    kDMax = speciesPars[spID == iSpID]$kDMax
    
    shadeMean = speciesPars[spID == iSpID]$kLy
    shadeCond = shadeF(LAI, shadeMean = shadeMean)
    
    envMean = speciesPars[spID == iSpID]$kDDMin/1100
    envCond = envF(envMean = envMean, env = env)
    
    # Cohort-specific size-dependent mortality
    gPSize = 0.1 * (D / kDMax)^kAlpha                    #annual mortality rate at kDMax = 10%
    
    mortP = min(c((1-shadeCond) + (1-envCond) + gPSize + tDist, 1))
    #print(paste0("sp=",speciesPars[spID == iSpID]$species,"    shadeMean=",shadeMean, "    envMean=", envMean, "   mortP=", mortP))
    if(cohorts[[i]]$nTrs > 0){
      nTrsDead = sum(rbinom(cohorts[[i]]$nTrs,1,mortP))
      cohorts[[i]]$nTrs = max(c(0, cohorts[[i]]$nTrs-nTrsDead))
    }
  }
  
  return(cohorts)
}


################ New implementation ###############


## speciesPars
# matrix, row = Species, columns = traits 
speciesPars = matrix(rnorm(5*11), 5, 11) 

## cohortMat
# matrix, row = Cohorte, columns = c(ID, outpus (e.g. dbh))
cohortMat = cbind(ID = 1:10, 
                  dbh = runif(10, 30, 100), 
                  nTree = rpois(10, 10), 
                  Species = sample.int(5, 5, replace = TRUE))

## envM
# matrix, row = Zeit, columns = Predictors
envM = matrix(rnorm(500*2), 500, 2)

compMat = matrix(0, nrow = 10L, ncol = 1L)

## envMfunctions
# list, with n = ncol(envM) functions
# ~E1*T2 + E2*T2 + E1*T3
# ET = abind::abind(list(outer(envM[,1], speciesPars[,2]), outer(envM[,2],speciesPars[,2]), outer(envM[,1],speciesPars[,3])), along = 0L)
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
  # HeightClass = function(cohortMat) {
  #   height = stateF$height(cohortMat)
  #   heightIntervals = createHeightIntervals(height)
  #   return(sapply(height, function(i) sum(i > heightIntervals)))
  # },
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

parGlobal = matrix(runif(5), 5, 1) 

parReg = matrix(runif(5*3), 5, 3)
regF = function(cohortMat, timestep, parReg) {
  AL = compF(cohortMat, height = 0)
  regP = 1*(AL >abs( parReg[,1] ))
  environment = envMfunctions$f2( par = parReg[,2:3], env = envM[timestep,,drop=FALSE])
  regeneration = rpois(nrow(speciesPars), exp(regP + environment))
  return(regeneration)
}


parMort = matrix(runif(5*4),5,4)
parMort[,4] = runif(5, 250, 350)
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

parGrowth = matrix(runif(5*4), 5, 4)
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


runModel = function(speciesPars, 
                     envM, 
                     envMfunctions = NULL, 
                     cohortMat, 
                     compF = NULL, 
                     stateF = NULL,
                     timesteps = 100) {
  
  old = cohortMat
  cohortMat = old
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
  }
  
}
