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
lai_f <- function(cohorts){
  lai = 0
  for(iCohort in cohorts){
    iSpID = iCohort[["spID"]]
    kA1 = species_pars[spID == iSpID,]$kA1
    kA2 = species_pars[spID == iSpID,]$kA2
    kC1 = species_pars[spID == iSpID,]$kC1
    kC2 = species_pars[spID == iSpID,]$kC2
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

shadeF = function(LAI, shadeMean, shadeSD = 0.1){
  AL = exp(-0.25*LAI)
  minP = pnorm(0, mean = shadeMean, sd = shadeSD, lower.tail = F)
  shadeCond = pnorm(AL, mean = shadeMean, sd = shadeSD, lower.tail = F)
  shadeCond = (minP-shadeCond)/minP
  return(shadeCond)
}
#x = seq(0,50,0.1)
#plot(x, sapply(x, function(x1) shadeF(x1, shadeMean = 0.01, shadeSD = 0.3)), ylim = c(0,1))

weather = function() runif(1,0,1)

regeneration_f = function(LAI, env){
  for(i in 1:length(pars[["species"]])){
    iSpID= pars[["species"]][i]
    shadeMean = species_pars[spID == iSpID]$kLy
    shadeCond = shadeF(LAI, shadeMean = shadeMean)
    
    envMean = species_pars[spID == iSpID]$kDDMin/1100
    envCond = envF(envMean = envMean, env = env)
    
    if(shadeCond*envCond != 0){
        regP[i] = min(c(shadeCond+envCond,1))
       
      }else(
        regP[i] = 0
      )
  }
  
  baseReg = 10
  if(sum(regP) == 0) regTrsMeans = regP else regTrsMeans = regP/sum(regP)*baseReg
  
  newID = cohorts[length(cohorts)][[1]][["cohortID"]]
  
  for(i in 1:length(pars[["species"]])){
    iSpID= pars[["species"]][i]
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

growth_f = function(cohorts){
  for (i in 1:length(cohorts)) {
    cohorts[[i]]$dbh = cohorts[[i]]$dbh + 5
  }
  return(cohorts)
}

runModel <- function(pars){
  cohorts = pars[["initPop"]]
  timesteps= pars[["timesteps"]]
  
  plot(0,sum(sapply(cohorts, function(x) x$nTrs)), xlim = c(0, pars$timesteps), ylim = c(0,1000))
  
  for(t in 1:timesteps){
    LAI = lai_f(cohorts)
    env = weather()
    cohorts = regeneration_f(LAI, env)
    cohorts = growth_f(cohorts)
    #cohorts = growth_f(cohorts, LAI[t], env[t])
    points(t,sum(sapply(cohorts, function(x) x$nTrs)))
    #cohorts = mortality_f(cohorts, LAI[t], env[t])
  }
}


