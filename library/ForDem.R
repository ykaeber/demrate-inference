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

weather = function() runif(1,0,1)

regeneration_f = function(cohorts, LAI, env, pars){
  speciesPars = pars$speciesPars
  regP = c()
  baseReg = pars$baseReg
  for(i in 1:length(pars[["species"]])){
    iSpID= pars[["species"]][i]
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

parabulaF <- function(x, b0, b1){
  (b0*x + b1*(x)^2)/(b0^2 + b1*(b0)^2)
}

growth_f = function(cohorts, LAI, env, pars){
  for (i in 1:length(cohorts)) {
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

mortality_f = function(cohorts, LAI, env, pars){
  for (i in 1:length(cohorts)) {
    nTrsDead = rbinom(1,1,0.1)
    cohorts[[i]]$nTrs = max(c(0, round(cohorts[[i]]$nTrs-nTrsDead)))
  }
  return(cohorts)
}


runModel <- function(pars){
  cohorts = pars[["initPop"]]
  timesteps = pars[["timesteps"]]
  out = data.table()
  
  if(!is.null(cohorts)){
    plot(0,sum(sapply(cohorts, function(x) x$nTrs)), xlim = c(0, pars$timesteps), ylim = c(0,1000))
  }else{
    plot(0,0, xlim = c(0, pars$timesteps), ylim = c(0,1000))
  }
  
  for(t in 1:timesteps){
    LAI = lai_f(cohorts, pars)
    env = weather()
    cohorts = regeneration_f(cohorts, LAI, env, pars)
    cohorts = growth_f(cohorts, LAI = LAI, env = env, pars)
    cohorts = mortality_f(cohorts, LAI = LAI, env = env, pars)
    #cohorts = growth_f(cohorts, LAI[t], env[t])
    points(t,sum(sapply(cohorts, function(x) x$nTrs)))
    #cohorts = mortality_f(cohorts, LAI[t], env[t], pars)
    cohorts_dt <- data.table(matrix(unlist(cohorts), ncol = length(cohorts[[1]]), byrow = T))
    colnames(cohorts_dt) <- names(cohorts[[1]])
    out = rbind(out, data.table(cohorts_dt, t = t))
  }
  return(out)
}

