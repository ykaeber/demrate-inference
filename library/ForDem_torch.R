library(data.table)
source("library/ForDem-species.R")
library(sjSDM)
torch = sjSDM:::pkg.env$torch
dt = torch$float32
as_ft = function(x) torch$tensor(x, dtype = dt)
as_ftg = function(x) torch$tensor(x, dtype = dt, requires_grad = TRUE)
to_r = function(r) reticulate::py_to_r(r$data$cpu()$numpy())
pyro = reticulate::import("pyro")

################################################################################
## define model default variables
################################################################################

# foliage function from ForClim
gFolA_f <- function(kA1, kA2, kC1, kC2, dbh, areaHectar = as_ft(0.01)) {
  #gFolW = kA1 * exp( kA2*log(dbh)) * kC1;
  gFolW = kA1$multiply( torch$exp(kA2$multiply(kA2, dbh$log())) )$multiply(kC1)
  #gFolA = kC2 * gFolW / kC1
  gFolA = kC2$multiply(gFolW)$div(kC1)$div(areaHectar)$div(as_ft(10000))$div(as_ft(2.0))
  #gFolA = gFolA/areaHectar/10000 # 0.01 is the patch size in ha
  return(gFolA)
}

# function for calculating the lai from the current cohorts
lai_f <- function(cohorts, pars){
  speciesPars = pars$speciesPars
  lai = as_ft(0.0)
  for(iCohort in cohorts){
    iSpID = iCohort[["spID"]]
    kA1 = speciesPars[speciesPars$spID == iSpID,]$kA1[[1]]
    kA2 = speciesPars[speciesPars$spID == iSpID,]$kA2[[1]]
    kC1 = speciesPars[speciesPars$spID == iSpID,]$kC1[[1]]
    kC2 = speciesPars[speciesPars$spID == iSpID,]$kC2[[1]]
    dbh = iCohort[["dbh"]]
    nTrs = iCohort[["nTrs"]]
    lai = lai$add(gFolA_f(kA1, kA2, kC1, kC2, dbh)$multiply(nTrs))
  }
  return(lai)
}

envF = function(envMean, envSD=as_ft(0.1), env){
  maxDens = torch$distributions$Normal( envMean,  envSD)$log_prob(envMean)$exp()
  pDens = torch$distributions$Normal( envMean,  envSD)$log_prob(env)$exp()
  return(pDens$div(maxDens))
}


AL_F = function(LAI) torch$exp(as_ft(-0.25)$multiply(LAI))

shadeF = function(LAI, shadeMean, shadeSD = as_ft(0.1)){
  AL = AL_F(LAI)
  #minP = pnorm(0, mean = shadeMean, sd = shadeSD, lower.tail = F)
  minP = as_ft(1.0)$sub(torch$distributions$Normal(shadeMean, shadeSD)$cdf(as_ft(0.0)))
  #shadeCond = pnorm(AL, mean = shadeMean, sd = shadeSD, lower.tail = F)
  shadeCond = as_ft(1.0)$sub(torch$distributions$Normal(shadeMean, shadeSD)$cdf(AL))
  shadeCond = minP$sub(shadeCond)$div(minP)
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
    shadeMean = speciesPars[speciesPars$spID == iSpID,]$kLy[[1]]
    shadeCond = shadeF(LAI, shadeMean = shadeMean)
    
    envMean = speciesPars[speciesPars$spID == iSpID,]$kDDMin[[1]]$div(as_ft(1100))
    envCond = envF(envMean = envMean, env = env)
    
    regP = append(regP, torch$min(torch$cat(c(as_ft(1.0)$view(1L), shadeCond$multiply(envCond)$view(1L))))$view(1L))
    # if(shadeCond*envCond != 0){
    #   regP[i] = min(c(shadeCond+envCond,1))
    #   
    # }else(
    #   regP[i] = 0
    # )
  }
  regP = torch$cat(regP)
  baseReg = as_ft(10)
  if(sum(to_r(regP)) == 0) regTrsMeans = regP else regTrsMeans = regP$div(regP$sum()$add(baseReg))
  
  if(!is.null(cohorts)) newID = cohorts[length(cohorts)][[1]][["cohortID"]] else newID = 0
  for(i in 1:length(pars[["actualSpecies"]])){
    iSpID = pars[["actualSpecies"]][i]
    regTrsMeans_tmp = torch$index_select(regTrsMeans, 0L, torch$tensor(as.integer(i-1L)))
    if(to_r(regTrsMeans_tmp) > 0){
      #newTrs = rpois(1, regTrsMeans)
      newTrs = torch$round(regTrsMeans_tmp$add(torch$distributions$Normal(0.0, 0.5)$rsample(list(1L)))$exp())
      if(to_r(newTrs) > 0){
        newID = newID + 1
        cohorts[[newID]] = list(
          cohortID = newID,
          spID = iSpID,
          nTrs = newTrs,
          dbh = as_ft(1)
        )
      }
    }
  }
  return(cohorts)
}

parabulaF <- function(x, b0, b1){
  #(b0*x + b1*(x)^2)/(b0^2 + b1*(b0)^2)
  b0$multiply(x)$add( b1$multiply(x$pow(2.0)) )$div( b0$pow(2.0)$add(b1$multiply(b0$pow(2.0))) )
}


growth_f = function(cohorts, LAI, env, pars){
  aliveCohorts = which(sapply(cohorts, function(x) to_r(x$nTrs)>0))
  for (i in aliveCohorts) {
    speciesPars = pars$speciesPars
    D = cohorts[[i]]$dbh
    iSpID = cohorts[[i]]$spID
    kLa = speciesPars[speciesPars$spID == iSpID,]$kLa[[1]]
    kHMax = speciesPars[speciesPars$spID == iSpID,]$kHMax[[1]]$multiply(as_ft(100))
    kG = speciesPars[speciesPars$spID == iSpID,]$kG[[1]]
    envMean = speciesPars[speciesPars$spID == iSpID,]$kDDMin[[1]]$div(as_ft(1100))
    
    kB1 = as_ft(137)
    #kE1 = 14 * (kLa / 3 + 3) + 13       #NH 2018       
    kE1 = as_ft(14)$multiply(kLa$div(as_ft(3.0))$add(as_ft(3.0)) )$add(as_ft(13))
    kSMin = as_ft(1.3)$multiply(kLa$div(as_ft(3))$add(as_ft(3.0)))$add( as_ft(39.5))  #NH 2018
    
    kSIn = kSMin$add(as_ft(0.75)$multiply( kE1 ))     #kS-value for height initialisation
    
    H = kB1$add(kHMax$sub(kB1))$multiply(as_ft(1)$sub(torch$exp(kSIn$neg()$multiply( D )$div(kHMax$sub( kB1)))))
    
    AL = AL_F(LAI)
    #Calculation of gS and related variables [LR]
    gS = kSMin$add( kE1$multiply (as_ft(1.0)$sub(AL)))
    gFun = torch$max(   torch$cat(c(as_ft(0)$view(1L),  gS$multiply(as_ft(1.0)$sub(H$sub(kB1))$div(kHMax$sub(kB1)))$view(1L) )  ))
    
    gGRF = envF(envMean = envMean, env = env)
    #New diameter growth rate [LR]
    #gRateD = max(c(0, gGRF*kG*D*(1-H/kHMax)/(2*H+gFun*D)))
    gRateD = torch$max( torch$cat(c(as_ft(0.0)$view(1L), gGRF$multiply(kG)$multiply(D)$multiply(as_ft(1)$sub(H$div(kHMax)))$div(as_ft(2.0)$multiply(H)$add(gFun$multiply(D)))$view(1L)      ) ))
    
    cohorts[[i]]$dbh = cohorts[[i]]$dbh$add( gRateD )
  }
  return(cohorts)
}

mortality_f = function(cohorts, LAI, env, tDist, pars){
  aliveCohorts = which(sapply(cohorts, function(x) to_r(x$nTrs)>0))
  for (i in aliveCohorts) {
    speciesPars = pars$speciesPars
    kAlpha = pars$bgMort
    D = cohorts[[i]]$dbh
    iSpID= cohorts[[i]]$spID
    kDMax = speciesPars[speciesPars$spID == iSpID,]$kDMax[[1]]
    
    shadeMean = speciesPars[speciesPars$spID == iSpID,]$kLy[[1]]
    shadeCond = shadeF(LAI, shadeMean = shadeMean)
    
    envMean = speciesPars[speciesPars$spID == iSpID,]$kDDMin[[1]]$div(as_ft(1100))
    envCond = envF(envMean = envMean, env = env)
    
    # Cohort-specific size-dependent mortality
    gPSize = as_ft(0.1)$multiply ( (D$div ( kDMax))$pow(kAlpha)         )           #annual mortality rate at kDMax = 10%
    
    #mortP = min(c((1-shadeCond) + (1-envCond) + gPSize + tDist, 1))
    mortP = torch$min(torch$cat(c(as_ft(2.0)$sub(shadeCond)$sub(envCond)$add( gPSize )$add(tDist)$view(1L), as_ft(1.0)$view(1L)   )))
    #print(paste0("sp=",speciesPars[spID == iSpID]$species,"    shadeMean=",shadeMean, "    envMean=", envMean, "   mortP=", mortP))
    if(to_r(cohorts[[i]]$nTrs) > 0){
      
      #nTrsDead = torch$round(torch$distributions$Normal(0.0, 1.0)$cdf(torch$distributions$Normal(as_ft(1.0)$sub(mortP), 0.2)$rsample(list(cohorts[[i]]$nTrs))))$sum()
      nTrsDead = torch$round(torch$nn$functional$relu(mortP$multiply(cohorts[[i]]$nTrs)$add(torch$distributions$Normal(0.0, 0.2)$rsample(list(1L)))))
      
      #nTrsDead = sum(rbinom(cohorts[[i]]$nTrs,1,mortP))
      #cohorts[[i]]$nTrs = max(c(0, cohorts[[i]]$nTrs-nTrsDead))
      cohorts[[i]]$nTrs = torch$max(torch$cat(c(as_ft(0)$view(1L), cohorts[[i]]$nTrs$sub(nTrsDead)$view(1L))))
    }
  }
  
  return(cohorts)
}


runModel <- function(pars){
  speciesPars = pars[["speciesPars"]]
  cohorts = pars[["initPop"]]
  timesteps = pars[["timesteps"]]
  out = data.table()
  
  for(t in 1:timesteps){
    LAI = lai_f(cohorts, pars)
    tDist = as_ft(rbinom(1,1,to_r(pars[["distP"]])))
    # env = weather()
    env = speciesPars[speciesPars$spID == 9,]$kDDMin[[1]]$div(as_ft(1100))
    cohorts = regeneration_f(cohorts, LAI, env, pars)
    cohorts = growth_f(cohorts, LAI = LAI, env = env, pars)
    cohorts = mortality_f(cohorts, LAI = LAI, env = env, tDist = tDist, pars)
    #cohorts = growth_f(cohorts, LAI[t], env[t])
    #points(t,sum(sapply(cohorts, function(x) x$nTrs)))
    #cohorts = mortality_f(cohorts, LAI[t], env[t], pars)
    cohorts_dt <- data.table(matrix(unlist(cohorts), ncol = length(cohorts[[1]]), byrow = T))
    colnames(cohorts_dt) <- names(cohorts[[1]])
    out = rbind(out, data.table(cohorts_dt, t = t))
  }
  return(out)
}

