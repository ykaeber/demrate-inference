library(data.table)
library(ggplot2)

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
      return(sum(pDens/maxDens))
    }
  },
  f3 = function(par, env) {
    minP = pnorm(0, mean = par, sd = 0.1, lower.tail = F)
    envCond = pnorm(env, mean = par, sd = 0.1, lower.tail = F)
    envCond = (minP-envCond)/minP
    return(envCond)
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
  cohortHeights = stateF$height(cohortMat, parGlobal[cohortMat[,4],1])
  BA_height = sapply(height, function(i) sum(BA[cohortHeights > i]))
  AL = 1-BA_height/minLight
  AL[AL < 0] = 0
  return(AL)
}

regF = function(cohortMat, timestep, parReg) {
  AL = compF(cohortMat, height = 0)
  regP = 1*(AL >abs( parReg[,1] ))
  environment = envMfunctions$f2( par = parReg[,2:3], env = envM[timestep,,drop=FALSE])
  regeneration = rpois(nrow(parReg), 0.9*exp(regP + environment))# TODO add model parameter for 0.9
  return(regeneration)
}

mortF = function(cohortMat, timestep, parMort) {
  # shade
  AL = compF(cohortMat, height = stateF$height(cohortMat, parGlobal[cohortMat[,4],1]))
  Shade = 1-envMfunctions$f3(par = parMort[cohortMat[,4],1], env = AL)
  # Umwelt
  environment = 1-envMfunctions$f2(env = envM[timestep,,drop=FALSE], par = parMort[cohortMat[,4],2:3])
  environment[environment<0] = 0
  # size
  gPSize = 0.1*(cohortMat[,2]/parMort[cohortMat[,4],4])^2.3
  mortP = sapply(gPSize+Shade+environment, function(x) min(c(x, 1)))
  mort = rbinom(nrow(cohortMat), cohortMat[,3] , mortP)
  return(mort)
}

# last parameter % of maximal Growth
growthF = function(cohortMat, timestep, parGrowth) {
  # Shade
  AL = compF(cohortMat, height = stateF$height(cohortMat, parGlobal[cohortMat[,4],1]))
  shade = envMfunctions$f3(env =  AL, par = parGrowth[cohortMat[,4], 1])
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
    if(nrow(cohortMat) > 0){
      g = growthF(cohortMat, tstep, parGrowth)
      cohortMat[,2] = cohortMat[,2] + g
      m = mortF(cohortMat, tstep, parMort)
      cohortMat[,3] = cohortMat[,3] - m
      cohortMat = cohortMat[cohortMat[,3] > 0,]
    }
    r = regF(cohortMat, tstep, parReg)
    if(length(r) > 0 & sum(r) > 0) {
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
# cohortMat = cbind(ID = 1:10, 
#                   dbh = runif(10, 30, 100), 
#                   nTree = rpois(10, 10), 
#                   Species = sample.int(5, 5, replace = TRUE))
# 
# ## envM
# # matrix, row = Zeit, columns = Predictors
# envM = matrix(rnorm(500*2), 500, 2)
# 
# # Parameter
# parGlobal = matrix(runif(5), 5, 1) 
# parReg = matrix(runif(5*3,0, 3), 5, 3)
# parMort = matrix(runif(5*4),5,4)
# parMort[,4] = runif(5, 250, 350)
# parGrowth = matrix(runif(5*4, 0, 7), 5, 4)
# 
# res = runModel(envM = envM, 
#                envMfunctions = envMfunctions, 
#                cohortMat = cohortMat,
#                compF = compF, 
#                stateF = stateF,
#                growthF = growthF,
#                mortF = mortF,
#                regF = regF,
#                parGlobal = parGlobal,
#                parReg = parReg,
#                parMort = parMort,
#                parGrowth = parGrowth, timesteps = 500L)
# 
# 
# 
# out_l = lapply(1:length(res[-1]), function(i) data.table(res[-1][[i]], timestep = i))
# out_dt = rbindlist(out_l[-1])
# circle_area <- function(d) pi * (d/100 / 2)^2
# out_dt[,ba := circle_area(dbh)*nTree,]
# ggplot(out_dt, aes(x = timestep, y = ba, color = factor(Species)))+
#   geom_line()

Nspecies = 4

# species1 Fagus sylvatica
# species2 Picea abies
# species3 Acer pseudoplatanus
# species4 Betula pendula
# species5 Abies alba

species_dt = data.table(
  names = c("Fagus sylvatica", "Picea abies", "Acer pseudoplatanus", "Betula pendula", "Abies alba"),
  dbh2height = c(0.6,0.6,0.6,0.6,0.6),
  Light = c(0.1,0.15,0.5,0.6,0.07),
  Temp = c(0.6,0.3,0.7,0.4,0.3),
  Prec = c(0.6,0.7,0.4,0.4,0.5),
  maxD = c(350,300,200,100,350),
  maxG = c(3,5,7,8,4)
)


# dbh height relationship
parGlobal = matrix(species_dt$dbh2height,ncol = 1)
parReg = as.matrix(species_dt[,.(Light,Temp,Prec)])
parMort = as.matrix(species_dt[,.(Light,Temp,Prec,maxD)])
parGrowth = as.matrix(species_dt[,.(Light,Temp,Prec,maxG)])

cohortMat = cbind(ID = 1:10, 
                  dbh = runif(10, 30, 100), 
                  nTree = rep(0, 10), 
                  Species = sample.int(5, 5, replace = TRUE))

out_dt_all  <- data.table()
env_out_all  <- data.table()
Nsites = 1000
Ntimesteps = 300
Npatches = 100
NsamplesPerEnv = ceiling(sqrt(Nsites))
EnvSiteAvg = seq(0,1,length.out = NsamplesPerEnv)
EnvSiteAvgM = expand.grid(list(EnvSiteAvg,EnvSiteAvg))
cat("start at")
Sys.time()
time = system.time(
  {
    for(i_site in 1:Nsites){
      envAvg1 = EnvSiteAvgM[i_site,1]
      envAvg2 = EnvSiteAvgM[i_site,2]
      rm(envM)
      envM = matrix(c(rnorm(Ntimesteps, envAvg1, 0.1),rnorm(Ntimesteps,envAvg2,0.1)), Ntimesteps, 2)
      envM[envM > 1] = 1
      envM[envM < 0] = 0
      env_out = data.table(
        timestep = 1:Ntimesteps,
        site = i_site,
        envAct1 = envM[,1],
        envAct2 = envM[,2],
        envAvg1 = envAvg1,
        envAvg2 = envAvg2
        )
      if(i_site == 1 & p == 1) fwrite(env_out, file = "evaluation3-env.csv", append = F) else fwrite(env_out, file = "evaluation3-env.csv", append = T)
      env_out_all = rbind(env_out_all, env_out)
      for(p in 1:Npatches){
        cat(paste0("\r running site ", i_site, "/", Nsites, " patch ", p, "/", Npatches, "          "))
        res2 = runModel(envM = envM, 
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
                        parGrowth = parGrowth, timesteps = Ntimesteps)
        
        out_l = lapply(1:length(res2[-1]), function(i) data.table(res2[-1][[i]], timestep = i))
        out_dt = rbindlist(out_l)
        out_dt$patch = p
        out_dt$site = i_site
        if(i_site == 1 & p == 1) fwrite(out_dt, file = "evaluation3-obs.csv", append = F) else fwrite(out_dt, file = "evaluation3-obs.csv", append = T)
        # out_dt_all = rbind(out_dt_all, out_dt)
      }
    }
  })
time

out_dt_all <- fread("evaluation2-obs.csv")
env_out_all <- fread("evaluation2-env.csv")

out_dt_all <- merge(out_dt_all, env_out_all, by = c("site", "timestep"))

ggplot(env_out_all, aes(x = factor(envAvg1), y = envAct1))+
  geom_boxplot()
ggplot(env_out_all, aes(x = factor(envAvg2), y = envAct2))+
  geom_boxplot()

# fwrite(out_dt_all, "evaluation-obs.csv")
# fwrite(env_out_all, "evaluation-env.csv")
# fwrite(envM, "evaluation-env.csv")

circle_area <- function(d) pi * (d/100 / 2)^2
out_dt_all[,ba := (circle_area(dbh)*nTree)/0.1,]
out_dt_all[,SpeciesName := factor(Species, levels = 1:5, labels = paste(1:5, species_dt$names)),]
out_dt_all2 <- out_dt_all[,.(
  ba = sum(ba),
  nTree = sum(nTree)
  ), by = .(site,timestep,SpeciesName,patch, envAvg1, envAvg2)]

p1 <- ggplot(out_dt_all2[,.(
  ba = mean(ba)
  ), by = .(site,timestep,SpeciesName,envAvg1,envAvg2)], aes(x = timestep, y = ba, color = SpeciesName))+
  geom_line()+
  facet_grid(round(envAvg1,2)~round(envAvg2,2))
p1
p2 <- ggplot(out_dt_all2[,.(
  nTree = sum(nTree)
  ), by = .(timestep,SpeciesName)], aes(x = timestep, y = nTree, color = SpeciesName))+
  geom_line()

ggpubr::ggarrange(p1,p2,ncol = 1)

ggplot(out_dt_all[patch == 1,.(timestep,SpeciesName,dbh)], aes(x = cut(timestep, breaks = 10), y = dbh, color = paste(SpeciesName)))+
  geom_boxplot()

View(out_dt_all[patch == 1,.(timestep,SpeciesName,dbh)])
