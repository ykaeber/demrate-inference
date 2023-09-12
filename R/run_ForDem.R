################################################################################
## Project: ForDem
## Script purpose: example for running ForDem
## Date:
## Author:
################################################################################
library("Rcpp")
library(data.table)
library(ggplot2)
sourceCpp("library/fordem2.cpp")

cohorts <- list(
  list(
    cohortID = 1,
    spID = 1,
    nTrs = 400,
    dbh = 30,
    H = -1
  ),
  list(
    cohortID = 2,
    spID = 17,
    nTrs = 200,
    dbh = 40,
    H = -1
  ),
  list(
    cohortID = 3,
    spID = 0,
    nTrs = 5,
    dbh = 10,
    H = -1
  )
)

mat <- matrix(data = unlist(cohorts), ncol = 5, byrow = T)

selected_species_pars <- 
  data.frame(
    stringsAsFactors = FALSE,
    spID = c(0L,1L,2L,3L,4L,5L,6L,7L,
             8L,9L,10L,11L,12L,13L,14L,15L,16L,17L,18L,19L,
             20L,21L,22L,23L,24L,25L,26L,27L,28L,29L),
    species = c("Abies alba","Larix decidua",
                "Picea abies","Pinus cembra","Pinus montana",
                "Pinus sylvestris","Taxus baccata","Acer campestre",
                "Acer platanoides","Acer pseudoplatanus","Alnus glutinosa",
                "Alnus incana","Alnus viridis","Betula pendula",
                "Carpinus betulus","Castanea sativa","Corylus avellana",
                "Fagus sylvatica","Fraxinus excelsior","Populus nigra",
                "Populus tremula","Quercus petraea","Quercus pubescens",
                "Quercus robur","Salix alba","Sorbus aria","Sorbus aucuparia",
                "Tilia cordata","Tilia platyphyllos","Ulmus glabra"),
    kHMax = c(60L,54L,63L,26L,23L,48L,
              22L,25L,35L,40L,40L,25L,6L,30L,35L,35L,15L,52L,
              42L,40L,42L,50L,25L,52L,35L,23L,27L,40L,40L,
              43L),
    kDMax = c(200L,250L,200L,200L,50L,
              150L,350L,150L,200L,250L,150L,100L,25L,100L,150L,
              350L,50L,250L,250L,250L,150L,350L,150L,350L,250L,
              100L,100L,350L,350L,300L),
    kG = c(296L,400L,342L,198L,239L,
           393L,175L,210L,360L,338L,380L,218L,476L,448L,360L,
           375L,245L,307L,363L,394L,390L,378L,226L,376L,
           403L,230L,205L,365L,365L,361L),
    kDDMin = c(641L,323L,385L,323L,436L,
               610L,1011L,1062L,1042L,898L,898L,610L,272L,610L,
               898L,1237L,898L,723L,980L,662L,610L,785L,1011L,
               1042L,1062L,898L,498L,1339L,1339L,1062L),
    kDrTol = c(0.23,0.25,0.15,
               0.3,0.37,0.37,0.23,0.33,0.25,0.25,0.08,0.08,
               0.16,0.16,0.25,0.33,0.33,0.25,0.16,0.08,
               0.25,0.33,0.33,0.25,0.08,0.33,0.33,0.33,0.25,
               0.25),
    kLy = c(0.03,0.5,0.05,0.075,0.5,
            0.4,0.03,0.2,0.075,0.05,0.2,0.2,0.3,0.5,0.075,0.1,
            0.2,0.03,0.075,0.3,0.3,0.2,0.3,0.2,0.3,0.3,
            0.1,0.075,0.075,0.075),
    kLa = c(1L,9L,5L,5L,9L,8L,1L,6L,
            4L,3L,6L,7L,7L,9L,4L,5L,6L,1L,5L,7L,7L,7L,
            8L,7L,7L,8L,6L,5L,5L,5L),
    kA1 = c(0.23,0.1,0.23,0.23,0.23,
            0.17,0.23,0.1,0.06,0.06,0.1,0.1,0.1,0.08,0.06,
            0.06,0.06,0.06,0.1,0.1,0.1,0.06,0.06,0.06,0.08,0.1,
            0.08,0.06,0.06,0.06),
    kA2 = c(1.56,1.43,1.56,1.56,1.56,
            1.4,1.56,1.43,1.7,1.7,1.43,1.43,1.43,1.43,1.7,
            1.7,1.7,1.7,1.43,1.43,1.43,1.7,1.7,1.7,1.43,1.43,
            1.43,1.7,1.7,1.7),
    kC1 = c(0.45,0.35,0.45,0.45,0.45,
            0.45,0.45,0.35,0.35,0.35,0.35,0.35,0.35,0.35,0.35,
            0.35,0.35,0.35,0.35,0.35,0.35,0.35,0.35,0.35,
            0.35,0.35,0.35,0.35,0.35,0.35),
    kC2 = c(6L,12L,6L,6L,6L,6L,6L,
            12L,12L,12L,12L,12L,12L,12L,12L,12L,12L,12L,12L,
            12L,12L,12L,12L,12L,12L,12L,12L,12L,12L,12L)
  )




selected_species_pars <- data.table(selected_species_pars)
selected_species_pars[, env := scale(kDDMin), ]
# selected_species_pars[, env := scale(-kDrTol), ]
# selected_species_pars[, env := scale(kDrTol)*scale(kDDMin), ]

selected_species_pars[, env := env+abs(min(env)), ]
selected_species_pars[, env := env/max(env), ]

selected_species_pars[spID %in% c(0,2,13,1,17,5, 27)]


sourceCpp("library/fordem2.cpp")


################################################################################
## proposal figures
################################################################################

p_dat_out <- data.table()
p_dat_allSP_out <- data.table()

for(i_patches in c(1,100)){
  parsModel <- list(
    timesteps = 250,
    sampleSteps = 1,
    outVars = c("spID", "lai", "nTrs", "dbh", "ba"),
    # actualSpecies = c(0,2,13,1,17,5, 27),
    actualSpecies = c(0,1,2,5,9,10,13,14,17,22,23,27),
    # actualSpecies = selected_species_pars$spID,
    baseReg = 100/0.1, # range: 1/0.1 to 5000/0.1
    baseRegP = 0.1, # (optional for inference) range 0.01 to 1
    regEnvEff = 0.2, # range 0 to 10
    regShadeEff = 1, # range 0 to 10
    mortEnvEff = 1, # range 0 to 10
    mortShadeEff = 3, # range 0 to 10
    bgMort = 2.3,
    distP = .05, # range 0 to 0.2
    distInt = 0.1, # (fraction of patches affected) range 0 to 1
    env = 0.42, # range 0 to 1
    nicheWidth = 3,
    # env = 0.7,
    patchesN = i_patches,
    areaHectar = 0.1,
    heightClassesN = 10,
    initPop = NULL,
    # initPop = cohortsIN,
    speciesPars = selected_species_pars
  )
  
  timeOut <- paste0(
    capture.output(
      outMat1 <- runModel(pars = parsModel)
    ), collapse = "\n")
  
  time_dt <- fread(text = timeOut)
  
  time_dt[,tot_runtime := sum[TimeMarker == "runModel"],]
  time_dt[,perc_runtime := sum/tot_runtime,]
  
  time_dt[order(perc_runtime)]
  
  outDT <- data.table(outMat1)
  names(outDT) <- c(parsModel$outVars,"t","p")
  out1 <- outDT
  out1 <- merge(out1, selected_species_pars[,c("spID","species")], by = "spID")
  
  p_dat <- out1[, .(
    lai = mean(lai),
    ba = mean(ba),
    dbh = mean(dbh),
    nTrs = mean(nTrs)
  ), by = .(t, species)]
  p_dat_allSP <- p_dat[, .(
    lai = sum(lai),
    ba = sum(ba),
    dbh = mean(dbh),
    nTrs = sum(nTrs),
    species = "all"
  ), by = .(t)]
  
  p_dat <- melt(p_dat, id.vars = c("t", "species"))
  p_dat_allSP <- melt(p_dat_allSP, id.vars = c("t", "species"))
  p_dat$Npatches <- i_patches
  p_dat_allSP$Npatches <- i_patches
  
  p_dat_out <- rbind(p_dat_out, p_dat)
  p_dat_allSP_out <- rbind(p_dat_allSP_out, p_dat_allSP)
  gc()
}


ggplot(
  p_dat_out[variable %in% c("ba", "nTrs")], aes(x = t, y = value))+
  geom_line(aes(color = factor(species)))+
  geom_line(data = p_dat_allSP_out[variable %in% c("ba", "nTrs")], linetype = 2)+
  facet_wrap(paste("N of patches =",Npatches)~variable, scales = "free_y")+
  xlab("time")+
  theme_bw()+
  guides(color = guide_legend(title = "species"))


#sourceCpp("library/fordem-benchemark.cpp")
parsModel <- list(
  timesteps = 500,
  sampleSteps = 1,
  outVars = c("spID", "lai", "nTrs", "dbh", "ba"),
  actualSpecies = c(0,2,13,1,17,5, 27),
  # actualSpecies = selected_species_pars$spID,
  baseReg = 50/0.1, # range: 1/0.1 to 5000/0.1
  baseRegP = 0.1, # (optional for inference) range 0.01 to 1
  regEnvEff = 1, # range 0 to 10
  regShadeEff = 1, # range 0 to 10
  mortEnvEff = 0.2, # range 0 to 10
  mortShadeEff = 1, # range 0 to 10
  bgMort = 2.3,
  distP = .3, # range 0 to 0.2
  distInt = 0.1, # (fraction of patches affected) range 0 to 1
  env = 0.5, # range 0 to 1
  nicheWidth = 4,
  # env = 0.7,
  patchesN = 1,
  areaHectar = 0.1,
  heightClassesN = 10,
  initPop = NULL,
  # initPop = cohortsIN,
  speciesPars = selected_species_pars
)


# system.time(
#   outMat1 <- runModel(pars = parsModel, speciesPars = selected_species_pars)
# )
# mat2 <- updateH(mat, pars = parsModel, speciesPars = selected_species_pars)
# 
# regeneration_f(mat2, 3, 0.4, pars = parsModel, speciesPars = selected_species_pars)
# 
# x = seq(0, 1 , 0.01)
# envCond = envF(envMean = x, envSD = 0.1, env = 0.4, widthPar = 4)
# plot(x, 1*envCond^3)
# 
# x = seq(0, 1 , 0.01)
# envCond = envFVec(envMean = x, envSD = rep(0.1, length(x)), env = 0.4, widthPar = 4)
# plot(x, envCond)

# runModel(pars = parsModel, speciesPars = selected_species_pars)
timeOut <- paste0(
  capture.output(
    outMat1 <- runModel(pars = parsModel)
  ), collapse = "\n")

time_dt <- fread(text = timeOut)

time_dt[,tot_runtime := sum[TimeMarker == "runModel"],]
time_dt[,perc_runtime := sum/tot_runtime,]

time_dt[order(perc_runtime)]

outDT <- data.table(outMat1)
names(outDT) <- c(parsModel$outVars,"t","p")
out1 <- outDT
out1 <- merge(out1, selected_species_pars[,c("spID","species")], by = "spID")

p_dat <- out1[, .(
  lai = mean(lai),
  ba = mean(ba),
  dbh = mean(dbh),
  nTrs = mean(nTrs)
), by = .(t, species)]
p_dat_allSP <- p_dat[, .(
  lai = sum(lai),
  ba = sum(ba),
  dbh = mean(dbh),
  nTrs = sum(nTrs),
  species = "all"
), by = .(t)]

p_dat <- melt(p_dat, id.vars = c("t", "species"))
p_dat_allSP <- melt(p_dat_allSP, id.vars = c("t", "species"))

ggplot(
  p_dat, aes(x = t, y = value))+
  geom_line(aes(color = factor(species)))+
  geom_line(data = p_dat_allSP, linetype = 2)+
  facet_wrap(~variable, scales = "free_y")+
  ggtitle(paste())
# tmp = data.frame(
# "distP" = T,
# "distInt" = T,
# "regEnvEff" = T,
# "regShadeEff" = T,
# "mortEnvEff" = T,
# "mortShadeEff" = T,
# "env" = T
# )
# 
# for(i in 1:100){
#   distP =          ifelse(tmp$distP,       runif(1, 0, 0.2), 0.01)
#   distInt =        ifelse(tmp$distInt,     runif(1, 0, 1.0), 0.01)
#   regEnvEff =      ifelse(tmp$regEnvEff,   runif(1, 0, 10), 1.0)
#   regShadeEff =    ifelse(tmp$regShadeEff, runif(1, 0, 10), 1.0)
#   mortEnvEff =     ifelse(tmp$mortEnvEff,  runif(1, 0, 10), 0.2)
#   mortShadeEff =   ifelse(tmp$mortShadeEff,runif(1, 0, 10), 1.00)
#   env =            ifelse(tmp$env,         runif(1, 0, 1.0), 0.4)
#   
#   parsModel <- list(
#     timesteps = 500,
#     sampleSteps = 1,
#     outVars = c("spID", "lai", "nTrs", "dbh", "ba"),
#     actualSpecies = c(0,2,13,1,17,5, 27),
#     baseReg = 50/0.1, # range: 1/0.1 to 5000/0.1
#     baseRegP = 0.1, # (optional for inference) range 0.01 to 1
#     regEnvEff = regEnvEff, # range 0 to 10
#     regShadeEff = regShadeEff, # range 0 to 10
#     mortEnvEff = mortEnvEff, # range 0 to 10
#     mortShadeEff = mortShadeEff, # range 0 to 10
#     bgMort = 2.3,
#     distP = distP, # range 0 to 0.2
#     distInt = distInt, # (fraction of patches affected) range 0 to 1
#     env = env, # range 0 to 1
#     nicheWidth = 4,
#     patchesN = 100,
#     areaHectar = 0.1,
#     heightClassesN = 50,
#     initPop = NULL,
#     speciesPars = selected_species_pars
#   )
#   outMat1 <- runModel(pars = parsModel, speciesPars = selected_species_pars)
# }
# runModel(pars = parsModel, speciesPars = selected_species_pars)

outDT <- data.table(outMat1)
names(outDT) <- c(parsModel$outVars,"t","p")
out1 <- outDT
out1 <- merge(out1, selected_species_pars[,c("spID","species")], by = "spID")

p_dat <- out1[, .(
  lai = mean(lai),
  ba = mean(ba),
  dbh = mean(dbh),
  nTrs = mean(nTrs)
), by = .(t, species)]
p_dat_allSP <- p_dat[, .(
  lai = sum(lai),
  ba = sum(ba),
  dbh = mean(dbh),
  nTrs = sum(nTrs),
  species = "all"
), by = .(t)]

p_dat <- melt(p_dat, id.vars = c("t", "species"))
p_dat_allSP <- melt(p_dat_allSP, id.vars = c("t", "species"))

p <- ggplot(
  p_dat, aes(x = t, y = value))+
  geom_line(aes(color = factor(species)))+
  geom_line(data = p_dat_allSP, linetype = 2)+
  facet_wrap(~variable, scales = "free_y")+
  ggtitle(paste())

# }
install.packages("rbenchmark")

library(rbenchmark)

parsModel <- list(
  timesteps = 500,
  sampleSteps = 1,
  outVars = c("spID", "lai", "nTrs", "dbh", "ba"),
  actualSpecies = c(0,2,13,1,17,5, 27),
  # actualSpecies = selected_species_pars$spID,
  baseReg = 50/0.1, # range: 1/0.1 to 5000/0.1
  baseRegP = 0.1, # (optional for inference) range 0.01 to 1
  regEnvEff = 1, # range 0 to 10
  regShadeEff = 1, # range 0 to 10
  mortEnvEff = 0.2, # range 0 to 10
  mortShadeEff = 0.7, # range 0 to 10
  bgMort = 2.3,
  distP = .01, # range 0 to 0.2
  distInt = 0.5, # (fraction of patches affected) range 0 to 1
  env = 0.5, # range 0 to 1
  nicheWidth = 4,
  # env = 0.7,
  patchesN = 10,
  areaHectar = 0.1,
  heightClassesN = 10,
  initPop = NULL,
  # initPop = cohortsIN,
  speciesPars = selected_species_pars
)
outMat1 <- runModel(pars = parsModel)
outDT <- data.table(outMat1)
names(outDT) <- c(parsModel$outVars,"t","p")
out1 <- outDT
out1 <- merge(out1, selected_species_pars[,c("spID","species")], by = "spID")
p_dat <- out1[, .(
  lai = mean(lai),
  ba = mean(ba),
  dbh = mean(dbh),
  nTrs = mean(nTrs)
), by = .(t, species)]
p_dat_allSP <- p_dat[, .(
  lai = sum(lai),
  ba = sum(ba),
  dbh = mean(dbh),
  nTrs = sum(nTrs),
  species = "all"
), by = .(t)]
true = p_dat$ba

fun <- function(pars){
  parsModel <- list(
    timesteps = 500,
    sampleSteps = 1,
    outVars = c("spID", "lai", "nTrs", "dbh", "ba"),
    actualSpecies = c(0,2,13,1,17,5, 27),
    # actualSpecies = selected_species_pars$spID,
    baseReg = 50/0.1, # range: 1/0.1 to 5000/0.1
    baseRegP = 0.1, # (optional for inference) range 0.01 to 1
    regEnvEff = 1, # range 0 to 10
    regShadeEff = 1, # range 0 to 10
    mortEnvEff = 0.2, # range 0 to 10
    mortShadeEff = pars, # range 0 to 10
    bgMort = 2.3,
    distP = .01, # range 0 to 0.2
    distInt = 0.5, # (fraction of patches affected) range 0 to 1
    env = 0.5, # range 0 to 1
    nicheWidth = 4,
    # env = 0.7,
    patchesN = 10,
    areaHectar = 0.1,
    heightClassesN = 10,
    initPop = NULL,
    # initPop = cohortsIN,
    speciesPars = selected_species_pars
  )
outMat1 <- runModel(pars = parsModel, speciesPars = selected_species_pars)
outDT <- data.table(outMat1)
names(outDT) <- c(parsModel$outVars,"t","p")
out1 <- outDT
out1 <- merge(out1, selected_species_pars[,c("spID","species")], by = "spID")
p_dat <- out1[, .(
  lai = mean(lai),
  ba = mean(ba),
  dbh = mean(dbh),
  nTrs = mean(nTrs)
), by = .(t, species)]
p_dat_allSP <- p_dat[, .(
  lai = sum(lai),
  ba = sum(ba),
  dbh = mean(dbh),
  nTrs = sum(nTrs),
  species = "all"
), by = .(t)]
sum(sqrt((p_dat$ba - true)^2))
}

out1 = optim(0.1,fun,lower = 0,upper = 1, method = "L-BFGS-B")

out1



