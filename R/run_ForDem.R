################################################################################
## Project: ForDem
## Script purpose: example for running ForDem
## Date:
## Author:
################################################################################
source("library/ForDem.R")

species_pars
cohorts = cohortsIN
cohortsIN <- list(
  list(
    cohortID = 1,
    spID = 1,
    nTrs = 3,
    dbh = 2
  ),
  list(
    cohortID = 2,
    spID = 17,
    nTrs = 40,
    dbh = 20
  )
)

# parsSpecies <- list(
#   spID = c(0,9,17,13),
#   shadeMean = runif(3),
#   shadeSD = runif(3,0,1),
#   envMean = runif(3),
#   envSD = runif(3,0,1)
# )

parsModel <- list(
  timesteps = 500,
  species = c(0,9,17,13),
  baseReg = 10,
  #intrinsicGrowth = 2,
  initPop = cohortsIN,
  #initPop = NULL,
  speciesPars = species_pars
)

cohorts  <- runModel(pars = parsModel)
library(ggplot2)
str(cohorts)
ggplot(cohorts, 
       aes(y = dbh, x = t))+
  geom_line(aes(color = factor(spID), group = cohortID))

ggplot(cohorts, 
       aes(y = nTrs, x = t))+
  geom_line(aes(color = factor(spID), group = cohortID))



ggplot(cohorts[, .(nTrs = sum(nTrs), dbh = mean(dbh)), by = .(spID, t)],
       aes(y = dbh, x = t))+
  geom_line(aes(color = factor(spID)))



ggplot(cohorts[, .(nTrs = sum(nTrs), dbh = mean(dbh)), by = .(spID, t)],
       aes(y = nTrs*dbh, x = t))+
  geom_line(aes(color = factor(spID)))

species_pars

