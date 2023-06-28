################################################################################
## Project: ForDem
## Script purpose: example for running ForDem
## Date:
## Author:
################################################################################
source("library/ForDem.R")

species_pars
cohorts <- list(
  list(
    cohortID = 1,
    spID = 1,
    nTrs = 3,
    dbh = 2
  ),
  list(
    cohortID = 2,
    spID = 17,
    nTrs = 1,
    dbh = 5
  )
)

Nspecies = 3
parsSpecies <- list(
  spID = c(0,9,17,13),
  shadeMean = runif(3),
  shadeSD = runif(3,0,1),
  envMean = runif(3),
  envSD = runif(3,0,1)
)

parsModel <- list(
  timesteps = 100,
  species = c(0,9,17,13),
  initPop = cohorts
)

pars = parsModel



