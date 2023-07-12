################################################################################
## Project: ForDem
## Script purpose: example for running ForDem
## Date:
## Author:
################################################################################
library("Rcpp")
library(data.table)
library(ggplot2)
sourceCpp("library/fordem.cpp")

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

sourceCpp("library/fordem.cpp")

# lai1 = seq(0,30, 0.2)
# 
# plot(lai1, sapply(lai1, function(x) shadeF(x, 0.3, 0.1)))
# plot(lai1, sapply(lai1, function(x) shadeF(x, 0.1, 0.1)))
# 
# x = seq(0,1,0.01)
# 
# gGRF = envF(envMean, 0.1, env, widthPar)
# gGRF2 = 1 - pow(1- gGRF,4)
# 
# plot(x, sapply(x, function(x) envF2(0.5, 0.7, x)), ylim = c(0,1))
# abline(v = 0.5)
# 
# f1 <- function(x1) 1 - (1- x1)^4
# plot(x, sapply(x, function(x) f1(envF2(0.5, 0.1, x))))
# abline(v = 0.5)
# 
# 
# plot(x, 0.2*(1-x)^3)

# sourceCpp("library/funtest.cpp")
# cohortsIN <- list(
#   # list(
#   #   cohortID = 1,
#   #   spID = 1,
#   #   nTrs = 1,
#   #   dbh = 1
#   # ),
#   list(
#     cohortID = 2,
#     spID = 17,
#     nTrs = 5,
#     dbh = 50
#   ),
#   list(
#     cohortID = 2,
#     spID = 17,
#     nTrs = 5,
#     dbh = 10
#   )
# )


selected_species_pars <- data.table(selected_species_pars)
selected_species_pars[, env := scale(kDDMin), ]
# selected_species_pars[, env := scale(-kDrTol), ]
# selected_species_pars[, env := scale(kDrTol)*scale(kDDMin), ]

selected_species_pars[, env := env+abs(min(env)), ]
selected_species_pars[, env := env/max(env), ]

sourceCpp("library/fordem.cpp")
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
    distP = 0.01, # range 0 to 0.2
    distInt = 0.01, # (fraction of patches affected) range 0 to 1
    env = 0.4, # range 0 to 1
    nicheWidth = 4,
    # env = 0.7,
    patchesN = 100,
    areaHectar = 0.1,
    heightClassesN = 50,
    initPop = NULL,
    # initPop = cohortsIN,
    speciesPars = selected_species_pars
  )

system.time(
  outMat1 <- runModel(pars = parsModel, speciesPars = selected_species_pars)
)
<<<<<<< HEAD
library(jointprof)
out_file <- tempfile("jointprof", fileext = ".out", tmpdir = ".")
start_profiler(out_file)
out <- capture.output(
  runModel(pars = parsModel, speciesPars = selected_species_pars)
)
stop_profiler()

out1 <- fread(paste0(out, collapse = "\n") )
=======
outDT <- data.table(outMat1)
names(outDT) <- c(parsModel$outVars,"t","p")
out1 <- outDT
>>>>>>> 5d1de56fd276b8f8ca9c9e9482506933100cebfe
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
  facet_wrap(~variable, scales = "free_y")

