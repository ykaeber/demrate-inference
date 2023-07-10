################################################################################
## Project: ForDem
## Script purpose: example for running ForDem
## Date:
## Author:
################################################################################
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


cohortsIN <- list(
  list(
    cohortID = 1,
    spID = 1,
    nTrs = 1,
    dbh = 1
  ),
  list(
    cohortID = 2,
    spID = 17,
    nTrs = 4,
    dbh = 30
  )
)

parsModel <- list(
  timesteps = 500,
  actualSpecies = c(0,13,1,17,5),
  baseReg = 10,
  bgMort = 2.3,
  distP = 0.01,
  env = 0.7,
  patchesN = 10,
  areaHectar = 0.01,
  initPop = cohortsIN,
  speciesPars = selected_species_pars
)

library("Rcpp")
sourceCpp("library/fordem.cpp")
system.time(
runModel(pars = parsModel, speciesPars = selected_species_pars)
)
system.time(
  out <- capture.output(
    runModel(pars = parsModel, speciesPars = selected_species_pars)
  )
)
library(jointprof)
out_file <- tempfile("jointprof", fileext = ".out", tmpdir = ".")
start_profiler(out_file)
out <- capture.output(
  runModel(pars = parsModel, speciesPars = selected_species_pars)
)
stop_profiler()

out1 <- fread(paste0(out, collapse = "\n") )
out1 <- merge(out1, selected_species_pars[,c("spID","species")], by = "spID")

allout1 <- out1[,.(baHectar = mean(baTot)*100), by = .(species, t)]
ggplot(allout1, aes(x = t, y = baHectar))+
  geom_line(aes(color = species))


hist(out1$baHectar)

ggplot(out1, aes(x = t, y = baHectar))+
  geom_line(aes(color = species, group = factor(p)), alpha = 0.5)

ggplot(out1, aes(x = t, y = nTrsSum*dbhMean))+
  geom_line(aes(color = species))



ggplot(out1, aes(x = t, y = 100*nTrsSum*pi*(((dbhMean/100)/2)^2)))+
  geom_line(aes(color = species))


