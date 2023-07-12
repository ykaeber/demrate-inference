#include <Rcpp.h>
#include <algorithm>
#include <chrono>
using namespace Rcpp;

// =============================================================================
// helper functions whith little or no relevance for the model formulation
// =============================================================================

double roundToDecimal(double x, int decimalPlaces) {
  double power = std::pow(10, decimalPlaces);
  return std::round(x * power) / power;
}

std::chrono::steady_clock::time_point startTime;

// Function to start the timer
// [[Rcpp::export]]
void startTimer() {
  startTime = std::chrono::steady_clock::now();
}

// Function to print the elapsed time
// [[Rcpp::export]]
void printElapsedTime(std::string additionalString = "") {
  // Stop the timer
  std::chrono::steady_clock::time_point endTime = std::chrono::steady_clock::now();
  
  // Calculate the elapsed time
  std::chrono::duration<double> elapsedSeconds = endTime - startTime;
  
  // Print the elapsed time with additional string
  Rcout << "Elapsed time" << (additionalString.empty() ? "" : " (" + additionalString + ")")
        << ": " << elapsedSeconds.count() << " seconds" << std::endl;
}


// [[Rcpp::export]]
NumericVector revCumSum(NumericVector vec) {
  int n = vec.size();
  NumericVector result(n);
  double sum = 0.0;
  
  for (int i = n - 1; i >= 0; i--) {
    sum += vec[i];
    result[i] = sum;
  }
  
  return result;
}


// [[Rcpp::export]]
NumericMatrix repMat(NumericMatrix mat, int N) {
  int numRows = mat.nrow();
  int numCols = mat.ncol();
  int newNumCols = numCols + 1;
  int newNumRows = numRows * N;
  
  NumericMatrix result(newNumRows, newNumCols);
  
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < numRows; j++) {
      for (int k = 0; k < numCols; k++) {
        result(i * numRows + j, k) = mat(j, k);
      }
      result(i * numRows + j, numCols) = i+1;
    }
  }
  
  return result;
}


// [[Rcpp::export]]
int findIntegerIndex(IntegerVector vec, int target) {
  for (int i = 0; i < vec.length(); i++) {
    if (vec[i] == target) {
      return i;  // Return the index of the matching element
    }
  }
  
  // If no match is found, return -1 or any other appropriate value
  return -1;
}


// [[Rcpp::export]]
int findCharacterIndex(CharacterVector vec, String target) {
  for (int i = 0; i < vec.length(); i++) {
    if (vec[i] == target) {
      return i;  // Return the index of the matching element
    }
  }
  // If no match is found, return -1 or any other appropriate value
  return -1;
}

// [[Rcpp::export]]
void printMatrix(const NumericMatrix& mat) {
  int nrows = mat.nrow();
  int ncols = mat.ncol();
  
  for (int i = 0; i < nrows; i++) {
    for (int j = 0; j < ncols; j++) {
      Rcpp::Rcout << mat(i, j) << " ";
    }
    Rcpp::Rcout << "\n";
  }
}

// [[Rcpp::export]]
NumericVector createHeightClasses(double Hmax, int heightClassesN) {
  NumericVector intervals(heightClassesN+1);
  
  double intervalSize = Hmax / heightClassesN;
  
  intervals[0] = 0;
  for (int i = 1; i < heightClassesN; i++) {
    intervals[i] = (i) * intervalSize;
  }
  
  intervals[intervals.size()-1] = 99999;
  
  return intervals;
}

// [[Rcpp::export]]
int getHeightClass(double H, NumericVector intervals) {
  H = H-0.00001;
  int intervalIndex = -1;
  int numIntervals = intervals.size() - 1;
  
  for (int i = 0; i < numIntervals; i++) {
    if (H >= intervals[i] && H < intervals[i + 1]) {
      intervalIndex = i;
      break;
    }
  }
  
  return intervalIndex;
}

// =============================================================================
// Functions used in Regeneration / Mortality / Growth and the StateClass
// =============================================================================

// foliage function from ForClim used for LAI calculation
double gFolA_f(double kA1, double kA2, double kC1, double kC2, double dbh, double areaHectar) {
  double gFolW = kA1 * exp(kA2 * log(dbh)) * kC1;
  double gFolA = kC2 * gFolW / kC1;
  gFolA = gFolA / areaHectar / 10000; // areaHectar is the patch size in ha
  return gFolA / 2;
}

// [[Rcpp::export]]
double envF(double envMean, double envSD, double env, double widthPar) {
  double maxDens = R::dnorm(envMean, envMean, envSD, 0);
  double pDens = R::dnorm(env, envMean, envSD, 0);
  double envCond = pDens / maxDens;
  envCond = 1 - pow(1-envCond, widthPar);
  return envCond;
}

// [[Rcpp::export]]
double envF2(double envMean, double envSD, double env, double widthPar) {
  double maxDens = R::dnorm(envMean, envMean, envSD, 0);
  double pDens = R::dnorm(env, envMean, envSD, 0);
  if(envMean < env){
    pDens = maxDens;
  }
  double envCond = pDens / maxDens;
  envCond = 1 - pow(1-envCond, widthPar);
  return pDens / maxDens;
}

// [[Rcpp::export]]
double AL_F(double LAI) {
  return exp(-0.25 * LAI);
}

// [[Rcpp::export]]
double shadeF(double LAI, double shadeMean, double shadeSD) {
  double AL = AL_F(LAI);
  double minP = R::pnorm(0, shadeMean, shadeSD, 0, 0);
  double shadeCond = R::pnorm(AL, shadeMean, shadeSD, 0, 0);
  shadeCond = (minP - shadeCond) / minP;
  return shadeCond;
}

// [[Rcpp::export]]
NumericMatrix updateH(NumericMatrix cohorts, List pars, List speciesPars) {
  NumericVector kLaVec = as<NumericVector>(speciesPars["kLa"]); // 1
  NumericVector kHMaxVec = as<NumericVector>(speciesPars["kHMax"]) * 100; // 1
  
  int numCohorts = cohorts.nrow();
  //int numVars = cohorts.ncol();
  
  for (int i = 0; i < numCohorts; i++) {
    double D = cohorts(i, 3);
    int ispID = cohorts(i, 1);
    
    double kLa = kLaVec[ispID];
    double kHMax = kHMaxVec[ispID];
    
    double kB1 = 137;
    double kE1 = 14 * (kLa / 3 + 3) + 13;
    double kSMin = 1.3 * (kLa / 3 + 3) + 39.5;
    double kSIn = kSMin + 0.75 * kE1;
    
    double H = kB1 + (kHMax - kB1) * (1 - exp(-kSIn * D / (kHMax - kB1)));
    
    cohorts(i, 2) = H;
  }
  
  return cohorts;
}

// =============================================================================
// StateClass refers to state variables that are relevant for the output Matrix
// and applying growth, regeneration, and mortality functions
// =============================================================================
class stateClass {
private:
  NumericMatrix data; // a matrix that stores the state of different variables for each species
  NumericVector laiVec;
  NumericVector heightClassesVec;
  
public:
  stateClass(int nrow, int ncol) : data(nrow, ncol) {
  }
  
  void updateLAI(NumericMatrix cohorts, List speciesPars, double areaHectar, IntegerVector actualSpecies, CharacterVector outVars, List pars) {
    NumericVector kA1Vec = as<NumericVector>(speciesPars["kA1"]);
    NumericVector kA2Vec = as<NumericVector>(speciesPars["kA2"]);
    NumericVector kC1Vec = as<NumericVector>(speciesPars["kC1"]);
    NumericVector kC2Vec = as<NumericVector>(speciesPars["kC2"]);
    
    int laiIDX = findCharacterIndex(outVars, "lai");
    int dbhIDX = findCharacterIndex(outVars, "dbh");
    int nTrsIDX = findCharacterIndex(outVars, "nTrs");
    int baIDX = findCharacterIndex(outVars, "ba");
    
    IntegerVector nTrsVec(actualSpecies.length());
    nTrsVec.fill(0);
    
    for (int i = 0; i < data.nrow(); i++) {
      data(i, laiIDX) = 0;
      data(i, dbhIDX) = 0;
      data(i, nTrsIDX) = 0;
      data(i, baIDX) = 0;
    }
    
    double maxH = 0;
    for (int i = 0; i < cohorts.nrow(); i++) {
      double H = cohorts(i, 2);
      if (H > maxH) maxH = H;
    }
    heightClassesVec = createHeightClasses(maxH, pars["heightClassesN"]);
    
    NumericVector laiVecTemp(heightClassesVec.size());
    for (int i = 0; i < cohorts.nrow(); i++) {
      int ispID = cohorts(i, 1);
      int spIDX = findIntegerIndex(actualSpecies, ispID);
      double dbh = cohorts(i, 3);
      int nTrs = cohorts(i, 0);
      
      double gFolA = gFolA_f(kA1Vec[ispID], kA2Vec[ispID], kC1Vec[ispID], kC2Vec[ispID], dbh, areaHectar);
      
      double H = cohorts(i, 2);
      int heightClass = getHeightClass(H, heightClassesVec);
      laiVecTemp[heightClass] += gFolA * nTrs;
      
      if (laiIDX != -1) data(spIDX, laiIDX) += gFolA * nTrs;
      if (dbhIDX != -1) data(spIDX, dbhIDX) += dbh * nTrs;
      if (baIDX != -1) data(spIDX, baIDX) += (nTrs * 3.14159 / 4 * dbh * dbh) / (areaHectar * 10000);
      
      if (nTrsIDX != -1) {
        data(spIDX, nTrsIDX) += nTrs;
      } else {
        nTrsVec[nTrsIDX] += nTrs;
      }
    }
    
    if (dbhIDX != -1) {
      for (int i = 0; i < data.nrow(); i++) {
        if (nTrsIDX != -1) {
          if (data(i, nTrsIDX) > 0) {
            data(i, dbhIDX) = data(i, dbhIDX) / data(i, nTrsIDX);
            data(i, nTrsIDX) = data(i, nTrsIDX) / areaHectar;
          } else {
            data(i, dbhIDX) = 0;
          }
        } else {
          if (data(i, nTrsIDX) > 0) {
            data(i, dbhIDX) = data(i, dbhIDX) / nTrsVec[i];
            nTrsVec[i] = nTrsVec[i] / areaHectar;
          } else {
            data(i, dbhIDX) = 0;
          }
        }
      }
    }
    laiVec = revCumSum(laiVecTemp);
  }
  
  NumericMatrix getMatrix() {
    return data;
  }
  
  double getLAI(int heightClass) {
    return laiVec[heightClass];
  }
  
  NumericVector getLAIVec() {
    return laiVec;
  }
  
  NumericVector getHeightClassesVec() {
    return heightClassesVec;
  }
};


// =============================================================================
// Regeneration
// =============================================================================

// [[Rcpp::export]]
NumericMatrix regeneration_f(NumericMatrix cohorts, double LAI, double env, List pars, List speciesPars) {
  IntegerVector actualSpecies = pars["actualSpecies"];
  NumericVector regP(actualSpecies.size());
  double areaHectar = pars["areaHectar"];
  double baseReg = pars["baseReg"];
  baseReg = baseReg * areaHectar;
  double baseRegP = pars["baseRegP"];
  double nicheWidth = pars["nicheWidth"];
  double regShadeEff = pars["regShadeEff"];
  double regEnvEff = pars["regEnvEff"];
  
  if (R::rbinom(1, baseRegP) == 1) {
    for (int i = 0; i < actualSpecies.size(); i++) {
      int ispID = actualSpecies[i];
      
      double shadeMean = as<NumericVector>(speciesPars["kLy"])[ispID];
      double shadeCond = shadeF(LAI, shadeMean, 0.1);
      double regShadeShape = 3;
      shadeCond = regShadeEff * pow(shadeCond, regShadeShape);
      
      double envMean = as<NumericVector>(speciesPars["env"])[ispID];
      double envCond = envF(envMean, 0.1, env, nicheWidth);
      double regEnvShape = 3;
      envCond = regEnvEff * pow(envCond, regEnvShape);
      
      double regCond = envCond * shadeCond;
      
      if (regCond != 0) {
        regP[i] = std::min(regCond, 1.0);
      } else {
        regP[i] = 0;
      }
    }
    
    NumericVector regTrsMeans;
    if (sum(regP) == 0) {
      regTrsMeans = regP;
    } else {
      regTrsMeans = regP / sum(regP) * baseReg;
    }
    int newID = 0;
    if (cohorts.nrow() > 0) {
      newID = cohorts(cohorts.nrow() - 1, 0);
    }
    
    int numActualSpecies = actualSpecies.size();
    NumericMatrix newCohorts(numActualSpecies, 4);
    
    int newRow = 0;
    for (int i = 0; i < numActualSpecies; i++) {
      int ispID = actualSpecies[i];
      if (regTrsMeans[i] > 0) {
        int newTrs = R::rpois(regTrsMeans[i]);
        if (newTrs > 0) {
          newID++;
          newCohorts(newRow, 0) = newID;
          newCohorts(newRow, 1) = ispID;
          newCohorts(newRow, 2) = newTrs;
          newCohorts(newRow, 3) = 1;
          newRow++;
        }
      }
    }
    
    if (newRow > 0) {
      NumericMatrix resizedCohorts(cohorts.nrow() + newRow, cohorts.ncol());
      for (int i = 0; i < cohorts.nrow(); i++) {
        for (int j = 0; j < cohorts.ncol(); j++) {
          resizedCohorts(i, j) = cohorts(i, j);
        }
      }
      for (int i = 0; i < newRow; i++) {
        for (int j = 0; j < cohorts.ncol(); j++) {
          resizedCohorts(cohorts.nrow() + i, j) = newCohorts(i, j);
        }
      }
      cohorts = resizedCohorts;
    }
  }
  
  return cohorts;
}


// =============================================================================
// Growth
// =============================================================================

// [[Rcpp::export]]
NumericMatrix growth_f(NumericMatrix cohorts, NumericVector laiVec, double env, List pars, List speciesPars, NumericVector heightClassesVec) {
  NumericVector kLaVec = as<NumericVector>(speciesPars["kLa"]);
  NumericVector kHMaxVec = as<NumericVector>(speciesPars["kHMax"])*100;
  NumericVector kGVec = as<NumericVector>(speciesPars["kG"]);
  NumericVector envMeanVec = as<NumericVector>(speciesPars["env"]);  
  
  int numCohorts = cohorts.nrow();
  //int numVars = cohorts.ncol();
  
  for (int i = 0; i < numCohorts; i++) {
    double D = cohorts(i, 3);
    int ispID = cohorts(i, 1);
    double H = cohorts(i, 2);
    int heightClass = getHeightClass(H, heightClassesVec);
    double LAI = laiVec[heightClass];
    
    double nicheWidth = pars["nicheWidth"];
    double kLa = kLaVec[ispID];
    double kHMax = kHMaxVec[ispID];
    double kG = kGVec[ispID];
    double envMean = envMeanVec[ispID];
    
    double kB1 = 137;
    double kE1 = 14 * (kLa / 3 + 3) + 13;
    double kSMin = 1.3 * (kLa / 3 + 3) + 39.5;
    
    double AL = AL_F(LAI);
    double gS = kSMin + kE1 * (1.0 - AL);
    double gFun = std::max(0.0, gS * (1 - (H - kB1) / (kHMax - kB1)));
    
    double gGRF = envF(envMean, 0.1, env, nicheWidth);
    double gGRF2 = 1 - pow(1- gGRF,4);
    double gRateD = std::max(0.0, gGRF2 * kG * D * (1 - H / kHMax) / (2 * H + gFun * D));
    
    cohorts(i, 3) = D + gRateD;
  }
  
  return cohorts;
}



// =============================================================================
// Mortality
// =============================================================================
// [[Rcpp::export]]
NumericMatrix mortality_f(NumericMatrix cohorts, NumericVector laiVec, double env, double tDist, List pars, List speciesPars, NumericVector heightClassesVec) {
  NumericVector kDMaxVec = as<NumericVector>(speciesPars["kDMax"]);
  NumericVector shadeMeanVec = as<NumericVector>(speciesPars["kLy"]);
  NumericVector envMeanVec = as<NumericVector>(speciesPars["env"]);
  
  int numCohorts = cohorts.nrow();
  int numVars = cohorts.ncol();
  
  NumericMatrix aliveCohorts(numCohorts, numVars);
  int aliveCount = 0;
  
  for (int i = 0; i < numCohorts; i++) {
    double D = cohorts(i, 3);
    int ispID = cohorts(i, 1);
    double H = cohorts(i, 2);
    int heightClass = getHeightClass(H, heightClassesVec);
    double LAI = laiVec[heightClass];
    
    double kAlpha = pars["bgMort"];
    double nicheWidth = pars["nicheWidth"];
    double mortEnvEff = pars["mortEnvEff"];
    double mortShadeEff = pars["mortShadeEff"];
    
    double kDMax = kDMaxVec[ispID];
    double shadeMean = shadeMeanVec[ispID];
    double envMean = envMeanVec[ispID];
    
    // shade dependent mortality
    double shadeCond = shadeF(LAI, shadeMean, 0.1);
    double mortShadeShape = 3;
    double gPShade = mortShadeEff * pow(1 - shadeCond, mortShadeShape);
    
    // environment dependent mortality
    double envCond = envF(envMean, 0.1, env, nicheWidth);
    double mortEnvShape = 3;
    double gPEnv = mortEnvEff * pow(1 - envCond, mortEnvShape);
    
    // size dependent mortality
    double gPSize = 0.1 * pow(D / kDMax, kAlpha);
    
    double mortP = std::min(gPShade + gPEnv + gPSize + tDist, 1.0);
    
    int nTrs = cohorts(i, 0);
    if (nTrs > 0) {
      int nTrsDead = Rcpp::sum(Rcpp::rbinom(nTrs, 1, mortP));
      cohorts(i, 0) = nTrs - nTrsDead;
      
      // Copy the cohort to aliveCohorts if it's still alive
      if (cohorts(i, 0) > 0) {
        for (int j = 0; j < numVars; j++) {
          aliveCohorts(aliveCount, j) = cohorts(i, j);
        }
        aliveCount++;
      }
    }
  }
  
  // Resize the aliveCohorts matrix to the actual number of alive cohorts
  if (aliveCount > 0) {
    aliveCohorts = aliveCohorts(Range(0, aliveCount - 1), Range(0, numVars - 1));
  } else {
    aliveCohorts = NumericMatrix(0, numVars);  // Return an empty matrix when no cohorts are alive
  }
  
  return aliveCohorts;
}


// =============================================================================
// Model Function
// =============================================================================

// [[Rcpp::export]]
NumericMatrix runModel(List pars, List speciesPars) {
  startTimer();
  int patchesN = pars["patchesN"];
  double env = pars["env"];
  int timesteps = pars["timesteps"];
  double areaHectar = pars["areaHectar"];
  double distInt = pars["distInt"];
  //Rcpp::Rcout << "start reading initPop worked\n";
  NumericMatrix initCohorts;
  if (!Rf_isNull(pars["initPop"])) {
    initCohorts = as<NumericMatrix>(pars["initPop"]);
  } else {
    // Empty matrix when initPop is NULL
    initCohorts = NumericMatrix(0, 4);
  }
  //Rcpp::Rcout << "finished reading initPop worked\n";
  IntegerVector actualSpeciesVec = pars["actualSpecies"];
  DataFrame out_df;
  
  CharacterVector outVars = pars["outVars"];
  stateClass stateVars(actualSpeciesVec.length(), outVars.length());
  stateClass stateVarsTot(actualSpeciesVec.length(), outVars.length());
  NumericMatrix outMat(actualSpeciesVec.length(), outVars.length());
  
  // Add values to the specified column
  for (int i = 0; i < stateVars.getMatrix().nrow(); i++) {
    stateVars.getMatrix()(i, 0) = actualSpeciesVec[i];
    stateVarsTot.getMatrix()(i, 0) = actualSpeciesVec[i];
    outMat(i, 0) = actualSpeciesVec[i];
  }
  
  outMat = repMat(outMat, timesteps);
  outMat = repMat(outMat, patchesN);
  
  IntegerMatrix distMat(timesteps, patchesN);
  for(int t = 1; t <= timesteps; t++){
    int tDist = R::rbinom(1, pars["distP"]);
    int distPatchesN = roundToDecimal(patchesN*distInt,0);
    IntegerVector distPatches = sample(patchesN, distPatchesN); // given
    for(int p = 1; p <= patchesN; p++){
      if(findIntegerIndex(distPatches, p) != -1){
        distMat(t, p) = tDist;
      }else {
        distMat(t, p) = 0;
      }
    }
  }
  
  for (int p = 1; p <= patchesN; p++) {
    NumericMatrix cohorts = clone(initCohorts);
    
    for (int t = 1; t <= timesteps; t++) {
      // Rcpp::Rcout << "######################### t = "<<t  << "\n";
      // startTimer();
      cohorts = updateH(cohorts, pars, speciesPars);
      //Rcpp::Rcout << "updateH worked\n";
      stateVars.updateLAI(cohorts, speciesPars, areaHectar, actualSpeciesVec, outVars,  pars);
      //Rcpp::Rcout << "updateLAI worked\n";
      
      for (int i = 0; i < stateVars.getMatrix().nrow(); i++) {
        for(int j = 1; j < stateVars.getMatrix().ncol(); j++){
          outMat( i + (t-1)*actualSpeciesVec.length() + (p-1)*timesteps*actualSpeciesVec.length(), j) = stateVars.getMatrix()(i, j); // FIXME
        }
      }
      //Rcpp::Rcout << "statevarsupdate worked\n";
      
      //int tDist = R::rbinom(1, pars["distP"]);
      int tDist = distMat(t, p);
      cohorts = growth_f(cohorts, stateVars.getLAIVec(), env, pars, speciesPars, stateVars.getHeightClassesVec());
      //Rcpp::Rcout << "growth worked\n";
      cohorts = mortality_f(cohorts, stateVars.getLAIVec(), env, tDist, pars, speciesPars, stateVars.getHeightClassesVec());
      //Rcpp::Rcout << "mortality worked\n";
      cohorts = regeneration_f(cohorts, stateVars.getLAI(0), env, pars, speciesPars);
      //Rcpp::Rcout << "regeneration worked\n";
    }
  }
  printElapsedTime();
  return outMat;
}

