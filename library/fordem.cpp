#include <Rcpp.h>
#include <algorithm>
#include <chrono>
using namespace Rcpp;

// foliage function from ForClim
double gFolA_f(double kA1, double kA2, double kC1, double kC2, double dbh, double areaHectar) {
  double gFolW = kA1 * exp(kA2 * log(dbh)) * kC1;
  double gFolA = kC2 * gFolW / kC1;
  gFolA = gFolA / areaHectar / 10000; // areaHectar is the patch size in ha
  return gFolA / 2;
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

class stateClass {
private:
  NumericMatrix data;
  NumericVector rowSums;
  NumericVector colSums;
  double LAI;
  
public:
  stateClass(int nrow, int ncol) : data(nrow, ncol), rowSums(nrow), colSums(ncol) {
    calculateSums();
  }
  
  void update(int row, int col, double value) {
    data(row, col) = value;
    calculateSums();
  }

  void updateLAI(List cohorts, List speciesPars, double areaHectar, IntegerVector actualSpecies, CharacterVector outVars){
    NumericVector kA1Vec = as<NumericVector>(speciesPars["kA1"]);
    NumericVector kA2Vec = as<NumericVector>(speciesPars["kA2"]);
    NumericVector kC1Vec = as<NumericVector>(speciesPars["kC1"]);
    NumericVector kC2Vec = as<NumericVector>(speciesPars["kC2"]);
    
    int laiIDX = findCharacterIndex(outVars, "lai");
    int dbhIDX = findCharacterIndex(outVars, "dbh");
    int nTrsIDX = findCharacterIndex(outVars, "nTrs");
    int baIDX = findCharacterIndex(outVars, "ba");
    
    IntegerVector countVec(actualSpecies.length());
    countVec.fill(0);
    
    IntegerVector nTrsVec(actualSpecies.length());
    nTrsVec.fill(0);
    
    for (int i = 0; i < data.nrow(); i++) {
      data(i, laiIDX) = 0; 
      data(i, dbhIDX) = 0; 
      data(i, nTrsIDX) = 0; 
      data(i, baIDX) = 0; 
    }
    
    for (int i = 0; i < cohorts.size(); i++) {
      List iCohort = cohorts[i];
      int iSpID = iCohort["spID"];
      int spIDX = findIntegerIndex(actualSpecies, iSpID);
      double dbh = iCohort["dbh"];
      int nTrs = iCohort["nTrs"];
      double gFolA = gFolA_f(kA1Vec[iSpID], kA2Vec[iSpID], kC1Vec[iSpID], kC2Vec[iSpID], dbh, areaHectar);
      
      if(laiIDX != -1) data(spIDX,laiIDX) += gFolA * nTrs;
      if(dbhIDX != -1) data(spIDX,dbhIDX) += dbh*nTrs;
      if(baIDX != -1) data(spIDX,baIDX) += (nTrs*3.14159/4*dbh*dbh)/(areaHectar*10000);
      // if(baIDX != -1) data(spIDX,baIDX) += nTrs*3.14159*std::pow(((dbh/100)/2), 2) * 1 / (areaHectar*10000);
      
      if(nTrsIDX != -1) {
        data(spIDX, nTrsIDX) += nTrs;
        }else{
          nTrsVec[nTrsIDX] += nTrs;
          }
        }
    if(dbhIDX != -1){
      for (int i = 0; i < data.nrow(); i++) {
        if(nTrsIDX != -1) {
          if(data(i, nTrsIDX) > 0){
            data(i, dbhIDX) = data(i, dbhIDX)/data(i, nTrsIDX);
          }else{
            data(i, dbhIDX) = 0;
          }
        }else{
          if(data(i, nTrsIDX) > 0){
            data(i, dbhIDX) = data(i, dbhIDX)/nTrsVec[i]; 
          }else{
            data(i, dbhIDX) = 0; 
          }
        }
      }
      }
    calculateSums();
    LAI = colSums[laiIDX];
    // Rcpp::Rcout << "====== in updateLAI ======\n";
    // Rcpp::Rcout << "" << "LAI = " << LAI << "\n";
  }
  
  void calculateSums() {
    for (int i = 0; i < data.nrow(); i++) {
      rowSums[i] = sum(data(i, _));
    }
    
    for (int j = 0; j < data.ncol(); j++) {
      colSums[j] = sum(data(_, j));
    }
  }
  
  NumericMatrix getMatrix() {
    return data;
  }
  
  NumericVector getRowSums() {
    return rowSums;
  }
  
  NumericVector getColSums() {
    return colSums;
  }
  double getLAI() {
    return LAI;
  }
};

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




double envF(double envMean, double envSD, double env) {
  double maxDens = R::dnorm(envMean, envMean, envSD, 0);
  double pDens = R::dnorm(env, envMean, envSD, 0);
  return pDens / maxDens;
}

double envF2(double envMean, double envSD, double env) {
  double maxDens = R::dnorm(envMean, envMean, envSD, 0);
  double pDens = R::dnorm(env, envMean, envSD, 0);
  if(envMean < env){
    pDens = maxDens;
  }
  return pDens / maxDens;
}

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
List regeneration_f(List cohorts, double LAI, double env, List pars, List speciesPars) {
  IntegerVector actualSpecies = pars["actualSpecies"];
  NumericVector regP(actualSpecies.size());
  double baseReg = pars["baseReg"];
  double baseRegP = pars["baseRegP"];
  
  if(R::rbinom(1, baseRegP) == 1){
    for (int i = 0; i < actualSpecies.size(); i++) {
      int iSpID = actualSpecies[i];
      double shadeMean = as<NumericVector>(speciesPars["kLy"])[iSpID];
      double shadeCond = shadeF(LAI, shadeMean, 0.1);
      
      double envMean = as<NumericVector>(speciesPars["kDDMin"])[iSpID]/1100;
      double envCond = envF2(envMean, 0.1, env);
      
      if (shadeCond * envCond != 0) {
        regP[i] = std::min(shadeCond + envCond, 1.0);
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
    if (cohorts.size() > 0) {
      List lastCohort = cohorts[cohorts.size()-1];
      int lastID =  lastCohort["cohortID"];
      newID = lastID;
    }
    
    for (int i = 0; i < actualSpecies.size(); i++) {
      int iSpID = actualSpecies[i];
      if (regTrsMeans[i] > 0) {
        int newTrs = R::rpois(regTrsMeans[i]);
        if (newTrs > 0) {
          newID++;
          cohorts.push_back(List::create(
              _["cohortID"] = newID,
              _["spID"] = iSpID,
              _["nTrs"] = newTrs,
              _["dbh"] = 1
          ));
        }
      }
    }
  }

  return cohorts;
}


// [[Rcpp::export]]
List growth_f(List cohorts, double LAI, double env, List pars, List speciesPars) {
  NumericVector kLaVec = as<NumericVector>(speciesPars["kLa"]);
  NumericVector kHMaxVec = as<NumericVector>(speciesPars["kHMax"])*100;
  NumericVector kGVec = as<NumericVector>(speciesPars["kG"]);
  NumericVector envMeanVec = as<NumericVector>(speciesPars["kDDMin"]) / 1100;  
  
  cohorts = clone(cohorts);
  for (int i = 0; i < cohorts.size(); i++) {
    List iCohort = cohorts[i];
    double D = iCohort["dbh"];
    int iSpID = iCohort["spID"];
    
    double kLa = kLaVec[iSpID];
    double kHMax = kHMaxVec[iSpID];
    double kG = kGVec[iSpID];
    double envMean = envMeanVec[iSpID];
    
    double kB1 = 137;
    double kE1 = 14 * (kLa / 3 + 3) + 13;
    double kSMin = 1.3 * (kLa / 3 + 3) + 39.5;
    double kSIn = kSMin + 0.75 * kE1;
    
    double H = kB1 + (kHMax - kB1) * (1 - exp(-kSIn * D / (kHMax - kB1)));
    double AL = AL_F(LAI);
    double gS = kSMin + kE1 * (1.0 - AL);
    double gFun = std::max(0.0, gS * (1 - (H - kB1) / (kHMax - kB1)));
    
    double gGRF = envF2(envMean, 0.1, env);
    double gRateD = std::max(0.0, gGRF * kG * D * (1 - H / kHMax) / (2 * H + gFun * D));
    
    iCohort["dbh"] = D + gRateD;
    cohorts[i] = iCohort;
  }
  return cohorts;
}

// [[Rcpp::export]]
List mortality_f(List cohorts, double LAI, double env, double tDist, List pars, List speciesPars) {
  NumericVector kDMaxVec = as<NumericVector>(speciesPars["kDMax"]);
  NumericVector shadeMeanVec = as<NumericVector>(speciesPars["kLy"]);
  NumericVector envMeanVec = as<NumericVector>(speciesPars["kDDMin"])/1100;

  cohorts = clone(cohorts);
  for (int i = 0; i < cohorts.size(); i++) {
    List iCohort = cohorts[i];
    double D = iCohort["dbh"];
    int iSpID = iCohort["spID"];
    double kAlpha = pars["bgMort"];
    
    double kDMax = kDMaxVec[iSpID];
    double shadeMean = shadeMeanVec[iSpID];
    double envMean = envMeanVec[iSpID];
    
    double shadeCond = shadeF(LAI, shadeMean, 0.1);
    
    double envCond = envF2(envMean, 0.1, env);
    
    double gPSize = 0.1 * pow(D / kDMax, kAlpha);
    
    
    double mortP = std::min(1.0 - shadeCond + 1.0 - envCond + gPSize + tDist, 1.0);
    
    // Rcpp::Rcout << "====== in Mortfun ======\n";
    // Rcpp::Rcout << "" << "LAI in Mort = " << LAI << "\n";
    // Rcpp::Rcout << "" << "shadeCond = " << shadeCond << "\n";
    // Rcpp::Rcout << "" << "envCond = " << envCond << "\n";
    // Rcpp::Rcout << "" << "gPSize = " << gPSize << "\n";
    // Rcpp::Rcout << "" << "tDist = " << tDist << "\n";
    
    int nTrs = iCohort["nTrs"];
    if (nTrs > 0) {
      int nTrsDead = Rcpp::sum(Rcpp::rbinom(nTrs, 1, mortP));
      iCohort["nTrs"] = nTrs-nTrsDead;
      cohorts[i] = iCohort;
    }
  }
  
  int num_cohorts = cohorts.size();
  List aliveCohorts;
  for (int i = 0; i < num_cohorts; i++) {
    List iCohort = cohorts[i];
    int nTrs = iCohort["nTrs"];
    if (nTrs > 0) {
      aliveCohorts.push_back(cohorts[i]);
    }
  }

  return aliveCohorts;
}

// [[Rcpp::export]]
NumericMatrix runModel(List pars, List speciesPars) {
  startTimer();
  int patchesN = pars["patchesN"];
  double env = pars["env"];
  int timesteps = pars["timesteps"];
  double areaHectar = pars["areaHectar"];
  List initCohorts = pars["initPop"];
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
  
  // start printing
  // for (int i = 0; i < outVars.length(); i++) {
  //   Rcpp::Rcout << outVars[i] << " ";
  // }
  //   Rcpp::Rcout << "t p\n";

  for (int p = 1; p <= patchesN; p++) {
    List cohorts = clone(initCohorts);
    for (int t = 1; t <= timesteps; t++) {
      // startTimer();
      stateVars.updateLAI(cohorts, speciesPars, areaHectar, actualSpeciesVec, outVars);
      
      for (int i = 0; i < stateVars.getMatrix().nrow(); i++) {
        for(int j = 1; j < stateVars.getMatrix().ncol(); j++){
          outMat( i + (t-1)*actualSpeciesVec.length() + (p-1)*timesteps*actualSpeciesVec.length(), j) = stateVars.getMatrix()(i, j); // FIXME
        }
      }
      
      // printMatrix(stateVars.getMatrix());
      // Rcpp::Rcout << "\n\n\n";
      
      double LAI = stateVars.getLAI();
      // Rcpp::Rcout << "====== in model Iteration ======\n";
      // Rcpp::Rcout << "" << "LAI in Mort = " << LAI << "\n";
      // if(t==timesteps && p == patchesN) printElapsedTime("finished LAI");
      int tDist = R::rbinom(1, pars["distP"]);
      //double env = 0.5;
      
      cohorts = regeneration_f(cohorts, LAI, env, pars, speciesPars);
      // if(t==timesteps && p == patchesN) printElapsedTime("finished regeneration / start growth");
      // startTimer();
      cohorts = growth_f(cohorts, LAI, env, pars, speciesPars);
      // if(t==timesteps && p == patchesN) printElapsedTime("finished growth / start mortality");
      // startTimer();
      cohorts = mortality_f(cohorts, LAI, env, tDist, pars, speciesPars);
      // if(t==timesteps && p == patchesN) printElapsedTime("finished mortality / start calculateMeansAndSums");
      
      // startTimer();
      //out_df = calculateMeansAndSums(cohorts, t, p, pars);
      // if(t==timesteps && p == patchesN) printElapsedTime("finished calculateMeansAndSums / start printDataFrame");

      // Rcpp::print(stateVars.getColSums());
      // printMatrix(stateVars.getMatrix(), t, p);
      //startTimer();
      //printDataFrame(out_df, p, t);
      //printElapsedTime("finished printDataFrame");
      
    }
  }
  printElapsedTime();
  //printMatrix(outMat);
  return outMat;
}
