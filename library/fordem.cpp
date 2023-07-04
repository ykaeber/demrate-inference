#include <Rcpp.h>
#include <algorithm>
#include <chrono>
using namespace Rcpp;

// foliage function from ForClim
double gFolA_f(double kA1, double kA2, double kC1, double kC2, double dbh, double areaHectar) {
  double gFolW = kA1 * exp(kA2 * log(dbh)) * kC1;
  double gFolA = kC2 * gFolW / kC1;
  gFolA = gFolA / areaHectar / 10000; // 0.01 is the patch size in ha
  return gFolA / 2;
}


// [[Rcpp::export]]
int findValueIndex(IntegerVector vec, int target) {
  for (int i = 0; i < vec.length(); i++) {
    if (vec[i] == target) {
      return i;  // Return the index of the matching element
    }
  }
  
  // If no match is found, return -1 or any other appropriate value
  return -1;
}

// [[Rcpp::export]]
void printMatrix(const NumericMatrix& mat, int t, int p) {
  int nrows = mat.nrow();
  int ncols = mat.ncol();
  
  for (int i = 0; i < nrows; i++) {
    for (int j = 0; j < ncols; j++) {
      Rcpp::Rcout << mat(i, j) << " ";
    }
    Rcpp::Rcout << t << " " << p <<"\n";
  }
}


class stateClass {
private:
  NumericMatrix data;
  NumericVector rowSums;
  NumericVector colSums;
  
public:
  stateClass(int nrow, int ncol) : data(nrow, ncol), rowSums(nrow), colSums(ncol) {
    calculateSums();
  }
  
  void update(int row, int col, double value) {
    data(row, col) = value;
    calculateSums();
  }
  
  void updateLAI(List cohorts, List speciesPars, double areaHectar, IntegerVector actualSpecies){
    NumericVector kA1Vec = as<NumericVector>(speciesPars["kA1"]);
    NumericVector kA2Vec = as<NumericVector>(speciesPars["kA2"]);
    NumericVector kC1Vec = as<NumericVector>(speciesPars["kC1"]);
    NumericVector kC2Vec = as<NumericVector>(speciesPars["kC2"]);
    
    double lai = 0;
    for (int i = 0; i < data.nrow(); i++) {
      data(i, 1) = 0; 
    }
    for (int i = 0; i < cohorts.size(); i++) {
      List iCohort = cohorts[i];
      int iSpID = iCohort["spID"];
      double dbh = iCohort["dbh"];
      int nTrs = iCohort["nTrs"];
      double gFolA = gFolA_f(kA1Vec[iSpID], kA2Vec[iSpID], kC1Vec[iSpID], kC2Vec[iSpID], dbh, areaHectar);
      lai += gFolA * nTrs;
      data(findValueIndex(actualSpecies, iSpID),1) += lai;
    }
    calculateSums();
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
DataFrame calculateMeansAndSums(List cohorts, int t, int p, List pars) {
  if(cohorts.size() == 0){
    List cohorts = List::create(Named("spID") = 0,
                               Named("cohortID") = 0,
                               Named("nTrs") = 0,
                               Named("dbh") = 0);
  }
    
    double areaHectar = pars["areaHectar"];
  
    std::unordered_map<int, double> nTrsMean;
    std::unordered_map<int, double> nTrsSum;
    std::unordered_map<int, double> dbhMean;
    std::unordered_map<int, double> dbhSum;
    std::unordered_map<int, double> ba;
    std::unordered_map<int, int> count;
  
    int numCohorts = cohorts.size();
    for (int i = 0; i < numCohorts; i++) {
      List cohort = cohorts[i];
      int spID = cohort["spID"];
      double nTrs = cohort["nTrs"];
      double dbh = cohort["dbh"];
  
      nTrsMean[spID] += nTrs;
      nTrsSum[spID] += nTrs;
      dbhMean[spID] += dbh*nTrs;
      dbhSum[spID] += dbh;
      //((dbh/100)/2)^2
      ba[spID] += nTrs*3.14159*std::pow(((dbh/100)/2), 2);
      count[spID]++;
    }
  
    int numSpecies = nTrsMean.size();
    IntegerVector spIDVec(numSpecies);
    NumericVector nTrsMeanVec(numSpecies);
    NumericVector nTrsSumVec(numSpecies);
    NumericVector dbhMeanVec(numSpecies);
    NumericVector dbhSumVec(numSpecies);
    NumericVector baHectarVec(numSpecies);
    NumericVector baTotVec(numSpecies);
  
    std::unordered_map<int, double>::iterator it;
    int index = 0;
    for (it = nTrsMean.begin(); it != nTrsMean.end(); ++it) {
      int spID = it->first;
      double nTrsMeanValue = it->second / count[spID];
      double nTrsSumValue = nTrsSum[spID];
      double dbhMeanValue = dbhMean[spID] / nTrsSum[spID];
      double dbhSumValue = dbhSum[spID];
      double baHectar =  ba[spID] * 1/ areaHectar;
      double baTot = ba[spID];
  
      if (std::isnan(dbhMeanValue)) dbhMeanValue = 0;
  
      spIDVec[index] = spID;
      nTrsMeanVec[index] = roundToDecimal(nTrsMeanValue, 2);
      nTrsSumVec[index] = roundToDecimal(nTrsSumValue, 2);
      dbhMeanVec[index] = roundToDecimal(dbhMeanValue, 2);
      dbhSumVec[index] = roundToDecimal(dbhSumValue, 2);
      baHectarVec[index] = roundToDecimal(baHectar, 4);
      baTotVec[index] = roundToDecimal(baTot, 4);
  
      index++;
    }
  
    DataFrame result = DataFrame::create(Named("spID") = spIDVec,
                                         Named("nTrsMean") = nTrsMeanVec,
                                         Named("nTrsSum") = nTrsSumVec,
                                         Named("dbhMean") = dbhMeanVec,
                                         Named("dbhSum") = dbhSumVec,
                                         Named("baHectar") = baHectarVec,
                                         Named("baTot") = baTotVec,
                                         Named("t") = t,
                                         Named("p") = p
                                           );
    
    return result;
}


// [[Rcpp::export]]
void printDataFrame(DataFrame df, int t, int p) {
  // Get the number of rows and columns in the data frame
  int numRows = df.nrows();
  int numCols = df.size();
  
  // Convert column names to std::vector<std::string>
  std::vector<std::string> colNames = as<std::vector<std::string>>(df.names());
  
  // Print column names
  if (t == 1 && p == 1){
    for (int i = 0; i < numCols; i++) {
      Rcout << colNames[i] << "\t";
    }
    Rcout << std::endl;
  }  
  // Iterate over the rows
  for (int j = 0; j < numRows; j++) {
    // Iterate over the columns
    for (int i = 0; i < numCols; i++) {
      // Get the i-th column from the data frame
      SEXP col = df[i];
      
      // Get the j-th value in the column
      switch (TYPEOF(col)) {
      case INTSXP:
        Rcout << INTEGER(col)[j] << "\t";
        break;
      case REALSXP:
        Rcout << REAL(col)[j] << "\t";
        break;
      case STRSXP:
        Rcout << CHAR(STRING_ELT(col, j)) << "\t";
        break;
      default:
        Rcout << "NA\t";
      break;
      }
    }
    Rcout << std::endl;
  }
}

// [[Rcpp::export]]
void runModel(List pars, List speciesPars) {
  int patchesN = pars["patchesN"];
  double env = pars["env"];
  int timesteps = pars["timesteps"];
  double areaHectar = pars["areaHectar"];
  IntegerVector actualSpeciesVec = pars["actualSpecies"];
  DataFrame out_df;
  
  CharacterVector colNames = CharacterVector::create("SpID", "lai", "nTrs", "dbh");
  stateClass stateVars(actualSpeciesVec.length(), colNames.length());
  
  // Add values to the specified column
  for (int i = 0; i < stateVars.getMatrix().nrow(); i++) {
    stateVars.getMatrix()(i, 0) += actualSpeciesVec[i];
    for(int j = 1; j < stateVars.getMatrix().nrow(); j++){
      stateVars.getMatrix()(i, j) += 0;
    }
  }
  
  for (int p = 1; p <= patchesN; p++) {
  List cohorts = pars["initPop"];
  
    for (int t = 1; t <= timesteps; t++) {
      // startTimer();
      stateVars.updateLAI(cohorts, speciesPars, areaHectar, actualSpeciesVec);
      
      double LAI = stateVars.getColSums()[1];
      // if(t==timesteps && p == patchesN) printElapsedTime("finished LAI");
      int tDist = R::rbinom(1, pars["distP"]);
      //double env = 0.5;
      
      // startTimer();
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
      printMatrix(stateVars.getMatrix(), t, p);
      //startTimer();
      //printDataFrame(out_df, p, t);
      //printElapsedTime("finished printDataFrame");
    }
  }
}
