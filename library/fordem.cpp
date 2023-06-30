#include <Rcpp.h>
#include <algorithm>
using namespace Rcpp;

List createSpeciesPars(Rcpp::DataFrame df) {
  // Extract the columns from the data frame
  IntegerVector spID = df["spID"];
  CharacterVector species = df["species"];
  IntegerVector kHMax = df["kHMax"];
  IntegerVector kDMax = df["kDMax"];
  IntegerVector kG = df["kG"];
  IntegerVector kDDMin = df["kDDMin"];
  NumericVector kLy = df["kLy"];
  IntegerVector kLa = df["kLa"];
  NumericVector kA1 = df["kA1"];
  NumericVector kA2 = df["kA2"];
  NumericVector kC1 = df["kC1"];
  IntegerVector kC2 = df["kC2"];
  
  // Create the Rcpp object
  List output;
  
  output["spID"] = spID;
  output["species"] = species;
  output["kHMax"] = kHMax;
  output["kDMax"] = kDMax;
  output["kG"] = kG;
  output["kDDMin"] = kDDMin;
  output["kLy"] = kLy;
  output["kLa"] = kLa;
  output["kA1"] = kA1;
  output["kA2"] = kA2;
  output["kC1"] = kC1;
  output["kC2"] = kC2;
  
  return output;
}


// foliage function from ForClim
double gFolA_f(double kA1, double kA2, double kC1, double kC2, double dbh, double areaHectar = 0.01) {
  double gFolW = kA1 * exp(kA2 * log(dbh)) * kC1;
  double gFolA = kC2 * gFolW / kC1;
  gFolA = gFolA / areaHectar / 10000; // 0.01 is the patch size in ha
  return gFolA / 2;
}

// function for calculating the lai from the current cohorts
// [[Rcpp::export]]
double lai_f(List cohorts, List speciesPars) {
  double lai = 0;
  for (int i = 0; i < cohorts.size(); i++) {
    List iCohort = cohorts[i];
    int iSpID = iCohort["spID"];
    double kA1 = as<NumericVector>(speciesPars["kA1"])[iSpID];
    double kA2 = as<NumericVector>(speciesPars["kA2"])[iSpID];
    double kC1 = as<NumericVector>(speciesPars["kC1"])[iSpID];
    double kC2 = as<NumericVector>(speciesPars["kC2"])[iSpID];
    double dbh = iCohort["dbh"];
    int nTrs = iCohort["nTrs"];
    double gFolA = gFolA_f(kA1, kA2, kC1, kC2, dbh);
    lai += gFolA * nTrs;
  }
  return lai;
}

double envF(double envMean, double envSD, double env) {
  double maxDens = R::dnorm(envMean, envMean, envSD, 0);
  double pDens = R::dnorm(env, envMean, envSD, 0);
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
  
  for (int i = 0; i < actualSpecies.size(); i++) {
    int iSpID = actualSpecies[i];
    double shadeMean = as<NumericVector>(speciesPars["kLy"])[iSpID];
    double shadeCond = shadeF(LAI, shadeMean, 0.1);
    
    double envMean = as<NumericVector>(speciesPars["kDDMin"])[iSpID]/1100;
    double envCond = envF(envMean, 0.1, env);
    
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
  
  return cohorts;
}


// [[Rcpp::export]]
List growth_f(List cohorts, double LAI, double env, List pars, List speciesPars) {
  int num_cohorts = cohorts.size();
  List aliveCohorts;
  for (int i = 0; i < num_cohorts; i++) {
    List iCohort = cohorts[i];
    int nTrs = iCohort["nTrs"];
    if (nTrs > 0) {
      aliveCohorts.push_back(i);
    }
  }
  
  int num_alive = aliveCohorts.size();

  for (int i = 0; i < num_alive; i++) {
    int cohortID = aliveCohorts[i];
    List iCohort = cohorts[cohortID];
    double D = iCohort["dbh"];
    int iSpID = iCohort["spID"];
    
    double kLa = as<NumericVector>(speciesPars["kLa"])[iSpID];
    double kHMax = as<NumericVector>(speciesPars["kHMax"])[iSpID]*100;
    double kG = as<NumericVector>(speciesPars["kG"])[iSpID];
    double envMean = as<NumericVector>(speciesPars["kDDMin"])[iSpID] / 1100;
    
    double kB1 = 137;
    double kE1 = 14 * (kLa / 3 + 3) + 13;
    double kSMin = 1.3 * (kLa / 3 + 3) + 39.5;
    double kSIn = kSMin + 0.75 * kE1;
    
    double H = kB1 + (kHMax - kB1) * (1 - exp(-kSIn * D / (kHMax - kB1)));
    double AL = AL_F(LAI);
    double gS = kSMin + kE1 * (1.0 - AL);
    double gFun = std::max(0.0, gS * (1 - (H - kB1) / (kHMax - kB1)));
    
    double gGRF = envF(envMean, 0.1, env);
    double gRateD = std::max(0.0, gGRF * kG * D * (1 - H / kHMax) / (2 * H + gFun * D));
    
    iCohort["dbh"] = D + gRateD;
    cohorts[cohortID] = iCohort;
  }
  return cohorts;
}

// [[Rcpp::export]]
List mortality_f(List cohorts, double LAI, double env, double tDist, List pars, List speciesPars) {
  for (int i = 0; i < cohorts.size(); i++) {
    List iCohort = cohorts[i];
    double D = iCohort["dbh"];
    int iSpID = iCohort["spID"];
    double kAlpha = pars["bgMort"];
    
    double kDMax = as<NumericVector>(speciesPars["kDMax"])[iSpID];
    
    double shadeMean = as<NumericVector>(speciesPars["kLy"])[iSpID];
    double shadeCond = shadeF(LAI, shadeMean, 0.1);
    
    double envMean = as<NumericVector>(speciesPars["kDDMin"])[iSpID]/1100;
    double envCond = envF(envMean, 0.1, env);
    
    double gPSize = 0.1 * pow(D / kDMax, kAlpha);
    double mortP = std::min(1.0 - shadeCond + 1.0 - envCond + gPSize + tDist, 1.0);
    
    int nTrs = iCohort["nTrs"];
    if (nTrs > 0) {
      int nTrsDead = Rcpp::sum(Rcpp::rbinom(nTrs, 1, mortP));
      iCohort["nTrs"] = nTrs-nTrsDead;
    }
  }
  return cohorts;
}

// [[Rcpp::export]]
DataFrame calculateMeansAndSums(List cohorts, int t) {
  std::unordered_map<int, double> nTrsMean;
  std::unordered_map<int, double> nTrsSum;
  std::unordered_map<int, double> dbhMean;
  std::unordered_map<int, double> dbhSum;
  std::unordered_map<int, int> count;
  
  int numCohorts = cohorts.size();
  for (int i = 0; i < numCohorts; i++) {
    List cohort = cohorts[i];
    int spID = cohort["spID"];
    double nTrs = cohort["nTrs"];
    double dbh = cohort["dbh"];
    
    nTrsMean[spID] += nTrs;
    nTrsSum[spID] += nTrs;
    dbhMean[spID] += dbh;
    dbhSum[spID] += dbh;
    count[spID]++;
  }
  
  int numSpecies = nTrsMean.size();
  IntegerVector spIDVec(numSpecies);
  NumericVector nTrsMeanVec(numSpecies);
  NumericVector nTrsSumVec(numSpecies);
  NumericVector dbhMeanVec(numSpecies);
  NumericVector dbhSumVec(numSpecies);
  
  std::unordered_map<int, double>::iterator it;
  int index = 0;
  for (it = nTrsMean.begin(); it != nTrsMean.end(); ++it) {
    int spID = it->first;
    double nTrsMeanValue = it->second / count[spID];
    double nTrsSumValue = nTrsSum[spID];
    double dbhMeanValue = dbhMean[spID] / count[spID];
    double dbhSumValue = dbhSum[spID];
    
    spIDVec[index] = spID;
    nTrsMeanVec[index] = nTrsMeanValue;
    nTrsSumVec[index] = nTrsSumValue;
    dbhMeanVec[index] = dbhMeanValue;
    dbhSumVec[index] = dbhSumValue;
    
    index++;
  }
  
  DataFrame result = DataFrame::create(Named("spID") = spIDVec,
                                       Named("nTrsMean") = nTrsMeanVec,
                                       Named("nTrsSum") = nTrsSumVec,
                                       Named("dbhMean") = dbhMeanVec,
                                       Named("dbhSum") = dbhSumVec,
                                       Named("t") = t);
  
  return result;
}

// [[Rcpp::export]]
void printDataFrame(DataFrame df, int t) {
  // Get the number of rows and columns in the data frame
  int numRows = df.nrows();
  int numCols = df.size();
  
  // Convert column names to std::vector<std::string>
  std::vector<std::string> colNames = as<std::vector<std::string>>(df.names());
  
  // Print column names
  if (t == 1){
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
  List cohorts = pars["initPop"];
  int timesteps = pars["timesteps"];
  DataFrame out_df;
  
  for (int t = 1; t <= timesteps; t++) {
    double LAI = lai_f(cohorts, speciesPars);
    int tDist = R::rbinom(1, pars["distP"]);
    double env = 0.5;
    
    cohorts = regeneration_f(cohorts, LAI, env, pars, speciesPars);
    cohorts = growth_f(cohorts, LAI, env, pars, speciesPars);
    cohorts = mortality_f(cohorts, LAI, env, tDist, pars, speciesPars);
    
    out_df = calculateMeansAndSums(cohorts, t);
    printDataFrame(out_df, t);
  }
}
