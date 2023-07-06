#include <Rcpp.h>
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

class MyMatrix {
private:
  NumericMatrix data;
  NumericVector rowSums;
  NumericVector colSums;
  
public:
  MyMatrix(int nrow, int ncol) : data(nrow, ncol), rowSums(nrow), colSums(ncol) {
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


