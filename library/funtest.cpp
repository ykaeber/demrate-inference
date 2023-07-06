#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector createIntervals(double Hmax, int heightClassesN) {
  NumericVector intervals(heightClassesN+1);
  
  double intervalSize = Hmax / heightClassesN;
  
  intervals[0] = 0;
  for (int i = 1; i < heightClassesN+1; i++) {
    intervals[i] = (i) * intervalSize;
  }
  
  return intervals;
}


// [[Rcpp::export]]
int getIntervalIndex(double H, NumericVector intervals) {
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

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector reverseCumulativeSum(NumericVector vec) {
  int n = vec.size();
  NumericVector result(n);
  double sum = 0.0;
  
  for (int i = n - 1; i >= 0; i--) {
    sum += vec[i];
    result[i] = sum;
  }
  
  return result;
}
