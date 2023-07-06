#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
int testDist(int nTrs, double mortP){
  return Rcpp::sum(Rcpp::rbinom(nTrs, 1, mortP));
}
