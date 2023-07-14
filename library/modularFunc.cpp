#include <Rcpp.h>
using namespace Rcpp;

// Original functions
double originalFunc1(double x1, double x2) {
  // Perform calculations
  return x1 + x2;
}

double originalFunc2(double previousResult) {
  // Perform calculations based on previous result
  return previousResult * 2;
}

double originalFunc3(double previousResult) {
  // Perform final calculations based on previous result
  return previousResult + 10;
}

// Rcpp function using original or custom functions
// [[Rcpp::export]]
double performCalculations(double x1, double x2,
                           SEXP func1 = R_NilValue,
                           SEXP func2 = R_NilValue,
                           SEXP func3 = R_NilValue) {
  
  if (TYPEOF(func1) == STRSXP && Rf_length(func1) == 1 && std::string(CHAR(STRING_ELT(func1, 0))) == "default")
    func1 = R_NilValue;
  
  if (TYPEOF(func2) == STRSXP && Rf_length(func2) == 1 && std::string(CHAR(STRING_ELT(func2, 0))) == "default")
    func2 = R_NilValue;
  
  if (TYPEOF(func3) == STRSXP && Rf_length(func3) == 1 && std::string(CHAR(STRING_ELT(func3, 0))) == "default")
    func3 = R_NilValue;
  
  bool fun1Default = !Rf_isNull(func1) && Rf_isFunction(func1);
  bool fun2Default = !Rf_isNull(func1) && Rf_isFunction(func2);
  bool fun3Default = !Rf_isNull(func1) && Rf_isFunction(func3);
  
  double result;
  
  if (fun1Default)
    result = as<double>(Function(func1)(x1, x2));
  else
    result = originalFunc1(x1, x2);
  
  if (fun2Default)
    result = as<double>(Function(func2)(result));
  else
    result = originalFunc2(result);
  
  if (fun3Default)
    result = as<double>(Function(func3)(result));
  else
    result = originalFunc3(result);
  
  return result;
}
