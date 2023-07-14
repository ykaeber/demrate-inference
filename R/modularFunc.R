library(Rcpp)

sourceCpp("library/modularFunc.cpp")

performCalculations(1, 3, "defualt", "default", "default")

fun1 = function(x1,x2) x1-x1 + x2
fun3 = function(x) x + 20
performCalculations(1, 3, "defualt", "default", fun3)
performCalculations(1, 3, fun1, "default", fun3)

