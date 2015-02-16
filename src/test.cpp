#include <Rcpp.h>
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

// [[Rcpp::export]]
int timesTwo(int x) {
   return x * 2;
}

// [[Rcpp::export]]
int cpower(int x) {
  Rcout<<"verbose activates";
   return (x*x);
}


extern "C" SEXP twoWrapper (SEXP xs){
  int x = Rcpp::as<int>(xs);
  int two = timesTwo(x);
  return (Rcpp::wrap(two));
}
extern "C" SEXP cpowWrapper (SEXP xs){
  int x = Rcpp::as<int>(xs);
  int p = cpower(x);
  return (Rcpp::wrap(p));
}
extern "C" SEXP IRTppWrapper (SEXP model, SEXP dataset , SEXP dimensions , SEXP initvals , SEXP xs){
  int x = Rcpp::as<int>(xs);
  int p = cpower(x);
  return (Rcpp::wrap(p));
}