#include <RcppArmadillo.h>

#include "ts_checks.h"

// Include polyroot for invertibility
#include "polyroot.h"

// Complex tools
#include "complex_tools.h"
using namespace Rcpp;


// [[Rcpp::export]]
double minroot(const arma::cx_vec& x){
  return min(
    // sqrt(x^2 + y^2)
    Mod_cpp(
      // Return roots
      do_polyroot_arma(x)
    )
  );
}

//' @title Check Invertibility Conditions
//' @description Checks the invertiveness of series of coefficients.
//' @param x A \code{cx_vec} that has a 1 appended before the coefficents. (e.g. c(1, x))
//' @return True (if outside unit circle) || False (if inside unit circle)
// [[Rcpp::export]]
bool invert_check(const arma::cx_vec& x){
  return minroot(x) > 1; // Outside the unit circle
}