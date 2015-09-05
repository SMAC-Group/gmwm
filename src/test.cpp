#include <RcppArmadillo.h>
using namespace Rcpp;


// [[Rcpp::export]]
void stringTest (std::vector <std::vector <std::string> >  x) {
  std::cout << "Test";
}


/*** R
stringTest(list(c("Hi","hi"),c("Rawr!","Poo")))
*/
