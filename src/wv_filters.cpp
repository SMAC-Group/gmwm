/* Copyright (C) 2014 - 2015  James Balamuta, Stephane Guerrier, Roberto Molinari
 *
 * This file is part of GMWM R Methods Package
 *
 * The `gmwm` R package is free software: you can redistribute it and/or modify it
 * under the terms of the Q Public License included within the packages source
 * as the LICENSE file.
 *
 * The `gmwm` R package is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * You should have received a copy of the Q Public License
 * along with `gmwm`.  If not, see <http://www.smac-group.com/licensing/>.
 * 
 */

#include <RcppArmadillo.h>
#include "wv_filters.h"
// We use reverse_vec
#include "armadillo_manipulations.h"
using namespace Rcpp;

//' @title Quadrature Mirror Filter
//' @description Calculate the series quadrature mirror filter (QMF). Requires a series of an even length.
//' @usage qmf(g, inverse)
//' @param g A \code{vector} that contains the filter constants.
//' @param inverse A \code{bool} that indicates whether the inverse quadrature mirror filter is computed. 
//' By default, the inverse quadrature mirror is computed.
//' @return A \code{vector} that contains either the forward QMF (evalute in order) or the inverse QMF (reverse order). 
//' @author JJB
//' @keywords internal
//' @examples
//' # Haar values
//' g = rep(1/sqrt(2),2)
//' qmf(g)
// [[Rcpp::export]]
arma::vec qmf(arma::vec g, bool inverse = true) {
  
  unsigned int L = g.n_elem;
  
  arma::vec rev_g = reverse_vec(g);
    
  for(unsigned int i = 0; i < L; i++){
  
    if( (i+!inverse) % 2 != 0){
      rev_g(i) = rev_g(i)*-1;
    }
    
  }
  
  return rev_g;
}

//' @title Haar filter construction
//' @description Creates the haar filter
//' @usage haar_filter()
//' @return A \code{field<vec>} that contains:
//' \itemize{
//'  \item{"L"}{A \code{integer} specifying the length of the filter}
//'  \item{"h"}{A \code{vector} containing the coefficients for the wavelet filter}
//'  \item{"g"}{A \code{vector} containing the coefficients for the scaling filter}
//' }
//' @details
//' This template can be used to increase the amount of filters available for selection.
//' @author JJB
//' @keywords internal
//' @examples
//' haar_filter()
// [[Rcpp::export]]
arma::field<arma::vec> haar_filter() {
  
    arma::vec L(1);
    L(0) = 2.0;
    
    arma::vec g(2);
    g.fill(0.7071067811865475);
    
    arma::vec h = qmf(g);
    
    arma::field<arma::vec> out(3);
    
    out(0)=L;
    out(1)=h;
    out(2)=g;
    
    return out;
}

//' @title Select the Wavelet Filter
//' @description Constructs the wavelet filter to be used.
//' @usage select_filter(filter_name)
//' @param filter_name A \code{String} that must receive: \code{"haar"}.
//' @return info A \code{field<vec>} that contains:
//' \itemize{
//'  \item{"L"}{A \code{integer} specifying the length of the filter}
//'  \item{"h"}{A \code{vector} containing the coefficients for the wavelet filter}
//'  \item{"g"}{A \code{vector} containing the coefficients for the scaling filter}
//' }
//' @details 
//' The package is oriented toward using only the haar filter. If the package extends at a later time, then the supporting infrastructure is there.
//' @author JJB
//' @keywords internal
//' @examples
//' select_filter("haar")
// [[Rcpp::export]]
arma::field<arma::vec> select_filter(std::string filter_name = "haar")
{
  
  arma::field<arma::vec> info(3);
  if(filter_name == "haar"){  
      info = haar_filter();
  }else{
      stop("Wave Filter is not supported! See ?select_filter for supported types."); 
  }
  
  return info;
}