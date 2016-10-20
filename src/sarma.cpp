/* Copyright (C) 2016 James Balamuta, Stephane Guerrier, Roberto Molinari
 *
 * This file is part of GMWM R Methods Package
 *
 * The `gmwm` R package is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The `gmwm` R package is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *  
 */


#include <RcppArmadillo.h>
#include "sarma.h"


//' Create the ts.model obj.desc given split values
//' 
//' Computes the total phi and total theta vector length.
//' @param ar  A \code{vec} containing the non-seasonal phi parameters.
//' @param ma  A \code{vec} containing the non-seasonal theta parameters.
//' @param sar A \code{vec} containing the seasonal phi parameters.
//' @param sma A \code{vec} containing the seasonal theta parameters.
//' @param s   An \code{unsigned integer} containing the frequency of seasonality.
//' @param i   An \code{unsigned integer} containing the number of non-seasonal differences.
//' @param si  An \code{unsigned integer} containing the number of seasonal differences.
//' @return A \code{vec} with rows:
//' \describe{
//' \item{np}{Number of Non-Seasonal AR Terms}
//' \item{nq}{Number of Non-Seasonal MA Terms}
//' \item{nsp}{Number of Seasonal AR Terms}
//' \item{nsq}{Number of Seasonal MA Terms}
//' \item{nsigma}{Number of Variances (always 1)}
//' \item{s}{Season Value}
//' \item{i}{Number of non-seasonal differences}
//' \item{si}{Number of Seasonal Differences}
//' }
// [[Rcpp::export]]
arma::vec sarma_objdesc(const arma::vec& ar, const arma::vec& ma,
                        const arma::vec& sar, const arma::vec& sma,
                        int s, int i, int si){
  
  // Number of ARMA(p,q) parameters
  unsigned int np = ar.n_elem, nq = ma.n_elem;
  
  // Number of Seasonal (P,Q)
  unsigned int nsp = sar.n_elem, nsq = sma.n_elem;
  
  arma::vec objdesc(8);
  
  // p, q, P, Q, 1, s, i, si
  objdesc(0) = np;
  objdesc(1) = nq;
  objdesc(2) = nsp;
  objdesc(3) = nsq;
  objdesc(4) = 1;
  objdesc(5) = s;
  objdesc(6) = i;
  objdesc(7) = si;
  
  return objdesc;
}

//' Calculates Length of Seasonal Padding
//' 
//' Computes the total phi and total theta vector length.
//' @param np  An \code{unsigned int} containing the number of non-seasonal phi parameters.
//' @param nq  An \code{unsigned int} containing the number of non-seasonal theta parameters.
//' @param nsp An \code{unsigned int} containing the number of seasonal phi parameters.
//' @param nsq An \code{unsigned int} containing the number of seasonal theta parameters.
//' @seealso \code{\link{sarma_components}}
//' @return A \code{vec} with rows:
//' \describe{
//'  \item{p}{Number of phi parameters}
//'  \item{q}{Number of theta parameters}
//' }
//' @keywords internal
//' 
//' 
// [[Rcpp::export]]
arma::vec sarma_calculate_spadding(unsigned int np, unsigned int nq,
                                   unsigned int nsp, unsigned int nsq,
                                   unsigned int ns){
  
  arma::vec o(2);

  o(0) = np + ns * nsp;
  o(1) = nq + ns * nsq;
    
  return o;
}

//' Determine parameter expansion based upon objdesc
//' 
//' Calculates the necessary vec space needed to pad the vectors
//' for seasonal terms. 
//' @param objdesc A \code{vec} with the appropriate sarima object description
//' @return A \code{vec} with the structure:
//' \describe{
//' \item{np}{Number of Non-Seasonal AR Terms}
//' \item{nq}{Number of Non-Seasonal MA Terms}
//' \item{nsp}{Number of Seasonal AR Terms}
//' \item{nsq}{Number of Seasonal MA Terms}
//' \item{ns}{Number of Seasons (e.g. 12 is year)}
//' \item{p}{Total number of phi terms}
//' \item{q}{Total number of theta terms}
//' }
//' @keywords internal
// [[Rcpp::export]]
arma::vec sarma_components(const arma::vec& objdesc){
  // Number of ARMA(p,q) parameters
  unsigned int np = objdesc(0), nq = objdesc(1);
  
  // Number of Seasonal (P,Q) and Number of Seasons (ns)
  unsigned int nsp = objdesc(2), nsq = objdesc(3), ns = objdesc(5);
  
  // Find the total number of parameters to expand
  arma::vec nparams = sarma_calculate_spadding(np, nq, nsp, nsq, ns);

  // Create an output vector
  arma::vec o(7);
  o(0) = np, o(1) = nq, o(2) = nsp, o(3) = nsq, o(4) = ns, o(5) = nparams(0), o(6) = nparams(1);
  
  return o;
}

//' Efficient way to merge items together
//' @keywords internal
// [[Rcpp::export]]
arma::vec sarma_params_construct(const arma::vec& ar, const arma::vec& ma,
                                 const arma::vec& sar, const arma::vec& sma){
  
  unsigned int nar = ar.n_elem;
  unsigned int nma = ma.n_elem;
  unsigned int nsar = sar.n_elem;
  unsigned int nsma = sma.n_elem;
  
  // Make a common
  arma::vec params(nar + nma + nsar + nsma);
  
  // count parameters
  unsigned int count = nar;
  
  
  // Begin filling the param vector
  if(nar > 0){
    params.rows(0, count - 1) = ar;
  }
  
  if(nma > 0){
    params.rows(count, count + nma - 1) = ma;
    count += nma;
  }
  
  if(nsar > 0){
    params.rows(count, count + nsar - 1) = sar;
    count += nsar;
  }
  
  if(nsma > 0){
    params.rows(count, count + nsma - 1) = sma;
    count += nsma;
  }
  
  return params;
}



//' (Internal) Expand the SARMA Parameters
//' @param params  A \code{vec} containing the theta values of the parameters.
//' @inheritParams sarma_calculate_spadding
//' @param p An \code{unsigned int} that is the total size of the phi vector. 
//' @param q An \code{unsigned int} that is the total size of the theta vector. 
//' @return A \code{field<vec>} that contains the expansion. 
//' @keywords internal
// [[Rcpp::export]]
arma::field<arma::vec> sarma_expand_unguided(const arma::vec& params,
                                             unsigned int np, unsigned int nq,
                                             unsigned int nsp, unsigned int nsq,
                                             unsigned int ns,
                                             unsigned int p,
                                             unsigned int q){
  // Loop variables
  unsigned int i, j;
  
  // Fill with zeros in advance.
  arma::vec phi   = arma::zeros<arma::vec>(p);
  arma::vec theta = arma::zeros<arma::vec>(q);
  
  /* Fill the non-seasonal components */
  
  // Fill AR(p)
  for (i = 0; i < np; i++) phi(i) = params(i);
  // Fill MA(q)
  for (i = 0; i < nq; i++) theta(i) = params(i + np);
  
  
  /* Is a seasonal component present? */
  if (ns > 0) {
    
    /* Fill the seasonal components */
    
    // Fill the Seasonal AR(P)
    for (j = 0; j < nsp; j++) {
      phi((j + 1) * ns - 1) += params(j + np + nq);
      
      // Handle the interaction between AR and SAR terms
      for (i = 0; i < np; i++){
        phi((j + 1) * ns + i) -= params(i) * params(j + np + nq);
      }
    }
    
    // Fill the Seasonal MA(Q)
    for (j = 0; j < nsq; j++) {
      theta( (j + 1) * ns - 1 ) += params(j + np + nq + nsp);
      
      // Handle interaction between MA and SMA terms
      for (i = 0; i < nq; i++){
        theta( (j + 1) * ns + i ) += params(i + np) * params(j + np + nq + nsp);
      }
    }
  }
  
  // Release transformed parameter space.
  arma::field<arma::vec> out(2);
  
  out(0) = phi;
  out(1) = theta;
  
  return out;
}

//' Expand Parameters for an SARMA object
//' 
//' Creates an expanded PHI and THETA vector for use in other objects. 
//' @param params  A \code{vec} containing the theta values of the parameters.
//' @param objdesc A \code{vec} containing the model term information.
//' @return A \code{field<vec>} of size two as follows:
//' \itemize{
//'   \item AR    values
//'   \item THETA values
//' }
//' @details 
//' The \code{objdesc} is assumed to have the structure of:
//' \itemize{
//' \item AR(p)
//' \item MA(q)
//' \item SAR(P)
//' \item SMA(Q)
//' \item Seasons
//' }
//' @keywords internal
//' @examples
//' # p, q, P, Q, 1, s, i, si
//' m = sarma_expand(c(0.5,.2,0,.1,.92,.83,.42,.33,.12), c(2,2,2,3,1,12,0,0))
// [[Rcpp::export]]
arma::field<arma::vec> sarma_expand(const arma::vec& params, const arma::vec& objdesc) {

    // Get a breakdown of relevant obj desc values.
    arma::vec o = sarma_components(objdesc);
    // O has different breakdown than the objdesc! ()

    return sarma_expand_unguided(params,
                                 o(0), o(1),
                                 o(2), o(3),
                                 o(4),
                                 o(5),o(6));
}
