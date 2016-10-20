/* Copyright (C) 2014 - 2016  James Balamuta, Stephane Guerrier, Roberto Molinari
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
#include "wv_filters.h"
// We use reverse_vec
#include "armadillo_manipulations.h"


// Creates an index to call functions by memory address.
struct A{
  
  // @title Construct Filter Selection Map
  static std::map<std::string, arma::field<arma::vec> (*)()> create_map()
  {
    // Makes a map of function addresses indexed by a call string
    std::map<std::string, arma::field<arma::vec> (*)()> filterMap; 
    
    
    filterMap["haar"] = &haar_filter;
    
    // Change to appropriate function (make sure to add to "wv_filters.h")
    filterMap["d4"] = &d4_filter;
    filterMap["d6"] = &d6_filter;
    filterMap["d8"] = &d8_filter;
    filterMap["d16"] = &d16_filter;
    
    filterMap["fk4"] = &fk4_filter;
    filterMap["fk8"] = &fk8_filter;
    filterMap["fk14"] = &fk14_filter;
    filterMap["fk22"] = &fk22_filter;
    
    filterMap["bl14"] = &bl14_filter;
    filterMap["bl20"] = &bl20_filter;
    
    filterMap["la8"] = &la8_filter;
    filterMap["la16"] = &la16_filter;
    filterMap["la20"] = &la20_filter;
    
    filterMap["mb4"] = &mb4_filter;
    filterMap["mb8"] = &mb8_filter;
    filterMap["mb16"] = &mb16_filter;
    filterMap["mb24"] = &mb24_filter;

    return filterMap;
  }
  
  static const std::map<std::string, arma::field<arma::vec> (*)()> filterMap;
};

const std::map<std::string, arma::field<arma::vec> (*)()> A::filterMap =  A::create_map();

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


//' @title d4 filter construction
//' @description Creates the d4 filter
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
//' d4_filter()
// [[Rcpp::export]]
arma::field<arma::vec> d4_filter() {
  
  arma::vec L(1);
  L(0) = 4.0;
  
  arma::vec g(4);
  g(0) = 0.4829629131445341;
  g(1) = 0.8365163037378077;
  g(2) = 0.2241438680420134; 
  g(3) = -0.1294095225512603;
  
  arma::vec h = qmf(g);
  
  arma::field<arma::vec> out(3);
  
  out(0)=L;
  out(1)=h;
  out(2)=g;
  
  return out;
}



//' @title mb4 filter construction
//' @description Creates the mb4 filter
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
//' mb4_filter()
// [[Rcpp::export]]
arma::field<arma::vec> mb4_filter() {
  
  arma::vec L(1);
  L(0) = 4.0;
  
  arma::vec g(4);
  g(0) = 4.801755e-01;
  g(1) = 8.372545e-01;
  g(2) = 2.269312e-01; 
  g(3) = -1.301477e-01;
  
  arma::vec h = qmf(g);
  
  arma::field<arma::vec> out(3);
  
  out(0)=L;
  out(1)=h;
  out(2)=g;
  
  return out;
}

//' @title w4 filter construction
//' @description Creates the w4 filter
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
//' w4_filter()
// [[Rcpp::export]]
arma::field<arma::vec> w4_filter() {
  
  arma::vec L(1);
  L(0) = 4.0;
  
  arma::vec g(4);
  g(0) = -1/8;
  g(1) = 3/8;
  g(2) = 3/8;
  g(3) = -1/8;
  
  arma::vec h(4);
  g(0) = -1/8;
  g(1) = 3/8;
  g(2) = -3/8;
  g(3) = 1/8;
  
  arma::field<arma::vec> out(3);
  
  out(0)=L;
  out(1)=h;
  out(2)=g;
  
  return out;
}

//' @title fk4 filter construction
//' @description Creates the fk4 filter
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
//' fk4_filter()
// [[Rcpp::export]]
arma::field<arma::vec> fk4_filter() {
  
  arma::vec L(1);
  L(0) = 4.0;
  
  arma::vec g(4);
  g(0) = .6539275555697651; 
  g(1) = .7532724928394872;
  g(2) = .5317922877905981e-1; 
  g(3) = -.4616571481521770e-1;
  
  arma::vec h = qmf(g);
  
  arma::field<arma::vec> out(3);
  
  out(0)=L;
  out(1)=h;
  out(2)=g;
  
  return out;
}

//' @title d6 filter construction
//' @description Creates the d6 filter
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
//' d6_filter()
// [[Rcpp::export]]
arma::field<arma::vec> d6_filter() {
  
  arma::vec L(1);
  L(0) = 6.0;
  
  arma::vec g(6);
  g(0) = 0.3326705529500827; 
  g(1) = 0.8068915093110928;
  g(2) = 0.4598775021184915;
  g(3) = -0.1350110200102546;
  g(4) = -0.0854412738820267; 
  g(5) = 0.0352262918857096;
  
  arma::vec h = qmf(g);
  
  arma::field<arma::vec> out(3);
  
  out(0)=L;
  out(1)=h;
  out(2)=g;
  
  return out;
}


//' @title fk6 filter construction
//' @description Creates the fk6 filter
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
//' fk6_filter()
// [[Rcpp::export]]
arma::field<arma::vec> fk6_filter() {
  
  arma::vec L(1);
  L(0) = 6.0;
  
  arma::vec g(6);
  g(0) = 0.4279150324223103; 
  g(1) = 0.8129196431369074;
  g(2) = 0.3563695110701871;
  g(3) = -0.1464386812725773;
  g(4) = -0.7717775740697006e-1; 
  g(5) = 0.4062581442323794e-1;
  
  arma::vec h = qmf(g);
  
  arma::field<arma::vec> out(3);
  
  out(0)=L;
  out(1)=h;
  out(2)=g;
  
  return out;
}

//' @title d8 filter construction
//' @description Creates the d8 filter
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
//' d8_filter()
// [[Rcpp::export]]
arma::field<arma::vec> d8_filter() {
  
  arma::vec L(1);
  L(0) = 8.0;
  
  arma::vec g(8);
  g(0) = 0.2303778133074431;
  g(1) = 0.7148465705484058; 
  g(2) = 0.6308807679358788;
  g(3) = -0.0279837694166834;
  g(4) = -0.1870348117179132;
  g(5) = 0.0308413818353661;
  g(6) = 0.0328830116666778;
  g(7) = -0.0105974017850021;
  
  arma::vec h = qmf(g);
  
  arma::field<arma::vec> out(3);
  
  out(0)=L;
  out(1)=h;
  out(2)=g;
  
  return out;
}


//' @title fk8 filter construction
//' @description Creates the fk8 filter
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
//' fk8_filter()
// [[Rcpp::export]]
arma::field<arma::vec> fk8_filter() {
  
  arma::vec L(1);
  L(0) = 8.0;
  
  arma::vec g(8);
  g(0) = .3492381118637999;
  g(1) = .7826836203840648; 
  g(2) = .4752651350794712;
  g(3) = -.9968332845057319e-1;
  g(4) = -.1599780974340301;
  g(5) = .4310666810651625e-1;
  g(6) = .4258163167758178e-1;
  g(7) = -.1900017885373592e-1;
  
  arma::vec h = qmf(g);
  
  arma::field<arma::vec> out(3);
  
  out(0)=L;
  out(1)=h;
  out(2)=g;
  
  return out;
}

//' @title la8 filter construction
//' @description Creates the la8 filter
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
//' la8_filter()
// [[Rcpp::export]]
arma::field<arma::vec> la8_filter() {
  
  arma::vec L(1);
  L(0) = 8.0;
  
  arma::vec g(8);
  g(0) = -0.07576571478935668;  
  g(1) = -0.02963552764596039;
  g(2) = 0.49761866763256290; 
  g(3) =  0.80373875180538600;
  g(4) = 0.29785779560560505;
  g(5) = -0.09921954357695636; 
  g(6) = -0.01260396726226383;
  g(7) = 0.03222310060407815;
  
  arma::vec h = qmf(g);
  
  arma::field<arma::vec> out(3);
  
  out(0)=L;
  out(1)=h;
  out(2)=g;
  
  return out;
}

//' @title mb8 filter construction
//' @description Creates the mb8 filter
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
//' mb8_filter()
// [[Rcpp::export]]
arma::field<arma::vec> mb8_filter() {
  
  arma::vec L(1);
  L(0) = 8.0;
  
  arma::vec g(8);
  g(0) = 6.436345e-02;  
  g(1) = 7.106015e-03;
  g(2) = -1.108673e-01; 
  g(3) =  2.947855e-01;
  g(4) = 7.351331e-01;
  g(5) = 5.725771e-01; 
  g(6) = 1.847751e-02;
  g(7) = -1.673619e-01;
  
  arma::vec h = qmf(g);
  
  arma::field<arma::vec> out(3);
  
  out(0)=L;
  out(1)=h;
  out(2)=g;
  
  return out;
}

//' @title bl14 filter construction
//' @description Creates the bl14 filter
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
//' bl14_filter()
// [[Rcpp::export]]
arma::field<arma::vec> bl14_filter() {
  
  arma::vec L(1);
  L(0) = 14.0;
  
  arma::vec g(14);
  g(0) = 0.0120154192834842;
  g(1) = 0.0172133762994439; 
  g(2) = -0.0649080035533744;
  g(3) =   -0.0641312898189170; 
  g(4) =   0.3602184608985549;
  g(5) =   0.7819215932965554;
  g(6) =   0.4836109156937821;
  g(7) =   -0.0568044768822707; 
  g(8) =   -0.1010109208664125;
  g(9) =   0.0447423494687405;
  g(10) =   0.0204642075778225; 
  g(11) =   -0.0181266051311065;
  g(12) =   -0.0032832978473081;
  g(13) =   0.0022918339541009;
  
  arma::vec h = qmf(g);
  
  arma::field<arma::vec> out(3);
  
  out(0)=L;
  out(1)=h;
  out(2)=g;
  
  return out;
}

//' @title fk14 filter construction
//' @description Creates the fk14 filter
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
//' fk14_filter()
// [[Rcpp::export]]
arma::field<arma::vec> fk14_filter() {
  
  arma::vec L(1);
  L(0) = 14.0;
  
  arma::vec g(14);
  g(0) = .2603717692913964; 
  g(1) = .6868914772395985;
  g(2) = .6115546539595115;
  g(3) = .5142165414211914e-1;
  g(4) = -.2456139281621916;
  g(5) = -.4857533908585527e-1;
  g(6) = .1242825609215128;
  g(7) = .2222673962246313e-1;
  g(8) = -.6399737303914167e-1;
  g(9) = -.5074372549972850e-2;
  g(10) = .2977971159037902e-1;
  g(11) = -.3297479152708717e-2;
  g(12) = -.9270613374448239e-2; 
  g(13) = .3514100970435962e-2;
  
  arma::vec h = qmf(g);
  
  arma::field<arma::vec> out(3);
  
  out(0)=L;
  out(1)=h;
  out(2)=g;
  
  return out;
}

//' @title d16 filter construction
//' @description Creates the d16 filter
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
//' d16_filter()
// [[Rcpp::export]]
arma::field<arma::vec> d16_filter() {
  
  arma::vec L(1);
  L(0) = 16.0;
  
  arma::vec g(16) ;
  
  g(0) = 0.0544158422431049;
  g(1) = 0.3128715909143031;
  g(2) = 0.6756307362972904;
  g(3) = 0.5853546836541907; 
  g(4) = -0.0158291052563816; 
  g(5) = -0.2840155429615702;
  g(6) = 0.0004724845739124;
  g(7) = 0.1287474266204837;
  g(8) = -0.0173693010018083;
  g(9) = -0.0440882539307952; 
  g(10) = 0.0139810279173995; 
  g(11) = 0.0087460940474061;
  g(12) = -0.0048703529934518;
  g(13) = -0.0003917403733770;
  g(14) = 0.0006754494064506;
  g(15) = -0.0001174767841248;
  
  arma::vec h = qmf(g);
  
  arma::field<arma::vec> out(3);
  
  out(0)=L;
  out(1)=h;
  out(2)=g;
  
  return out;
}

//' @title la16 filter construction
//' @description Creates the la16 filter
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
//' la16_filter()
// [[Rcpp::export]]
arma::field<arma::vec> la16_filter() {
  
  arma::vec L(1);
  L(0) = 16.0;
  
  arma::vec g(16) ;
  
  g(0) = -0.0033824159513594;
  g(1) = -0.0005421323316355; 
  g(2) = 0.0316950878103452;
  g(3) = 0.0076074873252848; 
  g(4) = -0.1432942383510542; 
  g(5) = -0.0612733590679088;
  g(6) = 0.4813596512592012;
  g(7) = 0.7771857516997478;
  g(8) = 0.3644418948359564;
  g(9) = -0.0519458381078751;
  g(10) = -0.0272190299168137;
  g(11) = 0.0491371796734768;
  g(12) = 0.0038087520140601;
  g(13) = -0.0149522583367926;
  g(14) = -0.0003029205145516;
  g(15) = 0.0018899503329007;
  
  arma::vec h = qmf(g);
  
  arma::field<arma::vec> out(3);
  
  out(0)=L;
  out(1)=h;
  out(2)=g;
  
  return out;
}

//' @title mb16 filter construction
//' @description Creates the mb16 filter
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
//' mb16_filter()
// [[Rcpp::export]]
arma::field<arma::vec> mb16_filter() {
  
  arma::vec L(1);
  L(0) = 16.0;
  
  arma::vec g(16) ;
  
  g(0) = 5.765899e-03;
  g(1) = 9.620427e-03; 
  g(2) = -4.984698e-02;
  g(3) = -2.483876e-02; 
  g(4) = 5.474628e-02; 
  g(5) = -1.987986e-02;
  g(6) = -5.656657e-02;
  g(7) = 2.345342e-01;
  g(8) = 6.701646e-01;
  g(9) = 6.349228e-01;
  g(10) = 1.188725e-01;
  g(11) = -2.278359e-01;
  g(12) = -5.776570e-02;
  g(13) = 1.136116e-01;
  g(14) = 2.173677e-02;
  g(15) = -1.302770e-02;
  
  arma::vec h = qmf(g);
  
  arma::field<arma::vec> out(3);
  
  out(0)=L;
  out(1)=h;
  out(2)=g;
  
  return out;
}

//' @title la20 filter construction
//' @description Creates the la20 filter
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
//' la20_filter()
// [[Rcpp::export]]
arma::field<arma::vec> la20_filter() {
  
  arma::vec L(1);
  L(0) = 20.0;
  
  arma::vec g(20) ;
  
  g(0) = 0.0007701598091030; 
  g(1) = 0.0000956326707837;
  g(2) = -0.0086412992759401;
  g(3) = -0.0014653825833465;
  g(4) = 0.0459272392237649;
  g(5) = 0.0116098939129724;
  g(6) = -0.1594942788575307;
  g(7) = -0.0708805358108615;
  g(8) = 0.4716906668426588;
  g(9) = 0.7695100370143388;
  g(10) = 0.3838267612253823;
  g(11) = -0.0355367403054689;
  g(12) = -0.0319900568281631;
  g(13) = 0.0499949720791560; 
  g(14) = 0.0057649120455518;
  g(15) = -0.0203549398039460;
  g(16) = -0.0008043589345370;
  g(17) = 0.0045931735836703;
  g(18) = 0.0000570360843390;
  g(19) = -0.0004593294205481;
  
  arma::vec h = qmf(g);
  
  arma::field<arma::vec> out(3);
  
  out(0)=L;
  out(1)=h;
  out(2)=g;
  
  return out;
}

//' @title bl20 filter construction
//' @description Creates the bl20 filter
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
//' bl20_filter()
// [[Rcpp::export]]
arma::field<arma::vec> bl20_filter() {
  
  arma::vec L(1);
  L(0) = 20.0;
  
  arma::vec g(20) ;
  
  g(0) = 0.0008625782242896;
  g(1) = 0.0007154205305517;
  g(2) = -0.0070567640909701;
  g(3) = 0.0005956827305406;
  g(4) = 0.0496861265075979;
  g(5) = 0.0262403647054251;
  g(6) = -0.1215521061578162;
  g(7) = -0.0150192395413644;
  g(8) = 0.5137098728334054;
  g(9) = 0.7669548365010849;
  g(10) = 0.3402160135110789;
  g(11) = -0.0878787107378667;
  g(12) = -0.0670899071680668;
  g(13) = 0.0338423550064691;
  g(14) = -0.0008687519578684;
  g(15) = -0.0230054612862905;
  g(16) = -0.0011404297773324;
  g(17) = 0.0050716491945793;
  g(18) = 0.0003401492622332;
  g(19) = -0.0004101159165852;
  
  arma::vec h = qmf(g);
  
  arma::field<arma::vec> out(3);
  
  out(0)=L;
  out(1)=h;
  out(2)=g;
  
  return out;
}


//' @title fk22 filter construction
//' @description Creates the fk22 filter
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
//' fk22_filter()
// [[Rcpp::export]]
arma::field<arma::vec> fk22_filter() {
  
  arma::vec L(1);
  L(0) = 22.0;
  
  arma::vec g(22) ;
  
  g(0) = .1938961077599566;
  g(1) = .5894521909294277;
  g(2) = .6700849629420265;
  g(3) = .2156298491347700;
  g(4) = -.2280288557715772;
  g(5) = -.1644657152688429;
  g(6) = .1115491437220700;
  g(7) = .1101552649340661;
  g(8) = -.6608451679377920e-1;
  g(9) = -.7184168192312605e-1;
  g(10) = .4354236762555708e-1;
  g(11) = .4477521218440976e-1;
  g(12) = -.2974288074927414e-1;
  g(13) = -.2597087308902119e-1;
  g(14) = .2028448606667798e-1;
  g(15) = .1296424941108978e-1;
  g(16) = -.1288599056244363e-1;
  g(17) = -.4838432636440189e-2;
  g(18) = .7173803165271690e-2;
  g(19) = .3612855622194901e-3;
  g(20) = -.2676991638581043e-2;
  g(21) = .8805773686384639e-3;
  
  arma::vec h = qmf(g);
  
  arma::field<arma::vec> out(3);
  
  out(0)=L;
  out(1)=h;
  out(2)=g;
  
  return out;
}

//' @title mb24 filter construction
//' @description Creates the mb24 filter
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
//' mb24_filter()
// [[Rcpp::export]]
arma::field<arma::vec> mb24_filter() {
  
  arma::vec L(1);
  L(0) = 24.0;
  
  arma::vec g(24) ;
  
  g(23) = -2.132706e-05; 
  g(22) = 4.745736e-04;
  g(21) = 7.456041e-04;
  g(20) = -4.879053e-03;
  g(19) = -1.482995e-03;
  g(18) = 4.199576e-02;
  g(17) = -2.658282e-03;
  g(16) = -6.559513e-03;
  g(15) = 1.019512e-01;
  g(14) = 1.689456e-01;
  g(13) = 1.243531e-01;
  g(12) = 1.949147e-01;
  g(11) = 4.581101e-01;
  g(10) = 6.176385e-01;
  g(9) = 2.556731e-01;
  g(8) = -3.091111e-01;
  g(7) = -3.622424e-01;
  g(6) = -4.575448e-03;
  g(5) = 1.479342e-01;
  g(4) = 1.027154e-02;
  g(3) = -1.644859e-02;
  g(2) = -2.062335e-03;
  g(1) = 1.193006e-03;
  g(0) = 5.361301e-05;
  
  
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
  
  
  std::map<std::string,arma::field<arma::vec> (*)()>::const_iterator it = A::filterMap.find(filter_name);
  if(it != A::filterMap.end())
  {
    //element found;
    info = (*(it->second))();
  }else{
    Rcpp::stop("Wave Filter is not supported! See ?select_filter for supported types."); 
  }
  
  return info;
}