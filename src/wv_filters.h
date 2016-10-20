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

#ifndef WV_FILTERS
#define WV_FILTERS

arma::vec qmf(arma::vec g, bool inverse);

// Filters
arma::field<arma::vec> haar_filter();

// Daubechies 
arma::field<arma::vec> d4_filter() ;

arma::field<arma::vec> d6_filter() ;
arma::field<arma::vec> d8_filter() ;
arma::field<arma::vec> d16_filter() ;

arma::field<arma::vec> fk4_filter() ;
arma::field<arma::vec> fk8_filter() ;
arma::field<arma::vec> fk14_filter() ;
arma::field<arma::vec> fk22_filter() ;

arma::field<arma::vec> bl14_filter() ;
arma::field<arma::vec> bl20_filter() ;

arma::field<arma::vec> la8_filter() ;
arma::field<arma::vec> la16_filter() ;
arma::field<arma::vec> la20_filter() ;

arma::field<arma::vec> mb4_filter() ;
arma::field<arma::vec> mb8_filter() ;
arma::field<arma::vec> mb16_filter() ;
arma::field<arma::vec> mb24_filter() ;


arma::field<arma::vec> select_filter(std::string filter_name);

#endif
