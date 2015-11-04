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

#ifndef DWT
#define DWT

arma::field<arma::vec> dwt_cpp(arma::vec x, std::string filter_name, 
                                 unsigned int nlevels, std::string boundary);
                                 
arma::field<arma::vec> modwt_cpp(arma::vec x, std::string filter_name, 
                                   unsigned int nlevels, std::string boundary);
                                   
arma::field<arma::vec> brick_wall(arma::field<arma::vec> x,  arma::field<arma::vec> wave_filter, std::string method);

#endif
