/* Copyright (C) 2014 - 2015  James Balamuta, Stephane Guerrier, Roberto Molinari
 *
 * This file is part of GMWM R Methods Package
 *
 * The `gmwm` R package is free software: you can redistribute it and/or modify it
 * under the terms of the Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)
 * as the LICENSE file.
 *
 * The `gmwm` R package is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * You should have received a copy of the Attribution-NonCommercial-ShareAlike 4.0 International 
 * (CC BY-NC-SA 4.0) along with `gmwm`.  If not, see <http://www.smac-group.com/licensing/>.
 * 
 */

#ifndef DWT
#define DWT

arma::field<arma::vec> brick_wall(arma::field<arma::vec> x,  
                                  arma::field<arma::vec> wave_filter, 
                                  std::string method = "modwt") ;

arma::field<arma::vec> dwt_cpp(arma::vec x, std::string filter_name = "haar", 
                               unsigned int nlevels = 4, std::string boundary = "periodic", bool brickwall = true);
                                 
arma::field<arma::vec> modwt_cpp(arma::vec x, std::string filter_name = "haar", 
                                 unsigned int nlevels = 4, std::string boundary = "periodic", bool brickwall = true);
#endif
