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
#ifndef TS_MODEL_CPP_H
#define TS_MODEL_CPP_H

arma::field<arma::vec> model_objdesc(std::vector<std::string> desc);

arma::vec model_theta(std::vector<std::string> desc);

#endif