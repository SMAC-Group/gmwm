# Copyright (C) 2014 - 2015  James Balamuta, Stephane Guerrier, Roberto Molinari
#
# This file is part of GMWM R Methods Package
#
# The `gmwm` R package is free software: you can redistribute it and/or modify it
# under the terms of the Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)
# included within the packages source as the LICENSE file.
#
# The `gmwm` R package is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# You should have received a copy of the Attribution-NonCommercial-ShareAlike 4.0 International 
# (CC BY-NC-SA 4.0) along with `gmwm`.  If not, see <http://www.smac-group.com/licensing/>.

.onAttach <- function(...) {
  if (!interactive()) return()
  
  tips = c(
    "Need help? E-mail the SMAC-Group: contact@smac-group.com",
    "Find out what's changed in the gmwm package at http://github.com/SMAC-Group/gmwm",
    "Use suppressPackageStartupMessages() to eliminate package startup messages."
  )
  
  tip = sample(tips, 1)
  
  tip = c(tip,"Please note, the optimization algorithm is a work in progress!")
  packageStartupMessage(paste(strwrap(tip), collapse = "\n"))
}