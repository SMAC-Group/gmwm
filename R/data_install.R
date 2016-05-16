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

#' @title Install IMU Data Package
#' @description Downloads and Install the IMU Data Package to use with the GMWM Package
#' @param type A \code{string} with value \code{"ONL"} or \code{"LOCAL"}
#' @param loc  A \code{string} that contains the file location.
#' @details 
#' The IMU Data Package contains data from real life IMUs and is approximately a
#'  20 MB download.
#' @author JJB
#' @examples 
#' \dontrun{
#' # Online install
#' install_imudata()
#' 
#' # Local install
#' install_imudata("LOCAL","C:/Users/James/Documents/imudata_x.y.z.tar.gz")
#' }
install_imudata = function(type="ONL",loc=NULL){
  install_data_package("imudata", type, loc)
}

#' @title Install Data Package
#' @description Downloads and Install the Data Package to use with the GMWM Package
#' @param type A \code{string} with value \code{"ONL"} or \code{"LOCAL"}
#' @param loc  A \code{string} that contains the file location.
#' @author JJB
#' @examples 
#' \dontrun{
#' # Online install
#' install_datapkg()
#' 
#' # Local install
#' install_datapkg("LOCAL","C:/Users/James/Documents/datapkg_x.y.z.tar.gz")
#' }
install_imudata = function(type="ONL",loc=NULL){
  install_data_package("datapkg", type, loc)
}


#' @title Install a SMAC Group Data Package
#' @description Downloads and Install the a SMAC Group Data Package to use.
#' @param pkg  A \code{string} with the package name to install.
#' @param type A \code{string} with value \code{"ONL"} or \code{"LOCAL"}
#' @param loc  A \code{string} that contains the file location.
#' @details 
#' For a list of data packages please see the SMAC-Group data repo \url{http://smac-group.com/datarepo}.
#' @author JJB
#' @keywords internal
#' @examples 
#' \dontrun{
#' # Online install
#' install_data_package("imudata")
#' 
#' # Local install
#' install_data_package("LOCAL","C:/Users/James/Documents/imudata_x.y.z.tar.gz")
#' }
install_data_package = function(pkg = "imudata", type = "ONL", loc = NULL){
  type = toupper(type)
  
  if(type == "ONL"){
    install.packages(pkg, repos = "http://smac-group.com/datarepo")
  }else{
    if(!is.null(loc)){ 
      message("Installing the package using a local .tar file.")
      message("Note: You must have a compiler installed!")
      install.packages(loc, type="source")
    }else{ 
      message("To install locally, you must also set the `loc` paramter!")
    }
  }
  
}

