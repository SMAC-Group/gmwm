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
#'  13.8 MB download.
#' @author JJB
#' @examples 
#' \dontrun{
#' # Online install
#' install_imudata()
#' 
#' # Local install
#' install_imudata("LOCAL","C:/Users/James/Documents/gmwm_x.y.z.tar.gz")
#' }
install_imudata = function(type="ONL",loc=NULL){
  type = toupper(type)
  
  if(type == "ONL"){
    install_("SMAC-Group","imudata")
  }else{
    if(!is.null(loc)){ 
      message("Installing the package using a local .tar file.")
      message("Note: You must have a compiler installed!")
      message("The following packages are also required:")
      message("ggplot2, RcppArmadillo, and gridExtra")
      install.packages(loc, type="source")
    }else{ 
      message("To install locally, you must also set the `loc` paramter!")
    }
  }
  
}

##

#' @title Light weight package install
#' @description
#' Downloads an R package from GitHub via the zip archive option,
#' unzips the repository, builds the R package source, and 
#' install the package
#' @param user      A \code{string} containing the user name
#' @param pkg.name  A \code{string} that contains the repository name
#' @details
#' This is a lightweight standalone method to download and install R packages from GitHub.
#' It is much safer to use devtools, though, it is one additional "download".
#' @author JJB
#' @examples 
#' \dontrun{
#' # Online install
#' install_("SMAC-Group","imudata")
#' }
install_ <- function (user,pkg.name) {
  # Step 1: Create a temp file
  tdir = tempdir()
  temp = file.path(tdir,paste0(pkg.name,".zip"))
  
  # Step 2: Download into temp file
  dl.file(paste0("https://github.com/",user,"/",pkg.name,"/archive/master.zip"),
                destfile = temp, mode = "wb")
  
  # Step 3: Extract
  extract_to = tempfile()
  unzip(temp,exdir=extract_to)
  
  # Step 4: Rename extract folder to just be package name.
  file.rename(file.path(extract_to,paste0(pkg.name,"-master")), file.path(extract_to,pkg.name))

  # Step 5: Build the .tar.gz source
  system(paste0('"',file.path(R.home("bin"), "R"),'" CMD build "', file.path(extract_to,pkg.name),'"'))
  
  # Step 6: Install the package from source
  install.packages(file.path(extract_to,pkg.name),repos=NULL,type="source")
  
  # Step 7: Unlink / delete on quit.
  on.exit(unlink(extract_to), add = TRUE)
  on.exit(unlink(tdir),add = TRUE)
}

#' @title Better file downloader
#' @description Replaces download.file's auto detect mechanism
#' @param url A \code{string} that contains a url
#' @param ... Additional parameters to be passed to \code{\link[utils]{download.file}}
#' @return A downloaded file...
#' @details 
#' See \url{https://support.rstudio.com/hc/en-us/articles/206827897-Secure-Package-Downloads-for-R}
#' for information on downloading via HTTPS. Also, read through the changes in \code{\link[utils]{download.file}}.
#' @seealso \code{\link[utils]{download.file}}
#' @keywords internal
dl.file = function(url, ...){
  if(grepl('^http:', url)){
    download.file(url, ...)
  }else{ 
    if(Sys.info()[['sysname']] != "Windows") {
      
      if(capabilities("libcurl")){         # Shipped lib
        method = "libcurl"
      }else if(nzchar(Sys.which("curl"))){ # Curl lib
        method = "curl"
      }else if(nzchar(Sys.which("wget"))){ # Historical, last resort
        method = "wget"
      }else{
        stop("Unable to download over HTTPS.\nPlease install `curl` or `wget` on your system.")
      }# end if
      
    }else{
      method = "auto"
    } # Let windows choose, since package requires >= 3.2, wininet will be automatically selected
    
    download.file(url, method = method, ...)
  }
}