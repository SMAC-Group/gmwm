#' @title Create a GMWM TS Object
#' @description Setups a time series oriented object that works well with graphing and summary utilities
#' @method gts
#' @param data A one column matrix, data.frame, or numeric vector.
#' @param name A string that provides an identifier to the data. If not supplied, the name of the first column will be used. Otherwise, it will be registered as "Dataset X"
#' @param freq A number that provides the rate of samples.
#' @param unit A string that contains the unit expression of the frequency. 
#' @return An S3 object with called gts with the following structure:
#' \itemize{
#'  \item{data}{Data set information}
#'  \item{name}{Name of the Dataset}
#'  \item{freq}{Numeric representation of frequency}
#'  \item{Unit}{String representation of the Unit}
#' }
#' @author JJB
#' @examples
#' m = data.frame(rnorm(5))
#' gts(5)
#' 
#' x = gen.ts(WN(sigma2=1), 50)
#' x = gts(x)
gts = function(data, name = "", freq = 1, unit = "sec"){

  # Force data.frame to matrix  
  if (is.data.frame(data)){ 
    data = data.matrix(data)
  }
  
  # Check if the data is in matrix form
  if (is.matrix(data)) {
    ndata = nrow(data)
  } else {
    ndata = length(data)
  }
  if (ndata == 0) {
    stop("Not a valid data object! Please supply a data set with one column that is in either a data.frame, matrix,  or numeric object.")
  }
  
  if(!is(freq,"numeric")){
    stop("freq must be numeric")
  }
  
  if(!is(name,"character") || !is(unit,"character")){
    stop("Name or unit is not a character")
  }

  x = structure(list(data = data, name = name, freq = freq, unit = unit), class = "gts")
  x
}