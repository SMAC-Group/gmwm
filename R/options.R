#' @title Place Legend
#' @description This function decides where the legend should be put (top left or bottom left)
#' @method placeLegend
#' @param wv_1 A \code{double} that indicates the first value of \code{wv.empir}
#' @param low_n A \code{doble} that indicates the last value of \code{ci.low}
#' @param high_n A \code{dobule} that indicates the last value of \code{ci.high}
#' @return A numeric vector containing 4 elements. The first two elements indicate legend.justification, the last two elements indicate legend.position (see \code{?theme}).
placeLegend = function(wv_1, low_n, high_n){
  if(log10(wv_1) > ( log10(low_n) + log10(high_n) )/2 ){
    # legend should be placed in bottom left
    legend.justification = c(0,0)
    legend.position = c(0,0)
    #x = wv_1[1]/xlim_length
    #y = high_n/ylim_length
  }
  else{
    # legend should be placed in top left
    legend.justification = c(0,1)
    legend.position = c(0,1)
    #x = wv_1[1]/xlim_length
    #y = low_n/ylim_length
  }
  return( c(legend.justification, legend.position) )
  
}