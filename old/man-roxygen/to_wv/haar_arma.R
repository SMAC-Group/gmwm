#' @section Process Haar Wavelet Variance Formula:
#' The Autoregressive Order \eqn{p} and Moving Average Order \eqn{q} (ARMA(\eqn{p},\eqn{q})) process has a Haar Wavelet Variance given by:
#' \deqn{\frac{{{\tau _j}\left[ {1 - \rho \left( {\frac{{{\tau _j}}}{2}} \right)} \right] + 2\sum\limits_{i = 1}^{\frac{{{\tau _j}}}{2} - 1} {i\left[ {2\rho \left( {\frac{{{\tau _j}}}{2} - i} \right) - \rho \left( i \right) - \rho \left( {{\tau _j} - i} \right)} \right]} }}{{\tau _j^2}}\sigma _X^2}{(tau[j]*(1-rho(tau[j]/2)) + 2*sum(i*(2*rho(tau[j]/2 - i) + rho(i) - rho(tau[j] - i))))/tau[j]^2 * sigma[x]^2}
#'  where \eqn{\sigma _X^2}{sigma[X]^2} is given by the variance of the ARMA process. 
#' Furthermore, this assumes that stationarity has been achieved as it directly 
