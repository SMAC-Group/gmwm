#' @section Process Haar Wavelet Variance Formula:
#' The Autoregressive Order \eqn{1} and Moving Average Order \eqn{1} (ARMA(\eqn{1},\eqn{1})) process has a Haar Wavelet Variance given by:
#' \deqn{\nu _j^2\left( {\phi ,\theta ,{\sigma ^2}} \right) =  - \frac{{2{\sigma ^2}\left( { - \frac{1}{2}{{(\theta  + 1)}^2}\left( {{\phi ^2} - 1} \right){\tau _j} - (\theta  + \phi )(\theta \phi  + 1)\left( {{\phi ^{{\tau _j}}} - 4{\phi ^{\frac{{{\tau _j}}}{2}}} + 3} \right)} \right)}}{{{{(\phi  - 1)}^3}(\phi  + 1)\tau _j^2}}}{nu[j]^2 (phi, theta, sigma2) = (-2*sigma2*((-(theta + phi))*(1 + theta*phi)*(3 - 4*phi^(tau[j]/2) + phi^tau[j]) - 0.5*(1 + theta)^2*(-1 + phi^2)*tau[j])) / ((-1 + phi)^3*(1 + phi)*tau[j]^2)}
#'  
