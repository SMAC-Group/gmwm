#' @section Process Haar WV Second Derivative:
#' Taking the derivative with respect to \eqn{\theta}{theta} yields:
#' \deqn{\frac{{{\partial ^2}}}{{\partial {\theta ^2}}}\nu _j^2\left( {\theta ,{\sigma ^2}} \right) = \frac{{2{\sigma ^2}}}{{{\tau _j}}}}{d^2/dtheta^2 nu[j]^2 (theta, sigma2) = (2*sigma2)/tau[j]}
#'
#' Taking the derivative with respect to \eqn{\sigma^2}{sigma^2} yields:
#' \deqn{\frac{{{\partial ^2}}}{{\partial {\sigma ^4}}}\nu _j^2\left( {\theta ,{\sigma ^2}} \right) = 0}{d^2/dsigma2^2 nu[j]^2 (theta, sigma2) = 0}
#'
#' Taking the first derivative with respect to \eqn{\theta}{theta} and \eqn{\sigma^2}{sigma^2} yields:
#' \deqn{\frac{\partial }{{\partial \theta }}\frac{\partial }{{\partial {\sigma ^2}}}\nu _j^2\left( {\theta ,{\sigma ^2}} \right) = \frac{{{\sigma ^2}\left( {2\left( {\theta  + 1} \right){\tau _j} - 6} \right)\left( {{{\left( {\theta  + 1} \right)}^2}{\tau _j} - 6\theta } \right)}}{{\tau _j^4}}}{d/dtheta * d/dsigma2 nu[j]^2 (theta, sigma2) = (sigma2*(2*(theta+1)*tau[j]-6)*((theta+1)^2*tau[j]-6*theta))/(tau[j]^4)}