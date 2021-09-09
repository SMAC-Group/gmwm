#' @section Process Haar WV Second Derivative:
#' Taking the second derivative with respect to \eqn{\phi}{phi} yields:
#' \deqn{\frac{{{\partial ^2}}}{{\partial {\phi ^2}}}\nu _j^2\left( \phi, \sigma ^2  \right) = \frac{2 \sigma ^2 \left(\left(\phi ^2-1\right) \tau _j \left(2 (\phi  (7 \phi +4)+1) \phi ^{\frac{\tau _j}{2}-1}-(\phi  (7 \phi +4)+1) \phi ^{\tau _j-1}+3 (\phi +1)^2\right)+\left(\phi ^2-1\right)^2 \tau _j^2 \left(\phi ^{\frac{\tau _j}{2}}-1\right) \phi ^{\frac{\tau _j}{2}-1}+4 (3 \phi +1) \left(\phi ^2+\phi +1\right) \left(\phi ^{\tau _j}-4 \phi ^{\frac{\tau _j}{2}}+3\right)\right)}{(\phi -1)^5 (\phi +1)^3 \tau _j^2} }{d^2/dphi^2 nu[j]^2(phi, sigma2) = 2*sigma2*(4*(1 + 3*phi)*(1 + phi + phi^2)*
#' (3 - 4*phi^(tau[j]/2) + phi^tau[j]) + (-1 + phi^2)*
#'   (3*(1 + phi)^2 + 2*phi^(tau[j]/2 - 1)*(1 + phi*(4 + 7*phi)) - 
#'      phi^(tau[j] - 1)*(1 + phi*(4 + 7*phi)))*
#'  tau[j] + phi^(tau[j]/2 - 1)*(-1 + phi^2)^2*(-1 + phi^(tau[j]/2))*tau[j]^2)/((-1 + phi)^5*(1 + phi)^3*tau[j]^2)}
#' 
#' Taking the second derivative with respect to \eqn{\sigma^2}{sigma^2} yields:
#' \deqn{\frac{{{\partial ^2}}}{{\partial {\sigma ^4}}}\nu _j^2\left( \sigma ^2  \right) = 0 }{d^2/dsigma2^2 nu[j]^2(phi, sigma2) = 0}
#' 
#' Taking the derivative with respect to \eqn{\phi}{phi} and \eqn{\sigma ^2}{sigma2} yields:
#' \deqn{\frac{{{\partial ^2}}}{{\partial {\phi } \partial {\sigma ^2}}}\nu _j^2\left( \phi, \sigma ^2  \right) = \frac{2 \left(\left(\phi ^2-1\right) \tau _j \left(\phi ^{\tau _j}-2 \phi ^{\frac{\tau _j}{2}}-\phi -1\right)-(\phi  (3 \phi +2)+1) \left(\phi ^{\tau _j}-4 \phi ^{\frac{\tau _j}{2}}+3\right)\right)}{(\phi -1)^4 (\phi +1)^2 \tau _j^2} }{d/dsigma * d/dphi nu[j]^2(phi, sigma2) = (2*((-(3 - 4*phi^(tau[j]/2) + phi^tau[j]))*(1 + phi*(2 + 3*phi)) + (-1 + phi^2)*(-1 - phi - 2*phi^(tau[j]/2) + phi^tau[j])*tau[j]))/((-1 + phi)^4*(1 + phi)^2*tau[j]^2)}
