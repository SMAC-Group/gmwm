#' @section Process Haar WV Second Derivative:
#' Taking the second derivative with respect to \eqn{\phi}{phi} yields:
#' \deqn{\frac{{{\partial ^2}}}{{\partial {\phi ^2}}}\nu _j^2\left( {\phi ,\theta ,{\sigma ^2}} \right) = \frac{{2{\sigma ^2}}}{{{{(\phi  - 1)}^5}{{(\phi  + 1)}^3}\tau _j^2}}\left( \begin{array}{cc}
#' &{(\phi  - 1)^2}\left( {{{(\phi  + 1)}^2}\left( {{\theta ^2}\phi  + \theta {\phi ^2} + \theta  + \phi } \right)\tau _j^2\left( {{\phi ^{\frac{{{\tau _j}}}{2}}} - 1} \right){\phi ^{\frac{{{\tau _j}}}{2} - 2}} + \left( {{\phi ^2} - 1} \right)\left( {{\theta ^2}( - \phi ) + \theta \left( {{\phi ^2} + 4\phi  + 1} \right) - \phi } \right){\tau _j}\left( {{\phi ^{\frac{{{\tau _j}}}{2}}} - 2} \right){\phi ^{\frac{{{\tau _j}}}{2} - 2}} - 2{{(\theta  - 1)}^2}\left( {{\phi ^{{\tau _j}}} - 4{\phi ^{\frac{{{\tau _j}}}{2}}} + 3} \right)} \right) \\
#' &- 12{(\phi  + 1)^2}\left( { - \frac{1}{2}{{(\theta  + 1)}^2}\left( {{\phi ^2} - 1} \right){\tau _j} - (\theta  + \phi )(\theta \phi  + 1)\left( {{\phi ^{{\tau _j}}} - 4{\phi ^{\frac{{{\tau _j}}}{2}}} + 3} \right)} \right) \\
#' &+ 6(\phi  + 1)(\phi  - 1)\left( {\frac{1}{2}{{(\theta  + 1)}^2}\left( {{\phi ^2} - 1} \right){\tau _j} + (\theta  + \phi )(\theta \phi  + 1)\left( {{\phi ^{{\tau _j}}} - 4{\phi ^{\frac{{{\tau _j}}}{2}}} + 3} \right) + (\phi  + 1)\left( { - (\theta  + \phi )(\theta \phi  + 1){\tau _j}\left( {{\phi ^{\frac{{{\tau _j}}}{2}}} - 2} \right){\phi ^{\frac{{{\tau _j}}}{2} - 1}} - \theta (\theta  + \phi )\left( {{\phi ^{{\tau _j}}} - 4{\phi ^{\frac{{{\tau _j}}}{2}}} + 3} \right) - (\theta \phi  + 1)\left( {{\phi ^{{\tau _j}}} - 4{\phi ^{\frac{{{\tau _j}}}{2}}} + 3} \right) - {{(\theta  + 1)}^2}\phi {\tau _j}} \right)} \right) \\ 
#' \end{array}  \right)}{d^2/dphi^2 nu[j]^2(phi, theta, sigma2) = (1/((phi - 1)^5*(phi + 1)^3*tau[j]^2))*(2*sigma2*((phi - 1)^2*
#'((phi + 1)^2*(theta^2*phi + theta*phi^2 + theta + phi)*tau[j]^2*
#'   (phi^(tau[j]/2) - 1)*phi^(tau[j]/2 - 2) + (phi^2 - 1)*(theta^2*(-phi) + theta*(phi^2 + 4*phi + 1) - phi)*tau[j]*(phi^(tau[j]/2) - 2)*phi^(tau[j]/2 - 2) - 
#'   2*(theta - 1)^2*(phi^tau[j] - 4*phi^(tau[j]/2) + 3))
#'- 12*(phi + 1)^2*(-((1/2)*(theta + 1)^2*(phi^2 - 1)*tau[j]) - (theta + phi)*(theta*phi + 1)*
#'                    (phi^tau[j] - 4*phi^(tau[j]/2) + 3)) + 6*(phi + 1)*(phi - 1)*((1/2)*(theta + 1)^2*(phi^2 - 1)*tau[j] + (theta + phi)*(theta*phi + 1)*(phi^tau[j] - 4*phi^(tau[j]/2) + 3) + 
#'                                                                                    (phi + 1)*(-((theta + phi)*(theta*phi + 1)*tau[j]*(phi^(tau[j]/2) - 2)*phi^(tau[j]/2 - 1)) - theta*(theta + phi)*(phi^tau[j] - 4*phi^(tau[j]/2) + 3) - 
#'                                                                                                 (theta*phi + 1)*(phi^tau[j] - 4*phi^(tau[j]/2) + 3) - (theta + 1)^2*phi*tau[j]))))}
#'
#' Taking the second derivative with respect to \eqn{\theta}{theta} yields:
#' \deqn{\frac{{{\partial ^2}}}{{\partial {\theta ^2}}}\nu _j^2\left( {\phi ,\theta ,{\sigma ^2}} \right) = \frac{{2{\sigma ^2}\left( {\left( {{\phi ^2} - 1} \right){\tau _j} + 2\phi \left( {{\phi ^{{\tau _j}}} - 4{\phi ^{\frac{{{\tau _j}}}{2}}} + 3} \right)} \right)}}{{{{(\phi  - 1)}^3}(\phi  + 1)\tau _j^2}} }{d^2/dtheta^2 nu[j]^2(phi, theta, sigma2) = (2*sigma2*(2*phi*(3 - 4*phi^(tau[j]/2) + phi^tau[j]) + (-1 + phi^2)*tau[j]))/((-1 + phi)^3*(1 + phi)*tau[j]^2)}
#' 
#' Taking the second derivative with respect to \eqn{\sigma ^2}{sigma^2} yields:
#' \deqn{\frac{{{\partial ^2}}}{{\partial {\sigma ^4}}}\nu _j^2\left( {\phi ,\theta ,{\sigma ^2}} \right) = 0 }{d^2/dsigma2^2 nu[j]^2(phi, theta, sigma2) = 0}
#'
#' Taking the derivative with respect to \eqn{\sigma^2}{sigma^2} and \eqn{\theta}{theta} yields:
#' 
#' \deqn{\frac{\partial }{{\partial \theta }}\frac{\partial }{{\partial {\sigma ^2}}}\nu _j^2\left( {\phi ,\theta ,{\sigma ^2}} \right) = \frac{2}{{{{(\phi  - 1)}^3}(\phi  + 1)\tau _j^2}}\left( {(\theta  + 1)\left( {{\phi ^2} - 1} \right){\tau _j} + \left( {2\theta \phi  + {\phi ^2} + 1} \right)\left( {{\phi ^{{\tau _j}}} - 4{\phi ^{\frac{{{\tau _j}}}{2}}} + 3} \right)} \right)}{(2*((theta + 1)*(phi^2 - 1)*tau + (2*theta*phi + phi^2 + 1)*
#' (phi^tau - 4*phi^(tau/2) + 3)))/((phi - 1)^3*(phi + 1)*tau^2)}
#'
#' Taking the derivative with respect to \eqn{\sigma^2}{sigma^2} and \eqn{\phi}{phi} yields:
#' \deqn{\frac{\partial }{{\partial \phi }}\frac{\partial }{{\partial {\sigma ^2}}}\nu _j^2\left( {\phi ,\theta ,{\sigma ^2}} \right) = \frac{2}{{{{(\phi  - 1)}^4}{{(\phi  + 1)}^2}\tau _j^2}}\left( \begin{array}{ll}
#' &- (\phi  - 1)(\phi  + 1)\left( \begin{array}{ll}
#'                                &- (\theta  + \phi )(\theta \phi  + 1){\tau _j}\left( {{\phi ^{\frac{{{\tau _j}}}{2}}} - 2} \right){\phi ^{\frac{{{\tau _j}}}{2} - 1}} \\
#'                                &- \theta (\theta  + \phi )\left( {{\phi ^{{\tau _j}}} - 4{\phi ^{\frac{{{\tau _j}}}{2}}} + 3} \right) \\
#'                                &- (\theta \phi  + 1)\left( {{\phi ^{{\tau _j}}} - 4{\phi ^{\frac{{{\tau _j}}}{2}}} + 3} \right) \\
#'                                &- {(\theta  + 1)^2}\phi {\tau _j} \\ 
#'                                \end{array}  \right) \\
#' &+ (\phi  - 1)\left( { - \frac{1}{2}{{(\theta  + 1)}^2}\left( {{\phi ^2} - 1} \right){\tau _j} - (\theta  + \phi )(\theta \phi  + 1)\left( {{\phi ^{{\tau _j}}} - 4{\phi ^{\frac{{{\tau _j}}}{2}}} + 3} \right)} \right) \\
#' &+ 3(\phi  + 1)\left( { - \frac{1}{2}{{(\theta  + 1)}^2}\left( {{\phi ^2} - 1} \right){\tau _j} - (\theta  + \phi )(\theta \phi  + 1)\left( {{\phi ^{{\tau _j}}} - 4{\phi ^{\frac{{{\tau _j}}}{2}}} + 3} \right)} \right) \\ 
#' \end{array}  \right)}{(1/((-1 + phi)^4*(1 + phi)^2*tau[j]^2))*
#' (2*((-1 + phi)*((-(theta + phi))*(1 + theta*phi)*(3 - 4*phi^(tau[j]/2) + phi^tau[j]) - (1/2)*(1 + theta)^2*(-1 + phi^2)*tau[j]) +
#'       3*(1 + phi)*((-(theta + phi))*(1 + theta*phi)*(3 - 4*phi^(tau[j]/2) + phi^tau[j]) - (1/2)*(1 + theta)^2*(-1 + phi^2)*tau[j]) - 
#'       (-1 + phi)*(1 + phi)*((-theta)*(theta + phi)*(3 - 4*phi^(tau[j]/2) + phi^tau[j]) -
#'                               (1 + theta*phi)*(3 - 4*phi^(tau[j]/2) + phi^tau[j]) - 
#'                               (1 + theta)^2*phi*tau[j] -
#'                               phi^(-1 + tau[j]/2)*(theta + phi)*(1 + theta*phi)*(-2 + phi^(tau[j]/2))*tau[j])))}
#' 
#' Taking the derivative with respect to \eqn{\phi}{phi} and \eqn{\theta}{theta} yields:
#' \deqn{\frac{\partial }{{\partial \theta }}\frac{\partial }{{\partial \phi }}\nu _j^2\left( {\phi ,\theta ,{\sigma ^2}} \right) =  - \frac{{2{\sigma ^2}}}{{{{(\phi  - 1)}^4}{{(\phi  + 1)}^2}\tau _j^2}}\left( \begin{array}{cc}
#' &{\tau _j}\left( \begin{array}{cc}
#'                 &2(\theta  + 1)(\phi  - 1){(\phi  + 1)^2} \\
#'                 &+ 2\left( {{\phi ^2} - 1} \right)\left( {2\theta \phi  + {\phi ^2} + 1} \right){\phi ^{\frac{{{\tau _j}}}{2} - 1}} \\
#'                 &- \left( {{\phi ^2} - 1} \right)\left( {2\theta \phi  + {\phi ^2} + 1} \right){\phi ^{{\tau _j} - 1}} \\ 
#'                 \end{array}  \right) \\
#' &+ 2\left( {\theta (\phi (3\phi  + 2) + 1) + \phi \left( {{\phi ^2} + \phi  + 3} \right) + 1} \right)\left( {{\phi ^{{\tau _j}}} - 4{\phi ^{\frac{{{\tau _j}}}{2}}} + 3} \right) \\ 
#' \end{array}  \right)}{d/dphi * d/dtheta nu[j]^2(phi, theta, sigma2) =(-(1/((-1 + phi)^4*(1 + phi)^2*tau[j]^2)))*2*sigma2*(2*(3 - 4*phi^(tau[j]/2) + phi^tau[j])*
#' (1 + phi*(3 + phi + phi^2) + 
#'   theta*(1 + phi*(2 + 3*phi))) + 
#'   (2*(1 + theta)*(-1 + phi)*(1 + phi)^2 + 
#'      2*phi^(tau[j]/2 - 1)*(-1 + phi^2)*
#'      (1 + 2*theta*phi + phi^2) - 
#'     phi^(tau[j] - 1)*(-1 + phi^2)*
#'     (1 + 2*theta*phi + phi^2))*tau[j])}
