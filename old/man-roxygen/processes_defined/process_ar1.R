#' @section Process Definition:
#' The Autoregressive order 1 (AR1) process with non-zero parameter \eqn{\phi \in (-1,+1)}{phi in (-1,1)} and \eqn{\sigma^2 \in {\rm I\!R}^{2}}{sigma^2 in R^{+}}.
#' This process is defined as: 
#' \deqn{{X_t} = {\phi _1}{X_{t - 1}} + {\varepsilon_t} }{X[t] = phi[1]X[t-1]  + W[t]},
#'  where \deqn{{\varepsilon_t}\mathop  \sim \limits^{iid} N\left( {0,\sigma^2} \right)}{W[t] ~ N(0,sigma^2) iid}
#' AR(1) processes are sometimes used as an approximation for Bias Instability noises.
