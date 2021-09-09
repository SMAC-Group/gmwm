#' @section Process Definition:
#' The Autoregressive order \eqn{p} and Moving Average order \eqn{q} (ARMA(\eqn{p},\eqn{q})) process with non-zero parameters \eqn{\phi_i \in (-1,+1)}{phi[i] in (-1,1)} for the AR components,
#'  \eqn{\theta_j \in (-1,+1)}{theta[j] in (-1,+1)} for the MA components, and \eqn{\sigma^2 \in {\rm I\!R}^{+}}{sigma^2 in R^{+}}.
#' This process is defined as:
#' 
#' \deqn{{X_t} = \sum\limits_{i = 1}^p {{\phi _i}{X_{t - i}}}  + \sum\limits_{i = 1}^q {{\theta _i}{\varepsilon _{t - i}}}  + {\varepsilon _t}}{X[t] = sum(phi[p]*X[t-1]) + sum(theta[q]*W[t-1]) + W[t]}  
#'  where
#'   \deqn{{\varepsilon_t}\mathop  \sim \limits^{iid} N\left( {0,\sigma^2} \right)}{W[t] ~ N(0,sigma^2) iid}
