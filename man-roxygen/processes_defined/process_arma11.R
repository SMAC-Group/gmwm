#' @section Process Definition:
#' The Autoregressive order 1 and Moving Average order 1 (ARMA(1,1)) process with non-zero parameters \eqn{\phi \in (-1,+1)}{phi in (-1,1)} for the AR component,
#'  \eqn{\theta \in (-1,+1)}{theta in (-1,+1)} for the MA component, and \eqn{\sigma^2 \in {\rm I\!R}^{+}}{sigma^2 in R^{+}}.
#' This process is defined as:
#'  \deqn{{X_t} = {\phi _1}{X_{t - 1}} + {\theta _1}{\varepsilon_{t - 1}} + {\varepsilon_t}}{W[t] = phi*W[t-1] + theta*W[t-1] + W[t]},
#'  where
#'   \deqn{{\varepsilon_t}\mathop  \sim \limits^{iid} N\left( {0,\sigma^2} \right)}{W[t] ~ N(0,sigma^2) iid}
