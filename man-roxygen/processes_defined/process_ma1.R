#' @section Process Definition:
#' The Moving Average order 1 (MA(1)) process with non-zero parameter \eqn{\theta \in (-1,+1)}{theta in (-1,+1)} 
#' and \eqn{\sigma^2 \in {\rm I\!R}^{+}}{sigma^2 in R^{+}}. This process is defined as:
#' \deqn{{x_t} = {\varepsilon_t} + {\theta _1}{\varepsilon_{t - 1}}}{x[t] = W[t] + theta*W[t-1]},
#'  where 
#'  \deqn{{\varepsilon_t}\mathop  \sim \limits^{iid} N\left( {0,\sigma^2} \right)}{W[t] ~ N(0,sigma^2) iid}
