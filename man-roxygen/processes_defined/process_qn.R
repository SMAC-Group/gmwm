#' @section Process Definition: 
#' Quantization Noise (QN) with parameter \eqn{Q^2 \in R^{+}}{Q^2 in R^{+}}. 
#' With i.i.d \eqn{Y_t \sim U(0,1)} (i.e. a standard uniform variable), this process is
#' defined as:
#' 
#' \deqn{X_t = \sqrt{12Q^2}(Y_{t}-Y_{t-1})}{X_t = sqrt(12*Q^2)*(Y[t]-Y[t-1])}
