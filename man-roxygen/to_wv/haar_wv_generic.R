## Requires

## #' @templateVar mname Gaussian White Noise
## #' @templateVar sname WN
## #' @templateVar deqn_latex \nu _j^2\left( {{\sigma ^2}} \right) = \frac{{{\sigma ^2}}}{{\tau _j^2}}
## #' @templateVar deqn_ascii nu[j]^2 (sigma2) = (sigma2)/(tau[j]^2)
## #' @template to_wv/haar_wv_generic
## #' @template misc/haar_wv_formulae_link
## #' @backref src/process_to_wv.cpp
## #' @backref src/process_to_wv.h

#' @section Process Haar Wavelet Variance Formula:
<%= paste0("#' ", mname, " (", sname ,") process has a Haar Wavelet Variance given by:") %>
<%= paste0("#' \\deqn{", deqn_latex, "}{",deqn_ascii,"}") %>
