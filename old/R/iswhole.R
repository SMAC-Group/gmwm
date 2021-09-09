#' @title Integer Check
#' @description Checks whether the submitted value is an integer
#' @keywords internal
#' @param x A \code{numeric} value.
#' @return A \code{boolean} value indicating whether the value is an integer or not.
#' @author James Balamuta
#' @export
#' @examples
#' is.whole(2.3)
#' is.whole(4)
#' is.whole(c(1,2,3))
#' is.whole(c(.4,.5,.6))
#' is.whole(c(7,.8,9))
is.whole = function(x){ is.numeric(x) && all(floor(x)==x) } 