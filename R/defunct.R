#' @title Defunct function(s) in the gmwm package
#' @description These functions have been removed from the gmwm package.
#' @rdname gmwm-defunct
#' @name gmwm-defunct
#' @param ... Grab old parameters
#' @aliases compare.gmwm compare.models compare.wvar gen.gts gen.lts gmwm.imu  auto.imu rank.models 
#' @section Details:
#' \tabular{rl}{
#'   \code{compare.gmwm} \tab has been removed in favor of \code{\link{compare_models}}\cr
#'   \code{compare.models} \tab has been removed in favor of \code{\link{compare_models}}\cr
#'   \code{compare.wvar} \tab has been removed in favor of \code{\link{compare_wvar}}\cr
#'   \code{gen.gts} \tab has been removed in favor of \code{\link{gen_gts}}\cr
#'   \code{gen.lts} \tab has been removed in favor of \code{\link{gen_lts}}\cr
#'   \code{demo.lts} \tab has been removed in favor of \code{\link{demo_lts}}\cr
#'   \code{gmwm.imu} \tab has been removed in favor of \code{\link{gmwm_imu}}\cr
#'   \code{rank.models} \tab has been removed in favor of \code{\link{rank_models}}\cr
#'   \code{auto.imu} \tab has been removed in favor of \code{\link{auto_imu}}\cr
#' }
compare.gmwm = function(...){
  .Defunct("compare_models", package="gmwm")
}

#' @rdname gmwm-defunct
compare.wvar = function(...){
  .Defunct("compare_wvar", package="gmwm")
}

#' @rdname gmwm-defunct
compare.models = function(...){
  .Defunct("compare_models", package="gmwm")
}

#' @rdname gmwm-defunct
gen.gts = function(...){
  .Defunct("gen_gts", package="gmwm")
}

#' @rdname gmwm-defunct
gen.lts = function(...){
  .Defunct("gen_lts", package="gmwm")
}

#' @rdname gmwm-defunct
demo.lts = function(...){
  .Defunct("demo_lts", package="gmwm")
}

#' @rdname gmwm-defunct
gmwm.imu = function(...){
  .Defunct("gmwm_imu", package="gmwm")
}

#' @rdname gmwm-defunct
rank.models = function(...){
  .Defunct("rank_models", package="gmwm")
}

#' @rdname gmwm-defunct
auto.imu = function(...){
  .Defunct("auto_imu", package="gmwm")
}
