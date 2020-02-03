#' Carlson elliptic integral RC
#' @description Evaluate the Carlson elliptic integral RC.
#'
#' @param x,y real or complex numbers, with \code{y} different from \code{0}
#' @param minerror bound of relative error
#'
#' @return A complex number.
#' @export
#'
#' @examples Carlson_RC(5, 2)
#' gsl::ellint_RC(5, 2)
Carlson_RC <- function(x, y, minerror = 2*.Machine$double.eps){
  if(y == 0) stop("`y` cannot be 0.")
  Carlson_RF(x, y, y, minerror)
}
