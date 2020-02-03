#' Carlson elliptic integral RC
#' @description Evaluate the Carlson elliptic integral RC.
#'
#' @param x,y real or complex numbers; at least...
#' @param minerror bound of relative error
#'
#' @return A complex number.
#' @export
#'
#' @examples Carlson_RC(5, 2)
#' gsl::ellint_RC(5, 2)
Carlson_RC <- function(x, y, minerror = .Machine$double.eps){
  Carlson_RF(x, y, y, minerror)
}
