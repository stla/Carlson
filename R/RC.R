#' Carlson elliptic integral RC
#' @description Evaluate the Carlson elliptic integral RC.
#'
#' @param x,y real or complex numbers, with \code{y} different from \code{0}
#' @param minerror bound on the relative error passed to \code{\link{Carlson_RF}}
#'
#' @return A complex number, the value of the Carlson elliptic integral
#'   \ifelse{html}{\out{R<sub>C</sub>(x,y)}}{\eqn{R_C(x,y)}{RC(x,y)}}.
#' @export
#'
#' @note The function returns a value when \code{x} or \code{y}
#'   are negative real numbers, but this value is not the one of the
#'   Carlson integral.
#'
#' @examples Carlson_RC(5, 2)
#' gsl::ellint_RC(5, 2)
Carlson_RC <- function(x, y, minerror = 1e-15){
  if(y == 0) stop("`y` cannot be 0.")
  Carlson_RF_(x, y, y, minerror)
}
