#' Carlson elliptic integral RF
#' @description Evaluate the Carlson elliptic integral RF.
#'
#' @param x,y,z real or complex numbers; at most one can be 0
#' @param minerror bound on relative error
#'
#' @return A complex number, the value of the Carlson elliptic integral
#'   \ifelse{html}{\out{R<sub>F</sub>(x,y,z)}}{\eqn{R_F(x,y,z)}{RF(x,y,z)}}.
#' @export
#'
#' @note The function returns a value when \code{x}, \code{y} or \code{z}
#'   are negative real numbers, but this value is not the one of the
#'   Carlson integral.
#'
#' @examples Carlson_RF(5, 2, 3)
#' gsl::ellint_RF(5, 2, 3)
Carlson_RF <- function(x, y, z, minerror = 1e-15){
  stopifnot(minerror > 0)
  if(sum(c(x, y, z) == 0) > 1L){
    stop("At most one of `x`, `y`, `z` can be 0.")
  }
  x <- as.complex(x); y <- as.complex(y); z <- as.complex(z)
  Carlson_RF_(x, y, z, minerror)
}
