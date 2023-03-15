#' Carlson elliptic integral RD
#' @description Evaluate the Carlson elliptic integral RD.
#'
#' @param x,y,z real or complex numbers; at most one can be 0
#' @param minerror bound on the relative error
#'
#' @return A complex number, the value of the Carlson elliptic integral
#' \ifelse{html}{\out{R<sub>D</sub>(x,y,z)}}{\eqn{R_D(x,y,z)}{RD(x,y,z)}}.
#' @export
#'
#' @note The function returns a value when \code{x}, \code{y} or \code{z}
#' are negative real numbers, but this value is not the one of the
#' Carlson integral.
#'
#' @examples Carlson_RD(5, 2, 3)
#' gsl::ellint_RD(5, 2, 3)
Carlson_RD <- function(x, y, z, minerror = 1e-15){
  stopifnot(minerror > 0)
  if(sum(c(x, y, z) == 0) > 1L){
    stop("At most one of `x`, `y`, `z` can be 0.")
  }
  x <- as.complex(x); y <- as.complex(y); z <- as.complex(z)
  Carlson_RD_(x, y, z, minerror)
}
