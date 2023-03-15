#' Carlson elliptic integral RJ
#' @description Evaluate the Carlson elliptic integral RJ.
#'
#' @param x,y,z,p real or complex numbers; at most one can be 0
#' @param minerror bound on the relative error
#'
#' @return A complex number, the value of the Carlson elliptic integral
#'   \ifelse{html}{\out{R<sub>J</sub>(x,y,z,t)}}{\eqn{R_J(x,y,z,t)}{RJ(x,y,z,t)}}.
#'
#' @export
#'
#' @note The function returns a value when \code{x}, \code{y}, \code{z} or
#'   \code{p} are negative real numbers, but this value is not the one of the
#'   Carlson integral.
#'
#' @examples Carlson_RJ(5, 2, 3, 4)
#' gsl::ellint_RJ(5, 2, 3, 4)
Carlson_RJ <- function(x, y, z, p, minerror = 1e-15){
  stopifnot(minerror > 0)
  if(sum(c(x, y, z, p) == 0) > 1L){
    stop("At most one of `x`, `y`, `z`, `p` can be 0.")
  }
  x <- as.complex(x); y <- as.complex(y); z <- as.complex(z); p <- as.complex(p)
  Carlson_RJ_(x, y, z, p, minerror)
}


