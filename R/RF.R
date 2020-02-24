#' Carlson elliptic integral RF
#' @description Evaluate the Carlson elliptic integral RF.
#'
#' @param x,y,z real or complex numbers; at most one can be 0
#' @param minerror bound on relative error
#'
#' @return A complex number, the value of the Carlson elliptic integral
#' \ifelse{html}{\out{R<sub>F</sub>(x,y,z)}}{\eqn{R_F(x,y,z)}{RF(x,y,z)}}.
#' @export
#'
#' @note The function returns a value when \code{x}, \code{y} or \code{z}
#' are negative real numbers, but this value is not the one of the
#' Carlson integral.
#'
#' @examples Carlson_RF(5, 2, 3)
#' gsl::ellint_RF(5, 2, 3)
Carlson_RF <- function(x, y, z, minerror = 2*.Machine$double.eps){
  stopifnot(minerror > 0)
  if(sum(c(x,y,z)==0) > 1){
    stop("At most one of `x`, `y`, `z` can be 0.")
  }
  x <- as.complex(x); y <- as.complex(y); z <- as.complex(z)
  dx <- dy <- dz <- Inf
  while(max(dx,dy,dz) > minerror){
    lambda <- sqrt(x)*sqrt(y) + sqrt(y)*sqrt(z) + sqrt(z)*sqrt(x)
    x <- (x + lambda) / 4
    y <- (y + lambda) / 4
    z <- (z + lambda) / 4
    A <- (x+y+z) / 3
    dx <- Mod(1 - x/A)
    dy <- Mod(1 - y/A)
    dz <- Mod(1 - z/A)
  }
  E2 <- dx*dy + dy*dz + dz*dx
  E3 <- dy*dx*dz
  (1 - E2/10 + E3/14 + E2*E2/24 - 3*E2*E3/44 - 5*E2*E2*E2/208 +
      3*E3*E3/104 + E2*E2*E3/16) / sqrt(A)
}
