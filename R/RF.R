#' Carlson elliptic integral RF
#' @description Evaluate the Carlson elliptic integral RF.
#'
#' @param x,y,z real or complex numbers; at least...
#' @param minerror minimum relative error
#'
#' @return A complex number.
#' @export
#'
#' @examples Carlson_RF(5, 2, 3)
#' gsl::ellint_RF(5, 2, 3)
Carlson_RF <- function(x, y, z, minerror = .Machine$double.eps){
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
