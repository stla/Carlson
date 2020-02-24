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
Carlson_RD <- function(x, y, z, minerror = 2*.Machine$double.eps){
  stopifnot(minerror > 0)
  if(sum(c(x,y,z)==0) > 1){
    stop("At most one of `x`, `y`, `z` can be 0.")
  }
  x <- as.complex(x); y <- as.complex(y); z <- as.complex(z)
  dx <- dy <- dz <- Inf
  s <- 0
  fac <- 1
  while(max(dx,dy,dz) > minerror){
    lambda <- sqrt(x)*sqrt(y) + sqrt(y)*sqrt(z) + sqrt(z)*sqrt(x)
    s <- s + fac / (sqrt(z) * (z + lambda))
    fac <- fac/4
    x <- (x + lambda) / 4
    y <- (y + lambda) / 4
    z <- (z + lambda) / 4
    A <- (x + y + 3*z) / 5
    dx <- Mod(1 - x/A)
    dy <- Mod(1 - y/A)
    dz <- Mod(1 - z/A)
  }
  E2 <- dx*dy + dy*dz + 3*dz*dz + 2*dz*dx + dx*dz + 2*dy*dz
  E3 <- dz*dz*dz + dx*dz*dz + 3*dx*dy*dz + 2*dy*dz*dz + dy*dz*dz + 2*dx*dz*dz
  E4 <- dy*dz*dz*dz + dx*dz*dz*dz + dx*dy*dz*dz + 2*dx*dy*dz*dz
  E5 <- dx*dy*dz*dz*dz
  3*s +
    fac * (1 - 3*E2/14 + E3/6 + 9*E2*E2/88 - 3*E4/22 - 9*E2*E3/52 + 3*E5/26 -
             E2*E2*E2/16 + 3*E3*E3/40 + 3*E2*E4/20 + 45*E2*E2*E3/272 -
             9*(E3*E4 + E2*E5)/68) / A / sqrt(A)
}
