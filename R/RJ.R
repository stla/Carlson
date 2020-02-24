#' Carlson elliptic integral RJ
#' @description Evaluate the Carlson elliptic integral RJ.
#'
#' @param x,y,z,p real or complex numbers; at most one can be 0
#' @param minerror bound on the relative error
#'
#' @return A complex number, the value of the Carlson elliptic integral
#' \ifelse{html}{\out{R<sub>J</sub>(x,y,z,t)}}{\eqn{R_J(x,y,z,t)}{RJ(x,y,z,t)}}.
#'
#' @export
#'
#' @note The function returns a value when \code{x}, \code{y}, \code{z} or
#' \code{p} are negative real numbers, but this value is not the one of the
#' Carlson integral.
#'
#' @examples Carlson_RJ(5, 2, 3, 4)
#' gsl::ellint_RJ(5, 2, 3, 4)
Carlson_RJ <- function(x, y, z, p, minerror = 2*.Machine$double.eps){
  stopifnot(minerror > 0)
  if(sum(c(x,y,z,p)==0) > 1){
    stop("At most one of `x`, `y`, `z`, `p` can be 0.")
  }
  x <- as.complex(x); y <- as.complex(y); z <- as.complex(z); p <- as.complex(p)
  A0 <- A <- (x + y + z + p + p) / 5
  delta <- (p-x)*(p-y)*(p-z)
  f <- fac <- 1
  d <- e <- c()
  Q <- (4/minerror)^(1/6) * max(Mod(A-x), Mod(A-y), Mod(A-z), Mod(A-p))
  while(Mod(A) <= Q){
    dnew <- (sqrt(p)+sqrt(x))*(sqrt(p)+sqrt(y))*(sqrt(p)+sqrt(z))
    d <- c(d, dnew*f)
    e <- c(e, fac * delta / dnew / dnew)
    f <- f * 4
    fac <- fac / 64
    lambda <- sqrt(x)*sqrt(y) + sqrt(y)*sqrt(z) + sqrt(z)*sqrt(x)
    x <- (x + lambda) / 4
    y <- (y + lambda) / 4
    z <- (z + lambda) / 4
    p <- (p + lambda) / 4
    A <- (A + lambda) / 4
    Q <- Q / 4
  }
  X <- (A0-x) / f / A
  Y <- (A0-y) / f / A
  Z <- (A0-z) / f / A
  P <- -(X+Y+Z) / 2
  E2 <- X*Y + X*Z + Y*Z - 3*P*P
  E3 <- X*Y*Z + 2*E2*P + 4*P*P*P
  E4 <- P*(2*X*Y*Z + E2*P + 3*P*P*P)
  E5 <- X*Y*Z*P*P

  (1 - 3*E2/14 + E3/6 + 9*E2*E2/88 - 3*E4/22 - 9*E2*E3/52 + 3*E5/26) /
    f / A / sqrt(A) +
    ifelse(length(e),
           6*sum(ifelse(e == 0, 1, atan(sqrt(e))/sqrt(e)) / d), 0)
}


