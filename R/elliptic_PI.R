#' Incomplete elliptic integral of the third kind
#' @description Evaluate the incomplete elliptic integral of the third kind.
#'
#' @param phi amplitude, real or complex number
#' @param n,m parameters, real or complex numbers
#' @param minerror the bound on the relative error passed to
#' \code{\link{Carlson_RF}} and \code{\link{Carlson_RJ}}
#'
#' @return A complex number.
#' @export
#'
#' @examples elliptic_PI(1, 0.8, 0.2)
#' gsl::ellint_P(1, sqrt(0.2), -0.8)
elliptic_PI <- function(phi, n, m, minerror = 2*.Machine$double.eps){
  if(phi == 0){
    0
  }else if(phi == pi/2 && m == 1){
    ifelse(n >= 1, -Inf, Inf)
  }else if(phi == pi/2 && n == 1){
    complex(real = Inf, imaginary = -Inf)
  }else if(Re(phi) >= -pi/2 && Re(phi) <= pi/2){
    sine <- sin(phi)
    sine2 <- sine*sine
    cosine2 <- 1 - sine2
    oneminusmsine2 <- 1 - m*sine2
    sine * (Carlson_RF(cosine2, oneminusmsine2, 1, minerror) +
              n * sine2 * Carlson_RJ(cosine2, oneminusmsine2, 1, 1-n*sine2,
                                     minerror) / 3)
  }else if(Re(phi) > pi/2){
    k <- 0
    while(Re(phi) > pi/2){
      phi <- phi - pi
      k <- k + 1
    }
    2*k*elliptic_PI(pi/2, n, m, minerror) + elliptic_PI(phi, n, m, minerror)
  }else{
    k <- 0
    while(Re(phi) < -pi/2){
      phi <- phi + pi
      k <- k - 1
    }
    2*k*elliptic_PI(pi/2, n, m, minerror) + elliptic_PI(phi, n, m, minerror)
  }
}
