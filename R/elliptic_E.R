#' Incomplete elliptic integral of the second kind
#' @description Evaluate the incomplete elliptic integral of the second kind.
#'
#' @param phi amplitude, real or complex number
#' @param m parameter, real or complex number
#' @param minerror the bound on the relative error passed to
#' \code{\link{Carlson_RF}} and \code{\link{Carlson_RD}}
#'
#' @return A complex number.
#' @export
#'
#' @examples elliptic_E(1, 0.2)
#' gsl::ellint_E(1, sqrt(0.2))
elliptic_E <- function(phi, m, minerror = 2*.Machine$double.eps){
  if(phi == 0){
    0
  }else if(phi == pi/2 && m == 1){
    1
  }else if(Re(phi) >= -pi/2 && Re(phi) <= pi/2){
    sine <- sin(phi)
    sine2 <- sine*sine
    cosine2 <- 1 - sine2
    oneminusmsine2 <- 1 - m*sine2
    sine * (Carlson_RF(cosine2, oneminusmsine2, 1, minerror) -
              m * sine2 * Carlson_RD(cosine2, oneminusmsine2, 1, minerror) / 3)
  }else if(Re(phi) > pi/2){
    k <- 0
    while(Re(phi) > pi/2){
      phi <- phi - pi
      k <- k + 1
    }
    2*k*elliptic_E(pi/2, m, minerror) + elliptic_E(phi, m, minerror)
  }else{
    k <- 0
    while(Re(phi) < -pi/2){
      phi <- phi + pi
      k <- k - 1
    }
    2*k*elliptic_E(pi/2, m, minerror) + elliptic_E(phi, m, minerror)
  }
}
