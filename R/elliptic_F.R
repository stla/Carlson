#' Incomplete elliptic integral of the first kind
#' @description Evaluate the incomplete elliptic integral of the first kind.
#'
#' @param phi amplitude, real or complex number
#' @param m parameter, real or complex number
#' @param minerror the bound on the relative error passed to
#' \code{\link{Carlson_RF}}
#'
#' @return A complex number.
#' @export
#'
#' @examples elliptic_F(1, 0.2)
#' gsl::ellint_F(1, sqrt(0.2))
elliptic_F <- function(phi, m, minerror = .Machine$double.eps){
  if(phi == 0){
    0
  }else if(phi == pi/2 && m == 1){
    Inf
  }else if(Re(phi) >= -pi/2 && Re(phi) <= pi/2){
    sine <- sin(phi)
    sine2 <- sine*sine
    cosine2 <- 1 - sine2
    oneminusmsine2 <- 1 - m*sine2
    sine * Carlson_RF(cosine2, oneminusmsine2, 1, minerror)
  }else if(Re(phi) > pi/2){
    k <- 0
    while(Re(phi) > pi/2){
      phi <- phi - pi
      k <- k + 1
    }
    2*k*elliptic_F(pi/2, m, minerror) + elliptic_F(phi, m, minerror)
  }else{
    k <- 0
    while(Re(phi) < -pi/2){
      phi <- phi + pi
      k <- k - 1
    }
    2*k*elliptic_F(pi/2, m, minerror) + elliptic_F(phi, m, minerror)
  }
}
