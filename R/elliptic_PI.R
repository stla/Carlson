#' Incomplete elliptic integral of the third kind
#' @description Evaluate the incomplete elliptic integral of the third kind.
#'
#' @param phi amplitude, real or complex number
#' @param n characteristic, real or complex number
#' @param m parameter, real or complex number
#' @param minerror the bound on the relative error passed to
#'   \code{\link{Carlson_RF}} and \code{\link{Carlson_RJ}}
#'
#' @return A complex number, the value of the incomplete elliptic integral
#'   \ifelse{html}{\out{&Pi;(&phi;,n,m)}}{\eqn{\Pi(\phi,n,m)}{PI(phi,n,m)}}.
#' @export
#'
#' @examples elliptic_PI(1, 0.8, 0.2)
#' gsl::ellint_P(1, sqrt(0.2), -0.8)
elliptic_PI <- function(phi, n, m, minerror = 1e-15){
  if(phi == 0 || (is.infinite(Re(m)) && Im(m) == 0) ||
     (is.infinite(Re(n)) && Im(n) == 0)){
    as.complex(0)
  }else if(phi == pi/2 && m == 1 && Im(n) == 0 && n != 1){
    ifelse(Re(n) > 1, -Inf, Inf)
  }else if(phi == pi/2 && n == 1){
    NaN
  }else if(phi == pi/2 && m == 0){
    pi/2/sqrt(as.complex(1-n))
  }else if(phi == pi/2 && n == m){
    elliptic_E(pi/2, m, minerror) / (1-m)
  }else if(phi == pi/2 && n == 0){
    elliptic_F(pi/2, m, minerror)
  }else if(Re(phi) >= -pi/2 && Re(phi) <= pi/2){
    sine <- sin(phi)
    if(is.infinite(Re(sine)) || is.infinite(Im(sine))){
      stop("`sin(phi)` is not finite.")
    }
    sine2 <- sine*sine
    cosine2 <- 1 - sine2
    oneminusmsine2 <- 1 - m*sine2
    sine * (Carlson_RF(cosine2, oneminusmsine2, 1, minerror) +
              n * sine2 * Carlson_RJ(cosine2, oneminusmsine2, 1, 1-n*sine2,
                                     minerror) / 3)
  }else if(Re(phi) > pi/2){
    k <- ceiling(Re(phi)/pi - 0.5)
    phi <- phi - k*pi
    2*k*elliptic_PI(pi/2, n, m, minerror) + elliptic_PI(phi, n, m, minerror)
  }else{
    k <- -floor(0.5 - Re(phi)/pi)
    phi <- phi - k*pi
    2*k*elliptic_PI(pi/2, n, m, minerror) + elliptic_PI(phi, n, m, minerror)
  }
}
