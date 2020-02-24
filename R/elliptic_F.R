#' Incomplete elliptic integral of the first kind
#' @description Evaluate the incomplete elliptic integral of the first kind.
#'
#' @param phi amplitude, real or complex number
#' @param m parameter, real or complex number
#' @param minerror the bound on the relative error passed to
#' \code{\link{Carlson_RF}}
#'
#' @return A complex number, the value of the incomplete elliptic integral
#' \ifelse{html}{\out{F(&phi;,m)}}{\eqn{F(\phi,m)}{F(phi,m)}}.
#' @export
#'
#' @examples elliptic_F(1, 0.2)
#' gsl::ellint_F(1, sqrt(0.2))
elliptic_F <- function(phi, m, minerror = 2*.Machine$double.eps){
  if(phi == 0 || m == Inf || m == -Inf){
    as.complex(0)
  }else if(Re(phi) == 0 && is.infinite(Im(phi)) && Im(m) == 0 && Re(m) > 0 &&
           Re(m) < 1){
    sign(Im(phi)) *
      (elliptic_F(pi/2,m,minerror) - elliptic_F(pi/2,1/m,minerror)/sqrt(m))
  }else if(abs(Re(phi)) == pi/2 && m == 1){
    NaN
  }else if(Re(phi) >= -pi/2 && Re(phi) <= pi/2){
    if(m == 1 && abs(Re(phi)) < pi/2){
      as.complex(asinh(tan(phi))) # or atanh(sin(phi))
    }else if(m == 0){
      as.complex(phi)
    }else{
      sine <- sin(phi) # sin(999i) = 0+Infi => pb sine2
      if(is.infinite(Re(sine)) || is.infinite(Im(sine))){
        stop("`sin(phi)` is not finite.")
      }
      sine2 <- sine*sine
      cosine2 <- 1 - sine2
      oneminusmsine2 <- 1 - m*sine2
      sine * Carlson_RF(cosine2, oneminusmsine2, 1, minerror)
    }
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
