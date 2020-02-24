#' Incomplete elliptic integral of the second kind
#' @description Evaluate the incomplete elliptic integral of the second kind.
#'
#' @param phi amplitude, real or complex number
#' @param m parameter, real or complex number
#' @param minerror the bound on the relative error passed to
#' \code{\link{Carlson_RF}} and \code{\link{Carlson_RD}}
#'
#' @return A complex number, the value of the incomplete elliptic integral
#' \ifelse{html}{\out{E(&phi;,m)}}{\eqn{E(\phi,m)}{E(phi,m)}}.
#' @export
#'
#' @examples elliptic_E(1, 0.2)
#' gsl::ellint_E(1, sqrt(0.2))
elliptic_E <- function(phi, m, minerror = 2*.Machine$double.eps){
  if(phi == 0){
    as.complex(0)
  }else if(is.infinite(Re(m)) && Im(m) == 0){
    NaN
  }else if(Re(phi) >= -pi/2 && Re(phi) <= pi/2){
    if(m == 0){
      as.complex(phi)
    }else if(m == 1){
      sin(as.complex(phi))
    }else{
      sine <- sin(phi)
      if(is.infinite(Re(sine)) || is.infinite(Im(sine))){
        stop("`sin(phi)` is not finite.")
      }
      sine2 <- sine*sine
      cosine2 <- 1 - sine2
      oneminusmsine2 <- 1 - m*sine2
      sine * (Carlson_RF(cosine2, oneminusmsine2, 1, minerror) -
                m * sine2 * Carlson_RD(cosine2, oneminusmsine2, 1, minerror) / 3)
    }
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
