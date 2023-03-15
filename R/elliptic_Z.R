#' Jacobi zeta function
#' @description Evaluate the Jacobi zeta function.
#'
#' @param phi amplitude, real or complex number
#' @param m parameter, real or complex number
#' @param minerror bound on relative error passed to \code{\link{elliptic_E}}
#'   and \code{\link{elliptic_F}}
#'
#' @return A complex number, the value of the Jacobi zeta function
#'   \ifelse{html}{\out{Z(&phi;,m)}}{\eqn{Z(\phi,m)}{Z(phi,m)}}.
#' @export
elliptic_Z <- function(phi, m, minerror = 1e-15){
  if(is.infinite(Re(m)) && Im(m) == 0){
    NaN
  }else if(m == 1){
    if(abs(Re(phi)) <= pi/2){
      sin(as.complex(phi))
    }else if(Re(phi) > pi/2){
      k <- ceiling(Re(phi)/pi - 0.5)
      phi <- phi - k*pi
      sin(as.complex(phi))
    }else{
      k <- -floor(0.5 - Re(phi)/pi)
      phi <- phi - k*pi
      sin(as.complex(phi))
    }
  }else{
    elliptic_E(phi, m, minerror) -
      elliptic_E(pi/2, m, minerror)/elliptic_F(pi/2, m, minerror) *
      elliptic_F(phi, m, minerror)
  }
}
