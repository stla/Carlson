#' Jacobi zeta function
#' @description Evaluate the Jacobi zeta function.
#'
#' @param phi amplitude, real or complex number
#' @param m parameter, real or complex number
#' @param minerror bound on relative error passed to \code{\link{EllipticE}}
#' and \code{\link{EllipticF}}
#'
#' @return A complex number.
#' @export
elliptic_Z <- function(phi, m, minerror = 2*.Machine$double.eps){
  if(m == 1){
    if(abs(Re(phi)) <= pi/2){
      sin(as.complex(phi))
    }else if(Re(phi) > pi/2){
      while(Re(phi) > pi/2){
        phi <- phi - pi
      }
      sin(as.complex(phi))
    }else{
      while(Re(phi) < pi/2){
        phi <- phi + pi
      }
      sin(as.complex(phi))
    }
  }else{
    elliptic_E(phi, m, minerror) -
      elliptic_E(pi/2, m, minerror)/elliptic_F(pi/2, m, minerror) *
      elliptic_F(phi, m, minerror)
  }
}
