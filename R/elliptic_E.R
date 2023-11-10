#' Incomplete elliptic integral of the second kind
#' @description Evaluate the incomplete elliptic integral of the second kind.
#'
#' @param phi amplitude, real or complex number/vector
#' @param m parameter, real or complex number/vector
#' @param minerror the bound on the relative error passed to
#'   \code{\link{Carlson_RF}} and \code{\link{Carlson_RD}}
#'
#' @return A complex number or vector, the value(s) of the incomplete elliptic
#'   integral \ifelse{html}{\out{E(&phi;,m)}}{\eqn{E(\phi,m)}{E(phi,m)}}.
#' @export
#'
#' @examples elliptic_E(1, 0.2)
#' gsl::ellint_E(1, sqrt(0.2))
elliptic_E <- function(phi, m, minerror = 1e-15){
  phi <- as.complex(phi)
  m   <- as.complex(m)
  if(length(phi) == 1L) {
    phi <- rep(phi, length(m))
  } else if(length(m) == 1L) {
    m <- rep(m, length(phi))
  } else if(length(phi) != length(m)) {
    stop("Incompatible lengths of `phi` and `m`.")
  }
  ellEcpp(phi, m, minerror)
  # if(phi == 0){
  #   as.complex(0)
  # }else if(is.infinite(Re(m)) && Im(m) == 0){
  #   NaN
  # }else if(Re(phi) >= -pi/2 && Re(phi) <= pi/2){
  #   if(m == 0){
  #     as.complex(phi)
  #   }else if(m == 1){
  #     sin(as.complex(phi))
  #   }else{
  #     sine <- sin(phi)
  #     if(is.infinite(Re(sine)) || is.infinite(Im(sine))){
  #       stop("`sin(phi)` is not finite.")
  #     }
  #     sine2 <- sine*sine
  #     cosine2 <- 1 - sine2
  #     oneminusmsine2 <- 1 - m*sine2
  #     sine * (Carlson_RF(cosine2, oneminusmsine2, 1, minerror) -
  #               m * sine2 * Carlson_RD(cosine2, oneminusmsine2, 1, minerror) / 3)
  #   }
  # }else if(Re(phi) > pi/2){
  #   k <- ceiling(Re(phi)/pi - 0.5)
  #   phi <- phi - k*pi
  #   2*k*elliptic_E(pi/2, m, minerror) + elliptic_E(phi, m, minerror)
  # }else{
  #   k <- -floor(0.5 - Re(phi)/pi)
  #   phi <- phi - k*pi
  #   2*k*elliptic_E(pi/2, m, minerror) + elliptic_E(phi, m, minerror)
  # }
}
