#' @title Heuman Lambda function
#' @description Evaluates the Heuman Lambda function.
#'
#' @param phi Jacobi amplitude, a complex number/vector
#' @param m parameter, a complex number/vector
#' @param minerror the bound on the relative error passed to
#'   \code{\link{elliptic_F}} and \code{\link{elliptic_Z}}
#'
#' @return A complex number or vector.
#' @export
Lambda0 <- function(phi, m, minerror = 1e-14) {
  if(m == 0) {
    elliptic_Z(phi, 1, minerror = minerror)
  } else {
    elliptic_F(phi, 1-m, minerror = minerror) /
      elliptic_F(pi/2, 1-m, minerror = minerror) +
      2/pi * elliptic_F(pi/2, m, minerror = minerror) *
      elliptic_Z(phi, 1-m, minerror = minerror)
  }
}
