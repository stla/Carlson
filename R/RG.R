#' Carlson elliptic integral RG
#' @description Evaluate the Carlson elliptic integral RG.
#'
#' @param x,y,z real or complex numbers; at most one can be zero
#' @param minerror bound on the relative error passed to
#' \code{\link{Carlson_RF}} and \code{\link{Carlson_RD}}
#'
#' @return A complex number.
#' @export
Carlson_RG <- function(x, y, z, minerror = 2*.Machine$double.eps){
  if(sum(c(x,y,z)==0) > 1){
    stop("At most one of `x`, `y`, `z` can be 0.")
  }
  x <- as.complex(x); y <- as.complex(y); z <- as.complex(z)
  (z*Carlson_RF(x, y, z, minerror) -
     (x-z)*(y-z)*Carlson_RD(x, y, z, minerror)/3 +
      sqrt(x)*sqrt(y)/sqrt(z)) / 2
}
