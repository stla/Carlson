#' Carlson elliptic integral RG
#' @description Evaluate the Carlson elliptic integral RG.
#'
#' @param x,y,z real or complex numbers; at least...
#' @param minerror bound of relative error
#'
#' @return A complex number.
#' @export
Carlson_RG <- function(x, y, z, minerror = 2*.Machine$double.eps){
  x <- as.complex(x); y <- as.complex(y); z <- as.complex(z)
  (z*Carlson_RF(x, y, z, minerror) -
     (x-z)*(y-z)*Carlson_RD(x, y, z, minerror)/3 +
      sqrt(x)*sqrt(y)/sqrt(z)) / 2
}
