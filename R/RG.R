#' Carlson elliptic integral RG
#' @description Evaluate the Carlson elliptic integral RG.
#'
#' @param x,y,z real or complex numbers; they can be zero
#' @param minerror bound on the relative error passed to
#' \code{\link{Carlson_RF}} and \code{\link{Carlson_RD}}
#'
#' @return A complex number, the value of the Carlson elliptic integral
#' \ifelse{html}{\out{R<sub>G</sub>(x,y,z)}}{\eqn{R_G(x,y,z)}{RG(x,y,z)}}.
#' @export
Carlson_RG <- function(x, y, z, minerror = 2*.Machine$double.eps){
  zeros <- sum(c(x,y,z)==0)
  if(zeros == 3L){
    return(0i)
  }
  if(zeros == 2L){
    nonzero <- which(c(x,y,z) != 0)
    return(sqrt(as.complex(c(x,y,z)[nonzero]))/2)
  }
  if(z == 0){
    return(Carlson_RG(y, z, x, minerror))
  }
  x <- as.complex(x); y <- as.complex(y); z <- as.complex(z)
  (z*Carlson_RF(x, y, z, minerror) -
     (x-z)*(y-z)*Carlson_RD(x, y, z, minerror)/3 +
      sqrt(x)*sqrt(y)/sqrt(z)) / 2
}
