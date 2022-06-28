library(Carlson)
library(rgl)
library(Rvcg)

mesh <- vcgSphere(subdivision = 8)

pow <- function(x, n) x^n
color <- apply(mesh$vb[-4L, ], 2L, function(xyz){
  if(sum(xyz == 0) >= 2){
    z <- as.matrix(NA_complex_)
  }else{
    a <- xyz[1]
    b <- xyz[2]
    c <- xyz[3]
    z <- as.matrix(Carlson_RJ(a, b, c, 1i, 1e-5))
  }
  return(RcppColors::colorMap1(z))
  r <- 38.0*r + 2748.0*pow(r,2.0) - 9193.0*pow(r,3.0) + 15152.0*pow(r,4.0)- 11011.0*pow(r,5.0) + 824.0*pow(r,6.0) + 1690.0*pow(r,7.0)
  g <- 246.0*g - 2744.0*pow(g,2.0) + 15047.0*pow(g,3.0) - 39169.0*pow(g,4.0) + 53844.0*pow(g,5.0) - 36715.0*pow(g,6.0) + 9746.0*pow(g,7.0)
  b <- 5.0 + 206.0*b + 7444.0*pow(b,2.0) - 53422.0*pow(b,3.0) + 157704.0*pow(b,4.0) - 238709.0*pow(b,5.0) + 179985.0*pow(b,6.0) - 53051.0*pow(b,7.0)
  rgb(r, g, b, maxColorValue = 255)
})
mesh$material <- list(color = color)

open3d(windowRect = c(50, 50, 562, 562), zoom=0.75)
bg3d("whitesmoke")
shade3d(mesh)

# # -- if you want an animation
M <- par3d("userMatrix")
movie3d(
  par3dinterp(
    time = seq(0, 1, len = 9),
    userMatrix = list(
      M,
      rotate3d(M, pi, 1, 0, 0),
      rotate3d(M, pi, 1, 1, 0),
      rotate3d(M, pi, 1, 1, 1),
      rotate3d(M, pi, 0, 1, 1),
      rotate3d(M, pi, 0, 1, 0),
      rotate3d(M, pi, 1, 0, 1),
      rotate3d(M, pi, 0, 0, 1),
      M
    )
  ),
  fps = 120,
  duration = 1,
  dir = ".",
  movie = "zzpic",
  convert = FALSE,
  webshot = FALSE
)

command <- "gifski --fps=9 --frames=zzpic*.png -o CarlsonBall.gif"
system(command)

pngfiles <- list.files(pattern = "^zzpic?.*png$")
file.remove(pngfiles)


