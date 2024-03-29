% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/elliptic_PI.R
\name{elliptic_PI}
\alias{elliptic_PI}
\title{Incomplete elliptic integral of the third kind}
\usage{
elliptic_PI(phi, n, m, minerror = 1e-15)
}
\arguments{
\item{phi}{amplitude, real or complex number/vector}

\item{n}{characteristic, real or complex number/vector}

\item{m}{parameter, real or complex number/vector}

\item{minerror}{the bound on the relative error passed to
\code{\link{Carlson_RF}} and \code{\link{Carlson_RJ}}}
}
\value{
A complex number or vector, the value(s) of the incomplete elliptic
  integral
  \ifelse{html}{\out{&Pi;(&phi;,n,m)}}{\eqn{\Pi(\phi,n,m)}{PI(phi,n,m)}}.
}
\description{
Evaluate the incomplete elliptic integral of the third kind.
}
\examples{
elliptic_PI(1, 0.8, 0.2)
gsl::ellint_P(1, sqrt(0.2), -0.8)
}
