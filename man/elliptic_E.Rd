% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/elliptic_E.R
\name{elliptic_E}
\alias{elliptic_E}
\title{Incomplete elliptic integral of the second kind}
\usage{
elliptic_E(phi, m, minerror = 1e-15)
}
\arguments{
\item{phi}{amplitude, real or complex number/vector}

\item{m}{parameter, real or complex number/vector}

\item{minerror}{the bound on the relative error passed to
\code{\link{Carlson_RF}} and \code{\link{Carlson_RD}}}
}
\value{
A complex number or vector, the value(s) of the incomplete elliptic
  integral \ifelse{html}{\out{E(&phi;,m)}}{\eqn{E(\phi,m)}{E(phi,m)}}.
}
\description{
Evaluate the incomplete elliptic integral of the second kind.
}
\examples{
elliptic_E(1, 0.2)
gsl::ellint_E(1, sqrt(0.2))
}
