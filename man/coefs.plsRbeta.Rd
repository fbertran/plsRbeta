% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coefs.plsRbeta.R
\name{coefs.plsRbeta}
\alias{coefs.plsRbeta}
\title{Coefficients function for bootstrap techniques}
\usage{
coefs.plsRbeta(
  dataset,
  ind,
  nt,
  modele,
  family = NULL,
  method = "logistic",
  link = NULL,
  link.phi = NULL,
  type = "ML",
  verbose = TRUE
)
}
\arguments{
\item{dataset}{dataset to resample}

\item{ind}{indices for resampling}

\item{nt}{number of components to use}

\item{modele}{name of the PLS glm or PLS beta model to be fitted
(\code{"pls"}, \code{"pls-glm-Gamma"}, \code{"pls-glm-gaussian"},
\code{"pls-glm-inverse.gaussian"}, \code{"pls-glm-logistic"},
\code{"pls-glm-poisson"}, \code{"pls-glm-polr"}, \code{"pls-beta"}). Use
\code{"modele=pls-glm-family"} to enable the \code{family} option.}

\item{family}{family to use if GLM model, see \link{plsRbeta}}

\item{method}{method for beta regression}

\item{link}{link for beta regression}

\item{link.phi}{link.phi for beta regression}

\item{type}{type of estimates}

\item{verbose}{should info messages be displayed ?}
}
\value{
Coefficients' Estimates on a sample.
}
\description{
Returns the coefficients of a \code{"plsRbeta"} model.
}
\examples{
\donttest{
data("GasolineYield",package="betareg")
modpls <- coefs.plsRbeta(GasolineYield[,-6],1:32,nt=3,modele="pls-beta")
}

}
\references{
Frédéric Bertrand, Nicolas Meyer,
Michèle Beau-Faller, Karim El Bayed, Izzie-Jacques Namer,
Myriam Maumy-Bertrand (2013). Régression Bêta
PLS. \emph{Journal de la Société Française de Statistique},
\bold{154}(3):143-159.
\url{http://publications-sfds.math.cnrs.fr/index.php/J-SFdS/article/view/215}
}
\seealso{
See also \code{\link{bootplsbeta}}.
}
\author{
Frédéric Bertrand\cr
\email{frederic.bertrand@math.unistra.fr}\cr
\url{http://www-irma.u-strasbg.fr/~fbertran/}
}
\keyword{models}
