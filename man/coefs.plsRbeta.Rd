\name{coefs.plsRbeta}
\alias{coefs.plsRbeta}
\title{Coefficients for bootstrap computations}
\description{
  ~~ A (1-5 lines) description of what the function does. ~~
}
\usage{
coefs.plsRbeta(dataset, ind, nt, modele, family=NULL, method="logistic", 
link=NULL, link.phi=NULL, type="ML")
}
\arguments{
  \item{dataset}{ ~~Explain \code{dataset} here~~ }
  \item{ind}{ ~~Explain \code{ind} here~~ }
  \item{nt}{ ~~Explain \code{nt} here~~ }
  \item{modele}{name of the PLS glm or PLS beta model to be fitted (\code{"pls"}, \code{"pls-glm-Gamma"}, \code{"pls-glm-gaussian"}, \code{"pls-glm-inverse.gaussian"}, \code{"pls-glm-logistic"}, \code{"pls-glm-poisson"}, \code{"pls-glm-polr"}, \code{"pls-beta"}). Use \code{"modele=pls-glm-family"} to enable the \code{family} option.}
  \item{family}{ ~~Explain \code{family} here~~ }
  \item{method}{ ~~Explain \code{method} here~~ }
  \item{link}{ ~~Explain \code{link} here~~ }
  \item{link.phi}{ ~~Explain \code{link.phi} here~~ }
  \item{type}{ ~~Explain \code{type} here~~ }
}
\details{
  ~~ More details than the description above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{
\enc{Frédéric}{Fr\'ed\'eric} Bertrand, Nicolas Meyer, \enc{Michèle}{Mich\`ele} Beau-Faller, Karim El Bayed, Izzie-Jacques Namer, Myriam Maumy-Bertrand (2013). \enc{Régression}{R\'egression} \enc{Bêta}{B\^eta} PLS. \emph{Journal de la \enc{Société}{Soci\'et'e} Française de Statistique}, \bold{154}(3):143-159.
\url{http://smf4.emath.fr/Publications/JSFdS/154_3/html/}
}
\author{\enc{Frédéric}{Fr\'ed\'eric} Bertrand\cr
\email{frederic.bertrand@math.unistra.fr}\cr
\url{http://www-irma.u-strasbg.fr/~fbertran/}
}
\note{ ~~some notes~~ 
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
data("GasolineYield",package="betareg")
}
\keyword{models}
