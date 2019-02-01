\name{permcoefs.plsRbeta}
\alias{permcoefs.plsRbeta}
\title{Coefficients function for permutation bootstrap techniques}
\description{
A function passed to \code{boot} to perform bootstrap.
}
\usage{
permcoefs.plsRbeta(dataset, ind, nt, modele, family=NULL, method="logistic",
link="logit",link.phi=NULL,type="ML",verbose=TRUE)
}
\arguments{
  \item{dataset}{dataset to resample}
  \item{ind}{indices for resampling}
  \item{nt}{number of components to use}
  \item{modele}{name of the PLS glm or PLS beta model to be fitted (\code{"pls"}, \code{"pls-glm-Gamma"}, \code{"pls-glm-gaussian"}, \code{"pls-glm-inverse.gaussian"}, \code{"pls-glm-logistic"}, \code{"pls-glm-poisson"}, \code{"pls-glm-polr"}, \code{"pls-beta"}). Use \code{"modele=pls-glm-family"} to enable the \code{family} option.}
  \item{family}{family to use if GLM model, see \link{plsRbeta}}
  \item{method}{method for beta regression}
  \item{link}{link for beta regression}
  \item{link.phi}{link.phi for beta regression}
  \item{type}{type of estimates}
  \item{verbose}{should info messages be displayed ?}
}
\value{Estimates on a bootstrap sample.}
\references{
\enc{Frédéric}{Fr\'ed\'eric} Bertrand, Nicolas Meyer, \enc{Michèle}{Mich\`ele} Beau-Faller, Karim El Bayed, Izzie-Jacques Namer, Myriam Maumy-Bertrand (2013). \enc{Régression}{R\'egression} \enc{Bêta}{B\^eta} PLS. \emph{Journal de la \enc{Société}{Soci\'et'e} Française de Statistique}, \bold{154}(3):143-159.
\url{http://publications-sfds.math.cnrs.fr/index.php/J-SFdS/article/view/215}
}
\author{\enc{Frédéric}{Fr\'ed\'eric} Bertrand\cr
\email{frederic.bertrand@math.unistra.fr}\cr
\url{http://www-irma.u-strasbg.fr/~fbertran/}
}
\seealso{See also \code{\link{bootplsbeta}}.}
\examples{\donttest{
data("GasolineYield",package="betareg")

GazYield.boot <- bootplsbeta(plsRbeta(yield~.,data=GasolineYield,nt=3,
modele="pls-beta", verbose=FALSE), sim="ordinary", stype="i", R=250, statistic=permcoefs.plsRbeta)
}
}\keyword{models}
