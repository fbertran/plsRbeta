\name{print.plsRbetamodel}
\alias{print.plsRbetamodel}
\title{Print method for plsRbeta models}
\description{
This function provides a print method for the class \code{"plsRbetamodel"}
}
\usage{
\method{print}{plsRbetamodel}(x, \dots)
}
\arguments{
  \item{x}{an object of the class \code{"plsRbetamodel"}}
  \item{\dots}{not used}
}
%\details{}
\value{
\code{NULL}
}
\references{
\enc{Frédéric}{Fr\'ed\'eric} Bertrand, Nicolas Meyer, \enc{Michèle}{Mich\`ele} Beau-Faller, Karim El Bayed, Izzie-Jacques Namer, Myriam Maumy-Bertrand (2013). \enc{Régression}{R\'egression} \enc{Bêta}{B\^eta} PLS. \emph{Journal de la \enc{Société}{Soci\'et'e} Française de Statistique}, \bold{154}(3):143-159.
\url{http://smf4.emath.fr/Publications/JSFdS/154_3/html/}
}
\author{\enc{Frédéric}{Fr\'ed\'eric} Bertrand\cr
\email{frederic.bertrand@math.unistra.fr}\cr
\url{http://www-irma.u-strasbg.fr/~fbertran/}
}
\seealso{\code{\link{print}}}
\examples{
data("GasolineYield",package="betareg")
modpls <- plsRbeta(yield~.,data=GasolineYield,nt=3,modele="pls-beta")
print(modpls)
}
\keyword{methods}
\keyword{print}
