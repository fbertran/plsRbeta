\name{kfolds2Chisq}
\alias{kfolds2Chisq}
\title{Computes Predicted Chisquare for kfold cross validated partial least squares beta regression models.}
\description{This function computes Predicted Chisquare for kfold cross validated partial least squares beta regression models.}
\usage{
kfolds2Chisq(pls_kfolds)
}
\arguments{
  \item{pls_kfolds}{a kfold cross validated partial least squares regression glm model}
}
%\details{}
\value{
  \item{list}{Total Predicted Chisquare vs number of components for the first group partition}
  \item{\dots}{\dots}
  \item{list}{Total Predicted Chisquare vs number of components for the last group partition}
}
\references{
\enc{Frédéric}{Fr\'ed\'eric} Bertrand, Nicolas Meyer, \enc{Michèle}{Mich\`ele} Beau-Faller, Karim El Bayed, Izzie-Jacques Namer, Myriam Maumy-Bertrand (2013). \enc{Régression}{R\'egression} \enc{Bêta}{B\^eta} PLS. \emph{Journal de la \enc{Société}{Soci\'et'e} Française de Statistique}, \bold{154}(3):143-159.
\url{http://publications-sfds.math.cnrs.fr/index.php/J-SFdS/article/view/215}
}
\author{\enc{Frédéric}{Fr\'ed\'eric} Bertrand\cr
\email{frederic.bertrand@math.unistra.fr}\cr
\url{http://www-irma.u-strasbg.fr/~fbertran/}
}
\note{Use \code{\link{PLS_beta_kfoldcv}} to create kfold cross validated partial least squares regression glm and beta models.}
\seealso{\code{\link[plsRglm]{kfolds2coeff}}, \code{\link[plsRglm]{kfolds2Press}}, \code{\link[plsRglm]{kfolds2Pressind}}, \code{\link{kfolds2Chisqind}}, \code{\link[plsRglm]{kfolds2Mclassedind}} and \code{\link[plsRglm]{kfolds2Mclassed}} to extract and transforms results from kfold cross validation.}
\examples{
\dontrun{
data("GasolineYield",package="betareg")
yGasolineYield <- GasolineYield$yield
XGasolineYield <- GasolineYield[,2:5]
bbb <- PLS_beta_kfoldcv(yGasolineYield,XGasolineYield,nt=3,modele="pls-beta")
kfolds2Chisq(bbb)
}
}
\keyword{models}
\keyword{regression}
