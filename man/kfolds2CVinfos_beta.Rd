\name{kfolds2CVinfos_beta}
\alias{kfolds2CVinfos_beta}
\title{Extracts and computes information criteria and fits statistics for kfold cross validated partial least squares beta regression models}
\description{
This function extracts and computes information criteria and fits statistics for kfold cross validated partial least squares beta regression models for both formula or classic specifications of the model.
}
\usage{
kfolds2CVinfos_beta(pls_kfolds, MClassed = FALSE)
}
\arguments{
  \item{pls_kfolds}{an object computed using \code{\link{PLS_beta_kfoldcv}}}
  \item{MClassed}{should number of miss classed be computed}
}
\details{
The Mclassed option should only set to \code{TRUE} if the response is binary. 
}
\value{
  \item{list}{table of fit statistics for first group partition}
  \item{\dots}{\dots}
  \item{list}{table of fit statistics for last group partition}
}
\references{
\enc{Frédéric}{Fr\'ed\'eric} Bertrand, Nicolas Meyer, \enc{Michèle}{Mich\`ele} Beau-Faller, Karim El Bayed, Izzie-Jacques Namer, Myriam Maumy-Bertrand (2013). \enc{Régression}{R\'egression} \enc{Bêta}{B\^eta} PLS. \emph{Journal de la \enc{Société}{Soci\'et'e} Française de Statistique}, \bold{154}(3):143-159.
\url{http://publications-sfds.math.cnrs.fr/index.php/J-SFdS/article/view/215}
}
\author{\enc{Frédéric}{Fr\'ed\'eric} Bertrand\cr
\email{frederic.bertrand@math.unistra.fr}\cr
\url{http://www-irma.u-strasbg.fr/~fbertran/}
}
\seealso{\code{\link[plsRglm]{kfolds2coeff}}, \code{\link[plsRglm]{kfolds2Pressind}}, \code{\link[plsRglm]{kfolds2Press}}, \code{\link[plsRglm]{kfolds2Mclassedind}} and \code{\link[plsRglm]{kfolds2Mclassed}} to extract and transforms results from kfold cross validation.}
\examples{
\dontrun{
data("GasolineYield",package="betareg")
bbb <- PLS_beta_kfoldcv_formula(yield~.,data=GasolineYield,nt=3,modele="pls-beta")
kfolds2CVinfos_beta(bbb)
}
}
\keyword{models}
\keyword{regression}
