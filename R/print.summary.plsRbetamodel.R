#' Print method for summaries of plsRbeta models
#' 
#' This function provides a print method for the class
#' \code{"summary.plsRbetamodel"}
#' 
#' 
#' @param x an object of the class \code{"summary.plsRbetamodel"}
#' @param \dots not used
#' @return \item{language}{call of the model}
#' @author Frédéric Bertrand\cr
#' \email{frederic.bertrand@@math.unistra.fr}\cr
#' \url{http://www-irma.u-strasbg.fr/~fbertran/}
#' @seealso \code{\link{print}} and \code{\link{summary}}
#' @references Frédéric Bertrand, Nicolas Meyer,
#' Michèle Beau-Faller, Karim El Bayed, Izzie-Jacques Namer,
#' Myriam Maumy-Bertrand (2013). Régression Bêta
#' PLS. \emph{Journal de la Société Française de Statistique},
#' \bold{154}(3):143-159.
#' \url{http://publications-sfds.math.cnrs.fr/index.php/J-SFdS/article/view/215}
#' @keywords methods print
#' @examples
#' 
#' data("GasolineYield",package="betareg")
#' modpls <- plsRbeta(yield~.,data=GasolineYield,nt=3,modele="pls-beta")
#' print(summary(modpls))
#' 
#' @export print.summary.plsRbetamodel
print.summary.plsRbetamodel <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
}
