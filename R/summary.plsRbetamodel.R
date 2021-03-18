#' Summary method for plsRbeta models
#' 
#' This function provides a summary method for the class \code{"plsRbetamodel"}
#' 
#' 
#' @param object an object of the class \code{"plsRbetamodel"}
#' @param \dots further arguments to be passed to or from methods.
#' @return \item{call }{function call of plsR beta models}
#' @author Frédéric Bertrand\cr
#' \email{frederic.bertrand@@math.unistra.fr}\cr
#' \url{http://www-irma.u-strasbg.fr/~fbertran/}
#' @seealso \code{\link{summary}}
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
#' summary(modpls)
#' 
summary.plsRbetamodel <- function(object, ...)
{
res <- list(call=object$call)
class(res) <- "summary.plsRbetamodel"
res
}
