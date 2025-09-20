#' Print method for plsRbeta models
#' 
#' This function provides a print method for the class \code{"plsRbetamodel"}
#' 
#' 
#' @param x an object of the class \code{"plsRbetamodel"}
#' @param \dots not used
#' @return \code{NULL}
#' @author Frédéric Bertrand\cr
#' \email{frederic.bertrand@@lecnam.net}\cr
#' \url{https://fbertran.github.io/homepage/}
#' @seealso \code{\link{print}}
#' @references Frédéric Bertrand, Nicolas Meyer,
#' Michèle Beau-Faller, Karim El Bayed, Izzie-Jacques Namer,
#' Myriam Maumy-Bertrand (2013). Régression Bêta
#' PLS. \emph{Journal de la Société Française de Statistique},
#' \bold{154}(3):143-159.
#' \url{https://ojs-test.apps.ocp.math.cnrs.fr/index.php/J-SFdS/article/view/215}
#' @keywords methods print
#' @examples
#' 
#' data("GasolineYield",package="betareg")
#' modpls <- plsRbeta(yield~.,data=GasolineYield,nt=3,modele="pls-beta")
#' print(modpls)
#' 
print.plsRbetamodel <- function(x,...)
{
  cat("Number of required components:\n")
  print(x$nt)
  cat("Number of successfully computed components:\n")
  print(x$computed_nt)
  cat("Coefficients:\n")
  print(x$Coeffs)
  cat("Information criteria and Fit statistics:\n")
  print(x$InfCrit)
  if (!is.null(x$family))
  {
    cat("Model with all the required components:\n")
    print(x$FinalModel)
  }
}
