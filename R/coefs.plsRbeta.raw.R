#' Raw coefficients function for bootstrap techniques
#' 
#' Returns the coefficients of a \code{"plsRbeta"} model.
#' 
#' 
#' @param dataset dataset to resample
#' @param ind indices for resampling
#' @param nt number of components to use
#' @param modele name of the PLS glm or PLS beta model to be fitted
#' (\code{"pls"}, \code{"pls-glm-Gamma"}, \code{"pls-glm-gaussian"},
#' \code{"pls-glm-inverse.gaussian"}, \code{"pls-glm-logistic"},
#' \code{"pls-glm-poisson"}, \code{"pls-glm-polr"}, \code{"pls-beta"}). Use
#' \code{"modele=pls-glm-family"} to enable the \code{family} option.
#' @param family family to use if GLM model, see \link{plsRbeta}
#' @param method method for beta regression
#' @param link link for beta regression
#' @param link.phi link.phi for beta regression
#' @param type type of estimates
#' @param verbose should info messages be displayed ?
#' @return Coefficients' Estimates on a sample.
#' @author Frédéric Bertrand\cr
#' \email{frederic.bertrand@@utt.fr}\cr
#' \url{https://fbertran.github.io/homepage/}
#' @seealso See also \code{\link{bootplsbeta}}.
#' @references Frédéric Bertrand, Nicolas Meyer,
#' Michèle Beau-Faller, Karim El Bayed, Izzie-Jacques Namer,
#' Myriam Maumy-Bertrand (2013). Régression Bêta
#' PLS. \emph{Journal de la Société Française de Statistique},
#' \bold{154}(3):143-159.
#' \url{http://publications-sfds.math.cnrs.fr/index.php/J-SFdS/article/view/215}
#' @keywords models
#' @examples
#' \donttest{
#' data("GasolineYield",package="betareg")
#' modpls <- coefs.plsRbeta.raw(GasolineYield[,-6],1:32,nt=3,modele="pls-beta")
#' }
#' 
coefs.plsRbeta.raw <- function(dataset, ind, nt, modele, family=NULL, method="logistic", link=NULL, link.phi=NULL, type="ML",verbose=TRUE) 
{
    tempcoefs <- try(PLS_beta_wvc(dataY = dataset[ind, 1], dataX = dataset[ind, 
        -1], nt = nt, modele = modele, family=family, method=method, link=link, keepcoeffs = TRUE, link.phi=link.phi, type=type, verbose=verbose)$coeffs, silent=TRUE)
    if (is.numeric(tempcoefs)) {
        return(tempcoefs)
    }
    else {
        return(as.matrix(as.numeric(ifelse(any(class(dataset[,1])=="factor"),rep(NA, ncol(dataset)+nlevels(dataset[,1])-1),rep(NA, ncol(dataset))))))
    }
}
