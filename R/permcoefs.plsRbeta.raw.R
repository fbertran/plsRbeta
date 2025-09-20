#' Raw coefficients function for permutation bootstrap techniques
#' 
#' A function passed to \code{boot} to perform bootstrap.
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
#' @param maxcoefvalues maximum values allowed for the estimates of the
#' coefficients to discard those coming from singular bootstrap samples
#' @param ifbootfail value to return if the estimation fails on a
#' @param verbose should info messages be displayed ?
#' @return Estimates on a bootstrap sample.
#' @author Frédéric Bertrand\cr
#' \email{frederic.bertrand@@lecnam.net}\cr
#' \url{https://fbertran.github.io/homepage/}
#' @seealso See also \code{\link{bootplsbeta}}.
#' @references Frédéric Bertrand, Nicolas Meyer,
#' Michèle Beau-Faller, Karim El Bayed, Izzie-Jacques Namer,
#' Myriam Maumy-Bertrand (2013). Régression Bêta
#' PLS. \emph{Journal de la Société Française de Statistique},
#' \bold{154}(3):143-159.
#' \url{https://ojs-test.apps.ocp.math.cnrs.fr/index.php/J-SFdS/article/view/215}
#' @keywords models
#' @examples
#' \donttest{
#' data("GasolineYield",package="betareg")
#' modplsbeta <- plsRbeta(yield~.,data=GasolineYield,nt=3, modele="pls-beta")
#' GazYield.boot.raw <- bootplsbeta(modplsbeta, sim="permutation", stype="i", 
#' R=250, statistic=coefs.plsRbeta.raw)
#' }
#' 
permcoefs.plsRbeta.raw <- function(dataset,ind,nt,modele,family=NULL,method="logistic",link="logit",link.phi=NULL,type="ML", maxcoefvalues,ifbootfail,verbose=TRUE){
  tempcoefs <- try(PLS_beta_wvc(dataY =dataset[,1], dataX=dataset[ind,-1], nt=nt, modele=modele, family=family, keepcoeffs=TRUE, method=method, link=link, link.phi=link.phi, type=type,verbose=verbose)$coeffs, silent=TRUE)
  Cond <- FALSE
  try(Cond<-is.numeric(tempcoefs)&all(abs(tempcoefs)<maxcoefvalues),silent=TRUE)
  if (Cond) {
    return(tempcoefs)
  }
  else {
    return(ifbootfail)
  }
}



