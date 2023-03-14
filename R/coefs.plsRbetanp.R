#' Coefficients for bootstrap computations of PLSBeta models
#' 
#' A function passed to \code{boot} to perform bootstrap.
#' 
#' 
#' @param dataRepYtt components' coordinates to bootstrap
#' @param ind indices for resampling
#' @param nt number of components to use
#' @param modele type of modele to use, see \link{plsRbeta}
#' @param family glm family to use, see \link{plsRbeta}
#' @param method method for beta regression
#' @param link link for beta regression
#' @param link.phi link.phi for beta regression
#' @param type type of estimates
#' @param verbose should info messages be displayed ?
#' @param maxcoefvalues maximum values allowed for the estimates of the
#' coefficients to discard those coming from singular bootstrap samples
#' @param wwetoile values of the Wstar matrix in the original fit
#' @param ifbootfail value to return if the estimation fails on a bootstrap
#' sample
#' @return estimates on a bootstrap sample or \code{ifbootfail} value if the
#' bootstrap computation fails.
#' @note ~~some notes~~
#' @author Frédéric Bertrand\cr
#' \email{frederic.bertrand@@utt.fr}\cr
#' \url{https://fbertran.github.io/homepage/}
#' @seealso See also \code{\link{bootplsbeta}}
#' @keywords models
#' @examples
#' \donttest{
#' data("GasolineYield",package="betareg")
#' bootplsbeta(plsRbeta(yield~.,data=GasolineYield,nt=3, modele="pls-beta"), typeboot="fmodel_np", 
#' R=250, statistic=coefs.plsRbetanp)
#' }
#' @export coefs.plsRbetanp
coefs.plsRbetanp <- function(dataRepYtt, ind, nt, modele, family = NULL, method="logistic", link=NULL, link.phi=NULL, type="ML",verbose=TRUE,maxcoefvalues, wwetoile, ifbootfail) 
{
dataRepYb=dataRepYtt[ind,1]
Tb=dataRepYtt[ind,-1]
tempCb=try(solve(t(Tb)%*%Tb)%*%t(Tb)%*%dataRepYb,silent=TRUE)
tempcoefs <- wwetoile%*%tempCb
    Cond <- FALSE
    try(Cond<-is.numeric(tempcoefs)&all(abs(tempcoefs)<maxcoefvalues),silent=TRUE)
    if (Cond) {
        return(tempcoefs)
    }
    else {
        return(ifbootfail)
    }
}
