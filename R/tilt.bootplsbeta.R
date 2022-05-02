#' Non-parametric tilted bootstrap for PLS beta regression models
#' 
#' Provides a wrapper for the bootstrap function \code{tilt.boot} from the
#' \code{boot} R package.\cr Implements non-parametric tilted bootstrap for PLS
#' beta regression models by case resampling : the \code{tilt.boot} function
#' will run an initial bootstrap with equal resampling probabilities (if
#' required) and will use the output of the initial run to find resampling
#' probabilities which put the value of the statistic at required values. It
#' then runs an importance resampling bootstrap using the calculated
#' probabilities as the resampling distribution.
#' 
#' 
#' @param object An object of class \code{plsRbetamodel} to bootstrap
#' @param typeboot The type of bootstrap. Either (Y,X) boostrap
#' (\code{typeboot="plsmodel"}) or (Y,T) bootstrap
#' (\code{typeboot="fmodel_np"}). Defaults to (Y,T) resampling.
#' @param statistic A function which when applied to data returns a vector
#' containing the statistic(s) of interest. \code{statistic} must take at least
#' two arguments. The first argument passed will always be the original data.
#' The second will be a vector of indices, frequencies or weights which define
#' the bootstrap sample. Further, if predictions are required, then a third
#' argument is required which would be a vector of the random indices used to
#' generate the bootstrap predictions. Any further arguments can be passed to
#' statistic through the \code{...} argument.
#' @param R The number of bootstrap replicates. Usually this will be a single
#' positive integer. For importance resampling, some resamples may use one set
#' of weights and others use a different set of weights. In this case \code{R}
#' would be a vector of integers where each component gives the number of
#' resamples from each of the rows of weights.
#' @param alpha The alpha level to which tilting is required. This parameter is
#' ignored if \code{R[1]} is 0 or if \code{theta} is supplied, otherwise it is
#' used to find the values of \code{theta} as quantiles of the initial uniform
#' bootstrap. In this case \code{R[1]} should be large enough that
#' \code{min(c(alpha, 1-alpha))*R[1] > 5}, if this is not the case then a
#' warning is generated to the effect that the \code{theta} are extreme values
#' and so the tilted output may be unreliable.
#' @param sim A character string indicating the type of simulation required.
#' Possible values are \code{"ordinary"} (the default), \code{"balanced"},
#' \code{"permutation"}, or \code{"antithetic"}.
#' @param stype A character string indicating what the second argument of
#' \code{statistic} represents. Possible values of stype are \code{"i"}
#' (indices - the default), \code{"f"} (frequencies), or \code{"w"} (weights).
#' @param index The index of the statistic of interest in the output from
#' \code{statistic}. By default the first element of the output of
#' \code{statistic} is used.
#' @return An object of class "boot".
#' @author Frédéric Bertrand\cr
#' \email{frederic.bertrand@@utt.fr}\cr
#' \url{https://fbertran.github.io/homepage/}
#' @seealso \code{\link[boot:boot]{tilt.boot}}
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
#' 
#' GazYield.tilt.boot <- tilt.bootplsbeta(plsRbeta(yield~.,data=GasolineYield,nt=3,
#' modele="pls-beta"), statistic=coefs.plsRbeta, R=c(499, 100, 100), 
#' alpha=c(0.025, 0.975), sim="balanced", stype="i", index=1)
#' boxplots.bootpls(GazYield.tilt.boot,1:2)
#' 
#' }
#' 
tilt.bootplsbeta <- function(object, typeboot="plsmodel", statistic=coefs.plsRbeta, R=c(499, 250, 250), alpha=c(0.025, 0.975), sim="ordinary", stype="i", index=1){
callplsRbeta <- object$call
dataset <- cbind(y = eval(callplsRbeta$dataY),eval(callplsRbeta$dataX))
nt <- eval(callplsRbeta$nt)
if(!is.null(callplsRbeta$modele)){modele <- eval(callplsRbeta$modele)} else {modele <- "pls"}
if(!is.null(callplsRbeta$family)){family <- eval(callplsRbeta$family)} else {family <- NULL}
if(!is.null(callplsRbeta$link)){link <- eval(callplsRbeta$link)} else {link <- "logit"}
if(!is.null(callplsRbeta$link.phi)){link.phi <- eval(callplsRbeta$link.phi)} else {link.phi <- NULL}
if(!is.null(callplsRbeta$type)){type <- eval(callplsRbeta$type)} else {type <- "ML"}
if(!is.null(callplsRbeta$verbose)){verbose <- eval(callplsRbeta$verbose)} else {verbose <- TRUE}
if(typeboot=="plsmodel"){
return(tilt.boot(data=dataset, statistic=if(!(sim=="permutation")){coefs.plsRbeta} else {permcoefs.plsRbeta}, sim=sim, stype=stype, R=R, nt=nt, modele=modele, family=family, link=link, link.phi=link.phi, type=type, verbose=verbose))
}
if(typeboot=="fmodel_np"){
return(tilt.boot(data=dataset, statistic=if(!(sim=="permutation")){coefs.plsRbeta} else {permcoefs.plsRbeta}, sim=sim, stype=stype, R=R, nt=nt, modele=modele, family=family, link=link, link.phi=link.phi, type=type, verbose=verbose))
}
if(typeboot=="fmodel_par"){
return(tilt.boot(data=dataset, statistic=if(!(sim=="permutation")){coefs.plsRbeta} else {permcoefs.plsRbeta}, sim=sim, stype=stype, R=R, nt=nt, modele=modele, family=family, link=link, link.phi=link.phi, type=type, verbose=verbose))
}
}
