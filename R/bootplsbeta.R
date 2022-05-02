#' Non-parametric Bootstrap for PLS beta regression models
#' 
#' Provides a wrapper for the bootstrap function \code{boot} from the
#' \code{boot} R package.\cr Implements non-parametric bootstrap for PLS beta
#' regression models by case resampling.
#' 
#' More details on bootstrap techniques are available in the help of the
#' \code{\link[boot:boot]{boot}} function.
#' 
#' @param object An object of class \code{plsRbetamodel} to bootstrap
#' @param typeboot The type of bootstrap. Either (Y,X) boostrap
#' (\code{typeboot="plsmodel"}) or (Y,T) bootstrap
#' (\code{typeboot="fmodel_np"}). Defaults to (Y,T) resampling.
#' @param R The number of bootstrap replicates. Usually this will be a single
#' positive integer. For importance resampling, some resamples may use one set
#' of weights and others use a different set of weights. In this case \code{R}
#' would be a vector of integers where each component gives the number of
#' resamples from each of the rows of weights.
#' @param statistic A function which when applied to data returns a vector
#' containing the statistic(s) of interest. \code{statistic} must take at least
#' two arguments. The first argument passed will always be the original data.
#' The second will be a vector of indices, frequencies or weights which define
#' the bootstrap sample. Further, if predictions are required, then a third
#' argument is required which would be a vector of the random indices used to
#' generate the bootstrap predictions. Any further arguments can be passed to
#' statistic through the \code{...} argument.
#' @param sim A character string indicating the type of simulation required.
#' Possible values are \code{"ordinary"} (the default), \code{"balanced"},
#' \code{"permutation"}, or \code{"antithetic"}.
#' @param stype A character string indicating what the second argument of
#' \code{statistic} represents. Possible values of stype are \code{"i"}
#' (indices - the default), \code{"f"} (frequencies), or \code{"w"} (weights).
#' @param stabvalue A value to hard threshold bootstrap estimates computed from
#' atypical resamplings.
#' @param \dots Other named arguments for \code{statistic} which are passed
#' unchanged each time it is called. Any such arguments to \code{statistic}
#' should follow the arguments which \code{statistic} is required to have for
#' the simulation. Beware of partial matching to arguments of \code{boot}
#' listed above.
#' @return An object of class \code{"boot"}. See the Value part of the help of
#' the function \code{\link[boot:boot]{boot}}.
#' @author Frédéric Bertrand\cr
#' \email{frederic.bertrand@@utt.fr}\cr
#' \url{https://fbertran.github.io/homepage/}
#' @seealso \code{\link[boot:boot]{boot}}
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
#' # Std coefficients
#' modplsbeta <- plsRbeta(yield~.,data=GasolineYield,nt=3, modele="pls-beta")
#' GazYield.boot <- bootplsbeta(modplsbeta, sim="ordinary", stype="i", R=250)
#' boot::boot.ci(GazYield.boot, conf = c(0.90,0.95), 
#' type = c("norm","basic","perc","bca"), index=1)
#' boot::boot.ci(GazYield.boot, conf = c(0.90,0.95), 
#' type = c("norm","basic","perc","bca"), index=2)
#' boot::boot.ci(GazYield.boot, conf = c(0.90,0.95), 
#' type = c("norm","basic","perc","bca"), index=3)
#' boot::boot.ci(GazYield.boot, conf = c(0.90,0.95), 
#' type = c("norm","basic","perc","bca"), index=4)
#' boot::boot.ci(GazYield.boot, conf = c(0.90,0.95), 
#' type = c("norm","basic","perc","bca"), index=5)
#' boot::boot.ci(GazYield.boot, conf = c(0.90,0.95), 
#' type = c("norm","basic","perc","bca"), index=6)
#' 
#' plsRglm::boxplots.bootpls(GazYield.boot)
#' plsRglm::confints.bootpls(GazYield.boot)
#' plsRglm::plots.confints.bootpls(plsRglm::confints.bootpls(GazYield.boot))
#' 
#' #Raw coefficients
#' modplsbeta <- plsRbeta(yield~.,data=GasolineYield,nt=3, modele="pls-beta")
#' GazYield.boot.raw <- bootplsbeta(modplsbeta, sim="ordinary", stype="i", 
#' R=250, statistic=coefs.plsRbeta.raw)
#' boot::boot.ci(GazYield.boot.raw, conf = c(0.90,0.95), 
#' type = c("norm","basic","perc","bca"), index=1)
#' boot::boot.ci(GazYield.boot.raw, conf = c(0.90,0.95), 
#' type = c("norm","basic","perc","bca"), index=2)
#' boot::boot.ci(GazYield.boot.raw, conf = c(0.90,0.95), 
#' type = c("norm","basic","perc","bca"), index=3)
#' boot::boot.ci(GazYield.boot.raw, conf = c(0.90,0.95), 
#' type = c("norm","basic","perc","bca"), index=4)
#' boot::boot.ci(GazYield.boot.raw, conf = c(0.90,0.95), 
#' type = c("norm","basic","perc","bca"), index=5)
#' boot::boot.ci(GazYield.boot.raw, conf = c(0.90,0.95), 
#' type = c("norm","basic","perc","bca"), index=6)
#' 
#' plsRglm::boxplots.bootpls(GazYield.boot)
#' plsRglm::confints.bootpls(GazYield.boot)
#' plsRglm::plots.confints.bootpls(plsRglm::confints.bootpls(GazYield.boot))
#' 
#' 
#' plot(GazYield.boot,index=2)
#' boot::jack.after.boot(GazYield.boot, index=2, useJ=TRUE, nt=3)
#' plot(GazYield.boot, index=2,jack=TRUE)
#' 
#' # PLS bootstrap balanced
#' 
#' GazYield.boot <- bootplsbeta(plsRbeta(yield~.,data=GasolineYield,nt=3,
#' modele="pls-beta"), sim="balanced", stype="i", R=250)
#' boot::boot.ci(GazYield.boot, conf = c(0.90,0.95), 
#' type = c("norm","basic","perc","bca"), index=1)
#' boot::boot.ci(GazYield.boot, conf = c(0.90,0.95), 
#' type = c("norm","basic","perc","bca"), index=2)
#' boot::boot.ci(GazYield.boot, conf = c(0.90,0.95), 
#' type = c("norm","basic","perc","bca"), index=3)
#' boot::boot.ci(GazYield.boot, conf = c(0.90,0.95), 
#' type = c("norm","basic","perc","bca"), index=4)
#' boot::boot.ci(GazYield.boot, conf = c(0.90,0.95), 
#' type = c("norm","basic","perc","bca"), index=5)
#' boot::boot.ci(GazYield.boot, conf = c(0.90,0.95), 
#' type = c("norm","basic","perc","bca"), index=6)
#' 
#' 
#' plsRglm::boxplots.bootpls(GazYield.boot)
#' plsRglm::confints.bootpls(GazYield.boot)
#' plsRglm::plots.confints.bootpls(plsRglm::confints.bootpls(GazYield.boot))
#' 
#' 
#' 
#' plot(GazYield.boot)
#' boot::jack.after.boot(GazYield.boot, index=1, useJ=TRUE, nt=3)
#' plot(GazYield.boot,jack=TRUE)
#' 
#' 
#' # PLS permutation bootstrap
#' 
#' GazYield.boot <- bootplsbeta(plsRbeta(yield~.,data=GasolineYield,nt=3,
#' modele="pls-beta"), sim="permutation", stype="i", R=250)
#' boot::boot.ci(GazYield.boot, conf = c(0.90,0.95), 
#' type = c("norm","basic","perc"), index=1)
#' boot::boot.ci(GazYield.boot, conf = c(0.90,0.95), 
#' type = c("norm","basic","perc"), index=2)
#' boot::boot.ci(GazYield.boot, conf = c(0.90,0.95), 
#' type = c("norm","basic","perc"), index=3)
#' boot::boot.ci(GazYield.boot, conf = c(0.90,0.95), 
#' type = c("norm","basic","perc"), index=4)
#' boot::boot.ci(GazYield.boot, conf = c(0.90,0.95), 
#' type = c("norm","basic","perc"), index=5)
#' boot::boot.ci(GazYield.boot, conf = c(0.90,0.95), 
#' type = c("norm","basic","perc"), index=6)
#' plsRglm::boxplots.bootpls(GazYield.boot)
#' plot(GazYield.boot)
#' }
#' 
bootplsbeta <- function(object, typeboot="plsmodel", R=250, statistic=NULL, sim="ordinary", stype="i", stabvalue=1e6, ...){
callplsRbeta <- object$call
maxcoefvalues <- stabvalue*abs(object$Coeffs)
#dataset <- cbind(y = eval(callplsRbeta$dataY),eval(callplsRbeta$dataX))
dataset <- cbind(y = object$dataY,object$dataX)
nt <- eval(callplsRbeta$nt)
ifbootfail <- as.matrix(as.numeric(ifelse(any(class(dataset[,1])=="factor"),rep(NA, ncol(dataset)+nlevels(dataset[,1])-1),rep(NA, ncol(dataset)))))
if(!is.null(callplsRbeta$modele)){modele <- eval(callplsRbeta$modele)} else {modele <- "pls"}
if(!is.null(callplsRbeta$family)){family <- eval(callplsRbeta$family)} else {family <- NULL}
if(!is.null(callplsRbeta$method)){method <- eval(callplsRbeta$method)} else {method <- "logistic"}
if(!is.null(callplsRbeta$link)){link <- eval(callplsRbeta$link)} else {link <- "logit"}
if(!is.null(callplsRbeta$link.phi)){link.phi <- eval(callplsRbeta$link.phi)} else {link.phi <- NULL}
if(!is.null(callplsRbeta$type)){type <- eval(callplsRbeta$type)} else {type <- "ML"}
if(!is.null(callplsRbeta$verbose)){verbose <- eval(callplsRbeta$verbose)} else {verbose <- TRUE}
if(typeboot=="plsmodel"){
temp.bootplsRbeta <- if(!(sim=="permutation")){if(is.null(statistic)){statistic=coefs.plsRbeta};boot(data=dataset, statistic=statistic, sim=sim, stype=stype, R=R, nt=nt, modele=modele, family=family, method=method, link=link, link.phi=link.phi,type=type, verbose=verbose,...)} else {
  if(is.null(statistic)){statistic=permcoefs.plsRbeta};boot(data=dataset, statistic=statistic, sim=sim, stype=stype, R=R, nt=nt, modele=modele, family=family, method=method, link=link, link.phi=link.phi, type=type, verbose=verbose)}
indices.temp.bootplsRbeta <- !is.na(temp.bootplsRbeta$t[,1])
temp.bootplsRbeta$t=temp.bootplsRbeta$t[indices.temp.bootplsRbeta,]
temp.bootplsRbeta$R=sum(indices.temp.bootplsRbeta)
temp.bootplsRbeta$call$R<-sum(indices.temp.bootplsRbeta)
return(temp.bootplsRbeta)
}
if(typeboot=="fmodel_np"){
  dataRepYtt <- cbind(y = object$RepY,object$tt)
  wwetoile <- object$wwetoile
  temp.bootplsRbeta <- if(!(sim=="permutation")){if(is.null(statistic)){statistic=coefs.plsRbetanp};boot(data=dataRepYtt, statistic=statistic, sim=sim, stype=stype, R=R, nt=nt, modele=modele, family=family, maxcoefvalues = maxcoefvalues[-(1:(length(object$Coeffs)-ncol(object$dataX)))], wwetoile = wwetoile, ifbootfail=ifbootfail, ...)} else {
    if(is.null(statistic)){statistic=permcoefs.plsRbetanp};boot(data=dataset, statistic=statistic, sim=sim, stype=stype, R=R, nt=nt, modele=modele, family=family, maxcoefvalues = maxcoefvalues[-(1:(length(object$Coeffs)-ncol(object$dataX)))], wwetoile = wwetoile, ifbootfail=ifbootfail)}
indices.temp.bootplsRbeta <- !is.na(temp.bootplsRbeta$t[,1])
temp.bootplsRbeta$t=temp.bootplsRbeta$t[indices.temp.bootplsRbeta,]
temp.bootplsRbeta$R=sum(indices.temp.bootplsRbeta)
temp.bootplsRbeta$call$R<-sum(indices.temp.bootplsRbeta)
return(temp.bootplsRbeta)
}
if(typeboot=="fmodel_par"){
temp.bootplsRbeta <- if(!(sim=="permutation")){if(is.null(statistic)){statistic=coefs.plsRbeta};boot(data=dataset, statistic=statistic, sim=sim, stype=stype, R=R, nt=nt, modele=modele, family=family, method=method, link=link, link.phi=link.phi,type=type, verbose=verbose,...)} else {
  if(is.null(statistic)){statistic=permcoefs.plsRbeta};boot(data=dataset, statistic=statistic, sim=sim, stype=stype, R=R, nt=nt, modele=modele, family=family, method=method, link=link, link.phi=link.phi,type=type, verbose=verbose)}
indices.temp.bootplsRbeta <- !is.na(temp.bootplsRbeta$t[,1])
temp.bootplsRbeta$t=temp.bootplsRbeta$t[indices.temp.bootplsRbeta,]
temp.bootplsRbeta$R=sum(indices.temp.bootplsRbeta)
temp.bootplsRbeta$call$R<-sum(indices.temp.bootplsRbeta)
return(temp.bootplsRbeta)
}
}
