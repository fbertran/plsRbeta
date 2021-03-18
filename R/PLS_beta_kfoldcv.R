#' Partial least squares regression beta models with kfold cross validation
#' 
#' This function implements kfold cross validation on complete or incomplete
#' datasets for partial least squares beta regression models
#' 
#' Predicts 1 group with the \code{K-1} other groups. Leave one out cross
#' validation is thus obtained for \code{K==nrow(dataX)}.
#' 
#' There are seven different predefined models with predefined link functions
#' available : \describe{ \item{list("\"pls\"")}{ordinary pls models}
#' \item{list("\"pls-glm-Gamma\"")}{glm gaussian with inverse link pls models}
#' \item{list("\"pls-glm-gaussian\"")}{glm gaussian with identity link pls
#' models} \item{list("\"pls-glm-inverse-gamma\"")}{glm binomial with square
#' inverse link pls models} \item{list("\"pls-glm-logistic\"")}{glm binomial
#' with logit link pls models} \item{list("\"pls-glm-poisson\"")}{glm poisson
#' with log link pls models} \item{list("\"pls-glm-polr\"")}{glm polr with
#' logit link pls models} } Using the \code{"family="} option and setting
#' \code{"modele=pls-glm-family"} allows changing the family and link function
#' the same way as for the \code{\link[stats]{glm}} function. As a consequence
#' user-specified families can also be used.  \describe{ \item{The }{accepts
#' the links (as names) \code{identity}, \code{log} and
#' \code{inverse}.}\item{list("gaussian")}{accepts the links (as names)
#' \code{identity}, \code{log} and \code{inverse}.}\item{ family}{accepts the
#' links (as names) \code{identity}, \code{log} and \code{inverse}.} \item{The
#' }{accepts the links \code{logit}, \code{probit}, \code{cauchit},
#' (corresponding to logistic, normal and Cauchy CDFs respectively) \code{log}
#' and \code{cloglog} (complementary log-log).}\item{list("binomial")}{accepts
#' the links \code{logit}, \code{probit}, \code{cauchit}, (corresponding to
#' logistic, normal and Cauchy CDFs respectively) \code{log} and \code{cloglog}
#' (complementary log-log).}\item{ family}{accepts the links \code{logit},
#' \code{probit}, \code{cauchit}, (corresponding to logistic, normal and Cauchy
#' CDFs respectively) \code{log} and \code{cloglog} (complementary log-log).}
#' \item{The }{accepts the links \code{inverse}, \code{identity} and
#' \code{log}.}\item{list("Gamma")}{accepts the links \code{inverse},
#' \code{identity} and \code{log}.}\item{ family}{accepts the links
#' \code{inverse}, \code{identity} and \code{log}.} \item{The }{accepts the
#' links \code{log}, \code{identity}, and
#' \code{sqrt}.}\item{list("poisson")}{accepts the links \code{log},
#' \code{identity}, and \code{sqrt}.}\item{ family}{accepts the links
#' \code{log}, \code{identity}, and \code{sqrt}.} \item{The }{accepts the links
#' \code{1/mu^2}, \code{inverse}, \code{identity} and
#' \code{log}.}\item{list("inverse.gaussian")}{accepts the links \code{1/mu^2},
#' \code{inverse}, \code{identity} and \code{log}.}\item{ family}{accepts the
#' links \code{1/mu^2}, \code{inverse}, \code{identity} and \code{log}.}
#' \item{The }{accepts the links \code{logit}, \code{probit}, \code{cloglog},
#' \code{identity}, \code{inverse}, \code{log}, \code{1/mu^2} and
#' \code{sqrt}.}\item{list("quasi")}{accepts the links \code{logit},
#' \code{probit}, \code{cloglog}, \code{identity}, \code{inverse}, \code{log},
#' \code{1/mu^2} and \code{sqrt}.}\item{ family}{accepts the links
#' \code{logit}, \code{probit}, \code{cloglog}, \code{identity},
#' \code{inverse}, \code{log}, \code{1/mu^2} and \code{sqrt}.} \item{The
#' function }{can be used to create a power link
#' function.}\item{list("power")}{can be used to create a power link function.}
#' }
#' 
#' Non-NULL weights can be used to indicate that different observations have
#' different dispersions (with the values in weights being inversely
#' proportional to the dispersions); or equivalently, when the elements of
#' weights are positive integers w_i, that each response y_i is the mean of w_i
#' unit-weight observations.
#' 
#' @param dataY response (training) dataset
#' @param dataX predictor(s) (training) dataset
#' @param nt number of components to be extracted
#' @param limQ2set limit value for the Q2
#' @param modele name of the PLS glm or PLS beta model to be fitted
#' (\code{"pls"}, \code{"pls-glm-Gamma"}, \code{"pls-glm-gaussian"},
#' \code{"pls-glm-inverse.gaussian"}, \code{"pls-glm-logistic"},
#' \code{"pls-glm-poisson"}, \code{"pls-glm-polr"}, \code{"pls-beta"}). Use
#' \code{"modele=pls-glm-family"} to enable the \code{family} option.
#' @param family a description of the error distribution and link function to
#' be used in the model. This can be a character string naming a family
#' function, a family function or the result of a call to a family function.
#' (See \code{\link[stats]{family}} for details of family functions.) To use
#' the family option, please set \code{modele="pls-glm-family"}. User defined
#' families can also be defined. See details.
#' @param K number of groups
#' @param NK number of times the group division is made
#' @param grouplist to specify the members of the \code{K} groups
#' @param random should the \code{K} groups be made randomly
#' @param scaleX scale the predictor(s) : must be set to TRUE for
#' \code{modele="pls"} and should be for glms pls.
#' @param scaleY scale the response : Yes/No. Ignored since non always possible
#' for glm responses.
#' @param keepcoeffs shall the coefficients for each model be returned
#' @param keepfolds shall the groups' composition be returned
#' @param keepdataY shall the observed value of the response for each one of
#' the predicted value be returned
#' @param keepMclassed shall the number of miss classed be returned
#' (unavailable)
#' @param tol_Xi minimal value for Norm2(Xi) and \eqn{\mathrm{det}(pp' \times
#' pp)}{det(pp'*pp)} if there is any missing value in the \code{dataX}. It
#' defaults to \eqn{10^{-12}}{10^{-12}}
#' @param weights an optional vector of 'prior weights' to be used in the
#' fitting process. Should be \code{NULL} or a numeric vector.
#' @param method logistic, probit, complementary log-log or cauchit
#' (corresponding to a Cauchy latent variable).
#' @param link character specification of the link function in the mean model
#' (mu). Currently, "\code{logit}", "\code{probit}", "\code{cloglog}",
#' "\code{cauchit}", "\code{log}", "\code{loglog}" are supported.
#' Alternatively, an object of class "\code{link-glm}" can be supplied.
#' @param link.phi character specification of the link function in the
#' precision model (phi). Currently, "\code{identity}", "\code{log}",
#' "\code{sqrt}" are supported. The default is "\code{log}" unless
#' \code{formula} is of type \code{y~x} where the default is "\code{identity}"
#' (for backward compatibility). Alternatively, an object of class
#' "\code{link-glm}" can be supplied.
#' @param type character specification of the type of estimator. Currently,
#' maximum likelihood ("\code{ML}"), ML with bias correction ("\code{BC}"), and
#' ML with bias reduction ("\code{BR}") are supported.
#' @param verbose should info messages be displayed ?
#' @return \item{results_kfolds}{list of \code{NK}. Each element of the list
#' sums up the results for a group division: \describe{ \item{}{list of
#' \code{K} matrices of size about \code{nrow(dataX)/K * nt} with the predicted
#' values for a growing number of components} \item{list()}{\dots{}}
#' \item{}{list of \code{K} matrices of size about \code{nrow(dataX)/K * nt}
#' with the predicted values for a growing number of components} }}
#' \item{folds}{list of \code{NK}. Each element of the list sums up the
#' informations for a group division: \describe{ \item{}{list of \code{K}
#' vectors of length about \code{nrow(dataX)} with the numbers of the rows of
#' \code{dataX} that were used as a training set} \item{list()}{\dots{}}
#' \item{}{list of \code{K} vectors of length about \code{nrow(dataX)} with the
#' numbers of the rows of \code{dataX} that were used as a training set} } }
#' \item{dataY_kfolds}{list of \code{NK}. Each element of the list sums up the
#' results for a group division: \describe{ \item{}{list of \code{K} matrices
#' of size about \code{nrow(dataX)/K * 1} with the observed values of the
#' response} \item{list()}{\dots{}} \item{}{list of \code{K} matrices of size
#' about \code{nrow(dataX)/K * 1} with the observed values of the response} } }
#' \item{call}{the call of the function}
#' @note Works for complete and incomplete datasets.
#' @author Frédéric Bertrand\cr
#' \email{frederic.bertrand@@math.unistra.fr}\cr
#' \url{http://www-irma.u-strasbg.fr/~fbertran/}
#' @seealso \code{\link[plsRglm]{kfolds2coeff}},
#' \code{\link[plsRglm]{kfolds2Pressind}}, \code{\link[plsRglm]{kfolds2Press}},
#' \code{\link[plsRglm]{kfolds2Mclassedind}},
#' \code{\link[plsRglm]{kfolds2Mclassed}} and
#' \code{\link[plsRbeta]{kfolds2CVinfos_beta}} to extract and transform results
#' from kfold cross validation.
#' @references Frédéric Bertrand, Nicolas Meyer,
#' Michèle Beau-Faller, Karim El Bayed, Izzie-Jacques Namer,
#' Myriam Maumy-Bertrand (2013). Régression Bêta
#' PLS. \emph{Journal de la Société Française de Statistique},
#' \bold{154}(3):143-159.
#' \url{http://publications-sfds.math.cnrs.fr/index.php/J-SFdS/article/view/215}
#' @keywords models regression
#' @examples
#' 
#' \dontrun{
#' data("GasolineYield",package="betareg")
#' yGasolineYield <- GasolineYield$yield
#' XGasolineYield <- GasolineYield[,2:5]
#' bbb <- PLS_beta_kfoldcv(yGasolineYield,XGasolineYield,nt=3,modele="pls-beta")
#' kfolds2CVinfos_beta(bbb)
#' }
#' 
PLS_beta_kfoldcv <- function(dataY,dataX,nt=2,limQ2set=.0975,modele="pls", family=NULL, K=nrow(dataX), NK=1, grouplist=NULL, random=FALSE, scaleX=TRUE, scaleY=NULL, keepcoeffs=FALSE, keepfolds=FALSE, keepdataY=TRUE, keepMclassed=FALSE, tol_Xi=10^(-12),weights,method,link=NULL,link.phi=NULL,type="ML",verbose=TRUE) {

    if (missing(weights)) {NoWeights <- TRUE} else {NoWeights <- FALSE}
    res <- NULL
    res$nr <- nrow(dataX)
        if (K > res$nr) {
          if(verbose){cat(paste("K cannot be > than nrow(dataX) =",res$nr,"\n"))}
          if(verbose){cat(paste("K is set to", nrow(dataX), "\n"))}
            K <- res$nr
            random = FALSE
        }
    call <- match.call(expand.dots=FALSE)
    nt <- eval(nt,parent.frame())
    if (is.null(modele) & !is.null(family)) {modele<-"pls-glm-family"}
    if (as.character(call["family"])=="NULL") {
        if (modele=="pls") {call$family<-NULL}
        if (modele=="pls-beta") {family<-NULL}
        if (modele=="pls-glm-Gamma") {call$family<-Gamma(link = "inverse")}
        if (modele=="pls-glm-gaussian") {call$family<-gaussian(link = "identity")}
        if (modele=="pls-glm-inverse.gaussian") {call$family<-inverse.gaussian(link = "1/mu^2")}
        if (modele=="pls-glm-logistic") {call$family<-binomial(link = "logit")}
        if (modele=="pls-glm-poisson") {call$family<-poisson(link = "log")}
        if (modele=="pls-glm-polr") {call$family<-NULL}
    }
    if (!is.null(call$family)) {
        if (is.character(call$family)) {call$family <- get(call$family, mode = "function", envir = parent.frame())}
        if (is.function(call$family)) {call$family <- call$family()}
        if (is.language(call$family)) {call$family <- eval(call$family)}
    }
    if (missing(method)){method<-"logistic"}
    if (is.null(link)){link<-"logit"} else {if(!(link %in% c("logit", "probit", "cloglog", "cauchit", "log", "loglog")) & !is(link,"link-glm")) {link<-"logit"}}
    if (modele=="pls") {if(verbose){cat("\nModel:", modele, "\n\n")}}
    if (modele %in% c("pls-glm-family","pls-glm-Gamma","pls-glm-gaussian","pls-glm-inverse.gaussian","pls-glm-logistic","pls-glm-poisson")) {if(verbose){if(verbose){cat(family,"\n\n")}}}
    if (modele %in% c("pls-glm-polr")) {if(verbose){cat("\nModel:", modele, "\n");cat("Method:", method, "\n\n")}}
    if (modele=="pls-beta") {if(verbose){cat("\nModel:", modele, "\n\n");cat("Link:", link, "\n\n");cat("Link.phi:", link.phi, "\n\n");cat("Type:", type, "\n\n")}}

    if (as.character(call["tol_Xi"])=="NULL") {call$tol_Xi <- 10^(-12)}
    if (as.character(call["modele"])=="NULL") {call$modele <- "pls"}
    if (as.character(call["limQ2set"])=="NULL") {call$limQ2set <- .0975}
    if (!is.data.frame(dataX)) {dataX <- data.frame(dataX)}
    folds_kfolds <-vector("list",NK)
    if (NK==1) {respls_kfolds <- list(vector("list", K))}
    else
    {
      if (NK>1)
      {
      respls_kfolds <-vector("list",NK)
        for (jj in 1:NK) {
          respls_kfolds[[jj]] <-vector("list",K)
        }
      }
    }
    if (modele=="pls-beta") {
    if (NK==1) {respls_kfolds_phi <- list(vector("list", K))}
    else
    {
      if (NK>1)
      {
      respls_kfolds_phi <-vector("list",NK)
        for (jj in 1:NK) {
          respls_kfolds_phi[[jj]] <-vector("list",K)
        }
      }
    }
    }
    if (keepdataY) {
    if (NK==1) {dataY_kfolds <- list(vector("list", K))}
    else
    {
      if (NK>1)
      {
      dataY_kfolds <-vector("list",NK)
        for (jj in 1:NK) {
          dataY_kfolds[[jj]] <-vector("list",K)
        }
      }
    }
    }
    if (keepcoeffs) {
    if (NK==1) {coeffs_kfolds <- list(vector("list", K))}
    else
    {
      if (NK>1)
      {
      coeffs_kfolds <-vector("list",NK)
        for (jj in 1:NK) {
          coeffs_kfolds[[jj]] <-vector("list",K)
        }
      }
    }
    }
    compl = function (part, set)
    {
        comp = c()
        for (z in set) {
            if (length(which(z == part)) == 0) {
                comp = c(comp, z)
            }
        }
        return(comp)
    }
    for (nnkk in 1:NK) {
      if(verbose){cat(paste("NK:", nnkk, "\n"))}
        if (K == res$nr) {
          if(verbose){cat("Leave One Out\n")}
            random = FALSE
        }
      if(verbose){cat(paste("Number of groups :", K, "\n"))}
        if (!is.list(grouplist)) {
            if (random == TRUE) {
                randsample = sample(1:res$nr, replace = FALSE)
                groups = suppressWarnings(split(randsample, as.factor(1:K)))
            }
            else {
                randsample = sample(1:res$nr, replace = FALSE)
                groups = suppressWarnings(split(randsample, as.factor(1:K)))
                be = 1
                en = 0
                for (z in 1:K) {
                  en = en + length(unlist(groups[z]))
                  groups[z] = list(z = c(be:en))
                  be = en + 1
                }
            }
        }
        else {
            nogroups = grouplist[[nnkk]]
            groups = c()
            for (i in 1:K) groups = c(groups, list(compl(as.vector(unlist(nogroups[i])),
                (1:res$nr))))
        }
        rnames = c()
        for (k in 1:K) rnames = c(rnames, rownames(dataX)[-as.vector(unlist(groups[k]))])
        if (K == 1) {rnames = rownames(dataX)}
        folds = c()
        for (ii in 1:K) {
            nofolds = as.vector(unlist(groups[ii]))
            if (K == 1) {
                folds = c(folds, list(nofolds))
                nofolds = NULL
            }
            else folds = c(folds, list(as.vector(unlist(groups[-ii]))))
            if (K == 1) {
                if(NoWeights){
                temptemp <- PLS_beta_wvc(dataY=dataY, dataX=dataX, nt=nt, dataPredictY=dataX, modele=modele,family=family,scaleX=scaleX,scaleY=scaleY,keepcoeffs=keepcoeffs,tol_Xi=tol_Xi,method=method,link=link,link.phi=link.phi,type=type,verbose=verbose)
                respls_kfolds[[nnkk]][[ii]] <- temptemp$valsPredict
                } else {
                temptemp <- PLS_beta_wvc(dataY=dataY, dataX=dataX, nt=nt, dataPredictY=dataX, modele=modele,family=family,scaleX=scaleX,scaleY=scaleY,keepcoeffs=keepcoeffs,tol_Xi=tol_Xi,,weights=weights,method=method,link=link,link.phi=link.phi,type=type,verbose=verbose)
                respls_kfolds[[nnkk]][[ii]] <- temptemp$valsPredict; attr(respls_kfolds[[nnkk]],"XWeights")=weights; attr(respls_kfolds[[nnkk]],"YWeights")=NULL}             
                if (keepdataY) {dataY_kfolds[[nnkk]][[ii]] = NULL}
                if (keepcoeffs) {coeffskfolds[[nnkk]][[ii]] = temptemp$coeffs}
                if (modele=="pls-beta") {respls_kfolds_phi[[nnkk]][[ii]] = temptemp$valsPredictPhis}
                }
            else {
              if(verbose){cat(paste(ii,"\n"))}
                  if(NoWeights){
                  temptemp <- PLS_beta_wvc(dataY=dataY[-nofolds], dataX=dataX[-nofolds,], nt=nt, dataPredictY=dataX[nofolds,], modele=modele,family=family,scaleX=scaleX,scaleY=scaleY,keepcoeffs=keepcoeffs,tol_Xi=tol_Xi,method=method,link=link,link.phi=link.phi,type=type,verbose=verbose)
                  respls_kfolds[[nnkk]][[ii]] <- temptemp$valsPredict
                  } else {
                  temptemp <- PLS_beta_wvc(dataY=dataY[-nofolds], dataX=dataX[-nofolds,], nt=nt, dataPredictY=dataX[nofolds,], modele=modele,family=family,scaleX=scaleX,scaleY=scaleY,keepcoeffs=keepcoeffs,tol_Xi=tol_Xi,weights=weights[-nofolds],method=method,link=link,link.phi=link.phi,type=type,verbose=verbose) 
                respls_kfolds[[nnkk]][[ii]] <- temptemp$valsPredict; attr(respls_kfolds[[nnkk]][[ii]],"XWeights")=weights[-nofolds]; attr(respls_kfolds[[nnkk]][[ii]],"YWeights")=weights[nofolds]}
                  if (keepdataY) {dataY_kfolds[[nnkk]][[ii]] = dataY[nofolds]}
                  if (keepcoeffs) {coeffs_kfolds[[nnkk]][[ii]] = temptemp$coeffs}
                  if (modele=="pls-beta") {respls_kfolds_phi[[nnkk]][[ii]] = temptemp$valsPredictPhis}
                  }
        }
        folds_kfolds[[nnkk]]<-folds
    }
results <- list(results_kfolds=respls_kfolds)
if (keepcoeffs) {results$coeffs_kfolds <- coeffs_kfolds}
if (keepfolds) {results$folds <- folds_kfolds}
if (keepdataY) {results$dataY_kfolds <- dataY_kfolds}
if (modele=="pls-beta") {results$results_kfolds_phi <- respls_kfolds_phi}
results$call <- call
results$call$nt <- nt
return(results)
}
