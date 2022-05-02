#' Partial least squares regression beta models with kfold cross validation
#' 
#' This function implements kfold cross validation on complete or incomplete
#' datasets for partial least squares beta regression models (formula
#' specification of the model).
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
#' A typical predictor has the form response ~ terms where response is the
#' (numeric) response vector and terms is a series of terms which specifies a
#' linear predictor for response. A terms specification of the form first +
#' second indicates all the terms in first together with all the terms in
#' second with any duplicates removed.
#' 
#' A specification of the form first:second indicates the the set of terms
#' obtained by taking the interactions of all terms in first with all terms in
#' second. The specification first*second indicates the cross of first and
#' second. This is the same as first + second + first:second.
#' 
#' The terms in the formula will be re-ordered so that main effects come first,
#' followed by the interactions, all second-order, all third-order and so on:
#' to avoid this pass a terms object as the formula.
#' 
#' Non-NULL weights can be used to indicate that different observations have
#' different dispersions (with the values in weights being inversely
#' proportional to the dispersions); or equivalently, when the elements of
#' weights are positive integers w_i, that each response y_i is the mean of w_i
#' unit-weight observations.
#' 
#' @param formula an object of class "\code{\link{formula}}" (or one that can
#' be coerced to that class): a symbolic description of the model to be fitted.
#' The details of model specification are given under 'Details'.
#' @param data an optional data frame, list or environment (or object coercible
#' by \code{\link{as.data.frame}} to a data frame) containing the variables in
#' the model. If not found in \code{data}, the variables are taken from
#' \code{environment(formula)}, typically the environment from which
#' \code{plsRglm} is called.
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
#' @param subset an optional vector specifying a subset of observations to be
#' used in the fitting process.
#' @param start starting values for the parameters in the linear predictor.
#' @param etastart starting values for the linear predictor.
#' @param mustart starting values for the vector of means.
#' @param offset this can be used to specify an \emph{a priori} known component
#' to be included in the linear predictor during fitting. This should be
#' \code{NULL} or a numeric vector of length equal to the number of cases. One
#' or more \code{\link{offset}} terms can be included in the formula instead or
#' as well, and if more than one is specified their sum is used. See
#' \code{\link{model.offset}}.
#' @param method \describe{ \item{for fitting glms with glm (}{the method to be
#' used in fitting the model. The default method \code{"glm.fit"} uses
#' iteratively reweighted least squares (IWLS). User-supplied fitting functions
#' can be supplied either as a function or a character string naming a
#' function, with a function which takes the same arguments as \code{glm.fit}.
#' If "model.frame", the model frame is
#' returned.}\item{list("\"pls-glm-Gamma\"")}{the method to be used in fitting
#' the model. The default method \code{"glm.fit"} uses iteratively reweighted
#' least squares (IWLS). User-supplied fitting functions can be supplied either
#' as a function or a character string naming a function, with a function which
#' takes the same arguments as \code{glm.fit}. If "model.frame", the model
#' frame is returned.}\item{, }{the method to be used in fitting the model. The
#' default method \code{"glm.fit"} uses iteratively reweighted least squares
#' (IWLS). User-supplied fitting functions can be supplied either as a function
#' or a character string naming a function, with a function which takes the
#' same arguments as \code{glm.fit}. If "model.frame", the model frame is
#' returned.}\item{list("\"pls-glm-gaussian\"")}{the method to be used in
#' fitting the model. The default method \code{"glm.fit"} uses iteratively
#' reweighted least squares (IWLS). User-supplied fitting functions can be
#' supplied either as a function or a character string naming a function, with
#' a function which takes the same arguments as \code{glm.fit}. If
#' "model.frame", the model frame is returned.}\item{, }{the method to be used
#' in fitting the model. The default method \code{"glm.fit"} uses iteratively
#' reweighted least squares (IWLS). User-supplied fitting functions can be
#' supplied either as a function or a character string naming a function, with
#' a function which takes the same arguments as \code{glm.fit}. If
#' "model.frame", the model frame is
#' returned.}\item{list("\"pls-glm-inverse.gaussian\"")}{the method to be used
#' in fitting the model. The default method \code{"glm.fit"} uses iteratively
#' reweighted least squares (IWLS). User-supplied fitting functions can be
#' supplied either as a function or a character string naming a function, with
#' a function which takes the same arguments as \code{glm.fit}. If
#' "model.frame", the model frame is returned.}\item{, }{the method to be used
#' in fitting the model. The default method \code{"glm.fit"} uses iteratively
#' reweighted least squares (IWLS). User-supplied fitting functions can be
#' supplied either as a function or a character string naming a function, with
#' a function which takes the same arguments as \code{glm.fit}. If
#' "model.frame", the model frame is
#' returned.}\item{list("\"pls-glm-logistic\"")}{the method to be used in
#' fitting the model. The default method \code{"glm.fit"} uses iteratively
#' reweighted least squares (IWLS). User-supplied fitting functions can be
#' supplied either as a function or a character string naming a function, with
#' a function which takes the same arguments as \code{glm.fit}. If
#' "model.frame", the model frame is returned.}\item{, }{the method to be used
#' in fitting the model. The default method \code{"glm.fit"} uses iteratively
#' reweighted least squares (IWLS). User-supplied fitting functions can be
#' supplied either as a function or a character string naming a function, with
#' a function which takes the same arguments as \code{glm.fit}. If
#' "model.frame", the model frame is
#' returned.}\item{list("\"pls-glm-poisson\"")}{the method to be used in
#' fitting the model. The default method \code{"glm.fit"} uses iteratively
#' reweighted least squares (IWLS). User-supplied fitting functions can be
#' supplied either as a function or a character string naming a function, with
#' a function which takes the same arguments as \code{glm.fit}. If
#' "model.frame", the model frame is returned.}\item{, }{the method to be used
#' in fitting the model. The default method \code{"glm.fit"} uses iteratively
#' reweighted least squares (IWLS). User-supplied fitting functions can be
#' supplied either as a function or a character string naming a function, with
#' a function which takes the same arguments as \code{glm.fit}. If
#' "model.frame", the model frame is
#' returned.}\item{list("\"modele=pls-glm-family\"")}{the method to be used in
#' fitting the model. The default method \code{"glm.fit"} uses iteratively
#' reweighted least squares (IWLS). User-supplied fitting functions can be
#' supplied either as a function or a character string naming a function, with
#' a function which takes the same arguments as \code{glm.fit}. If
#' "model.frame", the model frame is returned.}\item{)}{the method to be used
#' in fitting the model. The default method \code{"glm.fit"} uses iteratively
#' reweighted least squares (IWLS). User-supplied fitting functions can be
#' supplied either as a function or a character string naming a function, with
#' a function which takes the same arguments as \code{glm.fit}. If
#' "model.frame", the model frame is returned.}
#' \item{list("pls-glm-polr")}{logistic, probit, complementary log-log or
#' cauchit (corresponding to a Cauchy latent variable).}}
#' @param control a list of parameters for controlling the fitting process. For
#' \code{glm.fit} this is passed to \code{\link{glm.control}}.
#' @param contrasts an optional list. See the \code{contrasts.arg} of
#' \code{model.matrix.default}.
#' @param sparse should the coefficients of non-significant predictors
#' (<\code{alpha.pvals.expli}) be set to 0
#' @param sparseStop should component extraction stop when no significant
#' predictors (<\code{alpha.pvals.expli}) are found
#' @param naive Use the naive estimates for the Degrees of Freedom in plsR?
#' Default is \code{FALSE}.
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
#' sums up the results for a group division: \describe{ \item{list}{ of
#' \code{K} matrices of size about \code{nrow(dataX)/K * nt} with the predicted
#' values for a growing number of components} \item{list()}{\dots{}}
#' \item{list}{ of \code{K} matrices of size about \code{nrow(dataX)/K * nt}
#' with the predicted values for a growing number of components} }}
#' \item{folds}{list of \code{NK}. Each element of the list sums up the
#' informations for a group division: \describe{ \item{list}{ of \code{K}
#' vectors of length about \code{nrow(dataX)} with the numbers of the rows of
#' \code{dataX} that were used as a training set} \item{list()}{\dots{}}
#' \item{list}{ of \code{K} vectors of length about \code{nrow(dataX)} with the
#' numbers of the rows of \code{dataX} that were used as a training set} } }
#' \item{dataY_kfolds}{list of \code{NK}. Each element of the list sums up the
#' results for a group division: \describe{ \item{list}{ of \code{K} matrices
#' of size about \code{nrow(dataX)/K * 1} with the observed values of the
#' response} \item{list()}{\dots{}} \item{list}{ of \code{K} matrices of size
#' about \code{nrow(dataX)/K * 1} with the observed values of the response} } }
#' \item{call}{the call of the function}
#' @note Work for complete and incomplete datasets.
#' @author Frédéric Bertrand\cr
#' \email{frederic.bertrand@@utt.fr}\cr
#' \url{https://fbertran.github.io/homepage/}
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
#' bbb <- PLS_beta_kfoldcv_formula(yield~.,data=GasolineYield,nt=3,modele="pls-beta")
#' kfolds2CVinfos_beta(bbb)
#' }
#' 
PLS_beta_kfoldcv_formula <- function(formula,data=NULL,nt=2,limQ2set=.0975,modele="pls", family=NULL, K=nrow(dataX), NK=1, grouplist=NULL, random=FALSE, scaleX=TRUE, scaleY=NULL, keepcoeffs=FALSE, keepfolds=FALSE, keepdataY=TRUE, keepMclassed=FALSE, tol_Xi=10^(-12),weights,subset,start=NULL,etastart,mustart,offset,method,control=list(),contrasts=NULL,sparse=FALSE,sparseStop=TRUE,naive=FALSE,link=NULL,link.phi=NULL,type="ML",verbose=TRUE) {

    if (missing(weights)) {NoWeights <- TRUE} else {NoWeights <- FALSE}  
    if (missing(data)) {data <- environment(formula)}
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "etastart", "mustart", "offset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf$na.action <- na.pass    

if (modele %in% c("pls-beta")) {
oformula <- as.formula(formula)
    formula <- Formula::as.Formula(formula)
    if (length(formula)[2L] < 2L) {
        formula <- Formula::as.Formula(formula(formula), ~1)
        simple_formula <- TRUE
    }
    else {
        if (length(formula)[2L] > 2L) {
            formula <- Formula::Formula(formula(formula, rhs = 1:2))
            warning("formula must not have more than two RHS parts")
        }
        simple_formula <- FALSE
    }
    mf$formula <- formula
}

    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())

if (modele %in% c("pls-glm-family","pls-glm-Gamma","pls-glm-gaussian","pls-glm-inverse.gaussian","pls-glm-logistic","pls-glm-poisson")) {
if(match("method",names(call), 0L)==0L){method<-"glm.fit"}
}
if (modele %in% c("pls-glm-polr")) {
if(match("method",names(call), 0L)==0L){method<-"logistic"} else {if(!(call$method %in% c("logistic", "probit", "cloglog", "cauchit"))) {method<-"logistic"}}
}

if (modele %in% c("pls-glm-family","pls-glm-Gamma","pls-glm-gaussian","pls-glm-inverse.gaussian","pls-glm-logistic","pls-glm-poisson","pls","pls-glm-polr")) {
    mt <- attr(mf, "terms")
    attr(mt,"intercept")<-0L
    dataY <- model.response(mf, "any")
    if (length(dim(dataY)) == 1L) {
        nm <- rownames(dataY)
        dim(dataY) <- NULL
        if (!is.null(nm)) names(dataY) <- nm
        }
    dataX <- if (!is.empty.model(mt)) model.matrix(mt, mf, contrasts)
        else matrix(, NROW(dataY), 0L)
    weights <- as.vector(model.weights(mf))
    if (!is.null(weights) && !is.numeric(weights)) stop("'weights' must be a numeric vector")
    if (!is.null(weights) && any(weights < 0)) stop("negative weights not allowed")
    offset <- as.vector(model.offset(mf))
    if (!is.null(offset)) {
        if (length(offset) != NROW(dataY)) stop(gettextf("number of offsets is %d should equal %d (number of observations)", length(offset), NROW(dataY)), domain = NA)
        }
    }
if (modele %in% "pls-beta") {
mt <- terms(formula, data = data)
mtX <- terms(formula, data = data, rhs = 1L)
mtZ <- delete.response(terms(formula, data = data, rhs = 2L))
attr(mtX,"intercept")<-0L
attr(mtZ,"intercept")<-0L
dataY <- model.response(mf, "any")
dataX <- if (!is.empty.model(mtX)) model.matrix(mtX, mf) else matrix(, NROW(dataY), 0L)
#if (!is.empty.model(mtX)) model.matrix(mtX, mf, contrasts) else matrix(, NROW(dataY), 0L)
dataZ <- if (!is.empty.model(mtZ)) model.matrix(mtZ, mf) else matrix(, NROW(dataY), 0L)
#if (!is.empty.model(mtZ)) model.matrix(mtZ, mf, contrasts) else matrix(, NROW(dataY), 0L)
    if (length(dataY) < 1) 
        stop("empty model")
    if (!(min(dataY) > 0 & max(dataY) < 1)) 
        stop("invalid dependent variable, all observations must be in (0, 1)")
    n <- length(dataY)
    weights <- model.weights(mf)
    if (!is.null(weights) && !is.numeric(weights)) stop("'weights' must be a numeric vector")
    if (!is.null(weights) && any(weights < 0)) stop("negative weights not allowed")
    if (is.null(weights)) 
        weights <- 1
    if (length(weights) == 1) 
        weights <- rep.int(weights, n)
    weights <- as.vector(weights)
    names(weights) <- rownames(mf)
    expand_offset <- function(offset) {
        if (is.null(offset)) 
            offset <- 0
        if (length(offset) == 1) 
            offset <- rep.int(offset, n)
        as.vector(offset)
    }
    offsetX <- expand_offset(model.offset(model.part(formula, 
        data = mf, rhs = 1L, terms = TRUE)))
    offsetZ <- expand_offset(model.offset(model.part(formula, 
        data = mf, rhs = 2L, terms = TRUE)))
    if (!is.null(mf$offset)) 
        offsetX <- offsetX + expand_offset(mf[, "(offset)"])
    offset <- list(mean = offsetX, precision = offsetZ)
}
         
    if (!is.data.frame(dataX)) {dataX <- data.frame(dataX)}
    res <- NULL
    res$nr <- nrow(dataX)
        if (K > res$nr) {
          if(verbose){cat(paste("K cannot be > than nrow(dataX) =",res$nr,"\n"))}
          if(verbose){cat(paste("K is set to", nrow(dataX), "\n"))}
            K <- res$nr
            random = FALSE
        }
    call <- match.call(expand.dots=FALSE)
    
    
if (modele %in% c("pls-glm-family","pls-glm-Gamma","pls-glm-gaussian","pls-glm-inverse.gaussian","pls-glm-logistic","pls-glm-poisson")) {
if(match("method",names(call), 0L)==0L){method<-"glm.fit"}
}
if (modele %in% c("pls-glm-polr")) {
if(match("method",names(call), 0L)==0L){method<-"logistic"} else {if(!(call$method %in% c("logistic", "probit", "cloglog", "cauchit"))) {method<-"logistic"}}
}

    
    if (is.null(modele) & !is.null(family)) {modele<-"pls-glm-family"}
    if (!(modele %in% c("pls","pls-glm-logistic","pls-glm-family","pls-glm-Gamma","pls-glm-gaussian","pls-glm-inverse.gaussian","pls-glm-poisson","pls-glm-polr","pls-beta"))) {print(modele);stop("'modele' not recognized")}
    if (!(modele %in% "pls-glm-family") & !is.null(family)) {stop("Set 'modele=pls-glm-family' to use the family option")}
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
if (is.null(link)){link<-"logit"} else {if(!(link %in% c("logit", "probit", "cloglog", "cauchit", "log", "loglog")) & !is(link,"link-glm")) {link<-"logit"}}
    if (modele %in% c("pls-glm-family","pls-glm-Gamma","pls-glm-gaussian","pls-glm-inverse.gaussian","pls-glm-logistic","pls-glm-poisson")) {if(verbose){cat(family,"\n\n")}}
    if (modele %in% c("pls-glm-polr")) {if(verbose){cat("\nModel:", modele, "\n");cat("Method:", method, "\n\n")}}
    if (modele=="pls-beta") {if(verbose){cat("\nModel:", modele, "\n\n");cat("Link:", link, "\n\n");cat("Link.phi:", link.phi, "\n\n");cat("Type:", type, "\n\n")}}
    if (modele=="pls") {if(verbose){cat("\nModel:", modele, "\n\n")}}


    if (as.character(call["tol_Xi"])=="NULL") {call$tol_Xi <- 10^(-12)}
    if (as.character(call["modele"])=="NULL") {call$modele <- "pls"}
    if (as.character(call["limQ2set"])=="NULL") {call$limQ2set <- .0975}
    if (as.character(call["start"])=="NULL") {call$start <- NULL}
    if (as.character(call["method"])=="NULL") {call$method <- NULL}
    if (as.character(call["control"])=="NULL") {call$control <- list()}
    if (as.character(call["contrasts"])=="NULL") {call$contrasts <- NULL}
    if (as.character(call["sparse"])=="NULL") {call$sparse <- FALSE}
    if (as.character(call["sparseStop"])=="NULL") {call$sparseStop <- FALSE}
    if (as.character(call["naive"])=="NULL") {call$contrasts <- FALSE}
    if (as.character(call["link"])=="NULL") {call$link <- "logit"}
    if (as.character(call["link.phi"])=="NULL") {call$link.phi <- NULL}
    if (as.character(call["type"])=="NULL") {call$type <- "ML"}
    
    
    
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
                mf2 <- match.call(expand.dots = FALSE)
                m2 <- match(c("nt","modele","family","scaleX","scaleY","keepcoeffs","tol_Xi","weights","subset","start","etastart","mustart","offset","control","method","contrasts","sparse","sparseStop","naive","link","link.phi","type","verbose"), names(mf2), 0L)
                mf2 <- mf2[c(1L, m2)]
                mf2[[1L]] <- as.name("PLS_beta_wvc")
                mf2$family <- family
                mf2$weights <- weights
                mf2$dataY <- dataY
                mf2$dataX <- dataX
                mf2$nt <- eval(mf2$nt,parent.frame())
                mf2$dataPredictY <- dataX
if (modele %in% c("pls-glm-family","pls-glm-Gamma","pls-glm-gaussian","pls-glm-inverse.gaussian","pls-glm-logistic","pls-glm-poisson")) {
if(match("method",names(call), 0L)==0L){mf2$method<-"glm.fit"}
}
if (modele %in% c("pls-glm-polr")) {
if(match("method",names(call), 0L)==0L){mf2$method<-"logistic"} else {if(!(call$method %in% c("logistic", "probit", "cloglog", "cauchit"))) {mf2$method<-"logistic"}}
}
                temptemp <- eval(mf2, parent.frame())
                respls_kfolds[[nnkk]][[ii]] <- temptemp$valsPredict
                if(!NoWeights) {attr(respls_kfolds[[nnkk]],"XWeights")=weights; attr(respls_kfolds[[nnkk]],"YWeights")=NULL}
                if (keepcoeffs) {coeffskfolds[[nnkk]][[ii]] = temptemp$coeffs}
                if (keepdataY) {dataY_kfolds[[nnkk]][[ii]] = NULL}
                if (modele=="pls-beta") {respls_kfolds_phi[[nnkk]][[ii]] = temptemp$valsPredictPhis}
                }
            else {
              if(verbose){cat(paste(ii,"\n"))}
                  mf2 <- match.call(expand.dots = FALSE)
                  m2 <- match(c("nt","modele","family","scaleX","scaleY","keepcoeffs","tol_Xi","weights","subset","start","etastart","mustart","offset","control","method","contrasts","sparse","sparseStop","naive","link","link.phi","type","verbose"), names(mf2), 0L)
                  mf2 <- mf2[c(1L, m2)]
                  mf2[[1L]] <- as.name("PLS_beta_wvc")
                  mf2$family <- family
                  mf2$weights <- weights[-nofolds]
                  mf2$dataY <- dataY[-nofolds]
                  mf2$dataX <- dataX[-nofolds,]
                  mf2$nt <- eval(mf2$nt,parent.frame())
                  mf2$dataPredictY <- dataX[nofolds,]
                  temptemp <- eval(mf2, parent.frame())
                  respls_kfolds[[nnkk]][[ii]] <- temptemp$valsPredict
                  if(!NoWeights) {attr(respls_kfolds[[nnkk]][[ii]],"XWeights")=weights[-nofolds]; attr(respls_kfolds[[nnkk]][[ii]],"YWeights")=weights[nofolds]}
                  if (keepcoeffs) {coeffs_kfolds[[nnkk]][[ii]] = temptemp$coeffs}
                  if (keepdataY) {dataY_kfolds[[nnkk]][[ii]] = dataY[nofolds]}
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
results$call$nt <- mf2$nt
return(results)
}
