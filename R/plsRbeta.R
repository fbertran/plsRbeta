#' Partial least squares Regression beta regression models
#' 
#' This function implements Partial least squares Regression generalized linear
#' models complete or incomplete datasets.
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
#' The default estimator for Degrees of Freedom is the Kramer and Sugiyama's
#' one which only works for classical plsR models. For these models,
#' Information criteria are computed accordingly to these estimations. Naive
#' Degrees of Freedom and Information Criteria are also provided for comparison
#' purposes. For more details, see Kraemer, N., Sugiyama M. (2010). "The
#' Degrees of Freedom of Partial Least Squares Regression". preprint,
#' http://arxiv.org/abs/1002.4112.
#' 
#' @aliases plsRbeta plsRbetamodel.default plsRbetamodel.formula
#' @param object a response (training) dataset or an object of
#' class "\code{\link{formula}}" (or one that can be coerced to 
#' that class): a symbolic description of the model to be fitted.
#' @param dataX predictor(s) (training) dataset
#' @param formula 
#' The details of model specification are given under 'Details'.
#' @param data an optional data frame, list or environment (or object coercible
#' by \code{\link{as.data.frame}} to a data frame) containing the variables in
#' the model. If not found in \code{data}, the variables are taken from
#' \code{environment(formula)}, typically the environment from which
#' \code{plsRbeta} is called.
#' @param nt number of components to be extracted
#' @param limQ2set limit value for the Q2
#' @param dataPredictY predictor(s) (testing) dataset
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
#' @param typeVC type of leave one out cross validation. For back compatibility
#' purpose.  \describe{ \item{list("none")}{no cross validation}
#' \item{list("standard")}{no cross validation} \item{list("missingdata")}{no
#' cross validation} \item{list("adaptative")}{no cross validation} }
#' @param EstimXNA only for \code{modele="pls"}. Set whether the missing X
#' values have to be estimated.
#' @param scaleX scale the predictor(s) : must be set to TRUE for
#' \code{modele="pls"} and should be for glms pls.
#' @param scaleY scale the response : Yes/No. Ignored since non always possible
#' for glm responses.
#' @param pvals.expli should individual p-values be reported to tune model
#' selection ?
#' @param alpha.pvals.expli level of significance for predictors when
#' pvals.expli=TRUE
#' @param MClassed number of missclassified cases, should only be used for
#' binary responses
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
#' @param method the method to be used in fitting the model. The default method
#' \code{"glm.fit"} uses iteratively reweighted least squares (IWLS).
#' User-supplied fitting functions can be supplied either as a function or a
#' character string naming a function, with a function which takes the same
#' arguments as \code{glm.fit}.
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
#' @param \dots arguments to pass to \code{plsRmodel.default} or to
#' \code{plsRmodel.formula}
#' @return Depends on the model that was used to fit the model.
#' @note Use \code{plsRbeta} instead.
#' @author Frédéric Bertrand\cr
#' \email{frederic.bertrand@@utt.fr}\cr
#' \url{https://fbertran.github.io/homepage/}
#' @seealso \code{\link[plsRglm]{plsR}} and \code{\link[plsRglm]{plsRglm}}
#' @references Frédéric Bertrand, Nicolas Meyer,
#' Michèle Beau-Faller, Karim El Bayed, Izzie-Jacques Namer,
#' Myriam Maumy-Bertrand (2013). Régression Bêta
#' PLS. \emph{Journal de la Société Française de Statistique},
#' \bold{154}(3):143-159.
#' \url{http://publications-sfds.math.cnrs.fr/index.php/J-SFdS/article/view/215}
#' @keywords models regression
#' @examples
#' 
#' 
#' data("GasolineYield",package="betareg")
#' modpls <- plsRbeta(yield~.,data=GasolineYield,nt=3,modele="pls-beta")
#' modpls$pp
#' modpls$Coeffs
#' modpls$Std.Coeffs
#' modpls$InfCrit
#' modpls$PredictY[1,]
#' rm("modpls")
#' 
#' data("GasolineYield",package="betareg")
#' yGasolineYield <- GasolineYield$yield
#' XGasolineYield <- GasolineYield[,2:5]
#' modpls <- plsRbeta(yGasolineYield,XGasolineYield,nt=3,modele="pls-beta")
#' modpls$pp
#' modpls$Coeffs
#' modpls$Std.Coeffs
#' modpls$InfCrit
#' modpls$PredictY[1,]
#' rm("modpls")
#' 
#' 
#' @export plsRbeta
plsRbeta <- function(object, ...) UseMethod("plsRbetamodel")

#' @rdname plsRbeta
#' @aliases plsRbeta
#' @export plsRbetamodel
plsRbetamodel <- plsRbeta
