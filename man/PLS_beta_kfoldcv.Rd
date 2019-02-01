\name{PLS_beta_kfoldcv}
\alias{PLS_beta_kfoldcv}
\title{Partial least squares regression beta models with kfold cross validation}
\description{
This function implements kfold cross validation on complete or incomplete datasets for partial least squares beta regression models
}
\usage{
PLS_beta_kfoldcv(dataY, dataX, nt = 2, limQ2set = 0.0975, modele = "pls", 
family = NULL, K = nrow(dataX), NK = 1, grouplist = NULL, random = FALSE, 
scaleX = TRUE, scaleY = NULL, keepcoeffs = FALSE, keepfolds = FALSE, 
keepdataY = TRUE, keepMclassed=FALSE, tol_Xi = 10^(-12), weights, 
method,link=NULL,link.phi=NULL,type="ML",verbose=TRUE)
}
\arguments{
  \item{dataY}{response (training) dataset}
  \item{dataX}{predictor(s) (training) dataset}
  \item{nt}{number of components to be extracted}
  \item{limQ2set}{limit value for the Q2}
  \item{modele}{name of the PLS glm or PLS beta model to be fitted (\code{"pls"}, \code{"pls-glm-Gamma"}, \code{"pls-glm-gaussian"}, \code{"pls-glm-inverse.gaussian"}, \code{"pls-glm-logistic"}, \code{"pls-glm-poisson"}, \code{"pls-glm-polr"}, \code{"pls-beta"}). Use \code{"modele=pls-glm-family"} to enable the \code{family} option.}
  \item{family}{a description of the error distribution and link function to be used in the model. This can be a character string naming a family function, a family function or the result of a call to a family function. (See \code{\link[stats]{family}} for details of family functions.) To use the family option, please set \code{modele="pls-glm-family"}. User defined families can also be defined. See details.}
  \item{K}{number of groups}
  \item{NK}{number of times the group division is made}
  \item{grouplist}{to specify the members of the \code{K} groups}
  \item{random}{should the \code{K} groups be made randomly}
  \item{scaleX}{scale the predictor(s) : must be set to TRUE for \code{modele="pls"} and should be for glms pls.}
  \item{scaleY}{scale the response : Yes/No. Ignored since non always possible for glm responses.}
  \item{keepcoeffs}{shall the coefficients for each model be returned}
  \item{keepfolds}{shall the groups' composition be returned}
  \item{keepdataY}{shall the observed value of the response for each one of the predicted value be returned}
  \item{keepMclassed}{shall the number of miss classed be returned (unavailable)}
  \item{tol_Xi}{minimal value for Norm2(Xi) and \eqn{\mathrm{det}(pp' \times pp)}{det(pp'*pp)} if there is any missing value in the \code{dataX}. It defaults to \eqn{10^{-12}}{10^{-12}}}
  \item{weights}{an optional vector of 'prior weights' to be used in the fitting process. Should be \code{NULL} or a numeric vector.}
  \item{method}{logistic, probit, complementary log-log or cauchit (corresponding to a Cauchy latent variable).}
  \item{link}{character specification of the link function in the mean model (mu). Currently, "\code{logit}", "\code{probit}", "\code{cloglog}", "\code{cauchit}", "\code{log}", "\code{loglog}" are supported. Alternatively, an object of class "\code{link-glm}" can be supplied.}
  \item{link.phi}{character specification of the link function in the precision model (phi). Currently, "\code{identity}", "\code{log}", "\code{sqrt}" are supported. The default is "\code{log}" unless \code{formula} is of type \code{y~x} where the default is "\code{identity}" (for backward compatibility). Alternatively, an object of class "\code{link-glm}" can be supplied.}
  \item{type}{character specification of the type of estimator. Currently, maximum likelihood ("\code{ML}"), ML with bias correction ("\code{BC}"), and ML with bias reduction ("\code{BR}") are supported.}
  \item{verbose}{should info messages be displayed ?}
}
\details{
Predicts 1 group with the \code{K-1} other groups. Leave one out cross validation is thus obtained for \code{K==nrow(dataX)}.

There are seven different predefined models with predefined link functions available :
\describe{
\item{\code{"pls"}}{ordinary pls models}
\item{\code{"pls-glm-Gamma"}}{glm gaussian with inverse link pls models}
\item{\code{"pls-glm-gaussian"}}{glm gaussian with identity link pls models}
\item{\code{"pls-glm-inverse-gamma"}}{glm binomial with square inverse link pls models}
\item{\code{"pls-glm-logistic"}}{glm binomial with logit link pls models}
\item{\code{"pls-glm-poisson"}}{glm poisson with log link pls models}
\item{\code{"pls-glm-polr"}}{glm polr with logit link pls models}
}
Using the \code{"family="} option and setting \code{"modele=pls-glm-family"} allows changing the family and link function the same way as for the \code{\link[stats]{glm}} function. As a consequence user-specified families can also be used. 
\describe{
\item{The \code{gaussian} family}{accepts the links (as names) \code{identity}, \code{log} and \code{inverse}.}
\item{The \code{binomial} family}{accepts the links \code{logit}, \code{probit}, \code{cauchit}, (corresponding to logistic, normal and Cauchy CDFs respectively) \code{log} and \code{cloglog} (complementary log-log).} 
\item{The \code{Gamma} family}{accepts the links \code{inverse}, \code{identity} and \code{log}.}
\item{The \code{poisson} family}{accepts the links \code{log}, \code{identity}, and \code{sqrt}.}
\item{The \code{inverse.gaussian} family}{accepts the links \code{1/mu^2}, \code{inverse}, \code{identity} and \code{log}.}
\item{The \code{quasi} family}{accepts the links \code{logit}, \code{probit}, \code{cloglog}, \code{identity}, \code{inverse}, \code{log}, \code{1/mu^2} and \code{sqrt}.}
\item{The function \code{power}}{can be used to create a power link function.} 
}

Non-NULL weights can be used to indicate that different observations have different dispersions (with the values in weights being inversely proportional to the dispersions); or equivalently, when the elements of weights are positive integers w_i, that each response y_i is the mean of w_i unit-weight observations. 
}
\value{
  \item{results_kfolds}{list of \code{NK}. Each element of the list sums up the results for a group division:
  \describe{
    \item{}{list of \code{K} matrices of size about \code{nrow(dataX)/K * nt} with the predicted values for a growing number of components}
    \item{\dots}{\dots}
    \item{}{list of \code{K} matrices of size about \code{nrow(dataX)/K * nt} with the predicted values for a growing number of components}
    }}
  \item{folds}{list of \code{NK}. Each element of the list sums up the informations for a group division:
  \describe{
    \item{}{list of \code{K} vectors of length about \code{nrow(dataX)} with the numbers of the rows of \code{dataX} that were used as a training set}
    \item{\dots}{\dots}
    \item{}{list of \code{K} vectors of length about \code{nrow(dataX)} with the numbers of the rows of \code{dataX} that were used as a training set}
    }
  }
  \item{dataY_kfolds}{list of \code{NK}. Each element of the list sums up the results for a group division:
  \describe{
    \item{}{list of \code{K} matrices of size about \code{nrow(dataX)/K * 1} with the observed values of the response}
    \item{\dots}{\dots}
    \item{}{list of \code{K} matrices of size about \code{nrow(dataX)/K * 1} with the observed values of the response}
    }
  }
  \item{call}{the call of the function}
}
\references{
\enc{Frédéric}{Fr\'ed\'eric} Bertrand, Nicolas Meyer, \enc{Michèle}{Mich\`ele} Beau-Faller, Karim El Bayed, Izzie-Jacques Namer, Myriam Maumy-Bertrand (2013). \enc{Régression}{R\'egression} \enc{Bêta}{B\^eta} PLS. \emph{Journal de la \enc{Société}{Soci\'et'e} Française de Statistique}, \bold{154}(3):143-159.
\url{http://publications-sfds.math.cnrs.fr/index.php/J-SFdS/article/view/215}
}
\author{\enc{Frédéric}{Fr\'ed\'eric} Bertrand\cr
\email{frederic.bertrand@math.unistra.fr}\cr
\url{http://www-irma.u-strasbg.fr/~fbertran/}
}
\note{Works for complete and incomplete datasets.}
\seealso{\code{\link[plsRglm]{kfolds2coeff}}, \code{\link[plsRglm]{kfolds2Pressind}}, \code{\link[plsRglm]{kfolds2Press}}, \code{\link[plsRglm]{kfolds2Mclassedind}}, \code{\link[plsRglm]{kfolds2Mclassed}} and \code{\link[plsRbeta]{kfolds2CVinfos_beta}} to extract and transform results from kfold cross validation.}
\examples{
\dontrun{
data("GasolineYield",package="betareg")
yGasolineYield <- GasolineYield$yield
XGasolineYield <- GasolineYield[,2:5]
bbb <- PLS_beta_kfoldcv(yGasolineYield,XGasolineYield,nt=3,modele="pls-beta")
kfolds2CVinfos_beta(bbb)
}
}
\keyword{models}
\keyword{regression}
