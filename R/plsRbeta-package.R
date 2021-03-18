#' @title plsRbeta-package
#'
#' @description Provides Partial least squares Regression for (weighted) beta regression models (Bertrand 2013,  <http://journal-sfds.fr/article/view/215>) and k-fold cross-validation of such models using various criteria. It allows for missing data in the explanatory variables. Bootstrap confidence intervals constructions are also available.
#'
#' @docType package
#' @name plsRbeta-package
#' @references Partial least squares Regression for (weighted) beta regression models (Bertrand 2013,  <http://journal-sfds.fr/article/view/215>), \url{https://github.com/fbertran/plsRbeta/} et \url{https://fbertran.github.io/plsRbeta/}
#' 
#' @examples
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
NULL



