#' Computes Predicted Chisquare for kfold cross validated partial least squares
#' beta regression models.
#' 
#' This function computes Predicted Chisquare for kfold cross validated partial
#' least squares beta regression models.
#' 
#' 
#' @param pls_kfolds a kfold cross validated partial least squares regression
#' glm model
#' @return \item{list}{Total Predicted Chisquare vs number of components for
#' the first group partition} \item{list()}{\dots{}} \item{list}{Total
#' Predicted Chisquare vs number of components for the last group partition}
#' @note Use \code{\link{PLS_beta_kfoldcv}} to create kfold cross validated
#' partial least squares regression glm and beta models.
#' @author Frédéric Bertrand\cr
#' \email{frederic.bertrand@@lecnam.net}\cr
#' \url{https://fbertran.github.io/homepage/}
#' @seealso \code{\link[plsRglm]{kfolds2coeff}},
#' \code{\link[plsRglm]{kfolds2Press}}, \code{\link[plsRglm]{kfolds2Pressind}},
#' \code{\link{kfolds2Chisqind}}, \code{\link[plsRglm]{kfolds2Mclassedind}} and
#' \code{\link[plsRglm]{kfolds2Mclassed}} to extract and transforms results
#' from kfold cross validation.
#' @references Frédéric Bertrand, Nicolas Meyer,
#' Michèle Beau-Faller, Karim El Bayed, Izzie-Jacques Namer,
#' Myriam Maumy-Bertrand (2013). Régression Bêta
#' PLS. \emph{Journal de la Société Française de Statistique},
#' \bold{154}(3):143-159.
#' \url{https://ojs-test.apps.ocp.math.cnrs.fr/index.php/J-SFdS/article/view/215}
#' @keywords models regression
#' @examples
#' 
#' \dontrun{
#' data("GasolineYield",package="betareg")
#' yGasolineYield <- GasolineYield$yield
#' XGasolineYield <- GasolineYield[,2:5]
#' bbb <- PLS_beta_kfoldcv(yGasolineYield,XGasolineYield,nt=3,modele="pls-beta")
#' kfolds2Chisq(bbb)
#' }
#' 
kfolds2Chisq <- function(pls_kfolds) {
    if (!is.null(pls_kfolds$call$family)) {
        if (is.character(pls_kfolds$call$family)) {pls_kfolds$call$family <- get(pls_kfolds$call$family, mode = "function", envir = parent.frame())}
        if (is.function(pls_kfolds$call$family)) {pls_kfolds$call$family <- pls_kfolds$call$family()}
        if (is.language(pls_kfolds$call$family)) {pls_kfolds$call$family <- eval(pls_kfolds$call$family)}
        fam_var <- pls_kfolds$call$family$variance
        fam_name <- paste(pls_kfolds$call$family$family,pls_kfolds$call$family$link)
    } else {
        if (pls_kfolds$call$modele=="pls") {
            fam_var <- function(vals) {return(1)}
            fam_name <- "pls"
        }
        if (pls_kfolds$call$modele=="pls-glm-polr") {
            fam_name <- "pls-glm-polr"
            Varyy <- function(piVaryy) {
            diag(piVaryy[-length(piVaryy)])-piVaryy[-length(piVaryy)]%*%t(piVaryy[-length(piVaryy)])
            }
            Chisqcomp <- function(yichisq,pichisq) {
            t(yichisq[-length(yichisq)]-pichisq[-length(pichisq)])%*%MASS::ginv(Varyy(pichisq))%*%(yichisq[-length(yichisq)]-pichisq[-length(pichisq)])
            }
            Chiscompmatrix <- function(rowspi,rowsyi) {
            sum(mapply(Chisqcomp,rowsyi,rowspi))
            }
            Chiscompmatrixweight <- function(rowspi,rowsyi) {
            (mapply(Chisqcomp,rowsyi,rowspi))
            }
        }
        if (pls_kfolds$call$modele=="pls-beta") {            
            fam_beta <- function(vals,phis) {return(vals*(1-vals)/(1+phis))}
            fam_name <- "pls-beta"
        }
    }

            if (pls_kfolds$call$modele=="pls-beta") {
    max_nt <- rep(NA,length(pls_kfolds$results_kfolds))
    if (length(pls_kfolds$results_kfolds)==1) {
        max_nt[1] <- min(unlist(lapply(pls_kfolds$results_kfolds[[1]],ncol)))
        preChisq_kfolds <- list(rep(0, max_nt[1]))
    }
    else
    {
      if (length(pls_kfolds$results_kfolds)>1)
      {
      preChisq_kfolds <-vector("list",length(pls_kfolds$results_kfolds))
        for (jj in 1:length(pls_kfolds$results_kfolds))
        {
          max_nt[jj] <- min(unlist(lapply(pls_kfolds$results_kfolds[[jj]],ncol)))
          preChisq_kfolds[[jj]] <- rep(0,max_nt[jj])
        }
      rm(jj)
      }
    }

    for (nnkk in 1:length(pls_kfolds$results_kfolds))
    {
        for (ii in 1:length(pls_kfolds$results_kfolds[[1]]))
        {
            if (dim(pls_kfolds$results_kfolds[[nnkk]][[ii]])[1]==1)
            {
                if(is.null(attr(pls_kfolds$results_kfolds[[nnkk]][[ii]],"YWeights"))){
preChisq_kfolds[[nnkk]] <- preChisq_kfolds[[nnkk]]+(pls_kfolds$dataY_kfolds[[nnkk]][[ii]]-pls_kfolds$results_kfolds[[nnkk]][[ii]][1:max_nt[nnkk]])^2/fam_beta(pls_kfolds$results_kfolds[[nnkk]][[ii]][1:max_nt[nnkk]],pls_kfolds$results_kfolds_phi[[nnkk]][[ii]][1:max_nt[nnkk]])
                    } else {
preChisq_kfolds[[nnkk]] <- preChisq_kfolds[[nnkk]]+attr(pls_kfolds$results_kfolds[[nnkk]][[ii]],"YWeights")*(pls_kfolds$dataY_kfolds[[nnkk]][[ii]]-pls_kfolds$results_kfolds[[nnkk]][[ii]][1:max_nt[nnkk]])^2/fam_beta(pls_kfolds$results_kfolds[[nnkk]][[ii]][1:max_nt[nnkk]],pls_kfolds$results_kfolds_phi[[nnkk]][[ii]][1:max_nt[nnkk]])           
            }            
            }
            else
            {
                if(is.null(attr(pls_kfolds$results_kfolds[[nnkk]][[ii]],"YWeights"))){
                preChisq_kfolds[[nnkk]] <- preChisq_kfolds[[nnkk]]+colSums((apply(pls_kfolds$results_kfolds[[nnkk]][[ii]][,1:max_nt[nnkk],drop=FALSE],2,'-',pls_kfolds$dataY_kfolds[[nnkk]][[ii]]))^2/fam_beta(pls_kfolds$results_kfolds[[nnkk]][[ii]][,1:max_nt[nnkk],drop=FALSE],pls_kfolds$results_kfolds_phi[[nnkk]][[ii]][,1:max_nt[nnkk],drop=FALSE])) 
                    } else {
                preChisq_kfolds[[nnkk]] <- preChisq_kfolds[[nnkk]]+colSums(attr(pls_kfolds$results_kfolds[[nnkk]][[ii]],"YWeights")*(apply(pls_kfolds$results_kfolds[[nnkk]][[ii]][,1:max_nt[nnkk]],2,'-',pls_kfolds$dataY_kfolds[[nnkk]][[ii]]))^2/fam_beta(pls_kfolds$results_kfolds[[nnkk]][[ii]][,1:max_nt[nnkk]],pls_kfolds$results_kfolds_phi[[nnkk]][[ii]][1:max_nt[nnkk]]))             
             }            
            }
        }
    }
rm(ii)
rm(nnkk)
            } 
            
if (pls_kfolds$call$modele %in% c("pls-glm-family","pls-glm-Gamma","pls-glm-gaussian","pls-glm-inverse.gaussian","pls-glm-logistic","pls-glm-poisson")) {
    max_nt <- rep(NA,length(pls_kfolds$results_kfolds))
    if (length(pls_kfolds$results_kfolds)==1) {
        max_nt[1] <- min(unlist(lapply(pls_kfolds$results_kfolds[[1]],ncol)))
        preChisq_kfolds <- list(rep(0, max_nt[1]))
    }
    else
    {
      if (length(pls_kfolds$results_kfolds)>1)
      {
      preChisq_kfolds <-vector("list",length(pls_kfolds$results_kfolds))
        for (jj in 1:length(pls_kfolds$results_kfolds))
        {
          max_nt[jj] <- min(unlist(lapply(pls_kfolds$results_kfolds[[jj]],ncol)))
          preChisq_kfolds[[jj]] <- rep(0,max_nt[jj])
        }
      rm(jj)
      }
    }

    for (nnkk in 1:length(pls_kfolds$results_kfolds))
    {
        for (ii in 1:length(pls_kfolds$results_kfolds[[1]]))
        {
            if (dim(pls_kfolds$results_kfolds[[nnkk]][[ii]])[1]==1)
            {
                if(is.null(attr(pls_kfolds$results_kfolds[[nnkk]][[ii]],"YWeights"))){
preChisq_kfolds[[nnkk]] <- preChisq_kfolds[[nnkk]]+(pls_kfolds$dataY_kfolds[[nnkk]][[ii]]-pls_kfolds$results_kfolds[[nnkk]][[ii]][1:max_nt[nnkk]])^2/(fam_var(pls_kfolds$results_kfolds[[nnkk]][[ii]][1:max_nt[nnkk]]))
                    } else {
preChisq_kfolds[[nnkk]] <- preChisq_kfolds[[nnkk]]+attr(pls_kfolds$results_kfolds[[nnkk]][[ii]],"YWeights")*(pls_kfolds$dataY_kfolds[[nnkk]][[ii]]-pls_kfolds$results_kfolds[[nnkk]][[ii]][1:max_nt[nnkk]])^2/(fam_var(pls_kfolds$results_kfolds[[nnkk]][[ii]][1:max_nt[nnkk]]))            
            }            
            }
            else
            {
                if(is.null(attr(pls_kfolds$results_kfolds[[nnkk]][[ii]],"YWeights"))){
                preChisq_kfolds[[nnkk]] <- preChisq_kfolds[[nnkk]]+colSums((apply(pls_kfolds$results_kfolds[[nnkk]][[ii]][,1:max_nt[nnkk],drop=FALSE],2,'-',pls_kfolds$dataY_kfolds[[nnkk]][[ii]]))^2/(fam_var(pls_kfolds$results_kfolds[[nnkk]][[ii]][,1:max_nt[nnkk],drop=FALSE]))) 
                    } else {
                preChisq_kfolds[[nnkk]] <- preChisq_kfolds[[nnkk]]+colSums(attr(pls_kfolds$results_kfolds[[nnkk]][[ii]],"YWeights")*(apply(pls_kfolds$results_kfolds[[nnkk]][[ii]][,1:max_nt[nnkk]],2,'-',pls_kfolds$dataY_kfolds[[nnkk]][[ii]]))^2/(fam_var(pls_kfolds$results_kfolds[[nnkk]][[ii]][,1:max_nt[nnkk]])))             
             }            
            }
        }
    }
rm(ii)
rm(nnkk)
}

if (pls_kfolds$call$modele=="pls-glm-polr") {
    max_nt <- rep(NA,length(pls_kfolds$results_kfolds))
    if (length(pls_kfolds$results_kfolds)==1) {
        max_nt[1] <- min(unlist(lapply(pls_kfolds$results_kfolds[[1]],length)))
        preChisq_kfolds <- list(rep(0, max_nt[1]))
    }
    else
    {
      if (length(pls_kfolds$results_kfolds)>1)
      {
      preChisq_kfolds <-vector("list",length(pls_kfolds$results_kfolds))
        for (jj in 1:length(pls_kfolds$results_kfolds))
        {
          max_nt[jj] <- min(unlist(lapply(pls_kfolds$results_kfolds[[jj]],length)))
          preChisq_kfolds[[jj]] <- rep(0,max_nt[jj])
        }
      rm(jj)
      }
    }

    if (length(pls_kfolds$results_kfolds)==1) {preChisqind_kfolds <- list(vector("list", length(pls_kfolds$results_kfolds[[1]])))}
    else
    {
      if (length(pls_kfolds$results_kfolds)>1)
      {
      preChisqind_kfolds <-vector("list",length(pls_kfolds$results_kfolds))
        for (jj in 1:length(pls_kfolds$results_kfolds))
        {
          preChisqind_kfolds[[jj]] <-vector("list",length(pls_kfolds$results_kfolds[[jj]]))
        }
      rm(jj)
      }
    }



    for (nnkk in 1:length(pls_kfolds$results_kfolds))
    {
        for (ii in 1:length(pls_kfolds$results_kfolds[[1]]))
        {
                    fff <- ~pls_kfolds$dataY_kfolds[[nnkk]][[ii]]-1
                    m <- model.frame(fff, pls_kfolds$dataY_kfolds[[nnkk]][[ii]])
                    mat <- model.matrix(fff, model.frame(fff, pls_kfolds$dataY_kfolds[[nnkk]][[ii]]))
                    if(is.null(attr(pls_kfolds$results_kfolds[[nnkk]][[ii]],"YWeights"))){
                    preChisqind_kfolds[[nnkk]][[ii]] <- (unlist(lapply(lapply(pls_kfolds$results_kfolds[[nnkk]][[ii]],function(xxx) {as.list(as.data.frame(t(xxx)))}),
                    Chiscompmatrix,as.list(as.data.frame(t(mat))))))[1:max_nt[nnkk]]
                    } else {
                    preChisqind_kfolds[[nnkk]][[ii]] <- (unlist(lapply(lapply(lapply(lapply(pls_kfolds$results_kfolds[[nnkk]][[ii]],function(xxx) {as.list(as.data.frame(t(xxx)))}),
                    Chiscompmatrixweight,as.list(as.data.frame(t(mat)))),"*",attr(pls_kfolds$results_kfolds[[nnkk]][[ii]],"YWeights")),sum)))[1:max_nt[nnkk]]
                    }
                    rm(fff)
                    rm(m)
                    rm(mat)
        }
    }
    for (nnkk in 1:length(pls_kfolds$results_kfolds))
    {
                    preChisq_kfolds[[nnkk]] <- colSums(matrix(unlist(preChisqind_kfolds[[nnkk]]),ncol=max_nt[nnkk],byrow=TRUE))
    }
rm(ii)
rm(nnkk)
}

return(preChisq_kfolds)
}
