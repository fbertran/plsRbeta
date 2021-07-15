## ----setup, include = FALSE---------------------------------------------------
#file.edit(normalizePath("~/.Renviron"))
LOCAL <- identical(Sys.getenv("LOCAL"), "TRUE")
#LOCAL=TRUE
knitr::opts_chunk$set(purl = LOCAL)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.height = 3, 
  fig.width = 5
)

## -----------------------------------------------------------------------------
library(plsRglm)
library(plsRbeta)

## ---- eval=FALSE--------------------------------------------------------------
#  data(TxTum)
#  
#  TxTum.mod.bootBR6 <- bootplsbeta(plsRbeta(formula=I(CELTUMCO/100)~.,data=TxTum,nt=6,modele="pls-beta",type="BR"), sim="ordinary", stype="i", R=1000)
#  #save(TxTum.mod.bootBR6,file="TxTum.mod.bootBR6.Rdata")

## ---- warnings=FALSE----------------------------------------------------------
data("TxTum.mod.bootBR6")
temp.ci <- suppressWarnings(plsRglm::confints.bootpls(TxTum.mod.bootBR6))

## ---- fig.cap="Figure 1. Boxplot of bootstrap distributions, 6 components BR."----
boxplots.bootpls(TxTum.mod.bootBR6,indices = 2:nrow(temp.ci))

## ---- fig.cap="Figure 2. Bootstrap 95% confidence intervals, 6 components BR."----
plsRglm::plots.confints.bootpls(temp.ci,prednames=TRUE,indices = 2:nrow(temp.ci))

## -----------------------------------------------------------------------------
ind_BCa_nt6BR <- (temp.ci[,7]<0&temp.ci[,8]<0)|(temp.ci[,7]>0&temp.ci[,8]>0)
rownames(temp.ci)[ind_BCa_nt6BR]

## ---- fig.cap="Figure 3. Bootstrap 95% confidence intervals, BCa significant variables, 6 components BR."----
plsRglm::plots.confints.bootpls(temp.ci,prednames=TRUE,indices = (2:nrow(temp.ci))[ind_BCa_nt6BR[-1]],articlestyle=FALSE,legendpos="topright")

## -----------------------------------------------------------------------------
#TxTum <- read.table("MET_rev_edt.txt",sep="\t",header=TRUE)

data("TxTum.mod.bootBC1")

temp.ci <- plsRglm::confints.bootpls(TxTum.mod.bootBC1)

data("ind_BCa_nt1BC")
data("ind_BCa_nt2BC")
data("ind_BCa_nt3BC")
data("ind_BCa_nt4BC")
data("ind_BCa_nt5BC")
data("ind_BCa_nt6BC")
data("ind_BCa_nt1BR")
data("ind_BCa_nt2BR")
data("ind_BCa_nt3BR")
data("ind_BCa_nt4BR")
data("ind_BCa_nt5BR")
data("ind_BCa_nt6BR")

rownames(temp.ci)[ind_BCa_nt1BC]
(colnames(TxTum.mod.bootBC1$data)[-1])[]

indics <- cbind(ind_BCa_nt1BC,
ind_BCa_nt2BC,
ind_BCa_nt3BC,
ind_BCa_nt4BC,
ind_BCa_nt5BC,
ind_BCa_nt6BC,
ind_BCa_nt1BR,
ind_BCa_nt2BR,
ind_BCa_nt3BR,
ind_BCa_nt4BR,
ind_BCa_nt5BR,
ind_BCa_nt6BR)

nbpreds <- rbind(colSums(indics[,1:6]),colSums(indics[,7:12]))
colnames(nbpreds) <- c("1","2","3","4","5","6")
rownames(nbpreds) <- c("BC","BR")

## -----------------------------------------------------------------------------
nbpreds

## -----------------------------------------------------------------------------
inone <- indics[rowSums(indics)>0,c(1,7,2,8,3,9,4,10,5,11,6,12)]
inone <- t(apply(inone,1,as.numeric))
rownames(inone)
colnames(inone) <- c("BC1","BR1","BC2","BR2","BC3","BR3","BC4","BR4","BC5","BR5","BC6","BR6")
colnames(inone) <- c("1 BC","  BR","2 BC","  BR","3 BC","  BR","4 BC","  BR","5 BC","  BR","6 BC","  BR")

## -----------------------------------------------------------------------------
inone

## ---- fig.cap="Figure 4. Selected variables using the bootstrap BCa technique."----
library(bipartite)
bipartite::visweb(t(inone),type="None",labsize=2,square="b",box.col="grey25",pred.lablength=7)

## -----------------------------------------------------------------------------
data(colon)
orig <- colon

## ---- eval=FALSE--------------------------------------------------------------
#  modpls.boot3 <- bootplsbeta(plsRbeta(X..Cellules.tumorales~.,data=colon,nt=3,modele="pls-beta"), sim="ordinary", stype="i", R=250)
#  #save(modpls.boot3,file="modpls.boot_nt3.Rdata")

## -----------------------------------------------------------------------------
data("modpls.boot_nt3")
temp.ci <- suppressWarnings(plsRglm::confints.bootpls(modpls.boot3))
ind_BCa_nt3 <- (temp.ci[,7]<0&temp.ci[,8]<0)|(temp.ci[,7]>0&temp.ci[,8]>0)
#save(ind_BCa_nt3,file="ind_BCa_nt3.Rdata")

## ---- fig.cap="Figure 5. Boxplots of the bootstrap distribution of the coefficients of the predictors, 3 components BR"----
boxplots.bootpls(modpls.boot3,indices = 2:nrow(temp.ci))

## ---- fig.cap="Figure 6. Bootstrap 95% confidence intervals, 3 components BR"----
plsRglm::plots.confints.bootpls(temp.ci,prednames=TRUE,indices = 2:nrow(temp.ci))

## -----------------------------------------------------------------------------
data(file="ind_BCa_nt3")
rownames(temp.ci)[ind_BCa_nt3]

## ---- fig.cap="Figure 7. Bootstrap 95% confidence intervals, BCa significant variables, 3 components BR."----
plsRglm::plots.confints.bootpls(temp.ci,prednames=TRUE,indices = (2:nrow(temp.ci))[ind_BCa_nt3[-1]],articlestyle=FALSE)

## ---- fig.cap="Figure 8. BCa Bootstrap 95% confidence intervals, BCa significant variables, 3 components BR."----
plsRglm::plots.confints.bootpls(temp.ci,prednames=TRUE,indices = (2:nrow(temp.ci))[ind_BCa_nt3[-1]],articlestyle=FALSE,typeIC="BCa")

## -----------------------------------------------------------------------------
data("ind_BCa_nt3")
colon_sub4 <- colon[c(TRUE,ind_BCa_nt3[-1])]

## ---- eval=FALSE--------------------------------------------------------------
#  modpls_sub4 <- plsRbeta(X..Cellules.tumorales~.,data=colon_sub4,nt=10,modele="pls-beta")
#  #save(modpls_sub4,file="modpls_sub4.Rdata")

## -----------------------------------------------------------------------------
data("modpls_sub4")
modpls_sub4

## ---- fig.cap="Figure 9. Residuals index plot.", fig.keep='last'--------------
# Index plot
sfit <- modpls_sub4$FinalModel;
rd <- resid(sfit, type="sweighted2");
plot(seq(1,sfit$n), rd, xlab="Subject Index", ylab="Std Weigthed Resid 2",
     main="Index Plot");
abline(h=0, lty=3);

## -----------------------------------------------------------------------------
library(betareg)
# A Half-normal Plot with simulated envelop
phat <- predict(sfit, type="response");
phihat <- predict(sfit, type="precision");
tt <- modpls_sub4$tt
n <- sfit$n
absrds <- list();
for (i in 1:19)
  {
     tYwotNA <- rbeta(sfit$n, phat*phihat, (1-phat)*phihat);
     tsfit <- betareg(tYwotNA~tt, x=TRUE);
     trd <- residuals(tsfit, type="sweighted2");
     absrds[[i]]<-sort(abs(trd));
  }
lower <- upper <- middle <- rep(0, n);
for (j in 1:n) {
   min <- max <- sum <- absrds[[1]][j];
   min; max;
   for (i in 2:19) {
      if (min > absrds[[i]][j]) min <- absrds[[i]][j];
      if (max < absrds[[i]][j]) max <- absrds[[i]][j];
      sum <- sum + absrds[[i]][j];
   }
   lower[j] <- min;
   upper[j] <- max;
   middle[j] <- sum / 19;
}

qnval <- qnorm((1:n+n-1/8) / (2*n+1/2));

## ----simulated_envelope_submod4, fig.cap="Figure 10. Residuals and simulated envelop.", fig.keep='last'----
plot(qnval, sort(abs(rd)), ylim=range(0,3.5), xlab="Expected", ylab="Observed",
     main="Half-Normal Plot with Simulated Envelope");
lines(qnval, lower, lty=3);
lines(qnval, upper, lty=3);
lines(qnval, middle, lty=1, lwd=2);

## ----mod_sub4_nt2_12_beta, fig.cap="Figure 11. Representation of individuals on the plane formed by the first two components."----
plot(modpls_sub4$tt[,1],modpls_sub4$tt[,2],col=plotrix::color.scale(colon[,1],c(0,1,1),c(1,1,0),0),pch=16,ylab="",xlab="")

## ----plot_logit_xpred_yorig_submod4, fig.cap="Figure 12. Observed values versus predicted values on the logit scale.", fig.keep='last'----
plot(binomial()$linkfun(modpls_sub4$ValsPredictY),binomial()$linkfun(orig[,1]),xlab="Valeurs prédites",ylab="Valeurs observées")
abline(0,1,lwd=2,lty=2)

## ----plot_xpred_yorig_submod4, fig.cap="Figure 13. Observed values versus predicted values on the original scale.", fig.keep='last'----
plot(modpls_sub4$ValsPredictY,orig[,1],xlim=c(0,1),ylim=c(0,1))
abline(0,1,lwd=2,lty=2)

