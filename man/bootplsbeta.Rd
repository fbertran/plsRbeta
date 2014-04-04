\name{bootplsbeta}
\alias{bootplsbeta}
\title{Non-parametric Bootstrap for PLS generalized linear models}
\description{
Provides a wrapper for the bootstrap function \code{boot} from the \code{boot} R package.\cr Implements non-parametric bootstrap for PLS generalized linear models by case resampling.
}
\usage{
bootplsbeta(object, typeboot="plsmodel", R=250, statistic=coefs.plsRbeta, 
sim="ordinary", stype="i", ...)
}
\arguments{
  \item{object}{ ~~Explain \code{object} here~~ }
  \item{typeboot}{ ~~Explain \code{typeboot} here~~ }
  \item{R}{ ~~Explain \code{R} here~~ }
  \item{statistic}{ ~~Explain \code{statistic} here~~ }
  \item{sim}{ ~~Explain \code{sim} here~~ }
  \item{stype}{ ~~Explain \code{stype} here~~ }
  \item{\dots}{ ~~Explain \code{\dots} here~~ }
}
\details{
  ~~ More details than the description above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{
\enc{Frédéric}{Fr\'ed\'eric} Bertrand, Nicolas Meyer, \enc{Michèle}{Mich\`ele} Beau-Faller, Karim El Bayed, Izzie-Jacques Namer, Myriam Maumy-Bertrand (2013). \enc{Régression}{R\'egression} \enc{Bêta}{B\^eta} PLS. \emph{Journal de la \enc{Société}{Soci\'et'e} Française de Statistique}, \bold{154}(3):143-159.
\url{http://smf4.emath.fr/Publications/JSFdS/154_3/html/}
}
\author{\enc{Frédéric}{Fr\'ed\'eric} Bertrand\cr
\email{frederic.bertrand@math.unistra.fr}\cr
\url{http://www-irma.u-strasbg.fr/~fbertran/}
}
\note{ ~~some notes~~ 
}
\seealso{\code{\link[boot:boot]{boot}}}
\examples{\donttest{
data("GasolineYield",package="betareg")

GazYield.boot <- bootplsbeta(plsRbeta(yield~.,data=GasolineYield,nt=3,
modele="pls-beta"), sim="ordinary", stype="i", R=250)
boot::boot.ci(GazYield.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=1)
boot::boot.ci(GazYield.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=2)
boot::boot.ci(GazYield.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=3)
boot::boot.ci(GazYield.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=4)
boot::boot.ci(GazYield.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=5)
boot::boot.ci(GazYield.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=6)

boxplots.bootpls(GazYield.boot)
confints.bootpls(GazYield.boot)
plots.confints.bootpls(confints.bootpls(GazYield.boot))


plot(GazYield.boot,index=2)
jack.after.boot(GazYield.boot, index=2, useJ=TRUE, nt=3)
plot(GazYield.boot, index=2,jack=TRUE)

aze_compl.tilt.boot <- tilt.bootplsbeta(plsRbeta(yield~.,data=GasolineYield,
nt=3,modele="pls-beta"), statistic=coefs.plsR, R=c(499, 100, 100), 
alpha=c(0.025, 0.975), sim="ordinary", stype="i", index=1)

# PLS bootstrap balanced

GazYield.boot <- bootplsbeta(plsRbeta(yield~.,data=GasolineYield,nt=3,
modele="pls-beta"), sim="balanced", stype="i", R=250)
boot::boot.ci(GazYield.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=1)
boot::boot.ci(GazYield.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=2)
boot::boot.ci(GazYield.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=3)
boot::boot.ci(GazYield.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=4)
boot::boot.ci(GazYield.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=5)
boot::boot.ci(GazYield.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=6)


boxplots.bootpls(GazYield.boot)
confints.bootpls(GazYield.boot)
plots.confints.bootpls(confints.bootpls(GazYield.boot))



plot(GazYield.boot)
jack.after.boot(GazYield.boot, index=1, useJ=TRUE, nt=3)
plot(GazYield.boot,jack=TRUE)

aze_compl.tilt.boot <- tilt.bootplsbeta(plsRbeta(yield~.,data=GasolineYield,
nt=3,modele="pls-beta"), statistic=coefs.plsR, R=c(499, 100, 100), 
alpha=c(0.025, 0.975), sim="balanced", stype="i", index=1)


# PLS permutation bootstrap

GazYield.boot <- bootplsbeta(plsRbeta(yield~.,data=GasolineYield,nt=3,
modele="pls-beta"), sim="permutation", stype="i", R=250)
boot::boot.ci(GazYield.boot, conf = c(0.90,0.95), type = c("norm","basic","perc"), index=1)
boot::boot.ci(GazYield.boot, conf = c(0.90,0.95), type = c("norm","basic","perc"), index=2)
boot::boot.ci(GazYield.boot, conf = c(0.90,0.95), type = c("norm","basic","perc"), index=3)
boot::boot.ci(GazYield.boot, conf = c(0.90,0.95), type = c("norm","basic","perc"), index=4)
boot::boot.ci(GazYield.boot, conf = c(0.90,0.95), type = c("norm","basic","perc"), index=5)
boot::boot.ci(GazYield.boot, conf = c(0.90,0.95), type = c("norm","basic","perc"), index=6)
boxplots.bootpls(GazYield.boot)
plot(GazYield.boot)
}
}\keyword{models}
