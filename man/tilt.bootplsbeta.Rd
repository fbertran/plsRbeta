\name{tilt.bootplsbeta}
\alias{tilt.bootplsbeta}
\title{Tilted bootstrap for PLS models}
\description{
  ~~ A (1-5 lines) description of what the function does. ~~
}
\usage{
tilt.bootplsbeta(object, typeboot="plsmodel", statistic=coefs.plsRbeta, 
R=c(499, 250, 250), alpha=c(0.025, 0.975), sim="ordinary", stype="i", 
index=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ ~~Explain \code{object} here~~ }
  \item{typeboot}{ ~~Explain \code{typeboot} here~~ }
  \item{statistic}{ ~~Explain \code{statistic} here~~ }
  \item{R}{ ~~Explain \code{R} here~~ }
  \item{alpha}{ ~~Explain \code{alpha} here~~ }
  \item{sim}{ ~~Explain \code{sim} here~~ }
  \item{stype}{ ~~Explain \code{stype} here~~ }
  \item{index}{ ~~Explain \code{index} here~~ }
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
\url{http://publications-sfds.math.cnrs.fr/index.php/J-SFdS/article/view/215}
}
\author{\enc{Frédéric}{Fr\'ed\'eric} Bertrand\cr
\email{frederic.bertrand@math.unistra.fr}\cr
\url{http://www-irma.u-strasbg.fr/~fbertran/}
}
\note{ ~~some notes~~ 
}
\seealso{\code{\link[boot:boot]{tilt.boot}}}
\examples{\donttest{
data("GasolineYield",package="betareg")

GazYield.tilt.boot <- tilt.bootplsbeta(plsRbeta(yield~.,data=GasolineYield,nt=3,
modele="pls-beta"), statistic=coefs.plsRbeta, R=c(499, 100, 100), 
alpha=c(0.025, 0.975), sim="balanced", stype="i", index=1)
boxplots.bootpls(GazYield.tilt.boot,1:2)

}
}\keyword{models}
