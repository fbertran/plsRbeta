# Non-parametric tilted bootstrap for PLS beta regression models

Provides a wrapper for the bootstrap function `tilt.boot` from the
`boot` R package.  
Implements non-parametric tilted bootstrap for PLS beta regression
models by case resampling : the `tilt.boot` function will run an initial
bootstrap with equal resampling probabilities (if required) and will use
the output of the initial run to find resampling probabilities which put
the value of the statistic at required values. It then runs an
importance resampling bootstrap using the calculated probabilities as
the resampling distribution.

## Usage

``` r
tilt.bootplsbeta(
  object,
  typeboot = "plsmodel",
  statistic = coefs.plsRbeta,
  R = c(499, 250, 250),
  alpha = c(0.025, 0.975),
  sim = "ordinary",
  stype = "i",
  index = 1,
  stabvalue = 1e+06
)
```

## Arguments

- object:

  An object of class `plsRbetamodel` to bootstrap

- typeboot:

  The type of bootstrap. Either (Y,X) boostrap (`typeboot="plsmodel"`)
  or (Y,T) bootstrap (`typeboot="fmodel_np"`). Defaults to (Y,T)
  resampling.

- statistic:

  A function which when applied to data returns a vector containing the
  statistic(s) of interest. `statistic` must take at least two
  arguments. The first argument passed will always be the original data.
  The second will be a vector of indices, frequencies or weights which
  define the bootstrap sample. Further, if predictions are required,
  then a third argument is required which would be a vector of the
  random indices used to generate the bootstrap predictions. Any further
  arguments can be passed to statistic through the `...` argument.

- R:

  The number of bootstrap replicates. Usually this will be a single
  positive integer. For importance resampling, some resamples may use
  one set of weights and others use a different set of weights. In this
  case `R` would be a vector of integers where each component gives the
  number of resamples from each of the rows of weights.

- alpha:

  The alpha level to which tilting is required. This parameter is
  ignored if `R[1]` is 0 or if `theta` is supplied, otherwise it is used
  to find the values of `theta` as quantiles of the initial uniform
  bootstrap. In this case `R[1]` should be large enough that
  `min(c(alpha, 1-alpha))*R[1] > 5`, if this is not the case then a
  warning is generated to the effect that the `theta` are extreme values
  and so the tilted output may be unreliable.

- sim:

  A character string indicating the type of simulation required.
  Possible values are `"ordinary"` (the default), `"balanced"`,
  `"permutation"`, or `"antithetic"`.

- stype:

  A character string indicating what the second argument of `statistic`
  represents. Possible values of stype are `"i"` (indices - the
  default), `"f"` (frequencies), or `"w"` (weights).

- index:

  The index of the statistic of interest in the output from `statistic`.
  By default the first element of the output of `statistic` is used.

- stabvalue:

  A value to hard threshold bootstrap estimates computed from atypical
  resamplings.

## Value

An object of class "boot".

## References

Frédéric Bertrand, Nicolas Meyer, Michèle Beau-Faller, Karim El Bayed,
Izzie-Jacques Namer, Myriam Maumy-Bertrand (2013). Régression Bêta PLS.
*Journal de la Société Française de Statistique*, **154**(3):143-159.
<https://ojs-test.apps.ocp.math.cnrs.fr/index.php/J-SFdS/article/view/215>

## See also

[`tilt.boot`](https://rdrr.io/pkg/boot/man/boot.html)

## Author

Frédéric Bertrand  
<frederic.bertrand@lecnam.net>  
<https://fbertran.github.io/homepage/>

## Examples

``` r
# \donttest{
data("GasolineYield",package="betareg")
yGasolineYield <- GasolineYield$yield
XGasolineYield <- GasolineYield[,2:5]
modplsRbeta <- plsRbeta(yGasolineYield, XGasolineYield, nt=3, 
modele="pls-beta")
#> ____************************************************____
#> 
#> Model: pls-beta 
#> 
#> Link: logit 
#> 
#> Link.phi: 
#> 
#> Type: ML 
#> 
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Predicting X without NA neither in X nor in Y____
#> ****________________________________________________****
#> 
# GazYield.tilt.boot <- tilt.bootplsbeta(modplsRbeta,
# statistic=coefs.plsRbeta, R=c(499, 100, 100), 
# alpha=c(0.025, 0.975), sim="balanced", stype="i", index=1)
# boxplots.bootpls(GazYield.tilt.boot,1:2)
# }
```
