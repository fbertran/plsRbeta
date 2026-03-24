# Partial least squares beta regression models

This function implements Partial least squares beta regression models on
complete or incomplete datasets.

## Usage

``` r
PLS_beta(
  dataY,
  dataX,
  nt = 2,
  limQ2set = 0.0975,
  dataPredictY = dataX,
  modele = "pls",
  family = NULL,
  typeVC = "none",
  EstimXNA = FALSE,
  scaleX = TRUE,
  scaleY = NULL,
  pvals.expli = FALSE,
  alpha.pvals.expli = 0.05,
  MClassed = FALSE,
  tol_Xi = 10^(-12),
  weights,
  method,
  sparse = FALSE,
  sparseStop = TRUE,
  naive = FALSE,
  link = NULL,
  link.phi = NULL,
  type = "ML",
  verbose = TRUE
)
```

## Arguments

- dataY:

  response (training) dataset

- dataX:

  predictor(s) (training) dataset

- nt:

  number of components to be extracted

- limQ2set:

  limit value for the Q2

- dataPredictY:

  predictor(s) (testing) dataset

- modele:

  name of the PLS glm or PLS beta model to be fitted (`"pls"`,
  `"pls-glm-Gamma"`, `"pls-glm-gaussian"`, `"pls-glm-inverse.gaussian"`,
  `"pls-glm-logistic"`, `"pls-glm-poisson"`, `"pls-glm-polr"`,
  `"pls-beta"`). Use `"modele=pls-glm-family"` to enable the `family`
  option.

- family:

  a description of the error distribution and link function to be used
  in the model. This can be a character string naming a family function,
  a family function or the result of a call to a family function. (See
  [`family`](https://rdrr.io/r/stats/family.html) for details of family
  functions.) To use the family option, please set
  `modele="pls-glm-family"`. User defined families can also be defined.
  See details.

- typeVC:

  type of leave one out cross validation. For back compatibility
  purpose.

  list("none")

  :   no cross validation

  list("standard")

  :   no cross validation

  list("missingdata")

  :   no cross validation

  list("adaptative")

  :   no cross validation

- EstimXNA:

  only for `modele="pls"`. Set whether the missing X values have to be
  estimated.

- scaleX:

  scale the predictor(s) : must be set to TRUE for `modele="pls"` and
  should be for glms pls.

- scaleY:

  scale the response : Yes/No. Ignored since not always possible for glm
  responses.

- pvals.expli:

  should individual p-values be reported to tune model selection ?

- alpha.pvals.expli:

  level of significance for predictors when pvals.expli=TRUE

- MClassed:

  number of missclassified cases, should only be used for binary
  responses

- tol_Xi:

  minimal value for Norm2(Xi) and \\\mathrm{det}(pp' \times pp)\\ if
  there is any missing value in the `dataX`. It defaults to \\10^{-12}\\

- weights:

  an optional vector of 'prior weights' to be used in the fitting
  process. Should be `NULL` or a numeric vector.

- method:

  the link function for `pls-glm-polr`, logistic, probit, complementary
  log-log or cauchit (corresponding to a Cauchy latent variable).

- sparse:

  should the coefficients of non-significant predictors
  (\<`alpha.pvals.expli`) be set to 0

- sparseStop:

  should component extraction stop when no significant predictors
  (\<`alpha.pvals.expli`) are found

- naive:

  use the naive estimates for the Degrees of Freedom in plsR? Default is
  `FALSE`.

- link:

  character specification of the link function in the mean model (mu).
  Currently, "`logit`", "`probit`", "`cloglog`", "`cauchit`", "`log`",
  "`loglog`" are supported. Alternatively, an object of class
  "`link-glm`" can be supplied.

- link.phi:

  character specification of the link function in the precision model
  (phi). Currently, "`identity`", "`log`", "`sqrt`" are supported. The
  default is "`log`" unless `formula` is of type `y~x` where the default
  is "`identity`" (for backward compatibility). Alternatively, an object
  of class "`link-glm`" can be supplied.

- type:

  character specification of the type of estimator. Currently, maximum
  likelihood ("`ML`"), ML with bias correction ("`BC`"), and ML with
  bias reduction ("`BR`") are supported.

- verbose:

  should info messages be displayed ?

## Value

Depends on the model that was used to fit the model.

## Details

There are seven different predefined models with predefined link
functions available :

- list("\\pls\\"):

  ordinary pls models

- list("\\pls-glm-Gamma\\"):

  glm gaussian with inverse link pls models

- list("\\pls-glm-gaussian\\"):

  glm gaussian with identity link pls models

- list("\\pls-glm-inverse-gamma\\"):

  glm binomial with square inverse link pls models

- list("\\pls-glm-logistic\\"):

  glm binomial with logit link pls models

- list("\\pls-glm-poisson\\"):

  glm poisson with log link pls models

- list("\\pls-glm-polr\\"):

  glm polr with logit link pls models

Using the `"family="` option and setting `"modele=pls-glm-family"`
allows changing the family and link function the same way as for the
[`glm`](https://rdrr.io/r/stats/glm.html) function. As a consequence
user-specified families can also be used.

- The :

  accepts the links (as names) `identity`, `log` and `inverse`.

- list("gaussian"):

  accepts the links (as names) `identity`, `log` and `inverse`.

- family:

  accepts the links (as names) `identity`, `log` and `inverse`.

- The :

  accepts the links `logit`, `probit`, `cauchit`, (corresponding to
  logistic, normal and Cauchy CDFs respectively) `log` and `cloglog`
  (complementary log-log).

- list("binomial"):

  accepts the links `logit`, `probit`, `cauchit`, (corresponding to
  logistic, normal and Cauchy CDFs respectively) `log` and `cloglog`
  (complementary log-log).

- family:

  accepts the links `logit`, `probit`, `cauchit`, (corresponding to
  logistic, normal and Cauchy CDFs respectively) `log` and `cloglog`
  (complementary log-log).

- The :

  accepts the links `inverse`, `identity` and `log`.

- list("Gamma"):

  accepts the links `inverse`, `identity` and `log`.

- family:

  accepts the links `inverse`, `identity` and `log`.

- The :

  accepts the links `log`, `identity`, and `sqrt`.

- list("poisson"):

  accepts the links `log`, `identity`, and `sqrt`.

- family:

  accepts the links `log`, `identity`, and `sqrt`.

- The :

  accepts the links `1/mu^2`, `inverse`, `identity` and `log`.

- list("inverse.gaussian"):

  accepts the links `1/mu^2`, `inverse`, `identity` and `log`.

- family:

  accepts the links `1/mu^2`, `inverse`, `identity` and `log`.

- The :

  accepts the links `logit`, `probit`, `cloglog`, `identity`, `inverse`,
  `log`, `1/mu^2` and `sqrt`.

- list("quasi"):

  accepts the links `logit`, `probit`, `cloglog`, `identity`, `inverse`,
  `log`, `1/mu^2` and `sqrt`.

- family:

  accepts the links `logit`, `probit`, `cloglog`, `identity`, `inverse`,
  `log`, `1/mu^2` and `sqrt`.

- The function :

  can be used to create a power link function.

- list("power"):

  can be used to create a power link function.

The default estimator for Degrees of Freedom is the Kramer and
Sugiyama's one which only works for classical plsR models. For these
models, Information criteria are computed accordingly to these
estimations. Naive Degrees of Freedom and Information Criteria are also
provided for comparison purposes. For more details, see Kraemer, N.,
Sugiyama M. (2010). "The Degrees of Freedom of Partial Least Squares
Regression". preprint, http://arxiv.org/abs/1002.4112.

## Note

Use `plsRbeta` instead.

## References

Frédéric Bertrand, Nicolas Meyer, Michèle Beau-Faller, Karim El Bayed,
Izzie-Jacques Namer, Myriam Maumy-Bertrand (2013). Régression Bêta PLS.
*Journal de la Société Française de Statistique*, **154**(3):143-159.
<https://ojs-test.apps.ocp.math.cnrs.fr/index.php/J-SFdS/article/view/215>

## See also

[`PLS_beta_wvc`](https://fbertran.github.io/plsRbeta/reference/PLS_beta_wvc.md)
and
[`PLS_beta_kfoldcv`](https://fbertran.github.io/plsRbeta/reference/PLS_beta_kfoldcv.md)

## Author

Frédéric Bertrand  
<frederic.bertrand@lecnam.net>  
<https://fbertran.github.io/homepage/>

## Examples

``` r

data("GasolineYield",package="betareg")
yGasolineYield <- GasolineYield$yield
XGasolineYield <- GasolineYield[,2:5]
modpls <- PLS_beta(yGasolineYield,XGasolineYield,nt=3,modele="pls-beta")
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
modpls$pp
#>             Comp_ 1    Comp_ 2    Comp_ 3
#> gravity   0.4590380 -0.4538663 -2.5188256
#> pressure  0.6395524 -0.4733525  0.6488823
#> temp10   -0.5435643  0.5292108 -1.3295905
#> temp      0.5682795  0.5473174 -0.2156423
modpls$Coeffs
#>                   [,1]
#> Intercept -3.324462308
#> gravity    0.001577508
#> pressure   0.072027686
#> temp10    -0.008398771
#> temp       0.010365973
modpls$Std.Coeffs
#>                   [,1]
#> Intercept -1.547207763
#> gravity    0.008889933
#> pressure   0.188700277
#> temp10    -0.315301401
#> temp       0.723088389
modpls$InfCrit
#>                  AIC        BIC Chi2_Pearson_Y      RSS_Y pseudo_R2_Y      R2_Y
#> Nb_Comp_0  -52.77074  -49.83927       30.72004 0.35640772          NA        NA
#> Nb_Comp_1 -112.87383 -108.47662       30.57369 0.05211039   0.8498691 0.8537900
#> Nb_Comp_2 -136.43184 -130.56889       30.97370 0.02290022   0.9256771 0.9357471
#> Nb_Comp_3 -139.08440 -131.75572       31.08225 0.02022386   0.9385887 0.9432564
modpls$PredictY[1,]
#>   gravity  pressure    temp10      temp 
#>  2.049533  1.686655 -1.371820 -1.821977 
rm("modpls")

```
