# Partial least squares Regression beta regression models

This function implements Partial least squares Regression generalized
linear models complete or incomplete datasets.

## Usage

``` r
plsRbeta(object, ...)

plsRbetamodel(object, ...)

# Default S3 method
plsRbetamodel(
  object,
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
  verbose = TRUE,
  ...
)

# S3 method for class 'formula'
plsRbetamodel(
  object,
  data = NULL,
  nt = 2,
  limQ2set = 0.0975,
  dataPredictY,
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
  subset,
  start = NULL,
  etastart,
  mustart,
  offset,
  method = "glm.fit",
  control = list(),
  contrasts = NULL,
  sparse = FALSE,
  sparseStop = TRUE,
  naive = FALSE,
  link = NULL,
  link.phi = NULL,
  type = "ML",
  verbose = TRUE,
  ...
)
```

## Arguments

- object:

  a response (training) dataset or an object of class
  "[`formula`](https://rdrr.io/r/stats/formula.html)" (or one that can
  be coerced to that class): a symbolic description of the model to be
  fitted.

- ...:

  arguments to pass to `plsRmodel.default` or to `plsRmodel.formula`

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

  scale the response : Yes/No. Ignored since non always possible for glm
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

  the method to be used in fitting the model. The default method
  `"glm.fit"` uses iteratively reweighted least squares (IWLS).
  User-supplied fitting functions can be supplied either as a function
  or a character string naming a function, with a function which takes
  the same arguments as `glm.fit`.

- sparse:

  should the coefficients of non-significant predictors
  (\<`alpha.pvals.expli`) be set to 0

- sparseStop:

  should component extraction stop when no significant predictors
  (\<`alpha.pvals.expli`) are found

- naive:

  Use the naive estimates for the Degrees of Freedom in plsR? Default is
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

- data:

  an optional data frame, list or environment (or object coercible by
  [`as.data.frame`](https://rdrr.io/r/base/as.data.frame.html) to a data
  frame) containing the variables in the model. If not found in `data`,
  the variables are taken from `environment(formula)`, typically the
  environment from which `plsRbeta` is called.

- subset:

  an optional vector specifying a subset of observations to be used in
  the fitting process.

- start:

  starting values for the parameters in the linear predictor.

- etastart:

  starting values for the linear predictor.

- mustart:

  starting values for the vector of means.

- offset:

  this can be used to specify an *a priori* known component to be
  included in the linear predictor during fitting. This should be `NULL`
  or a numeric vector of length equal to the number of cases. One or
  more [`offset`](https://rdrr.io/r/stats/offset.html) terms can be
  included in the formula instead or as well, and if more than one is
  specified their sum is used. See
  [`model.offset`](https://rdrr.io/r/stats/model.extract.html).

- control:

  a list of parameters for controlling the fitting process. For
  `glm.fit` this is passed to
  [`glm.control`](https://rdrr.io/r/stats/glm.control.html).

- contrasts:

  an optional list. See the `contrasts.arg` of `model.matrix.default`.

- formula:

  The details of model specification are given under 'Details'.

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

A typical predictor has the form response ~ terms where response is the
(numeric) response vector and terms is a series of terms which specifies
a linear predictor for response. A terms specification of the form
first + second indicates all the terms in first together with all the
terms in second with any duplicates removed.

A specification of the form first:second indicates the the set of terms
obtained by taking the interactions of all terms in first with all terms
in second. The specification first\*second indicates the cross of first
and second. This is the same as first + second + first:second.

The terms in the formula will be re-ordered so that main effects come
first, followed by the interactions, all second-order, all third-order
and so on: to avoid this pass a terms object as the formula.

Non-NULL weights can be used to indicate that different observations
have different dispersions (with the values in weights being inversely
proportional to the dispersions); or equivalently, when the elements of
weights are positive integers w_i, that each response y_i is the mean of
w_i unit-weight observations.

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

[`plsR`](https://fbertran.github.io/plsRglm/reference/plsR.html) and
[`plsRglm`](https://fbertran.github.io/plsRglm/reference/plsRglm.html)

## Author

Frédéric Bertrand  
<frederic.bertrand@lecnam.net>  
<https://fbertran.github.io/homepage/>

## Examples

``` r

data("GasolineYield",package="betareg")
modpls <- plsRbeta(yield~.,data=GasolineYield,nt=3,modele="pls-beta")
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
#> ____Predicting X without NA neither in X or Y____
#> ****________________________________________________****
#> 
modpls$pp
#>              Comp_ 1     Comp_ 2     Comp_ 3
#> gravity   0.37895923 -0.42864981  0.50983922
#> pressure  0.61533000 -0.41618828 -0.01737302
#> temp10   -0.50627633  0.47379983 -0.47750566
#> temp      0.30248369  0.60751756  0.28239621
#> batch1    0.50274128 -0.30221156 -0.25801764
#> batch2   -0.14241033 -0.13859422  0.80068659
#> batch3   -0.04388172 -0.17303214  0.48564161
#> batch4    0.11299471 -0.08302689  0.04755182
#> batch5    0.23341035  0.08396326 -0.51238456
#> batch6    0.07974302  0.07209943 -0.30710455
#> batch7   -0.37365392 -0.02133356  0.81852001
#> batch8   -0.12891598  0.16967195 -0.06904725
#> batch9   -0.02230288  0.19425476 -0.57189134
#> batch10  -0.25409429  0.28587553 -0.61277072
modpls$Coeffs
#>                    [,1]
#> Intercept -4.1210566075
#> gravity    0.0157208676
#> pressure   0.0305159627
#> temp10    -0.0074167766
#> temp       0.0108057945
#> batch1     0.0910284844
#> batch2     0.1398537354
#> batch3     0.2287070464
#> batch4    -0.0008124326
#> batch5     0.1018679027
#> batch6     0.1147971957
#> batch7    -0.1005469609
#> batch8    -0.0447907428
#> batch9    -0.0706292318
#> batch10   -0.1984703429
modpls$Std.Coeffs
#>                    [,1]
#> Intercept -1.5526788976
#> gravity    0.0885938394
#> pressure   0.0799466277
#> temp10    -0.2784359925
#> temp       0.7537685874
#> batch1     0.0305865495
#> batch2     0.0414169259
#> batch3     0.0677303525
#> batch4    -0.0002729861
#> batch5     0.0301676274
#> batch6     0.0339965674
#> batch7    -0.0337848600
#> batch8    -0.0132645358
#> batch9    -0.0173701781
#> batch10   -0.0587759166
modpls$InfCrit
#>                  AIC        BIC Chi2_Pearson_Y      RSS_Y pseudo_R2_Y      R2_Y
#> Nb_Comp_0  -52.77074  -49.83927       30.72004 0.35640772          NA        NA
#> Nb_Comp_1  -87.96104  -83.56383       31.31448 0.11172576   0.6879757 0.6865226
#> Nb_Comp_2 -114.10269 -108.23975       33.06807 0.04650238   0.8671800 0.8695248
#> Nb_Comp_3 -152.71170 -145.38302       30.69727 0.01138837   0.9526757 0.9680468
modpls$PredictY[1,]
#>    gravity   pressure     temp10       temp     batch1     batch2     batch3 
#>  2.0495333  1.6866554 -1.3718198 -1.8219769  2.6040833 -0.3165683 -0.3165683 
#>     batch4     batch5     batch6     batch7     batch8     batch9    batch10 
#> -0.3720119 -0.3165683 -0.3165683 -0.3720119 -0.3165683 -0.2541325 -0.3165683 
rm("modpls")

data("GasolineYield",package="betareg")
yGasolineYield <- GasolineYield$yield
XGasolineYield <- GasolineYield[,2:5]
modpls <- plsRbeta(yGasolineYield,XGasolineYield,nt=3,modele="pls-beta")
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
