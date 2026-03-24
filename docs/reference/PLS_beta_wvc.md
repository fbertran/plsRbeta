# Light version of PLS_beta for cross validation purposes

Light version of `PLS_beta` for cross validation purposes either on
complete or incomplete datasets.

## Usage

``` r
PLS_beta_wvc(
  dataY,
  dataX,
  nt = 2,
  dataPredictY = dataX,
  modele = "pls",
  family = NULL,
  scaleX = TRUE,
  scaleY = NULL,
  keepcoeffs = FALSE,
  keepstd.coeffs = FALSE,
  tol_Xi = 10^(-12),
  weights,
  method = "logistic",
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

- scaleX:

  scale the predictor(s) : must be set to TRUE for `modele="pls"` and
  should be for glms pls.

- scaleY:

  scale the response : Yes/No. Ignored since non always possible for glm
  responses.

- keepcoeffs:

  whether the coefficients of the linear fit on link scale of
  unstandardized eXplanatory variables should be returned or not.

- keepstd.coeffs:

  whether the coefficients of the linear fit on link scale of
  standardized eXplanatory variables should be returned or not.

- tol_Xi:

  minimal value for Norm2(Xi) and \\\mathrm{det}(pp' \times pp)\\ if
  there is any missing value in the `dataX`. It defaults to \\10^{-12}\\

- weights:

  an optional vector of 'prior weights' to be used in the fitting
  process. Should be `NULL` or a numeric vector.

- method:

  logistic, probit, complementary log-log or cauchit (corresponding to a
  Cauchy latent variable).

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

- valsPredict:

  `nrow(dataPredictY) * nt` matrix of the predicted values

- list("coeffs"):

  If the coefficients of the eXplanatory variables were requested:  
  i.e. `keepcoeffs=TRUE`.  
  `ncol(dataX) * 1` matrix of the coefficients of the the eXplanatory
  variables

## Details

This function is called by
[`PLS_beta_kfoldcv_formula`](https://fbertran.github.io/plsRbeta/reference/PLS_beta_kfoldcv_formula.md)
in order to perform cross validation either on complete or incomplete
datasets.

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

Non-NULL weights can be used to indicate that different observations
have different dispersions (with the values in weights being inversely
proportional to the dispersions); or equivalently, when the elements of
weights are positive integers w_i, that each response y_i is the mean of
w_i unit-weight observations.

## References

Frédéric Bertrand, Nicolas Meyer, Michèle Beau-Faller, Karim El Bayed,
Izzie-Jacques Namer, Myriam Maumy-Bertrand (2013). Régression Bêta PLS.
*Journal de la Société Française de Statistique*, **154**(3):143-159.
<https://ojs-test.apps.ocp.math.cnrs.fr/index.php/J-SFdS/article/view/215>

## See also

[`PLS_beta`](https://fbertran.github.io/plsRbeta/reference/PLS_beta.md)
for more detailed results,
[`PLS_beta_kfoldcv`](https://fbertran.github.io/plsRbeta/reference/PLS_beta_kfoldcv.md)
for cross validating models and
[`PLS_lm_wvc`](https://fbertran.github.io/plsRglm/reference/PLS_lm_wvc.html)
for the same function dedicated to plsR models

## Author

Frédéric Bertrand  
<frederic.bertrand@lecnam.net>  
<https://fbertran.github.io/homepage/>

## Examples

``` r

data("GasolineYield",package="betareg")
yGasolineYield <- GasolineYield$yield
XGasolineYield <- GasolineYield[,2:5]
modpls <- PLS_beta_wvc(yGasolineYield,XGasolineYield,nt=3,modele="pls-beta")
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
#> ____Predicting X without NA neither in X nor in Y____
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ****________________________________________________****
#> 
modpls
#> $valsPredict
#>          [,1]       [,2]       [,3]
#> 1  0.18187561 0.11903943 0.10952240
#> 2  0.28144600 0.21696140 0.20261811
#> 3  0.40832135 0.36230655 0.34425424
#> 4  0.53264139 0.51765777 0.49957441
#> 5  0.09146549 0.06943426 0.07502378
#> 6  0.13577318 0.11596641 0.12544647
#> 7  0.22233799 0.21890398 0.23599703
#> 8  0.10602196 0.07790980 0.07966260
#> 9  0.16157484 0.13521897 0.13883717
#> 10 0.25042473 0.23903046 0.24599354
#> 11 0.11891955 0.09303271 0.09654273
#> 12 0.18591623 0.16653971 0.17329467
#> 13 0.27871743 0.28018396 0.29138506
#> 14 0.35738689 0.38179969 0.39599093
#> 15 0.17161846 0.15752793 0.15803548
#> 16 0.25185901 0.25707499 0.25904081
#> 17 0.29688709 0.31561682 0.31849155
#> 18 0.11068533 0.09850211 0.10778989
#> 19 0.20894958 0.22098448 0.24058365
#> 20 0.27062891 0.30383939 0.32869099
#> 21 0.06553811 0.05499223 0.05039617
#> 22 0.08837118 0.08064164 0.07436527
#> 23 0.15947909 0.17048346 0.15961238
#> 24 0.23275612 0.27146324 0.25733075
#> 25 0.08053348 0.07770608 0.07413182
#> 26 0.14334064 0.16067020 0.15503960
#> 27 0.24074422 0.30093955 0.29386264
#> 28 0.12150336 0.13720602 0.13615557
#> 29 0.19979036 0.25165193 0.25145287
#> 30 0.09173971 0.11343104 0.10618428
#> 31 0.10162505 0.12869524 0.12076561
#> 32 0.14395482 0.19625333 0.18584054
#> 
#> $valsPredictPhis
#>        [,1]     [,2]     [,3]
#> 1  89.08427 199.6073 231.1113
#> 2  89.08427 199.6073 231.1113
#> 3  89.08427 199.6073 231.1113
#> 4  89.08427 199.6073 231.1113
#> 5  89.08427 199.6073 231.1113
#> 6  89.08427 199.6073 231.1113
#> 7  89.08427 199.6073 231.1113
#> 8  89.08427 199.6073 231.1113
#> 9  89.08427 199.6073 231.1113
#> 10 89.08427 199.6073 231.1113
#> 11 89.08427 199.6073 231.1113
#> 12 89.08427 199.6073 231.1113
#> 13 89.08427 199.6073 231.1113
#> 14 89.08427 199.6073 231.1113
#> 15 89.08427 199.6073 231.1113
#> 16 89.08427 199.6073 231.1113
#> 17 89.08427 199.6073 231.1113
#> 18 89.08427 199.6073 231.1113
#> 19 89.08427 199.6073 231.1113
#> 20 89.08427 199.6073 231.1113
#> 21 89.08427 199.6073 231.1113
#> 22 89.08427 199.6073 231.1113
#> 23 89.08427 199.6073 231.1113
#> 24 89.08427 199.6073 231.1113
#> 25 89.08427 199.6073 231.1113
#> 26 89.08427 199.6073 231.1113
#> 27 89.08427 199.6073 231.1113
#> 28 89.08427 199.6073 231.1113
#> 29 89.08427 199.6073 231.1113
#> 30 89.08427 199.6073 231.1113
#> 31 89.08427 199.6073 231.1113
#> 32 89.08427 199.6073 231.1113
#> 
rm("modpls")

```
