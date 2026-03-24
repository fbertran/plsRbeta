# Partial least squares regression beta models with kfold cross validation

This function implements kfold cross validation on complete or
incomplete datasets for partial least squares beta regression models
(formula specification of the model).

## Usage

``` r
PLS_beta_kfoldcv_formula(
  formula,
  data = NULL,
  nt = 2,
  limQ2set = 0.0975,
  modele = "pls",
  family = NULL,
  K = nrow(dataX),
  NK = 1,
  grouplist = NULL,
  random = FALSE,
  scaleX = TRUE,
  scaleY = NULL,
  keepcoeffs = FALSE,
  keepfolds = FALSE,
  keepdataY = TRUE,
  keepMclassed = FALSE,
  tol_Xi = 10^(-12),
  weights,
  subset,
  start = NULL,
  etastart,
  mustart,
  offset,
  method,
  control = list(),
  contrasts = NULL,
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

- formula:

  an object of class "[`formula`](https://rdrr.io/r/stats/formula.html)"
  (or one that can be coerced to that class): a symbolic description of
  the model to be fitted. The details of model specification are given
  under 'Details'.

- data:

  an optional data frame, list or environment (or object coercible by
  [`as.data.frame`](https://rdrr.io/r/base/as.data.frame.html) to a data
  frame) containing the variables in the model. If not found in `data`,
  the variables are taken from `environment(formula)`, typically the
  environment from which `plsRglm` is called.

- nt:

  number of components to be extracted

- limQ2set:

  limit value for the Q2

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

- K:

  number of groups

- NK:

  number of times the group division is made

- grouplist:

  to specify the members of the `K` groups

- random:

  should the `K` groups be made randomly

- scaleX:

  scale the predictor(s) : must be set to TRUE for `modele="pls"` and
  should be for glms pls.

- scaleY:

  scale the response : Yes/No. Ignored since non always possible for glm
  responses.

- keepcoeffs:

  shall the coefficients for each model be returned

- keepfolds:

  shall the groups' composition be returned

- keepdataY:

  shall the observed value of the response for each one of the predicted
  value be returned

- keepMclassed:

  shall the number of miss classed be returned (unavailable)

- tol_Xi:

  minimal value for Norm2(Xi) and \\\mathrm{det}(pp' \times pp)\\ if
  there is any missing value in the `dataX`. It defaults to \\10^{-12}\\

- weights:

  an optional vector of 'prior weights' to be used in the fitting
  process. Should be `NULL` or a numeric vector.

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

- method:

  for fitting glms with glm (

  :   the method to be used in fitting the model. The default method
      `"glm.fit"` uses iteratively reweighted least squares (IWLS).
      User-supplied fitting functions can be supplied either as a
      function or a character string naming a function, with a function
      which takes the same arguments as `glm.fit`. If "model.frame", the
      model frame is returned.

  list("\\pls-glm-Gamma\\")

  :   the method to be used in fitting the model. The default method
      `"glm.fit"` uses iteratively reweighted least squares (IWLS).
      User-supplied fitting functions can be supplied either as a
      function or a character string naming a function, with a function
      which takes the same arguments as `glm.fit`. If "model.frame", the
      model frame is returned.

  ,

  :   the method to be used in fitting the model. The default method
      `"glm.fit"` uses iteratively reweighted least squares (IWLS).
      User-supplied fitting functions can be supplied either as a
      function or a character string naming a function, with a function
      which takes the same arguments as `glm.fit`. If "model.frame", the
      model frame is returned.

  list("\\pls-glm-gaussian\\")

  :   the method to be used in fitting the model. The default method
      `"glm.fit"` uses iteratively reweighted least squares (IWLS).
      User-supplied fitting functions can be supplied either as a
      function or a character string naming a function, with a function
      which takes the same arguments as `glm.fit`. If "model.frame", the
      model frame is returned.

  ,

  :   the method to be used in fitting the model. The default method
      `"glm.fit"` uses iteratively reweighted least squares (IWLS).
      User-supplied fitting functions can be supplied either as a
      function or a character string naming a function, with a function
      which takes the same arguments as `glm.fit`. If "model.frame", the
      model frame is returned.

  list("\\pls-glm-inverse.gaussian\\")

  :   the method to be used in fitting the model. The default method
      `"glm.fit"` uses iteratively reweighted least squares (IWLS).
      User-supplied fitting functions can be supplied either as a
      function or a character string naming a function, with a function
      which takes the same arguments as `glm.fit`. If "model.frame", the
      model frame is returned.

  ,

  :   the method to be used in fitting the model. The default method
      `"glm.fit"` uses iteratively reweighted least squares (IWLS).
      User-supplied fitting functions can be supplied either as a
      function or a character string naming a function, with a function
      which takes the same arguments as `glm.fit`. If "model.frame", the
      model frame is returned.

  list("\\pls-glm-logistic\\")

  :   the method to be used in fitting the model. The default method
      `"glm.fit"` uses iteratively reweighted least squares (IWLS).
      User-supplied fitting functions can be supplied either as a
      function or a character string naming a function, with a function
      which takes the same arguments as `glm.fit`. If "model.frame", the
      model frame is returned.

  ,

  :   the method to be used in fitting the model. The default method
      `"glm.fit"` uses iteratively reweighted least squares (IWLS).
      User-supplied fitting functions can be supplied either as a
      function or a character string naming a function, with a function
      which takes the same arguments as `glm.fit`. If "model.frame", the
      model frame is returned.

  list("\\pls-glm-poisson\\")

  :   the method to be used in fitting the model. The default method
      `"glm.fit"` uses iteratively reweighted least squares (IWLS).
      User-supplied fitting functions can be supplied either as a
      function or a character string naming a function, with a function
      which takes the same arguments as `glm.fit`. If "model.frame", the
      model frame is returned.

  ,

  :   the method to be used in fitting the model. The default method
      `"glm.fit"` uses iteratively reweighted least squares (IWLS).
      User-supplied fitting functions can be supplied either as a
      function or a character string naming a function, with a function
      which takes the same arguments as `glm.fit`. If "model.frame", the
      model frame is returned.

  list("\\modele=pls-glm-family\\")

  :   the method to be used in fitting the model. The default method
      `"glm.fit"` uses iteratively reweighted least squares (IWLS).
      User-supplied fitting functions can be supplied either as a
      function or a character string naming a function, with a function
      which takes the same arguments as `glm.fit`. If "model.frame", the
      model frame is returned.

  )

  :   the method to be used in fitting the model. The default method
      `"glm.fit"` uses iteratively reweighted least squares (IWLS).
      User-supplied fitting functions can be supplied either as a
      function or a character string naming a function, with a function
      which takes the same arguments as `glm.fit`. If "model.frame", the
      model frame is returned.

  list("pls-glm-polr")

  :   logistic, probit, complementary log-log or cauchit (corresponding
      to a Cauchy latent variable).

- control:

  a list of parameters for controlling the fitting process. For
  `glm.fit` this is passed to
  [`glm.control`](https://rdrr.io/r/stats/glm.control.html).

- contrasts:

  an optional list. See the `contrasts.arg` of `model.matrix.default`.

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

## Value

- results_kfolds:

  list of `NK`. Each element of the list sums up the results for a group
  division:

  list

  :   of `K` matrices of size about `nrow(dataX)/K * nt` with the
      predicted values for a growing number of components

  list()

  :   ...

  list

  :   of `K` matrices of size about `nrow(dataX)/K * nt` with the
      predicted values for a growing number of components

- folds:

  list of `NK`. Each element of the list sums up the informations for a
  group division:

  list

  :   of `K` vectors of length about `nrow(dataX)` with the numbers of
      the rows of `dataX` that were used as a training set

  list()

  :   ...

  list

  :   of `K` vectors of length about `nrow(dataX)` with the numbers of
      the rows of `dataX` that were used as a training set

- dataY_kfolds:

  list of `NK`. Each element of the list sums up the results for a group
  division:

  list

  :   of `K` matrices of size about `nrow(dataX)/K * 1` with the
      observed values of the response

  list()

  :   ...

  list

  :   of `K` matrices of size about `nrow(dataX)/K * 1` with the
      observed values of the response

- call:

  the call of the function

## Details

Predicts 1 group with the `K-1` other groups. Leave one out cross
validation is thus obtained for `K==nrow(dataX)`.

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

## Note

Work for complete and incomplete datasets.

## References

Frédéric Bertrand, Nicolas Meyer, Michèle Beau-Faller, Karim El Bayed,
Izzie-Jacques Namer, Myriam Maumy-Bertrand (2013). Régression Bêta PLS.
*Journal de la Société Française de Statistique*, **154**(3):143-159.
<https://ojs-test.apps.ocp.math.cnrs.fr/index.php/J-SFdS/article/view/215>

## See also

[`kfolds2coeff`](https://fbertran.github.io/plsRglm/reference/kfolds2coeff.html),
[`kfolds2Pressind`](https://fbertran.github.io/plsRglm/reference/kfolds2Pressind.html),
[`kfolds2Press`](https://fbertran.github.io/plsRglm/reference/kfolds2Press.html),
[`kfolds2Mclassedind`](https://fbertran.github.io/plsRglm/reference/kfolds2Mclassedind.html),
[`kfolds2Mclassed`](https://fbertran.github.io/plsRglm/reference/kfolds2Mclassed.html)
and
[`kfolds2CVinfos_beta`](https://fbertran.github.io/plsRbeta/reference/kfolds2CVinfos_beta.md)
to extract and transform results from kfold cross validation.

## Author

Frédéric Bertrand  
<frederic.bertrand@lecnam.net>  
<https://fbertran.github.io/homepage/>

## Examples

``` r
if (FALSE) { # \dontrun{
data("GasolineYield",package="betareg")
bbb <- PLS_beta_kfoldcv_formula(yield~.,data=GasolineYield,nt=3,modele="pls-beta")
kfolds2CVinfos_beta(bbb)
} # }
```
