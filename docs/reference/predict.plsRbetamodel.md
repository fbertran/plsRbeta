# Predict method for plsRbeta models

This function predicts responses or component scores from a fitted
`"plsRbetamodel"` object.

## Usage

``` r
predict.plsRbetamodel(
  object,
  newdata,
  comps = object$computed_nt,
  type = c("response", "scores"),
  methodNA = "adaptative",
  verbose = TRUE,
  ...
)
```

## Arguments

- object:

  an object of class `"plsRbetamodel"`

- newdata:

  optional new predictor data. When omitted, fitted values or training
  scores are returned.

- comps:

  number of components to use.

- type:

  type of prediction. `"response"` returns predicted response values and
  `"scores"` returns component scores.

- methodNA:

  method to use when `newdata` contains missing values. `"adaptative"`
  uses the standard score computation for complete rows and the
  missing-data score reconstruction otherwise. `"missingdata"` applies
  the missing-data reconstruction to every row.

- verbose:

  should info messages be displayed ?

- ...:

  additional arguments passed to the underlying `predict` method when a
  final model is available.

## Value

A numeric vector/matrix of predictions or a matrix of component scores.

## Examples

``` r
data("GasolineYield", package = "betareg")
modpls <- plsRbeta(yield ~ ., data = GasolineYield, nt = 2, modele = "pls-beta")
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
#> ____Predicting X without NA neither in X or Y____
#> ****________________________________________________****
#> 
head(predict(modpls))
#>          1          2          3          4          5          6 
#> 0.13822439 0.22621610 0.34762680 0.47557444 0.07015158 0.10786977 
head(predict(modpls, newdata = GasolineYield[1:3, -1, drop = FALSE]))
#>         1         2         3 
#> 0.1382244 0.2262161 0.3476268 
head(predict(modpls, type = "scores"))
#>      Comp_ 1    Comp_ 2
#> 1  0.8338521 -3.3585512
#> 2  1.5693096 -2.2635486
#> 3  2.3047672 -1.1685460
#> 4  2.9561724 -0.1986865
#> 5 -1.2822054 -1.9011017
#> 6 -0.7043458 -1.0407424
```
