# Coefficients for permutation bootstrap computations of PLSBeta models

A function passed to `boot` to perform bootstrap.

## Usage

``` r
permcoefs.plsRbetanp(
  dataRepYtt,
  ind,
  nt,
  modele,
  family = NULL,
  maxcoefvalues,
  wwetoile,
  ifbootfail
)
```

## Arguments

- dataRepYtt:

  components' coordinates to bootstrap

- ind:

  indices for resampling

- nt:

  number of components to use

- modele:

  type of modele to use, see
  [plsRbeta](https://fbertran.github.io/plsRbeta/reference/plsRbeta.md)

- family:

  glm family to use, see
  [plsRbeta](https://fbertran.github.io/plsRbeta/reference/plsRbeta.md)

- maxcoefvalues:

  maximum values allowed for the estimates of the coefficients to
  discard those coming from singular bootstrap samples

- wwetoile:

  values of the Wstar matrix in the original fit

- ifbootfail:

  value to return if the estimation fails on a bootstrap sample

## Value

estimates on a bootstrap sample or `ifbootfail` value if the bootstrap
computation fails.

## Note

\~~some notes\~~

## See also

See also
[`bootplsbeta`](https://fbertran.github.io/plsRbeta/reference/bootplsbeta.md)

## Author

Frédéric Bertrand  
<frederic.bertrand@lecnam.net>  
<https://fbertran.github.io/homepage/>

## Examples

``` r
# \donttest{
data("GasolineYield",package="betareg")
modplsbeta <- plsRbeta(yield~.,data=GasolineYield,nt=3, modele="pls-beta")
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
bootplsbeta(modplsbeta, R=250, statistic=permcoefs.plsRbetanp, typeboot="fmodel_np")
#> 
#> ORDINARY NONPARAMETRIC BOOTSTRAP
#> 
#> 
#> Call:
#> boot(data = dataRepYtt, statistic = statistic, R = 250, sim = sim, 
#>     stype = stype, nt = nt, modele = modele, family = family, 
#>     maxcoefvalues = maxcoefvalues[-(1:(length(object$Coeffs) - 
#>         ncol(object$dataX)))], wwetoile = wwetoile, ifbootfail = ifbootfail)
#> 
#> 
#> Bootstrap Statistics :
#>           original        bias     std. error
#> t1*   0.0139841003 -1.431575e-02 0.0038229899
#> t2*   0.0129719714 -1.382620e-02 0.0077454137
#> t3*  -0.0437412382  4.400794e-02 0.0134061619
#> t4*   0.1099552630 -1.100199e-01 0.0209100182
#> t5*   0.0043782576 -5.044095e-03 0.0080290774
#> t6*   0.0068269103 -6.554950e-03 0.0061983660
#> t7*   0.0105674939 -1.039483e-02 0.0051347876
#> t8*  -0.0000215689 -6.398954e-05 0.0009976748
#> t9*   0.0034693882 -3.817977e-03 0.0070556227
#> t10*  0.0044693385 -4.589373e-03 0.0030738554
#> t11* -0.0037252887  4.333924e-03 0.0108662479
#> t12* -0.0020821555  2.235699e-03 0.0014725934
#> t13* -0.0031646188  3.099113e-03 0.0032364535
#> t14* -0.0089486739  9.144411e-03 0.0019703566
# }
```
