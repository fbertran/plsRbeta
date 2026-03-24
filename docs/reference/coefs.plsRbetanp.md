# Coefficients for bootstrap computations of PLSBeta models

A function passed to `boot` to perform bootstrap.

## Usage

``` r
coefs.plsRbetanp(
  dataRepYtt,
  ind,
  nt,
  modele,
  family = NULL,
  method = "logistic",
  link = NULL,
  link.phi = NULL,
  type = "ML",
  verbose = TRUE,
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

- method:

  method for beta regression

- link:

  link for beta regression

- link.phi:

  link.phi for beta regression

- type:

  type of estimates

- verbose:

  should info messages be displayed ?

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
bootplsbeta(plsRbeta(yield~.,data=GasolineYield,nt=3, modele="pls-beta"), typeboot="fmodel_np", 
R=250, statistic=coefs.plsRbetanp)
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
#>           original        bias    std. error
#> t1*   0.0139841003 -6.917754e-04 0.008434712
#> t2*   0.0129719714 -1.791867e-03 0.015879635
#> t3*  -0.0437412382  5.654104e-04 0.030145370
#> t4*   0.1099552630  5.071299e-04 0.050028234
#> t5*   0.0043782576 -1.333676e-03 0.015709145
#> t6*   0.0068269103  5.293340e-04 0.012730336
#> t7*   0.0105674939  3.558548e-04 0.011115232
#> t8*  -0.0000215689 -1.773027e-04 0.002033377
#> t9*   0.0034693882 -6.139327e-04 0.013252135
#> t10*  0.0044693385 -1.775854e-04 0.005557751
#> t11* -0.0037252887  1.122405e-03 0.020672308
#> t12* -0.0020821555  3.158314e-04 0.002934637
#> t13* -0.0031646188 -9.781348e-05 0.006505060
#> t14* -0.0089486739  3.823206e-04 0.004190039
# }
```
