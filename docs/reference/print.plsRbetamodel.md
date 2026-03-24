# Print method for plsRbeta models

This function provides a print method for the class `"plsRbetamodel"`

## Usage

``` r
# S3 method for class 'plsRbetamodel'
print(x, ...)
```

## Arguments

- x:

  an object of the class `"plsRbetamodel"`

- ...:

  not used

## Value

`NULL`

## References

Frédéric Bertrand, Nicolas Meyer, Michèle Beau-Faller, Karim El Bayed,
Izzie-Jacques Namer, Myriam Maumy-Bertrand (2013). Régression Bêta PLS.
*Journal de la Société Française de Statistique*, **154**(3):143-159.
<https://ojs-test.apps.ocp.math.cnrs.fr/index.php/J-SFdS/article/view/215>

## See also

[`print`](https://rdrr.io/r/base/print.html)

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
print(modpls)
#> Number of required components:
#> [1] 3
#> Number of successfully computed components:
#> [1] 3
#> Coefficients:
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
#> Information criteria and Fit statistics:
#>                  AIC        BIC Chi2_Pearson_Y      RSS_Y pseudo_R2_Y      R2_Y
#> Nb_Comp_0  -52.77074  -49.83927       30.72004 0.35640772          NA        NA
#> Nb_Comp_1  -87.96104  -83.56383       31.31448 0.11172576   0.6879757 0.6865226
#> Nb_Comp_2 -114.10269 -108.23975       33.06807 0.04650238   0.8671800 0.8695248
#> Nb_Comp_3 -152.71170 -145.38302       30.69727 0.01138837   0.9526757 0.9680468
```
