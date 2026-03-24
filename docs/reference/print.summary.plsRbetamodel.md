# Print method for summaries of plsRbeta models

This function provides a print method for the class
`"summary.plsRbetamodel"`

## Usage

``` r
# S3 method for class 'summary.plsRbetamodel'
print(x, ...)
```

## Arguments

- x:

  an object of the class `"summary.plsRbetamodel"`

- ...:

  not used

## Value

- language:

  call of the model

## References

Frédéric Bertrand, Nicolas Meyer, Michèle Beau-Faller, Karim El Bayed,
Izzie-Jacques Namer, Myriam Maumy-Bertrand (2013). Régression Bêta PLS.
*Journal de la Société Française de Statistique*, **154**(3):143-159.
<https://ojs-test.apps.ocp.math.cnrs.fr/index.php/J-SFdS/article/view/215>

## See also

[`print`](https://rdrr.io/r/base/print.html) and
[`summary`](https://rdrr.io/r/base/summary.html)

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
print(summary(modpls))
#> Call:
#> plsRbetamodel.formula(object = yield ~ ., data = GasolineYield, 
#>     nt = 3, modele = "pls-beta")
```
