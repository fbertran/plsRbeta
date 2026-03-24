# Computes Predicted Chisquare for kfold cross validated partial least squares beta regression models.

This function computes Predicted Chisquare for kfold cross validated
partial least squares beta regression models.

## Usage

``` r
kfolds2Chisq(pls_kfolds)
```

## Arguments

- pls_kfolds:

  a kfold cross validated partial least squares regression glm model

## Value

- list:

  Total Predicted Chisquare vs number of components for the first group
  partition

- list():

  ...

- list:

  Total Predicted Chisquare vs number of components for the last group
  partition

## Note

Use
[`PLS_beta_kfoldcv`](https://fbertran.github.io/plsRbeta/reference/PLS_beta_kfoldcv.md)
to create kfold cross validated partial least squares regression glm and
beta models.

## References

Frédéric Bertrand, Nicolas Meyer, Michèle Beau-Faller, Karim El Bayed,
Izzie-Jacques Namer, Myriam Maumy-Bertrand (2013). Régression Bêta PLS.
*Journal de la Société Française de Statistique*, **154**(3):143-159.
<https://ojs-test.apps.ocp.math.cnrs.fr/index.php/J-SFdS/article/view/215>

## See also

[`kfolds2coeff`](https://fbertran.github.io/plsRglm/reference/kfolds2coeff.html),
[`kfolds2Press`](https://fbertran.github.io/plsRglm/reference/kfolds2Press.html),
[`kfolds2Pressind`](https://fbertran.github.io/plsRglm/reference/kfolds2Pressind.html),
[`kfolds2Chisqind`](https://fbertran.github.io/plsRbeta/reference/kfolds2Chisqind.md),
[`kfolds2Mclassedind`](https://fbertran.github.io/plsRglm/reference/kfolds2Mclassedind.html)
and
[`kfolds2Mclassed`](https://fbertran.github.io/plsRglm/reference/kfolds2Mclassed.html)
to extract and transforms results from kfold cross validation.

## Author

Frédéric Bertrand  
<frederic.bertrand@lecnam.net>  
<https://fbertran.github.io/homepage/>

## Examples

``` r
if (FALSE) { # \dontrun{
data("GasolineYield",package="betareg")
yGasolineYield <- GasolineYield$yield
XGasolineYield <- GasolineYield[,2:5]
bbb <- PLS_beta_kfoldcv(yGasolineYield,XGasolineYield,nt=3,modele="pls-beta")
kfolds2Chisq(bbb)
} # }
```
