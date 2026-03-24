# Extracts and computes information criteria and fits statistics for kfold cross validated partial least squares beta regression models

This function extracts and computes information criteria and fits
statistics for kfold cross validated partial least squares beta
regression models for both formula or classic specifications of the
model.

## Usage

``` r
kfolds2CVinfos_beta(pls_kfolds, MClassed = FALSE)
```

## Arguments

- pls_kfolds:

  an object computed using
  [`PLS_beta_kfoldcv`](https://fbertran.github.io/plsRbeta/reference/PLS_beta_kfoldcv.md)

- MClassed:

  should number of miss classed be computed

## Value

- list:

  table of fit statistics for first group partition

- list():

  ...

- list:

  table of fit statistics for last group partition

## Details

The Mclassed option should only set to `TRUE` if the response is binary.

## References

Frédéric Bertrand, Nicolas Meyer, Michèle Beau-Faller, Karim El Bayed,
Izzie-Jacques Namer, Myriam Maumy-Bertrand (2013). Régression Bêta PLS.
*Journal de la Société Française de Statistique*, **154**(3):143-159.
<https://ojs-test.apps.ocp.math.cnrs.fr/index.php/J-SFdS/article/view/215>

## See also

[`kfolds2coeff`](https://fbertran.github.io/plsRglm/reference/kfolds2coeff.html),
[`kfolds2Pressind`](https://fbertran.github.io/plsRglm/reference/kfolds2Pressind.html),
[`kfolds2Press`](https://fbertran.github.io/plsRglm/reference/kfolds2Press.html),
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
bbb <- PLS_beta_kfoldcv_formula(yield~.,data=GasolineYield,nt=3,modele="pls-beta")
kfolds2CVinfos_beta(bbb)
} # }
```
