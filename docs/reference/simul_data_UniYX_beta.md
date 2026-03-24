# Data generating function for univariate beta plsR models

This function generates a single univariate rate response value \\Y\\
and a vector of explanatory variables \\(X_1,\ldots,X\_{totdim})\\ drawn
from a model with a given number of latent components.

## Usage

``` r
simul_data_UniYX_beta(
  totdim,
  ncomp,
  disp = 1,
  link = "logit",
  type = "a",
  phi0 = 20
)
```

## Arguments

- totdim:

  Number of columns of the X vector (from `ncomp` to hardware limits)

- ncomp:

  Number of latent components in the model (from 2 to 6)

- disp:

  Tune the shape of the beta distribution (defaults to 1)

- link:

  Character specification of the link function in the mean model (mu).
  Currently, "`logit`", "`probit`", "`cloglog`", "`cauchit`", "`log`",
  "`loglog`" are supported. Alternatively, an object of class "link-glm"
  can be supplied.

- type:

  Simulation scheme

- phi0:

  Simulation scheme "a" parameter

## Value

- vector:

  \\(Y,X_1,\ldots,X\_{totdim})\\

## Details

This function should be combined with the replicate function to give
rise to a larger dataset. The algorithm used is a modification of a port
of the one described in the article of Li which is a multivariate
generalization of the algorithm of Naes and Martens.

## References

Frédéric Bertrand, Nicolas Meyer, Michèle Beau-Faller, Karim El Bayed,
Izzie-Jacques Namer, Myriam Maumy-Bertrand (2013). Régression Bêta PLS.
*Journal de la Société Française de Statistique*, **154**(3):143-159.
<https://ojs-test.apps.ocp.math.cnrs.fr/index.php/J-SFdS/article/view/215>

T. Naes, H. Martens (1985). Comparison of prediction methods for
multicollinear data. *Commun. Stat., Simul.*, **14**:545-576.
\<doi:10.1080/03610918508812458\>

Baibing Li, Julian Morris, Elaine B. Martin (2002). Model selection for
partial least squares regression, *Chemometrics and Intelligent
Laboratory Systems*, **64**:79-89.
\<doi:110.1016/S0169-7439(02)00051-5\>

## See also

[`simul_data_UniYX`](https://fbertran.github.io/plsRglm/reference/simul_data_UniYX.html)

## Author

Frédéric Bertrand  
<frederic.bertrand@lecnam.net>  
<https://fbertran.github.io/homepage/>

## Examples

``` r
# \donttest{
# logit link
layout(matrix(1:4,nrow=2))
hist(t(replicate(100,simul_data_UniYX_beta(4,4)))[,1])
hist(t(replicate(100,simul_data_UniYX_beta(4,4,disp=3)))[,1])
hist(t(replicate(100,simul_data_UniYX_beta(4,4,disp=5)))[,1])
hist(t(replicate(100,simul_data_UniYX_beta(4,4,disp=15)))[,1])

layout(1)

# probit link
layout(matrix(1:4,nrow=2))
hist(t(replicate(100,simul_data_UniYX_beta(4,4,link="probit")))[,1])
hist(t(replicate(100,simul_data_UniYX_beta(4,4,disp=3,link="probit")))[,1])
hist(t(replicate(100,simul_data_UniYX_beta(4,4,disp=5,link="probit")))[,1])
hist(t(replicate(100,simul_data_UniYX_beta(4,4,disp=15,link="probit")))[,1])

layout(1)

# cloglog link
layout(matrix(1:4,nrow=2))
hist(t(replicate(100,simul_data_UniYX_beta(4,4,link="cloglog")))[,1])
hist(t(replicate(100,simul_data_UniYX_beta(4,4,disp=3,link="cloglog")))[,1])
hist(t(replicate(100,simul_data_UniYX_beta(4,4,disp=5,link="cloglog")))[,1])
hist(t(replicate(100,simul_data_UniYX_beta(4,4,disp=15,link="cloglog")))[,1])

layout(1)

# cauchit link
layout(matrix(1:4,nrow=2))
hist(t(replicate(100,simul_data_UniYX_beta(4,4,link="cauchit")))[,1])
hist(t(replicate(100,simul_data_UniYX_beta(4,4,disp=3,link="cauchit")))[,1])
hist(t(replicate(100,simul_data_UniYX_beta(4,4,disp=5,link="cauchit")))[,1])
hist(t(replicate(100,simul_data_UniYX_beta(4,4,disp=15,link="cauchit")))[,1])

layout(1)

# loglog link
layout(matrix(1:4,nrow=2))
hist(t(replicate(100,simul_data_UniYX_beta(4,4,link="loglog")))[,1])
hist(t(replicate(100,simul_data_UniYX_beta(4,4,disp=3,link="loglog")))[,1])
hist(t(replicate(100,simul_data_UniYX_beta(4,4,disp=5,link="loglog")))[,1])
hist(t(replicate(100,simul_data_UniYX_beta(4,4,disp=15,link="loglog")))[,1])

layout(1)

# log link
layout(matrix(1:4,nrow=2))
hist(t(replicate(100,simul_data_UniYX_beta(4,4,link="log")))[,1])
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
hist(t(replicate(100,simul_data_UniYX_beta(4,4,disp=3,link="log")))[,1])
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
hist(t(replicate(100,simul_data_UniYX_beta(4,4,disp=5,link="log")))[,1])
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
hist(t(replicate(100,simul_data_UniYX_beta(4,4,disp=15,link="log")))[,1])
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced
#> Warning: NAs produced

layout(1)
# }
```
