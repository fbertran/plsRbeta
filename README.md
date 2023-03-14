<!-- README.md is generated from README.Rmd. Please edit that file -->



# plsRbeta <img src="man/figures/logo.png" align="right" width="200"/>

# Partial Least Squares Regression for Beta Regression Models
## Frédéric Bertrand and Myriam Maumy-Bertrand

<!-- badges: start -->
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-green.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-CMD-check](https://github.com/fbertran/plsRbeta/workflows/R-CMD-check/badge.svg)](https://github.com/fbertran/plsRbeta/actions)
[![Codecov test coverage](https://codecov.io/gh/fbertran/plsRbeta/branch/master/graph/badge.svg)](https://app.codecov.io/gh/fbertran/plsRbeta?branch=master)
[![CRAN status](https://www.r-pkg.org/badges/version/plsRbeta)](https://cran.r-project.org/package=plsRbeta)
[![CRAN RStudio mirror downloads](https://cranlogs.r-pkg.org/badges/plsRbeta)](https://cran.r-project.org/package=plsRbeta)
[![GitHub Repo stars](https://img.shields.io/github/stars/fbertran/plsRbeta?style=social)](https://github.com/fbertran/plsRbeta)
[![DOI](https://zenodo.org/badge/18454088.svg)](https://zenodo.org/badge/latestdoi/18454088)
<!-- badges: end -->

The goal of plsRbeta is to provide Partial least squares Regression for (weighted) beta regression models (Bertrand 2013,  <http://journal-sfds.fr/article/view/215>) and k-fold cross-validation of such models using various criteria. It allows for missing data in the explanatory variables. Bootstrap confidence intervals constructions are also available.


The package was accepted for presentation at the the useR! 2021 international conference. A technical note for the package was created and published on the website of the conference. It can be accessed  here: [https://user2021.r-project.org/participation/technical_notes/t138/technote/](https://user2021.r-project.org/participation/technical_notes/t138/technote/). It is not only an english translation of most of the contents of the original article that was published in French but it also contains the R code reproduce the two examples that were presented in the article.


This website and these examples were created by F. Bertrand and M. Maumy-Bertrand.

## Installation

You can install the released version of plsRbeta from [CRAN](https://CRAN.R-project.org) with:


```r
install.packages("plsRbeta")
```

You can install the development version of plsRbeta from [github](https://github.com) with:


```r
devtools::install_github("fbertran/plsRbeta")
```

## Example

### Using a model matrix
Fit a plsRbeta model using a model matrix.


```r
data("GasolineYield",package="betareg")
yGasolineYield <- GasolineYield$yield
XGasolineYield <- GasolineYield[,2:5]
library(plsRbeta)
modpls <- plsRbeta(yGasolineYield,XGasolineYield,nt=3,modele="pls-beta")
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
#> ____Predicting X without NA neither in X nor in Y____
#> ****________________________________________________****
print(modpls)
#> Number of required components:
#> [1] 3
#> Number of successfully computed components:
#> [1] 3
#> Coefficients:
#>                   [,1]
#> Intercept -3.324462301
#> gravity    0.001577508
#> pressure   0.072027686
#> temp10    -0.008398771
#> temp       0.010365973
#> Information criteria and Fit statistics:
#>                  AIC        BIC Chi2_Pearson_Y
#> Nb_Comp_0  -52.77074  -49.83927       30.72004
#> Nb_Comp_1 -112.87383 -108.47662       30.57369
#> Nb_Comp_2 -136.43184 -130.56889       30.97370
#> Nb_Comp_3 -139.08440 -131.75572       31.08224
#>                RSS_Y pseudo_R2_Y      R2_Y
#> Nb_Comp_0 0.35640772          NA        NA
#> Nb_Comp_1 0.05211039   0.8498691 0.8537900
#> Nb_Comp_2 0.02290022   0.9256771 0.9357471
#> Nb_Comp_3 0.02022386   0.9385887 0.9432564
```

Additionnal values can be retrieved from the fitted model.

```r
modpls$pp
#>             Comp_ 1    Comp_ 2    Comp_ 3
#> gravity   0.4590380 -0.4538663 -2.5188256
#> pressure  0.6395524 -0.4733525  0.6488823
#> temp10   -0.5435643  0.5292108 -1.3295905
#> temp      0.5682795  0.5473174 -0.2156423
modpls$Coeffs
#>                   [,1]
#> Intercept -3.324462301
#> gravity    0.001577508
#> pressure   0.072027686
#> temp10    -0.008398771
#> temp       0.010365973
modpls$Std.Coeffs
#>                   [,1]
#> Intercept -1.547207760
#> gravity    0.008889933
#> pressure   0.188700277
#> temp10    -0.315301400
#> temp       0.723088387
modpls$InfCrit
#>                  AIC        BIC Chi2_Pearson_Y      RSS_Y
#> Nb_Comp_0  -52.77074  -49.83927       30.72004 0.35640772
#> Nb_Comp_1 -112.87383 -108.47662       30.57369 0.05211039
#> Nb_Comp_2 -136.43184 -130.56889       30.97370 0.02290022
#> Nb_Comp_3 -139.08440 -131.75572       31.08224 0.02022386
#>           pseudo_R2_Y      R2_Y
#> Nb_Comp_0          NA        NA
#> Nb_Comp_1   0.8498691 0.8537900
#> Nb_Comp_2   0.9256771 0.9357471
#> Nb_Comp_3   0.9385887 0.9432564
modpls$PredictY[1,]
#>   gravity  pressure    temp10      temp 
#>  2.049533  1.686655 -1.371820 -1.821977
rm("modpls")
```

###Formula support

Fit a plsRbeta model using formula support.

```r
data("GasolineYield",package="betareg")
modpls <- plsRbeta(yield~.,data=GasolineYield,nt=3,modele="pls-beta", verbose=FALSE)
print(modpls)
#> Number of required components:
#> [1] 3
#> Number of successfully computed components:
#> [1] 3
#> Coefficients:
#>                    [,1]
#> Intercept -4.1210566077
#> gravity    0.0157208676
#> pressure   0.0305159627
#> temp10    -0.0074167766
#> temp       0.0108057945
#> batch1     0.0910284843
#> batch2     0.1398537354
#> batch3     0.2287070465
#> batch4    -0.0008124326
#> batch5     0.1018679027
#> batch6     0.1147971957
#> batch7    -0.1005469609
#> batch8    -0.0447907428
#> batch9    -0.0706292318
#> batch10   -0.1984703429
#> Information criteria and Fit statistics:
#>                  AIC        BIC Chi2_Pearson_Y
#> Nb_Comp_0  -52.77074  -49.83927       30.72004
#> Nb_Comp_1  -87.96104  -83.56383       31.31448
#> Nb_Comp_2 -114.10269 -108.23975       33.06807
#> Nb_Comp_3 -152.71170 -145.38302       30.69727
#>                RSS_Y pseudo_R2_Y      R2_Y
#> Nb_Comp_0 0.35640772          NA        NA
#> Nb_Comp_1 0.11172576   0.6879757 0.6865226
#> Nb_Comp_2 0.04650238   0.8671800 0.8695248
#> Nb_Comp_3 0.01138837   0.9526757 0.9680468
```

Additionnal values can be retrieved from the fitted model.

```r
modpls$pp
#>              Comp_ 1     Comp_ 2     Comp_ 3
#> gravity   0.37895923 -0.42864981  0.50983922
#> pressure  0.61533000 -0.41618828 -0.01737302
#> temp10   -0.50627633  0.47379983 -0.47750566
#> temp      0.30248369  0.60751756  0.28239621
#> batch1    0.50274128 -0.30221156 -0.25801764
#> batch2   -0.14241033 -0.13859422  0.80068659
#> batch3   -0.04388172 -0.17303214  0.48564161
#> batch4    0.11299471 -0.08302689  0.04755182
#> batch5    0.23341035  0.08396326 -0.51238456
#> batch6    0.07974302  0.07209943 -0.30710455
#> batch7   -0.37365392 -0.02133356  0.81852001
#> batch8   -0.12891598  0.16967195 -0.06904725
#> batch9   -0.02230288  0.19425476 -0.57189134
#> batch10  -0.25409429  0.28587553 -0.61277072
modpls$Coeffs
#>                    [,1]
#> Intercept -4.1210566077
#> gravity    0.0157208676
#> pressure   0.0305159627
#> temp10    -0.0074167766
#> temp       0.0108057945
#> batch1     0.0910284843
#> batch2     0.1398537354
#> batch3     0.2287070465
#> batch4    -0.0008124326
#> batch5     0.1018679027
#> batch6     0.1147971957
#> batch7    -0.1005469609
#> batch8    -0.0447907428
#> batch9    -0.0706292318
#> batch10   -0.1984703429
modpls$Std.Coeffs
#>                    [,1]
#> Intercept -1.5526788976
#> gravity    0.0885938394
#> pressure   0.0799466278
#> temp10    -0.2784359925
#> temp       0.7537685874
#> batch1     0.0305865495
#> batch2     0.0414169259
#> batch3     0.0677303525
#> batch4    -0.0002729861
#> batch5     0.0301676274
#> batch6     0.0339965674
#> batch7    -0.0337848600
#> batch8    -0.0132645358
#> batch9    -0.0173701781
#> batch10   -0.0587759166
modpls$InfCrit
#>                  AIC        BIC Chi2_Pearson_Y      RSS_Y
#> Nb_Comp_0  -52.77074  -49.83927       30.72004 0.35640772
#> Nb_Comp_1  -87.96104  -83.56383       31.31448 0.11172576
#> Nb_Comp_2 -114.10269 -108.23975       33.06807 0.04650238
#> Nb_Comp_3 -152.71170 -145.38302       30.69727 0.01138837
#>           pseudo_R2_Y      R2_Y
#> Nb_Comp_0          NA        NA
#> Nb_Comp_1   0.6879757 0.6865226
#> Nb_Comp_2   0.8671800 0.8695248
#> Nb_Comp_3   0.9526757 0.9680468
modpls$PredictY[1,]
#>    gravity   pressure     temp10       temp     batch1 
#>  2.0495333  1.6866554 -1.3718198 -1.8219769  2.6040833 
#>     batch2     batch3     batch4     batch5     batch6 
#> -0.3165683 -0.3165683 -0.3720119 -0.3165683 -0.3165683 
#>     batch7     batch8     batch9    batch10 
#> -0.3720119 -0.3165683 -0.2541325 -0.3165683
```

###Information criteria and cross validation


```r
data("GasolineYield",package="betareg")
set.seed(1)
bbb <- PLS_beta_kfoldcv_formula(yield~.,data=GasolineYield,nt=3,modele="pls-beta",verbose=FALSE)
kfolds2CVinfos_beta(bbb)
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
#> NK: 1
#> [[1]]
#>                  AIC        BIC Q2Chisqcum_Y
#> Nb_Comp_0  -52.77074  -49.83927           NA
#> Nb_Comp_1  -87.96104  -83.56383    -1.121431
#> Nb_Comp_2 -114.10269 -108.23975    -5.291744
#> Nb_Comp_3 -152.71170 -145.38302   -11.583916
#>            limQ2 Q2Chisq_Y PREChi2_Pearson_Y
#> Nb_Comp_0     NA        NA                NA
#> Nb_Comp_1 0.0975 -1.121431          65.17044
#> Nb_Comp_2 0.0975 -1.965802          92.87255
#> Nb_Comp_3 0.0975 -1.000068          66.13838
#>           Chi2_Pearson_Y      RSS_Y pseudo_R2_Y
#> Nb_Comp_0       30.72004 0.35640772          NA
#> Nb_Comp_1       31.31448 0.11172576   0.6879757
#> Nb_Comp_2       33.06807 0.04650238   0.8671800
#> Nb_Comp_3       30.69727 0.01138837   0.9526757
#>                R2_Y
#> Nb_Comp_0        NA
#> Nb_Comp_1 0.6865226
#> Nb_Comp_2 0.8695248
#> Nb_Comp_3 0.9680468
```

###Bootstrap of the coefficients

Computing bootstrap distributions

```r
data("GasolineYield",package="betareg")
set.seed(1)
GazYield.boot <- bootplsbeta(modpls, sim="ordinary", stype="i", R=250)
```

Boxplots of the bootstrap distributions

```r
plsRglm::boxplots.bootpls(GazYield.boot)
```

<img src="man/figures/README-bootboxplots-1.png" title="plot of chunk bootboxplots" alt="plot of chunk bootboxplots" width="100%" />

Confidence intervals for the coefficients of the model based on the bootstrap distributions

```r
plsRglm::confints.bootpls(GazYield.boot)
#>                                                
#> Intercept -1.796887447 -1.298797470 -1.79109655
#> gravity    0.007803426  0.203529463 -0.03031919
#> pressure  -0.114413178  0.178241939 -0.10016933
#> temp10    -0.500300165 -0.196296503 -0.50450721
#> temp       0.634667387  0.964477695  0.64140043
#> batch1    -0.103808147  0.123669771 -0.09078670
#> batch2    -0.043844906  0.118181125 -0.05804124
#> batch3    -0.039650496  0.160223180 -0.02620071
#> batch4    -0.063189142  0.069329059 -0.05878901
#> batch5    -0.046868693  0.090317880 -0.04970864
#> batch6    -0.036189372  0.084497622 -0.04342852
#> batch7    -0.130445774  0.072421206 -0.10384760
#> batch8    -0.127087903  0.103619226 -0.09754607
#> batch9    -0.070998169  0.032240075 -0.06309787
#> batch10   -0.136043809  0.008565401 -0.14272981
#>                                                          
#> Intercept -1.32785762 -1.77750018 -1.31426124 -1.75724986
#> gravity    0.19625824 -0.01907056  0.20750687  0.01728695
#> pressure   0.23040737 -0.07051412  0.26006259 -0.22781373
#> temp10    -0.21483215 -0.34203983 -0.05236477 -0.40987882
#> temp       0.99074204  0.51679514  0.86613674  0.62994281
#> batch1     0.14234706 -0.08117396  0.15195980 -0.14041823
#> batch2     0.12705691 -0.04422306  0.14087509 -0.05179246
#> batch3     0.19773676 -0.06227605  0.16166141 -0.08981571
#> batch4     0.09310470 -0.09365068  0.05824304 -0.09749153
#> batch5     0.08458056 -0.02424531  0.11004389 -0.09315423
#> batch6     0.10265439 -0.03466126  0.11142165 -0.04818180
#> batch7     0.10180298 -0.16937270  0.03627788 -0.25198453
#> batch8     0.14968985 -0.17621892  0.07101700 -0.21517753
#> batch9     0.04674180 -0.08148215  0.02835751 -0.12674384
#> batch10    0.01478130 -0.13233313  0.02517798 -0.15107466
#>                       
#> Intercept -1.263641413
#> gravity    0.240794215
#> pressure   0.136939906
#> temp10    -0.175141922
#> temp       0.900503031
#> batch1     0.120479458
#> batch2     0.110789411
#> batch3     0.128573856
#> batch4     0.052650981
#> batch5     0.082446108
#> batch6     0.065003348
#> batch7     0.017661871
#> batch8     0.052435236
#> batch9     0.010888555
#> batch10    0.004957851
#> attr(,"typeBCa")
#> [1] TRUE
```

Plot of the confidence intervals for the coefficients of the model based on the bootstrap distributions

```r
plsRglm::plots.confints.bootpls(plsRglm::confints.bootpls(GazYield.boot))
```

<img src="man/figures/README-bootplotconfint-1.png" title="plot of chunk bootplotconfint" alt="plot of chunk bootplotconfint" width="100%" />
