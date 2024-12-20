# Bayesian mixed models, shrinkage, sparse and interaction kernel regression

[![Build
Status](https://travis-ci.org/MarniTausen/BayzR.svg?branch=master)](https://travis-ci.org/MarniTausen/BayzR)[![Coverage
Status
coveralls](https://coveralls.io/repos/github/MarniTausen/BayzR/badge.svg?branch=master)](https://coveralls.io/github/MarniTausen/BayzR?branch=master)

Sections in this page - [About R/bayz](#About%20R/bayz) - [Obtaining
R/bayz](#Obtaining%20R/bayz) - [A very short
tour](#A%20very%20short%20tour%20of%20R/bayz) - [Example](#examples) -
[Summarizing and using output](usingOutput.html)

# About R/bayz

R/bayz is in the basis a mixed model packages, but it can do much more
than what other common R mixed model packages offer. R/bayz has wide
support to insert covariance / similarity / relationship structures on
random effects, either through a supplied “kernel”, with modelled
structures such as spatial, unstructured and factor analytic covariance
structures, and variaous approaches to model heterogeneous variances.
For interactions of multiple random effects, R/bayz can build their
covariance structure as a Kronecker product, and these Kronecker
products are handled efficiently by avoiding to explicitly build the
product matrix. A common application of this in plant breeding is to
model Cultivar:Environment interactions with a Kronecker product of a
genomic and a enviromic similarity matrix. R/bayz can go beyond
interactions and Kronecker products of two terms, for instance, allowing
to add a Trait dimension as well, to fit large multi-trait models using
a Factor Analytic approach. R/bayz also supports using (multiple) large
sets of predictors that can be modelled with different shrinkage
regression options. Lastly, being Bayesian, R/bayz goes beyond standard
mixed model capabilities by allowing different shrinkage and sparse
shrinkage distributions on effects, as well as on kernels.

The R help ?bayz contains basic help information; full details are
available on \<ljanss.github.io/Rbayz\>.

# Obtaining R/bayz

R/bayz us currently not an officially (cran) released R package. It can
be downloaded and istalled with the options shown below. Not being on
cran implies that dependencies need to be installed manually, which is
also shown below.

### 1. Precompiled binary versions (Windows, Mac)

You can try your luck if you can install some of the following
precompiled binary packages for Windows and MacOS using the below
commands in the R terminal. If this fails, go to option 2 to install
source from Github.

Installing package dependencies:

``` r
install.packages("Rcpp")
install.packages("nlme")
install.packages("coda")
```

For Windows (compiled on Win11 64-bit)

``` r
install.packages("https://ljanss.github.io/Rbayz/Rbayz_0.9.0.zip",repos=NULL,type="win.binary")
```

For MacOS (compiled on Mac Silicon)

``` r
install.packages("https://ljanss.github.io/Rbayz/Rbayz_0.9.0.tgz",repos=NULL,type="mac.binary")
```

### 2. Install source from Github repository (all systems including linux)

This requires a development environment in R, which needs Rtools on
windows, or “command line tools” on MacOS, and the devtools package. On
linux the development tools may often be pre-installed. The below
commands run in the R terminal.

Installing devtools and package dependencies:

``` r
install.packages("devtools")
install.packages("Rcpp")
install.packages("nlme")
install.packages("coda")
```

Download and install/compile Rbayz using:

``` r
library(devtools)
devtools::install_github("ljanss/Rbayz")
```

# A very short tour of R/bayz

R/bayz has a single main function bayz() that accepts model formulas in
an extended R-formula syntax where all explanatory (right-hand-side)
variables are wrapped by a “function-like” to specify how to fit it in
the model. This may look like Yield ~ fx(Year) + rn(Variety) to fit
Yield with Year as a fixed factor and Variety as a random factor (fx()
and rn() imply the variable used should be (converted to) factor). For
random effects, variance structures can be added with a V= option inside
the rn(), for instance with a relationship / similarity matrix
rn(Variety,V=Gmat), which implies the variance-covariance structure
*G**m**a**t**σ*<sub>*g*</sub><sup>2</sup>.

Interactions between fators are specified using the colon, such as
fx(Year:Location) for fixed effects of Year-Location interactions.
R/bayz does not support automatic expansion with main effects, and main
effects, if desired, should be explicitly added in the model. For
interactions in random effects, the variance specification is expanded
to a product of terms (interpreted as Kronecker product), written as
rn(Variety:Environment, V=Gmat\*Emat), implying the variance-covariance
structure for the interaction effects as
*G**m**a**t* ⊗ *E**m**a**t**σ*<sub>*g**e*</sub><sup>2</sup>. Using
interactions and Kronecker products involving estimated unstructured
covariance matrices allows to fit multi-trait models by supplying data
in a “melted” form and specifying interactions with Trait. Interactions
and the matching variance structures can be specified to any order.
R/bayz can additionally fit large sets of predictors using the rr()
(ridge or random regression) model function, regular fixed (nested)
regressions using rg(), and random slope models using rs().

R/bayz is used as many other R model functions by capturing the output
object, and running methods on that output like summary(), plot(),
various methods to extract part of the estimates, and a contrast and
predict method. The summary() gives MCMC convergence diagnostics and
Highest Posterior Density regions for selected “traced” parameters,
while plot() produces trace plots and density plots.

Check the following manual pages for additional information on:

-   An extended tour of R/bayz with examples
-   Full syntax for the bayz() function call
-   Overview of all model-function terms and their options
-   Overview of all variance-model specifications
-   Specification of residual variance structures
-   Using the output with summary() and other methods
