
- <a href="#ddecompose-detailed-distributional-decompositions"
  id="toc-ddecompose-detailed-distributional-decompositions"><code>ddecompose</code>:
  Detailed Distributional Decompositions</a>
  - <a href="#overview" id="toc-overview">Overview</a>
  - <a href="#installation" id="toc-installation">Installation</a>
  - <a href="#background" id="toc-background">Background</a>
    - <a href="#oacaxa-blinder-decomposition"
      id="toc-oacaxa-blinder-decomposition">Oacaxa-Blinder Decomposition</a>
    - <a href="#reweighting-decomposition"
      id="toc-reweighting-decomposition">Reweighting Decomposition</a>
    - <a href="#sequential-decomposition"
      id="toc-sequential-decomposition">Sequential Decomposition</a>
    - <a href="#doubly-robust-oaxaca-blinder-decomposition"
      id="toc-doubly-robust-oaxaca-blinder-decomposition">‘Doubly Robust’
      Oaxaca-Blinder Decomposition</a>
    - <a href="#reweighted-rif-regression-decomposition"
      id="toc-reweighted-rif-regression-decomposition">Reweighted RIF
      Regression Decomposition</a>
    - <a href="#inference" id="toc-inference">Inference</a>
  - <a href="#examples" id="toc-examples">Examples</a>
    - <a href="#oaxaca-blinder-decomposition"
      id="toc-oaxaca-blinder-decomposition">Oaxaca-Blinder Decomposition</a>
    - <a href="#doubly-robust-decompostion"
      id="toc-doubly-robust-decompostion">‘Doubly Robust’ Decompostion</a>
    - <a href="#rif-regression-decomposition"
      id="toc-rif-regression-decomposition">RIF Regression Decomposition</a>
    - <a href="#reweighting-decomposition-1"
      id="toc-reweighting-decomposition-1">Reweighting Decomposition</a>
  - <a href="#replication-of-firpo-fortin-and-lemieux-2018"
    id="toc-replication-of-firpo-fortin-and-lemieux-2018">Replication of
    Firpo, Fortin, and Lemieux (2018)</a>
    - <a href="#loading-data" id="toc-loading-data">Loading Data</a>
    - <a href="#reweighted-rif-regression-decomposition-table-4"
      id="toc-reweighted-rif-regression-decomposition-table-4">Reweighted RIF
      Regression Decomposition (Table 4)</a>
    - <a href="#discussion-of-results"
      id="toc-discussion-of-results">Discussion of Results</a>
  - <a href="#credits" id="toc-credits">Credits</a>
  - <a href="#references" id="toc-references">References</a>

<!-- README.md is generated from README.Rmd. Please edit that file -->

# `ddecompose`: Detailed Distributional Decompositions

<!-- badges: start -->

[![R-CMD-check](https://github.com/samumei/ddecompose/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/samumei/ddecompose/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/samumei/ddecompose/branch/master/graph/badge.svg)](https://app.codecov.io/gh/samumei/ddecompose?branch=master)
[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![CRAN
status](https://www.r-pkg.org/badges/version/ddecompose)](https://CRAN.R-project.org/package=ddecompose)
[![License: GPL (\>=
3)](https://img.shields.io/badge/License-GPL%20%28%3E%3D%203%29-blue.svg)](https://choosealicense.com/licenses/gpl-3.0/)
<!-- badges: end -->

## Overview

The R package **ddecompose** implements the [Oaxaca
(1973)](https://www.jstor.org/stable/2525981)-[Blinder
(1973)](https://www.jstor.org/stable/144855) decomposition method and
generalizations of it that decompose differences in distributional
statistics beyond the mean.

`ob_decompose()` decomposes differences in the mean outcome between two
groups into one part explained by different covariates (composition
effect) and into another part due to differences in the way covariates
are linked to the outcome variable (structure effect). The function
further divides the two effects into the contribution of each covariate
and allows for weighted ‘doubly robust’ decompositions. For
distributional statistics beyond the mean, the function performs the RIF
decomposition proposed by [Firpo, Fortin, and Lemieux
(2018)](https://doi.org/10.3390/econometrics6020028).

`dfl_decompose()` divides differences in distributional statistics into
an composition effect and a structure effect using inverse probability
weighting as introduced by [DiNardo, Fortin, and Lemieux
(1996)](https://www.jstor.org/stable/2171954). The function also allows
to sequentially decompose the composition effect into the contribution
of single covariates.

The package contains generic summary, print and plot functions for the
results and computes standard errors. This documentation provides a
brief overview of the functions implemented in the package. For a more
detailed discussion of decomposition methods, including their respective
assumptions and limitations, refer to [Fortin, Lemieux, and Firpo
(2011)](https://economics.ubc.ca/wp-content/uploads/sites/38/2013/05/pdf_paper_nicole-fortin-decomposition-methods.pdf)
and Firpo et al. (2018).

## Installation

You can hopefully soon install the CRAN version of `ddecompose`

``` r
install.packages("ddecompose")
```

Until now, refer to the latest development version from GitHub:

``` r
devtools::install_github("samumei/ddecompose")
```

## Background

### Oacaxa-Blinder Decomposition

The original decomposition method introduced by Oaxaca (1973) and
Blinder (1973) divides the difference in the mean of an outcome variable
(e.g., hourly wages) between two groups $g = 0, 1$ into a part explained
by differences in the mean of the covariates (e.g., educational level or
experience) and into another part due to different linear regression
coefficients (e.g., returns to education) that link the covariates to
the outcome variable.

The method linearly models the relationship between the outcome $Y$ and
covariates $X$
$$Y_{g,i} = \beta_{g,0} + \sum^K_{k=1}X_{k,i}\beta_{g,k} + \varepsilon_{g,i},$$
where $\beta_{g,0}$ is the intercept and $\beta_{g,k}$ are the slope
coefficients of covariates $k = 1,\ldots, K$. Moreover, it is assumed
that the error term $\varepsilon$ is conditionally independent of $X$,
i.e., $E( \varepsilon_{g,i} | X_{1,i}, \ldots ,X_{k,i}) = 0$, and that
there is an overlap in observable characteristics across groups (‘common
support’).

The coefficients are estimated with OLS. Together with the sample means
of the covariates, one can derive a counterfactual mean  
$$\overline Y_C = \widehat \beta_{0,0} + \sum^K_{k=1}\overline X_{1,k} \widehat \beta_{0,k}$$
that would be observed if group $1$ had the same coefficients like group
$0$. By adding and subtracting the counterfactual, the observed
difference
$$\widehat\Delta^\mu_O = (\overline Y_1 - \overline Y_C) + (\overline Y_C - \overline Y_0) = \widehat\Delta^\mu_S + \widehat\Delta^\mu_C, $$

is divided into the aggregate structure effect
$$\widehat\Delta^\mu_S  = (\widehat \beta_{1,0} - \widehat \beta_{0,0}) + \sum^K_{k=1}\overline X_{1,k}(\widehat \beta_{1,k} - \widehat \beta_{0,k}),$$

that captures outcome differences due to different coefficients, and the
composition effect
$$\widehat\Delta^\mu_X = \sum^K_{k=1} (\overline X_{1,k} - \overline X_{0,k})\widehat \beta_{0,k}, $$

which accounts for mean differences of the covariates. Note that we
could also combine the coefficients from group 1 with the covariates of
group 0 to define a counterfactual. Such a change in the “reference”
group generally leads to different results.

$\widehat\Delta^\mu_S$ and $\widehat\Delta^\mu_X$ denote the aggregate
composition terms. Since the method hinges on an additive linear model,
the terms can be further divided into the contribution of each covariate
$k$ in a *detailed decomposition*. For instance, the contribution of
covariate $k = 1$ to the structure effect is
$\overline X_{1,1}(\widehat \beta_{1,1} - \widehat \beta_{0,1})$, while
$(\overline X_{1,1} - \overline X_{0,1})\widehat \beta_{0,1}$ is the
covariate’s contribution to the composition effect.

### Reweighting Decomposition

The Oaxaca-Blinder decomposition is limited to the mean. Moreover, the
decomposition terms are biased if the expectation of the outcome
conditional on the covariates is not linear (see [Barsky et al.,
2002](https://www.jstor.org/stable/3085702)). DiNardo, Fortin, and
Lemieux (1996), DFL hereafter, propose an alternative approach that
overcomes both shortcomings. Instead of modelling the conditional mean,
the method uses inverse probability weighting to define a counterfactual
outcome distribution that combines the conditional outcome distribution
of one group and the covariates distribution of the other group. For
instance, if we are interested in the outcomes of group 0 with
covariates of group 1, we would reweight the outcome distribution of
group 0 such that its covariates distribution matches that of group 1
$$F_{Y_C}(y) = \int F_{Y_0}(y|x)dF_{X_1} (x)= \int F_{Y_0}(y|x)\Psi_X(x)dF_{X_0}(x).$$

By applying Bayes’ rule, the reweighting factor,
$$\Psi_X(x) = \frac{dF_{X_1}(x)}{dF_{X_0}(x)} = \frac{P(g=0)P(g=1|x)}{P(g=1)P(g=0|x)},$$

can be expressed in terms of $P(g)$ and $P(g|x)$, the (conditional)
probabilities of belonging to group $g$. This allows us to estimate the
reweighting factor using sample probabilities of each group and fitting
conditional probability models (e.g., logit) in the joint sample. The
estimated factors are then used to estimate weighted distributional
statistics of interest (e.g., mean, quantiles or Gini coefficient) in
the reference sample – group 0 in the present example. The resulting
counterfactual distributional statistic,
$\widehat\nu_C=\widehat\nu(F_{Y_C})$, is then contrasted with the
observed difference
$$\widehat\Delta_O^{\nu} = (\widehat\nu_1 - \widehat\nu_C) + (\widehat\nu_C - \widehat\nu_0) = \widehat\Delta_S^\nu + \widehat\Delta_X^\nu,$$

which yields again an aggregate structure effect and aggregate
composition effect.

The two decomposition terms account for the contribution of the
covariates and of the conditional outcome distribution, respectively,
assuming common support and ignorability. The latter condition asserts
that the distribution of unobserved covariates $\varepsilon$ conditional
on observed covariates $X$ is independent of group $g$.

### Sequential Decomposition

In contrast to the Oaxaca-Blinder decomposition, where the contributions
of each covariate simply add up, detailed decompositions are not
straightforward in the reweighting framework. However, DFL show that we
can sequentially alter the covariates distributions to decompose the
composition effect into the contribution of single covariates. For
instance, assume we want to distinguish the effect of covariate $X_1$
(e.g., union status) from that of covariate $X_2$ (e.g., industry). We
begin again with the counterfactual distribution based on the
conditional outcome distribution of group 0 and the covariates
distribution of group 1
$$F_{Y_{C}}(y) = \iint F_{Y_0}(y|x_1,x_2)dF_{X_{1,1}}(x_1|x_2)dF_{X_{1,2}}(x_2)$$

and introduce a second counterfactual where we combine the conditional
outcome distribution of group 0 as well as the conditional covariate
distribution of $X_1$ given $X_2$ (e.g., union coverage by industry) of
group 0 with the covariates distribution $X_2$ (e.g., employment by
industry) of group 1
$$F_{Y_{C,X_2}}(y) = \iint F_{Y_0}(y|x_1,x_2)dF_{X_{0,1}}(x_1|x_2)dF_{X_{1,2}}(x_2) $$

which can be expressed as the outcome distribution 0
$$F_{Y_{C,X_2}}(y) = \iint F_{Y_0}(y|x_1,x_2)dF_{X_{0,1}}(x_1|x_2)\Psi_{X_2}(x_2)dF_{X_{0,2}}(x_2),$$
reweighted by the factor
$$\Psi_{X_2}(x_2) = \frac{dF_{X_{1,2}}(x_2)}{dF_{X_{0,2}}(x_2)} =  \frac{P(g=0)P(g=1|x_2)}{P(g=1)P(g=0|x_2)}.$$

With the distributional statistics of the additional counterfactual, we
can divide the aggregate decomposition effect into the contribution of
each covariate
$$\widehat \Delta_X^{\nu} =  (\widehat \nu_C - \widehat \nu_{C,X_2}) + (\widehat \nu_{C,X_2} - \widehat \nu_0) = \widehat \Delta_{X_1}^\nu + \widehat \Delta_{X_2}^\nu.$$

Sequential decompositions are path dependent because the detailed
composition effects attributed to single covariates depend on the order
of which we include the variables into the sequence. For instance, it
matters if we reweight union coverage by industry $F(x_1|x_2)$ or the
industry employment given union coverage $F(x_2|x_1)$.

Moreover, we get different results if we do not combine the marginal
covariate distribution $X_2$ (e.g., industry employment) of group 1 with
the conditional distribution of $X_1$ given $X_2$ (e.g., union density
by industry) of group 0 but rather, combine the marginal of group 0 with
the conditional distribution of group 1 to derive the counterfactual,
e.g.,
$$F_{Y_{C,X_1}}(y) = \iint F_{Y_0}(y|x_1,x_2)dF_{X_{1,1}}(x_1|x_2)dF_{X_0,2}(x_2).$$
where we would reweight group 0 with a slightly different factor
$$\Psi_{X_1}(x_1,x_2) = \frac{dF_{X_{1,1}}(x_1|x_2)}{dF_{X_{0,1}}(x_1|x_2)} =  \frac{P(g=0|x_2)P(g=1|x_1,x_2)}{P(g=1|x_2)P(g=0|x_1,x_2)}.$$

### ‘Doubly Robust’ Oaxaca-Blinder Decomposition

A robust and path independent alternative for detailed decompositions at
the mean is to combine DFL reweighting with the linear Oaxaca-Blinder
method (see Fortin et al., 2011: 48-51). This approach has the valuable
side effect of accounting for potential errors introduced by an
incomplete inverse probability weighting and the linear model
specification, respectively.

The estimation proceeds in two steps. First, the reweighting function
$\widehat\Psi_X(x)$, which matches the characteristics of group $0$ (if
group 0 is the reference group) to those of group $1$, is derived.
Second, the linear model and the covariate means are estimated in the
two observed samples as well as in the reweighted sample $0$, yielding
$\widehat \beta_{C,0}$, $\widehat \beta_{C,k}$, and $\overline X_{C,k}$.

The mean of the reweighted sample builds the main counterfactual to
derive the aggregate structure and composition effects, respectively.
The detailed decomposition terms are also evaluated with respect to the
statistics of the reweighted sample, i.e.,

$$\widehat\Delta^\mu_{O,R} = (\widehat \beta_{1,0} - \widehat \beta_{C,0}) + \sum^K_{k=1} (\overline X_{1,k}\widehat \beta_{1,k} - \overline X_{C,k}\widehat \beta_{C,k}) + \sum^K_{k=1} (\overline X_{C,k}\beta_{C,k} - \overline X_{0,k}\widehat \beta_{0,k}) = \widehat\Delta^\mu_{S,R}  + \widehat\Delta^\mu_{X,R}.$$
These decomposition terms can be further decomposed into a structure and
into a composition effect, respectively. This separates the errors from
reweighting $\widehat\Delta^\mu_{S,e}$ and the linear specification
$\widehat\Delta^\mu_{X,e}$, respectively, and yields a “pure”
composition effect $\widehat\Delta^\mu_{X,p}$ and a “pure” structure
effect $\widehat\Delta^\mu_{S,p}$ for each covariate.

Specifically, the structure effect can be written as

$$\widehat\Delta^\mu_{S,R} = (\widehat \beta_{1,0} - \widehat \beta_{C,0}) + \sum^K_{k=1}\overline X_{1,k}(\widehat \beta_{1,k} - \widehat \beta_{C,k}) + \sum^K_{k=1} (\overline X_{1,k} - \overline X_{C,k})\widehat \beta_{C,k} = \widehat\Delta^\mu_{S,p} + \widehat\Delta^\mu_{S,e}.$$

By comparing $\beta_{1,k}$ to the coefficients of the reweighted group
$0$, $\widehat \beta_{C,k}$, the true underlying structure effect can be
identified. Moreover, the reweighting error, $\widehat\Delta^\mu_{S,e}$,
indicates how accurately the reweighting function $\widehat\Psi(X_i)$
reweights the covariates distribution of group $0$ to that of group $1$.
Thus, identifying the reweighting error allows to measure the “pure”
structure effect, $\widehat\Delta^\mu_{S,p}$. If the composition of
group $1$ is equal to the reweighted group $0$, the reweighting error is
zero.

Similarly, the additional decomposition of the composition effect reads
as

$$\widehat\Delta^\mu_{X,R} = \sum^K_{k=1} (\overline X_{C,k} - \overline X_{0,k})\widehat \beta_{0,k} + (\widehat \beta_{C,0} - \widehat \beta_{0,0}) + \sum^K_{k=1}\overline X_{C,k}(\widehat \beta_{C,k} - \widehat \beta_{0,k}) = \widehat\Delta^\mu_{X,p} + \widehat\Delta^\mu_{X,e}.$$

The specification error, $\widehat\Delta^\mu_{X,e}$, measures the extent
to which the estimated coefficients change due to a different
distribution of covariates if the (wage) structure remains the same.
Thereby, the “pure” effect, $\overline X_{C,k} - \overline X_{0,k}$, for
each covariate can be estimated, summing to the aggregate pure structure
effect, $\widehat\Delta^\mu_{X,p}$. If the model is truly linear, i.e.,
correctly specified, $\beta_C$ of the reweighted group $0$ will be
identical to $\beta_0$ and the specification error will be zero (Fortin
et al. 2011, p. 49-50).

The reweighted OB decomposition is “doubly robust” as it yields
consistent estimates even if either the linear model or the reweighting
estimator is misspecified. In contrast to the simple reweighting or
linear approach it does not hinge on a single correctly specified model.
While the reweighted OB decomposition is doubly robust and path
independent, it is limited to the mean.

### Reweighted RIF Regression Decomposition

A path independent method that goes beyond the mean is the RIF
decomposition of Firpo, Fortin, and Lemieux (2018). The approach
approximates the expected value of the ‘recentered influence function’
(RIF) of the distributional statistic (e.g., quantile, variance, or Gini
coefficient) of an outcome variable conditional on covariates with
linear regressions. RIF regression coefficients can be consistent
estimates of the marginal effect of a small change in the expected value
of a covariate to the distributional statistics of an outcome variable
(see [Firpo et al.,
2009a](https://www.econometricsociety.org/publications/econometrica/2009/05/01/unconditional-quantile-regressions)
and the documentation of the companion package
[`rifreg`](https://github.com/samumei/rifreg)). Thus, they can be used
to decompose between-group difference in distributional statistics.
Firpo et al. (2018) combine the RIF regressions again with the
reweighting estimator to avoid specification errors.

First, the approach computes the recentered influence function (RIF) of
the outcome variable $Y$ and a distributional statistic of interest
$\nu$. Then, an OLS regression of the transformed outcome variable $Y$
is run on the explanatory variables $X$. Thereafter, the decomposition
method is analogous to the original OB method:

$$\widehat\Delta^\nu_{O,R} =  \widehat\Delta^\nu_{S,p} + \widehat\Delta^\nu_{S,e} + \widehat\Delta^\nu_{X,p} +  \widehat\Delta^\nu_{X,e}.$$
With $\widehat\Delta^\nu_{S,p}$ and $\widehat\Delta^\nu_{S,e}$
estimating the pure structure effect and the reweighting error as

$$\widehat\Delta^\nu_{S,R} =  \widehat\Delta^\nu_{S,p} + \widehat\Delta^\nu_{S,e} = (\widehat \beta_{1,0} - \widehat \beta_{C,0}) + \sum^K_{k=1}\overline X_{1,k}(\widehat \beta_{1,k} - \widehat \beta_{C,k}) + \sum^K_{k=1} (\overline X_{1,k} - \overline X_{C,k})\widehat \beta_{C,k}$$

and $\widehat\Delta^\nu_{X,p}$ and $\widehat\Delta^\nu_{X,e}$ estimating
the pure coefficient effect and the specification error:

$$\widehat\Delta^\nu_{X,R} = \widehat\Delta^\nu_{X,p} + \widehat\Delta^\nu_{X,e} = \sum^K_{k=1} (\overline X_{C,k} - \overline X_{0,k})\widehat \beta_{0,k} + (\widehat \beta_{C,0} - \widehat \beta_{0,0}) + \sum^K_{k=1}\overline X_{C,k}(\widehat \beta_{C,k} - \widehat \beta_{0,k})$$
with the RIF regression coefficients $\widehat \beta$ and the covariate
means $\overline X$ of groups $0$, $1$, and $C$, the reweighted
reference group 0, respectively.

Again, the specification error increases if the conditional expectation
of the RIF is not well approximated by the linear model. Moreover, as
the RIF regression coefficients capture the marginal effects of small
location shifts in the covariates distribution on the distributional
statistic of the outcome variable, the specification error can be large
if the location shifts are substantial or if the between-group
differences in the covariates distribution relate to higher moments or
the dependence structure (see also [Rothe, 2012:
16-19](https://docs.iza.org/dp6397.pdf) or [Rothe, 2015:
328](https://doi.org/10.1080/07350015.2014.948959)).

### Inference

Analytical standard errors are straightforward for the original OB
decomposition under the assumption of independence between the groups
(see Jann, [2005](https://boris.unibe.ch/69506/1/oaxaca_se_handout.pdf)
and
[2008](https://journals.sagepub.com/doi/abs/10.1177/1536867X0800800401)).
[Firpo
(2007)](https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1468-0262.2007.00738.x)
and [Firpo and Pinto (2016)](https://www.jstor.org/stable/26609621)
develop asymptotic theory for the reweighting estimator. [Firpo, Fortin,
and Lemieux
(2009b)](https://www.econometricsociety.org/publications/econometrica/2009/05/01/unconditional-quantile-regressions/supp/6822_extensions_0.pdf)
and Firpo et al. (2018) do the same for the RIF estimator and for the
RIF decomposition terms, respectively. The authors propose to bootstrap
standard errors. `ddecompose` allows bootstrapping standard errors in
both `ob_decompose()` and `dfl_decompose()`. For the classical OB
decomposition, analytical standard errors are available, too.

## Examples

The following examples illustrate the operation of the main
decomposition functions in `ddecompose`. We use a sample from the
National Longitudinal Survey (NLSY) 79 containing wage data from the
year 2000 for workers who aged 35 to 43. The data are from O’Neill and
O’Neill (2006) and were used to illustrate the Oxaca-Blinder mean
decomposition in Fortin, Lemieux, and Firpo (2011). The data contains
2655 male and 2654 female observations, respectively.

``` r
library(ddecompose)
#> Lade nötiges Paket: ggplot2
#> Warning: Paket 'ggplot2' wurde unter R Version 4.3.2 erstellt
data("nlys00")
```

### Oaxaca-Blinder Decomposition

We decompose the gender wage gap into composition and structure effect,
using the standard Oaxaca-Blinder decomposition without reweighting.
First, we specify a wage structure model and then run the estimation.

``` r

model <- log(wage) ~  age + region + education + 
                        years_worked_civilian + years_worked_military + 
                        part_time + family_responsibility + industry 

gender_gap_decomposition <- ob_decompose(formula = model,
                                        data = nlys00,
                                        group = female)
```

In the default settings, heteroscedasticity-consistent standard errors
are computed and factor variables are not normalized. Group 0 is defined
as the reference group (the group with the lower factor level) and is
subtracted from group 1. With `summary()`, we can display the
decomposition formula and the estimation results.

``` r
summary(gender_gap_decomposition)
#> 
#> 
#> Oaxaca-Blinder decomposition of mean difference
#> between female == 'no' (group 1) and female == 'yes' (group 0). 
#> The reference group is 'yes'. 
#> 
#> Group 0: female == 'yes' (2654 observations) 
#> Group 1: female == 'no' (2655 observations) 
#> 
#> Composition Effect: (X1 - X0) * b0 
#>   Structure Effect: X1 * (b1 - b0) 
#> 
#> Aggregate decomposition:
#> 
#>                      Estimate Std. Error    CI [Low     High]
#> Observed difference 0.2333003 0.01466550 0.20455644 0.2620441
#> Composition effect  0.1323814 0.01420535 0.10453942 0.1602234
#> Structure effect    0.1009189 0.01655736 0.06846707 0.1333707
#> 
#> 
#> Observed difference:
#> 
#>                          Estimate  Std. Error     CI [Low        High]
#> (Intercept)            0.09813187 0.213742888 -0.32079649  0.517060232
#> age                   -0.19459845 0.231908532 -0.64913082  0.259933920
#> region                 0.03970035 0.027067240 -0.01335047  0.092751163
#> education              0.05379461 0.033738154 -0.01233096  0.119920179
#> years_worked_civilian  0.27105219 0.061730110  0.15006340  0.392040980
#> years_worked_military  0.01683904 0.003284003  0.01040251  0.023275567
#> part_time             -0.01470307 0.008433550 -0.03123253  0.001826383
#> family_responsibility  0.03558673 0.011343387  0.01335410  0.057819361
#> industry              -0.07250298 0.031630225 -0.13449708 -0.010508877
#> 
#> 
#> Structure effect:
#> 
#>                            Estimate  Std. Error      CI [Low        High]
#> (Intercept)            0.0981318689 0.213742888 -0.320796494  0.517060232
#> age                   -0.2011988641 0.231292303 -0.654523448  0.252125720
#> region                 0.0352173874 0.026838967 -0.017386022  0.087820796
#> education              0.0667616317 0.032733435  0.002605277  0.130917986
#> years_worked_civilian  0.2155020831 0.063469186  0.091104764  0.339899403
#> years_worked_military -0.0085748963 0.008335066 -0.024911325  0.007761532
#> part_time             -0.0298335540 0.005684065 -0.040974117 -0.018692991
#> family_responsibility  0.0008774603 0.004232174 -0.007417448  0.009172369
#> industry              -0.0759642206 0.029614806 -0.134008173 -0.017920268
#> 
#> 
#> Composition effect:
#> 
#>                           Estimate  Std. Error       CI [Low        High]
#> (Intercept)            0.000000000 0.000000000  0.0000000000  0.000000000
#> age                    0.006600415 0.001918589  0.0028400506  0.010360779
#> region                 0.004482959 0.002448367 -0.0003157525  0.009281671
#> education             -0.012967020 0.005805598 -0.0243457828 -0.001588258
#> years_worked_civilian  0.055550105 0.005518193  0.0447346451  0.066365565
#> years_worked_military  0.025413935 0.007335229  0.0110371494  0.039790721
#> part_time              0.015130482 0.004237526  0.0068250830  0.023435882
#> family_responsibility  0.034709271 0.008300282  0.0184410164  0.050977526
#> industry               0.003461243 0.005102148 -0.0065387825  0.013461268
```

`ddecompose` comes with a handy plotting function. To plot the overall
composition and structure effects, we need to set
`detailed_effects = FALSE`.

``` r
plot(gender_gap_decomposition, detailed_effects = FALSE)
```

<img src="man/figures/README-unnamed-chunk-7-1.png" width="100%" />

To change the reference group, we have to set `reference_0 = FALSE`.
Note that this does not change the subtraction terms in the
decomposition, i.e., (X1 - X0) and (b1 - b0). To change the direction of
the subtraction, the parameter `subtract_1_from_0` must be set to
`TRUE`. This parameter only changes the sign of the decomposition
results. We can specify the weights column of the `data` in the
`weights` parameter and normalize the factor variables by setting
`normalize_factors = TRUE`.

### ‘Doubly Robust’ Decompostion

To compute a more robust decomposition, we can add reweighting to the
decompostion. Thereby, we also estimate the specification and
reweighting error of the results. By default, the same function model is
used for reweighting as for the decompostion. However, it is advisable
to add a more flexible reweighting model, taking interactions into
account. The reweighting formula is added to the decomposition formula,
separated by a `|`. For this decompostion method, standard errors need
to be computed using bootstrap. By increasing the number of `cores` the
bootstrap replications are computed on, we can reduce computation time.
Per default, 100 bootstrap replications are calculated.

``` r

model_w_reweighting <- log(wage) ~  
                           age + region + education + 
                           years_worked_civilian + years_worked_military + 
                           part_time + family_responsibility + industry |
                           age + region + education + 
                           years_worked_civilian + years_worked_military + 
                           part_time + family_responsibility + industry +
                           education:region + age:education 



gender_gap_decomposition_w_reweighting <- 
    ob_decompose(formula = model_w_reweighting,
                 data = nlys00,
                 group = female, 
                 reweighting = TRUE, 
                 bootstrap = TRUE, 
                 cores = 4)
#> 
#> Bootstrapping standard errors...
```

The default method for fitting and predicting conditional probabilities,
used to derive the reweighting factor, is a logit model. However, you
can also use `reweighting_method = "fastglm"` for a faster logit model
computation, or `random_forest` for a different reweighting method.
Setting `trimming = TRUE` will trim observations with dominant
reweighting factor values. By default, reweighting factor values are
trimmed according to the rule of Huber, Lechner, and Wunsch (2013).
Thereby, the `trimming_threshold`, i.e., the maximum accepted relative
weight of the reweighting factor value (inverse probability weight) of a
single observation, is set to `sqrt(N)/N`, where `N` is the number of
observations in the reference group. The trimming threshold can also be
manually set to a numeric value.

We can aggregate the detailed effects displayed in the `summary()` and
`plot()` function. For example, if we want to separate personal and
contextual factors that explain the wage gap, we can aggregate these
variables in a list.

``` r

custom_aggregation <- list(`Personal Factors` = c("age", 
                                                  "education<10 yrs", 
                                               "educationHS grad (diploma)",
                                               "educationHS grad (GED)",
                                               "educationSome college",
                                               "educationBA or equiv. degree",
                                               "educationMA or equiv. degree",
                                               "educationPh.D or prof. degree",
                                               "part_time",
                                               "family_responsibility"), 
                           `Experience` = c("years_worked_civilian",
                                            "years_worked_military"), 
                           `Contextual Factors` = c("regionNorth-central",
                                                 "regionSouth",
                                                 "regionWest",
                                                 "industryManufacturing",
                                                 "industryEducation, Health, Public Admin.",
                                                 "industryOther services"))

plot(gender_gap_decomposition_w_reweighting, 
     custom_aggregation = custom_aggregation)
```

<img src="man/figures/README-unnamed-chunk-9-1.png" width="100%" />

If we add reweighting, the `plot()` and `summary()` functions will also
display the specification and reweighting error. To reduce the
reweighting error, a more flexible reweighting model is recommended.

### RIF Regression Decomposition

To decompose group differences beyond the mean with `ob_deco` we use RIF
regressions. In the following examples, we will analyze the changes in
wage inequality between 1983/85 and 2003/05 and assess which covariates
contribute to explaining these changes. First, we look at the changes in
the variance. Then, we decompose the wage gap at each decile. The data
used are from the Merged Outgoing Rotation Group of the Current
Population Survey, used in Fortin, Lemieux & Firpo (2011).

``` r

data("men8305")

model_rifreg <- log(wage) ~ union*(education + experience) + education*experience 

# Variance
variance_decomposition <- ob_decompose(formula = model_rifreg,
                 data = men8305,
                 group = year, 
                 reweighting = TRUE, 
                 rifreg_statistic = "variance",
                 bootstrap = TRUE, 
                 cores = 4)
#> 
#> 
#> The same model specification is used for decomposition and to compute reweighting factors.
#> Bootstrapping standard errors...

# Deciles
deciles_decomposition <- ob_decompose(formula = model_rifreg,
                 data = men8305,
                 group = year, 
                 reweighting = TRUE, 
                 rifreg_statistic = "quantiles",
                 rifreg_probs = c(1:9)/10,
                 bootstrap = TRUE, 
                 cores = 4)
#> 
#> 
#> The same model specification is used for decomposition and to compute reweighting factors.
#> Warning in density.default(x = dep_var, weights = weights/sum(weights, na.rm =
#> TRUE), : Selecting bandwidth *not* using 'weights'

#> Warning in density.default(x = dep_var, weights = weights/sum(weights, na.rm =
#> TRUE), : Selecting bandwidth *not* using 'weights'

#> Warning in density.default(x = dep_var, weights = weights/sum(weights, na.rm =
#> TRUE), : Selecting bandwidth *not* using 'weights'

#> Warning in density.default(x = dep_var, weights = weights/sum(weights, na.rm =
#> TRUE), : Selecting bandwidth *not* using 'weights'

#> Warning in density.default(x = dep_var, weights = weights/sum(weights, na.rm =
#> TRUE), : Selecting bandwidth *not* using 'weights'

#> Warning in density.default(x = dep_var, weights = weights/sum(weights, na.rm =
#> TRUE), : Selecting bandwidth *not* using 'weights'

#> Warning in density.default(x = dep_var, weights = weights/sum(weights, na.rm =
#> TRUE), : Selecting bandwidth *not* using 'weights'

#> Warning in density.default(x = dep_var, weights = weights/sum(weights, na.rm =
#> TRUE), : Selecting bandwidth *not* using 'weights'

#> Warning in density.default(x = dep_var, weights = weights/sum(weights, na.rm =
#> TRUE), : Selecting bandwidth *not* using 'weights'
#> 
#> Bootstrapping standard errors...

# Plotting the deciles
plot(deciles_decomposition)
```

<img src="man/figures/README-unnamed-chunk-10-1.png" width="100%" />

The RIF functions for the following statistics are currently
implemented: `"quantiles"`, `"mean"`, `"variance"`, `"gini"`,
`"interquantile_range"`, and `"interquantile_ratio"`. However,
`ob_decompose` also allows you to pass a custom RIF function for the
decomposition, by setting `rifreg_statistic = "custom"` and passing the
custom function to `custom_rif_function`. [Cowell and Flachaire
(2007)](https://doi.org/10.1016/j.jeconom.2007.01.001), [Essama-Nssah &
Lambert (2012)](https://doi.org/10.1108/S1049-2585(2012)0000020009), and
[Rios-Avila (2020)](https://doi.org/10.1177/1536867X20909690) derive the
influence functions for an array of distributional statistics. More
information about RIF regressions can be found in the documentation of
the companion package [`rifreg`](https://github.com/samumei/rifreg).

Custom RIF functions must specify a `dep_var` parameter for the outcome
variable $Y$, a `weights` for potential sample weights, and `probs`. If
they are not needed, they must be set to `NULL` in the function
definition (e.g. `probs = NULL`).

The following example shows how to write the RIF for the top 10 percent
income share and, then, to estimate the RIF regression decomposition
using this custom function. The formula for this specific RIF can be
found in Essam-Nssah & Lambert (2012) or Rios-Avila (2020).

``` r

# custom RIF function for top 10% percent income share
custom_top_inc_share <- function(dep_var,
                                 weights,
                                 probs = NULL,
                                 top_share = 0.1) {
  top_share <- 1 - top_share
  weighted_mean <- weighted.mean(
    x = dep_var,
    w = weights
  )
  weighted_quantile <- Hmisc::wtd.quantile(
    x = dep_var,
    weights = weights,
    probs = top_share
  )
  lorenz_ordinate <- sum(dep_var[which(dep_var <= weighted_quantile)] *
    weights[which(dep_var <= weighted_quantile)]) /
    sum(dep_var * weights)
  if_lorenz_ordinate <- -(dep_var / weighted_mean) * lorenz_ordinate +
    ifelse(dep_var < weighted_quantile,
      dep_var - (1 - top_share) * weighted_quantile,
      top_share * weighted_quantile
    ) / weighted_mean
  rif_top_income_share <- (1 - lorenz_ordinate) - if_lorenz_ordinate
  rif <- data.frame(rif_top_income_share, weights)
  names(rif) <- c("rif_top_income_share", "weights")
  return(rif)
}

custom_decomposition <- ob_decompose(formula = model_rifreg,
                 data = men8305,
                 group = year, 
                 reweighting = TRUE, 
                 rifreg_statistic = "custom",
                 custom_rif_function = custom_top_inc_share,
                 bootstrap = TRUE, 
                 cores = 4)
#> 
#> 
#> The same model specification is used for decomposition and to compute reweighting factors.
#> Bootstrapping standard errors...

plot(custom_decomposition)
```

<img src="man/figures/README-unnamed-chunk-11-1.png" width="100%" />

### Reweighting Decomposition

Now, we present the use of the other main function: `dfl_decompose()`.
Fortin, Lemieux, and Firpo (FLF, 2011, p. 79-88) decompose the increase
in U.S. male wage inequality between the early 1980s and the early 2000s
using the CPS data. In this example, we replicate their results and
treat the observations from 1983 to 1985 as the reference group, which
is then reweighted.

``` r
library(ddecompose)
data("men8305")

flf_model <- log(wage) ~ union*(education + experience) + education*experience
flf_male_inequality  <- dfl_decompose(flf_model, 
                                  data = men8305,
                                  weights = weights,
                                  group = year, 
                                  bootstrap = TRUE, 
                                  cores = 1)
#> Bootstrapping standard errors...
```

We can summarize the results:

``` r
summary(flf_male_inequality)
#> Decomposition of difference between year == '2003-2005' (group 1) and
#>  year == '1983-1985' (group 0)
#> 
#> Reweighted reference group: year == '1983-1985' (group 0) 
#>  
#> Composition effect accounts for between-group differences
#> in the distribution of the following covariates:
#> 
#> union, education, experience 
#> 
#> ---------------------------------------------------------------------------------
#> Decomposition of difference at conditional quantiles:
#> 
#> Observed difference: 
#> ---------------------------------------------------------------------------------
#>  Quantile  Estimate Std. Error Pointwise CI: [low   high] Uniform CI: [low
#>       0.1  0.062533   0.008742           0.045399 0.07967         0.029976
#>       0.2  0.031163   0.010749           0.010095 0.05223        -0.008869
#>       0.3  0.007079   0.010900          -0.014285 0.02844        -0.033516
#>       0.4 -0.002060   0.011053          -0.023724 0.01960        -0.043224
#>       0.5  0.004041   0.008321          -0.012268 0.02035        -0.026948
#>       0.6  0.022516   0.007955           0.006925 0.03811        -0.007110
#>       0.7  0.044445   0.011113           0.022664 0.06622         0.003059
#>       0.8  0.101837   0.010088           0.082065 0.12161         0.064267
#>       0.9  0.174073   0.011648           0.151243 0.19690         0.130694
#>    high]
#>  0.09509
#>  0.07120
#>  0.04767
#>  0.03910
#>  0.03503
#>  0.05214
#>  0.08583
#>  0.13941
#>  0.21745
#> 
#> Composition effect: 
#> ---------------------------------------------------------------------------------
#>  Quantile Estimate Std. Error Pointwise CI: [low   high] Uniform CI: [low
#>       0.1  0.03342   0.005798            0.02205 0.04478          0.01375
#>       0.2  0.05429   0.006767            0.04103 0.06756          0.03133
#>       0.3  0.07207   0.006848            0.05865 0.08549          0.04883
#>       0.4  0.07697   0.006376            0.06447 0.08947          0.05533
#>       0.5  0.09504   0.007078            0.08117 0.10891          0.07102
#>       0.6  0.08023   0.005606            0.06924 0.09122          0.06121
#>       0.7  0.08993   0.005379            0.07938 0.10047          0.07168
#>       0.8  0.11981   0.005915            0.10822 0.13140          0.09974
#>       0.9  0.12265   0.007971            0.10703 0.13828          0.09561
#>    high]
#>  0.05309
#>  0.07725
#>  0.09531
#>  0.09860
#>  0.11906
#>  0.09925
#>  0.10818
#>  0.13988
#>  0.14970
#> 
#> Structure effect: 
#> ---------------------------------------------------------------------------------
#>  Quantile Estimate Std. Error Pointwise CI: [low      high] Uniform CI: [low
#>       0.1  0.02911   0.006795            0.01580  0.0424315         0.001335
#>       0.2 -0.02313   0.011057           -0.04480 -0.0014571        -0.068333
#>       0.3 -0.06499   0.009420           -0.08346 -0.0465288        -0.103504
#>       0.4 -0.07903   0.009178           -0.09702 -0.0610402        -0.116553
#>       0.5 -0.09100   0.007778           -0.10624 -0.0757548        -0.122797
#>       0.6 -0.05771   0.006562           -0.07057 -0.0448517        -0.084541
#>       0.7 -0.04548   0.009258           -0.06363 -0.0273374        -0.083330
#>       0.8 -0.01797   0.009550           -0.03669  0.0007441        -0.057019
#>       0.9  0.05142   0.010937            0.02998  0.0728540         0.006705
#>      high]
#>   0.056893
#>   0.022076
#>  -0.026480
#>  -0.041506
#>  -0.059201
#>  -0.030885
#>  -0.007634
#>   0.021070
#>   0.096131
#> 
#> Decomposition of difference for other distributional statistics
#> 
#> Observed difference: 
#> ---------------------------------------------------------------------------------
#>                               Statistic Estimate Std. Error  CI [low    high]
#>                                    Mean  0.06205   0.006512  0.04929  0.07482
#>                                Variance  0.06211   0.004562  0.05317  0.07105
#>  Gini of untransformed Y (=exp(log(Y)))  0.04299   0.002310  0.03847  0.04752
#>             Interquantile range p90-p10  0.11154   0.013881  0.08433  0.13875
#>             Interquantile range p90-p50  0.17003   0.012453  0.14562  0.19444
#>             Interquantile range p50-p10 -0.05849   0.009937 -0.07797 -0.03902
#> 
#> Composition effect: 
#> ---------------------------------------------------------------------------------
#>                               Statistic Estimate Std. Error  CI [low    high]
#>                                    Mean 0.080587   0.004432 0.071902 0.089273
#>                                Variance 0.023498   0.002027 0.019526 0.027471
#>  Gini of untransformed Y (=exp(log(Y))) 0.007114   0.001102 0.004955 0.009274
#>             Interquantile range p90-p10 0.089236   0.008793 0.072001 0.106470
#>             Interquantile range p90-p50 0.027615   0.008781 0.010405 0.044825
#>             Interquantile range p50-p10 0.061621   0.006963 0.047973 0.075268
#> 
#> Structure effect: 
#> ---------------------------------------------------------------------------------
#>                               Statistic Estimate Std. Error    CI [low
#>                                    Mean -0.01853   0.005560 -0.0294307
#>                                Variance  0.03861   0.004586  0.0296252
#>  Gini of untransformed Y (=exp(log(Y)))  0.03588   0.002270  0.0314313
#>             Interquantile range p90-p10  0.02230   0.011890 -0.0009994
#>             Interquantile range p90-p50  0.14242   0.012948  0.1170390
#>             Interquantile range p50-p10 -0.12011   0.009017 -0.1377866
#>      high]
#>  -0.007634
#>   0.047601
#>   0.040330
#>   0.045607
#>   0.167795
#>  -0.102439
#> 
#> Summary statistics of reweighting factors
#> 
#> Number of trimmed observations (not included in statistics): 0 (0%)
#> 
#> Psi_X1: 
#> ---------------------------------------------------------------------------------
#>              Estimate Std. Error
#> Min.           0.1179         NA
#> 10%-quantile   0.3407    0.02310
#> 25%-quantile   0.6318    0.02295
#> 50%-quantile   0.7669    0.03478
#> 75%-quantile   1.2396    0.04051
#> 90%-quantile   1.7638    0.07983
#> Max.           3.2585         NA
```

Using `plot()`, we can illustrate the decomposition across different
quantiles.

``` r
plot(flf_male_inequality)
```

<img src="man/figures/README-unnamed-chunk-14-1.png" width="100%" />

## Replication of Firpo, Fortin, and Lemieux (2018)

To validate the functions and provide users with an additional example,
we replicate the reweighted RIF regression decomposition estimates in
Firpo et al. (2018, p. 30). In their empirical example, Firpo et
al. focus on changes in wage inequality in the US between 1988 and 2016.
Using a large sample of male log wages in 1988-1990 (268,494
observations) and in 2014-2016 (236,296 observations) based on the
Outgoing Rotation Group (ORG) supplement of the Current Population
Survey (CPS), the authors attribute changes in variance, Gini, and
interquantile ranges between the two time periods to changes in
composition (e.g., union coverage, education, and experience) and
changes in wage structure.

This replication follows the Stata replication code, that the authors
published alongside the paper and that is available on one author’s
[website](https://sites.google.com/view/nicole-m-fortin/data-and-programs)
and
[here](https://drive.google.com/file/d/1sab0RuBPRmch3DhreTraQj_3nUBwGMEJ/view).
The highest wages in the public CPS dataset are not available due to
privacy concerns. The Firpo et al. impute these wages using a random
uniform Pareto distribution. Since random numbers are generated
differently in R and Stata, even with an equivalent seed, we perform all
data preparation in Stata up to the estimation of the decomposition in
Stata, i.e., the “oaxaca” commands in stata. This ensures that changes
in the estimation results between our function and the original paper
are not due to different input data.

To reproduce the code below, be sure to download the entire dataset and
Stata code from the journal’s website, execute the code up to the oaxaca
command, save the data, and load it into your R environment.

### Loading Data

``` r
# Make sure you execute the Stata Code until the "oaxaca" commands 
# and save the data in the appropriate folder.
men8816_t4 <- readstata13::read.dta13("data-raw/usmen8816_t4.dta")

# Removing redundant observations - we replicate the reweighting within the function
men8816_t4 <- men8816_t4[men8816_t4$time <= 1,] 
```

### Reweighted RIF Regression Decomposition (Table 4)

The model is specified as in the Stata files, containing weights,
computing bootstrapped standard errors with 100 iterations, and setting
a fixed bandwidth of 0.06 for the Epanechnikov kernel density
estimation.

``` r

set.seed(987421)

ffl_model_with_reweighting <- as.formula(
    paste("lwage2 ~ covered + nonwhite + nmarr + 
          ed0 + ed1 + ed3 + ed4 + ed5 + ",
          paste(grep(paste0("^ex(", paste(c(1:4, 6:9), collapse = "|"), ")$"), 
                     names(men8816_t4), value = T), collapse = " + "), " + ",
          paste(grep(paste0("^occd(", paste(c(11:60, 80:91), collapse = "|"), ")$"),
                     names(men8816_t4), value = T), collapse = " + "), " + ",                     
          paste(grep(paste0("^indd(", paste(c(1, 3:14), collapse = "|"), ")$"), 
                     names(men8816_t4), value = T), collapse = " + "), " + pub | ",
          paste("covered + nonwhite +",
                paste(grep("^marr", 
                           names(men8816_t4), value = TRUE), collapse = " + "), "+",
                paste(c("ed0", "ed1", "ed3", "ed4", "ed5"), collapse = " + "), "+",
                paste(grep("^ex[1-4]|^ex[6-9]", 
                           names(men8816_t4), value = TRUE), collapse = " + "), "+",
                paste(grep("^uned", 
                           names(men8816_t4), value = TRUE), collapse = " + "), "+",
                paste(grep("^unex", 
                           names(men8816_t4), value = TRUE), collapse = " + "), "+",
                paste(grep("^ex[1-9]ed", 
                           names(men8816_t4), value = TRUE), collapse = " + "), "+",
                paste(grep("^pub", 
                           names(men8816_t4), value = TRUE), collapse = " + "), "+",
                paste(grep("^indd(1|1e|[3-9]|10|11|13|14)(?!2)", 
                           names(men8816_t4), perl = TRUE, value = TRUE), collapse = " + "), "+",
                paste(grep("^occd", names(men8816_t4), value = TRUE), collapse = " + "))))


# Interquantile Ratio 90-10 
decompose_90_10  <- ob_decompose(formula = ffl_model_with_reweighting,
                                data = men8816_t4,
                                weights = eweight,
                                group = time,
                                reference_0 = TRUE,
                                rifreg_statistic = "interquantile_range",
                                rifreg_probs = c(0.9, 0.1),
                                bw = 0.065,
                                kernel = "epanechnikov",
                                reweighting = TRUE,
                              reweighting_method = "fastglm",
                                trimming = TRUE,
                                trimming_threshold = 100,
                                bootstrap = TRUE,
                                bootstrap_iterations = 100)
## Variance
set.seed(23904875)
decompose_variance  <- ddecompose::ob_decompose(formula = ffl_model_with_reweighting,
                                 data = men8816_t4,
                                 weights = eweight,
                                 group = time,
                                 reference_0 = TRUE,
                                 rifreg_statistic = "variance",
                                 reweighting = TRUE,
                              reweighting_method = "fastglm",
                                 trimming = TRUE,
                                 trimming_threshold = 100,
                                 bootstrap = TRUE,
                                 bootstrap_iterations = 100)

## Gini 

# Updating the model
gini_model_raw <- update(ffl_model_with_reweighting, exp(.) ~ .)
gini_model_character <- as.character(gini_model_raw)
gini_model_split <- strsplit(gini_model_character, "~")
ffl_model_with_reweighting_gini <- 
  as.formula(paste(gini_model_split[[2]], "~",
                   gsub("\\(|\\)", "", gini_model_split[[3]][1])))

set.seed(130234976)
decompose_gini  <- ddecompose::ob_decompose(formula = ffl_model_with_reweighting_gini,
                               data = men8816_t4,
                               weights = eweight,
                               group = time,
                               rifreg_statistic = "gini",
                               reweighting = TRUE,
                              reweighting_method = "fastglm",
                               trimming = TRUE,
                               trimming_threshold = 100,
                               bootstrap = TRUE,
                               bootstrap_iterations = 100, 
                              cores = 1)
```

### Discussion of Results

``` r

# Presenting the results

variables <- decompose_variance[["variance"]][["decomposition_terms"]][["Variable"]]


ffl_aggregation <- list(`Union` = "covered", 
                        `Other` = c("nonwhite", "nmarr", 
                                    grep("ex", variables, value = TRUE)),
                        `Education` = grep("ed[0-9]", variables, value = TRUE), 
                        `Occupation` = grep("occd", variables, value = TRUE), 
                        `Industry` = grep("indd", variables, value = TRUE))

summary(decompose_90_10, custom_aggregation = ffl_aggregation)
#> 
#> 
#> Reweighted RIF regression decomposition of difference in interquantile_range
#> between time == '1' (group 1) and time == '0' (group 0). 
#> The reference group is '0'. 
#> 
#> Group 0: time == '0' (268492 observations) 
#> Group 1: time == '1' (236287 observations) 
#> Group C: time == '0' (reference group) reweighted
#>          to match the characteristics of the other group (268492 observations).
#> 
#> Pure Composition Effect: (XC - X0) * b0 
#>   Pure Structure Effect: XC * (bC - b0) 
#>     Specification Error: (X1 - XC) * bC 
#>       Reweighting Error: X1 * (b1 - bC) 
#> 
#> Aggregate decomposition:
#> 
#>                          Estimate  Std. Error      CI [Low         High]
#> Observed difference  0.1256140541 0.004701494  0.116399295  0.1348288131
#> Composition effect   0.0909463043 0.002916652  0.085229771  0.0966628379
#> Structure effect     0.0381813178 0.005043205  0.028296818  0.0480658174
#> Specification error -0.0005160126 0.003691050 -0.007750338  0.0067183127
#> Reweighting error   -0.0029975554 0.001399051 -0.005739645 -0.0002554653
#> 
#> 
#> Observed difference:
#> 
#>                      Estimate  Std. Error     CI [Low        High]
#> Union              0.02977302 0.001658415  0.02652259  0.033023458
#> Other             -0.02891464 0.010296904 -0.04909620 -0.008733084
#> Education          0.02497305 0.004515913  0.01612202  0.033824073
#> Occupation         0.06949078 0.007782856  0.05423666  0.084744897
#> Industry          -0.04775471 0.010027084 -0.06740744 -0.028101989
#> (Other variables)  0.07804656 0.014766600  0.04910456  0.106988565
#> 
#> 
#> Pure structure effect:
#> 
#>                       Estimate  Std. Error      CI [Low       High]
#> Union              0.012603651 0.001468132  0.009726165  0.01548114
#> Other             -0.046989402 0.012405103 -0.071302956 -0.02267585
#> Education          0.058884186 0.006324962  0.046487489  0.07128088
#> Occupation         0.009694848 0.009540193 -0.009003587  0.02839328
#> Industry          -0.079592423 0.011883520 -0.102883695 -0.05630115
#> (Other variables)  0.083580458 0.017346762  0.049581429  0.11757949
#> 
#> 
#> Pure composition effect:
#> 
#>                      Estimate   Std. Error     CI [Low       High]
#> Union             0.016196176 0.0006300468 0.014961307 0.017431045
#> Other             0.019067346 0.0012139205 0.016688106 0.021446587
#> Education         0.008089822 0.0015352464 0.005080794 0.011098850
#> Occupation        0.021190773 0.0013468292 0.018551036 0.023830510
#> Industry          0.024753937 0.0013971957 0.022015483 0.027492390
#> (Other variables) 0.001648250 0.0002679935 0.001122993 0.002173508
#> 
#> 
#> Specification error:
#> 
#>                        Estimate    Std. Error       CI [Low        High]
#> Union              0.0008615798  0.0008615798 -0.0008270856  0.002550245
#> Other             -0.0010756721 -0.0010756721  0.0010326065 -0.003183951
#> Education         -0.0423315700 -0.0423315700  0.0406367826 -0.125299922
#> Occupation         0.0410382285  0.0410382285 -0.0393952214  0.121471678
#> Industry           0.0082166703  0.0082166703 -0.0078877076  0.024321048
#> (Other variables) -0.0072252492 -0.0072252492  0.0069359791 -0.021386478
#> 
#> 
#> Reweighting error:
#> 
#>                        Estimate    Std. Error       CI [Low         High]
#> Union              1.116175e-04  1.116175e-04 -1.071488e-04  0.0003303839
#> Other              8.308368e-05  8.308368e-05 -7.975734e-05  0.0002459247
#> Education          3.306077e-04  3.306077e-04 -3.173715e-04  0.0009785869
#> Occupation        -2.433069e-03 -2.433069e-03  2.335659e-03 -0.0072017978
#> Industry          -1.132897e-03 -1.132897e-03  1.087540e-03 -0.0033533329
#> (Other variables)  4.310162e-05  4.310162e-05 -4.137601e-05  0.0001275793
#> 
#> Summary statistics of reweighting factors
#> 
#> Number of trimmed observations (not included in statistics): 0 (0%)
#> 
#>                   Psi_X1
#> Min.          0.01240973
#> 10%-quantile  0.24044303
#> 25%-quantile  0.42188219
#> 50%-quantile  0.73438754
#> 75%-quantile  1.25697488
#> 90%-quantile  2.00825420
#> Max.         22.41294385
summary(decompose_variance, custom_aggregation = ffl_aggregation)
#> 
#> 
#> Reweighted RIF regression decomposition of difference in variance
#> between time == '1' (group 1) and time == '0' (group 0). 
#> The reference group is '0'. 
#> 
#> Group 0: time == '0' (268492 observations) 
#> Group 1: time == '1' (236287 observations) 
#> Group C: time == '0' (reference group) reweighted
#>          to match the characteristics of the other group (268492 observations).
#> 
#> Pure Composition Effect: (XC - X0) * b0 
#>   Pure Structure Effect: XC * (bC - b0) 
#>     Specification Error: (X1 - XC) * bC 
#>       Reweighting Error: X1 * (b1 - bC) 
#> 
#> Aggregate decomposition:
#> 
#>                         Estimate   Std. Error      CI [Low       High]
#> Observed difference  0.077751372 0.0017482219  0.074324920 0.081177823
#> Composition effect   0.042647761 0.0012056833  0.040284665 0.045010857
#> Structure effect     0.033735963 0.0020922369  0.029635254 0.037836672
#> Specification error  0.002737420 0.0008451794  0.001080899 0.004393941
#> Reweighting error   -0.001369773 0.0007649579 -0.002869063 0.000129517
#> 
#> 
#> Observed difference:
#> 
#>                        Estimate  Std. Error     CI [Low        High]
#> Union              1.150440e-02 0.000717005  0.01009910  0.012909707
#> Other             -4.375222e-06 0.005854285 -0.01147856  0.011469812
#> Education          2.149355e-02 0.002355699  0.01687646  0.026110634
#> Occupation         6.005131e-02 0.003420293  0.05334766  0.066754957
#> Industry          -2.046241e-02 0.005448059 -0.03114041 -0.009784408
#> (Other variables)  5.168896e-03 0.007777351 -0.01007443  0.020412223
#> 
#> 
#> Pure structure effect:
#> 
#>                       Estimate   Std. Error      CI [Low        High]
#> Union              0.003532492 0.0007291712  0.002103343  0.004961642
#> Other             -0.009084234 0.0065988197 -0.022017683  0.003849215
#> Education          0.024947096 0.0032184078  0.018639133  0.031255059
#> Occupation         0.025238396 0.0040670528  0.017267119  0.033209673
#> Industry          -0.034523565 0.0063051274 -0.046881387 -0.022165742
#> (Other variables)  0.023625777 0.0085958611  0.006778199  0.040473356
#> 
#> 
#> Pure composition effect:
#> 
#>                      Estimate   Std. Error     CI [Low       High]
#> Union             0.007084019 0.0002531646 0.006587826 0.007580213
#> Other             0.010060413 0.0005488449 0.008984697 0.011136129
#> Education         0.006374034 0.0005910651 0.005215567 0.007532500
#> Occupation        0.007480422 0.0005493005 0.006403812 0.008557031
#> Industry          0.010320617 0.0005818809 0.009180152 0.011461083
#> (Other variables) 0.001328256 0.0001487646 0.001036683 0.001619830
#> 
#> 
#> Specification error:
#> 
#>                        Estimate    Std. Error       CI [Low        High]
#> Union              0.0008415311  0.0008415311 -0.0008078395  0.002490902
#> Other             -0.0008730819 -0.0008730819  0.0008381272 -0.002584291
#> Education         -0.0100362859 -0.0100362859  0.0096344730 -0.029707045
#> Occupation         0.0284495154  0.0284495154 -0.0273105101  0.084209541
#> Industry           0.0041833718  0.0041833718 -0.0040158863  0.012382630
#> (Other variables) -0.0198276304 -0.0198276304  0.0190338111 -0.058689072
#> 
#> 
#> Reweighting error:
#> 
#>                        Estimate    Std. Error       CI [Low         High]
#> Union              4.636017e-05  4.636017e-05 -4.450409e-05  0.0001372244
#> Other             -1.074720e-04 -1.074720e-04  1.031693e-04 -0.0003181133
#> Education          2.087049e-04  2.087049e-04 -2.003492e-04  0.0006177590
#> Occupation        -1.117026e-03 -1.117026e-03  1.072305e-03 -0.0033063578
#> Industry          -4.428319e-04 -4.428319e-04  4.251026e-04 -0.0013107664
#> (Other variables)  4.249232e-05  4.249232e-05 -4.079110e-05  0.0001257757
#> 
#> Summary statistics of reweighting factors
#> 
#> Number of trimmed observations (not included in statistics): 0 (0%)
#> 
#>                   Psi_X1
#> Min.          0.01240973
#> 10%-quantile  0.24044303
#> 25%-quantile  0.42188219
#> 50%-quantile  0.73438754
#> 75%-quantile  1.25697488
#> 90%-quantile  2.00825420
#> Max.         22.41294385
summary(decompose_gini, custom_aggregation = ffl_aggregation)
#> 
#> 
#> Reweighted RIF regression decomposition of difference in gini
#> between time == '1' (group 1) and time == '0' (group 0). 
#> The reference group is '0'. 
#> 
#> Group 0: time == '0' (268492 observations) 
#> Group 1: time == '1' (236287 observations) 
#> Group C: time == '0' (reference group) reweighted
#>          to match the characteristics of the other group (268492 observations).
#> 
#> Pure Composition Effect: (XC - X0) * b0 
#>   Pure Structure Effect: XC * (bC - b0) 
#>     Specification Error: (X1 - XC) * bC 
#>       Reweighting Error: X1 * (b1 - bC) 
#> 
#> Aggregate decomposition:
#> 
#>                          Estimate   Std. Error      CI [Low        High]
#> Observed difference  0.0659908665 0.0015135722  0.063024320 6.895741e-02
#> Composition effect   0.0199765902 0.0009934684  0.018029428 2.192375e-02
#> Structure effect     0.0451291062 0.0017981373  0.041604822 4.865339e-02
#> Specification error  0.0014785630 0.0005137232  0.000471684 2.485442e-03
#> Reweighting error   -0.0005933929 0.0003444451 -0.001268493 8.170719e-05
#> 
#> 
#> Observed difference:
#> 
#>                       Estimate   Std. Error       CI [Low       High]
#> Union              0.009960677 0.0005996969  0.0087852925 0.011136061
#> Other              0.003106513 0.0054771363 -0.0076284769 0.013841503
#> Education          0.004133062 0.0017218577  0.0007582833 0.007507842
#> Occupation         0.023964506 0.0023517073  0.0193552446 0.028573768
#> Industry          -0.005860137 0.0036927292 -0.0130977536 0.001377479
#> (Other variables)  0.030686245 0.0062595714  0.0184177107 0.042954780
#> 
#> 
#> Pure structure effect:
#> 
#>                        Estimate   Std. Error      CI [Low        High]
#> Union              0.0022173824 0.0006105551  0.001020716  0.003414048
#> Other             -0.0007297231 0.0061460445 -0.012775749  0.011316303
#> Education          0.0126849132 0.0018750239  0.009009934  0.016359893
#> Occupation         0.0130636777 0.0026466873  0.007876266  0.018251090
#> Industry          -0.0121351390 0.0041780504 -0.020323967 -0.003946311
#> (Other variables)  0.0300279951 0.0069379047  0.016429952  0.043626039
#> 
#> 
#> Pure composition effect:
#> 
#>                       Estimate   Std. Error      CI [Low       High]
#> Union             0.0063399542 0.0002157312 0.0059171289 0.006762779
#> Other             0.0047845108 0.0005063463 0.0037920903 0.005776931
#> Education         0.0019375518 0.0004602056 0.0010355653 0.002839538
#> Occupation        0.0013849785 0.0004112917 0.0005788616 0.002191095
#> Industry          0.0047148753 0.0004569234 0.0038193219 0.005610429
#> (Other variables) 0.0008147197 0.0001167025 0.0005859871 0.001043452
#> 
#> 
#> Specification error:
#> 
#>                        Estimate    Std. Error       CI [Low         High]
#> Union              0.0013650891  0.0013650891 -0.0013104364  0.0040406147
#> Other             -0.0009370702 -0.0009370702  0.0008995537 -0.0027736942
#> Education         -0.0105565134 -0.0105565134  0.0101338727 -0.0312468996
#> Occupation         0.0099999030  0.0099999030 -0.0095995468  0.0295993528
#> Industry           0.0017927361  0.0017927361 -0.0017209621  0.0053064343
#> (Other variables) -0.0001855816 -0.0001855816  0.0001781517 -0.0005493149
#> 
#> 
#> Reweighting error:
#> 
#>                        Estimate    Std. Error       CI [Low         High]
#> Union              3.825111e-05  3.825111e-05 -3.671969e-05  1.132219e-04
#> Other             -1.120429e-05 -1.120429e-05  1.075572e-05 -3.316430e-05
#> Education          6.711095e-05  6.711095e-05 -6.442409e-05  1.986460e-04
#> Occupation        -4.840530e-04 -4.840530e-04  4.646734e-04 -1.432779e-03
#> Industry          -2.326096e-04 -2.326096e-04  2.232969e-04 -6.885162e-04
#> (Other variables)  2.911195e-05  2.911195e-05 -2.794643e-05  8.617033e-05
#> 
#> Summary statistics of reweighting factors
#> 
#> Number of trimmed observations (not included in statistics): 0 (0%)
#> 
#>                   Psi_X1
#> Min.          0.01240973
#> 10%-quantile  0.24044303
#> 25%-quantile  0.42188219
#> 50%-quantile  0.73438754
#> 75%-quantile  1.25697488
#> 90%-quantile  2.00825420
#> Max.         22.41294385
```

The results presented below are similar to those in Table 4 of Firpo et
al. (2018, p. 30). However, some of the coefficients calculated here
differ by several percentages from the results in the paper. There are
several reasons for these differences.

1.  Reweighting: An important difference is the reweighting factors.
    Firpo et al. use the entire dataset to compute the reweighting
    factors. However, for the decomposition they remove very high wages
    from the dataset. In `ob_decompose()`, the same (trimmed) dataset is
    used for reweighting and the decomposition estimation.

2.  Different decomposition formula: In the paper, the formula presented
    for the pure structure and reweighting errors is identical to the
    formula presented in the background chapter above. For instance, the
    pure wage structure effect is computed as
    $(\widehat \beta_{1,0} - \widehat \beta_{C,0}) + \sum^K_{k=1}\overline X_{1,k}(\widehat \beta_{1,k} - \widehat \beta_{C,k})$.
    However, the Stata code calculates a slightly different formula:
    $(\widehat \beta_{C,0} - \widehat \beta_{1,0}) + \sum^K_{k=1}\overline X_{C,k}(\widehat \beta_{C,k} - \widehat \beta_{1,k})$.
    Thus, in the Stata code, the results are multiplied by -1 so that
    the composition and and structure effects add up to the observed
    difference. When we calculate the results in Stata, using the
    formula presented in the paper, our function produces very similar
    results to the Stata output, with a deviations of less than 0.1 in
    most cases.

3.  Density estimation: In statistics where quantiles need to be
    estimated (e.g., interquartile range), the differences are generally
    larger. We attribute this to the different density estimations in
    Stata and R (even when using the same kernel function and
    bandwidth). Specifically, the way the grids that define the locus of
    the density estimates are set. These differences result in different
    RIF values and thus different regression coefficients. Even tough,
    if the differences are mostly within a few percentages, it indicates
    that the RIF regression method is less robust for quantiles than for
    other statistics, since changes in the density estimation have a
    relatively large impact on the results.

The bootstrapped standard errors are relatively similar to those
reported in the paper. With only 100 bootstrap replications and
different seeds, some variance in the terms is not surprising. In
addition, we also included the reweighting procedure in the bootstrap
estimation, while Fortin et al. only include the RIF regression
estimation.

In summary, our `ob_decompose()` function produces very similar results
than those calculated in Stata and presented in the original paper.
Using an identical formula, the deviations are mostly below 0.1 percent.
However, some values based on RIF estimations of quantiles have slightly
higher differences. We have also replicated the results of Table 1-3 in
Firpo et al. (2018, p. 21-29), where the differences are generally even
smaller. The replication files are available upon request. This
validation example illustrates that the `ddecompose` package works as
intended in computing reweighted RIF regression decompositions and
reliably produces the expected results.

## Credits

David Gallusser & Samuel Meier

## References

Barsky, Robert, John Bound, Kerwin Kofi Charles, and Joseph P. Lupton.
2002. “Accounting for the Black-White Wealth Gap: A Nonparametric
Approach.” *Journal of the American Statistical Association* 97: 663–73.

Cowell, Frank A., and Emmanuel Flachaire. 2007. “Income distribution and
inequality measurement: The problem of extreme values.” *Journal of
Econometrics* 141: 1044–1072.

Essama-Nssah, Boniface, and Peter J. Lambert. 2012. “Influence Functions
for Policy Impact Analysis.” In John A. Bishop and Rafael Salas, eds.,
*Inequality, Mobility and Segregation: Essays in Honor of Jacques
Silber*. Bigley, UK: Emerald Group Publishing Limited.

Firpo, Sergio. 2007. “Efficient Semiparametric Estimation of Quantile
Treatment Effects.” *Econometrica* 75(1): 259-276.

Firpo, Sergio, Nicole M. Fortin, and Thomas Lemieux. 2007a.
“Unconditional Quantile Regressions.” *Technical Working Paper 339,
National Bureau of Economic Research*. Cambridge, MA.

Firpo, Sergio, Nicole M. Fortin, and Thomas Lemieux. 2009a.
“Unconditional Quantile Regressions.” *Econometrica* 77(3): 953-973.

Firpo, Sergio, Nicole M. Fortin, and Thomas Lemieux. 2009b. “Supplement
to ‘Unconditional Quantile Regressions’.” *Econometrica Supplemental
Material*, 77.

Firpo, Sergio, Nicole M. Fortin, and Thomas Lemieux. 2018. “Decomposing
Wage Distributions Using Recentered Influence Function Regressions.”
*Econometrics* 6(2):28.

Fortin, Nicole M., Thomas Lemieux, and Sergio Firpo. 2011.
“Decomposition Methods in Economics.” *National Bureau of Economic
Research - Working Paper Series*, 16045.

Firpo, Sergio, and Pinto, Christine. 2016. “Identification and
Estimation of Distributional Impacts of Interventions Using Changes in
Inequality Measures.” *Journal of Applied Econometrics*, 31(3), 457–486.

Huber, Martin, Michael Lechner, and Conny Wunsch. 2013. “The performance
of estimators based on the propensity score.” *Journal of Econometrics*
175(1): 1-21.

Jann, Ben. 2005. “Standard errors for the Blinder-Oaxaca decomposition.”
*3rd German Stata Users’ Group Meeting 2005*. Available from
<https://boris.unibe.ch/69506/1/oaxaca_se_handout.pdf>.

Jann, Ben. 2008. “The Blinder–Oaxaca Decomposition for Linear Regression
Models”. *Stata Journal* 8: 435–479.

Kline, Pat. 2009. “Blinder-Oaxaca as a Reweighting Estimator”. *UC
Berkeley mimeo.*

Rios-Avila, Fernando. 2020. “Recentered influence functions (RIFs) in
Stata: RIF regression and RIF decomposition.” *The Stata Journal* 20(1):
51-94.

Rothe, Christoph. 2012. “Decomposing the Composition Effect. The Role of
Covariates in Determining Between-Group Differences in Economic
Outcomes.” *IZA Discussion Paper* No. 6397.

Rothe, Christoph. 2015. “Decomposing the Composition Effect. The Role of
Covariates in Determining Between-Group Differences in Economic
Outcomes.” *Journal of Business & Economic Statistics* 33(3): 323-337.
