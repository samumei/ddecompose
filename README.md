
- [`ddecompose`: Detailed Distributional Decompositions of Between-Group
  Differences in
  R](#ddecompose-detailed-distributional-decompositions-of-between-group-differences-in-r)
  - [Overview](#overview)
  - [Installation](#installation)
  - [Background](#background)
    - [Oacaxa-Blinder Decomposition](#oacaxa-blinder-decomposition)
    - [Reweighting Decomposition](#reweighting-decomposition)
    - [Sequential Decomposition](#sequential-decomposition)
    - [‘Doubly Robust’ Oaxaca-Blinder
      Decomposition](#doubly-robust-oaxaca-blinder-decomposition)
    - [Reweighted RIF Regression
      Decomposition](#reweighted-rif-regression-decomposition)
    - [Inference](#inference)
  - [Examples (unfinished)](#examples-unfinished)
    - [Oaxaca-Blinder Decomposition](#oaxaca-blinder-decomposition)
    - [Reweighting Decomposition](#reweighting-decomposition-1)
    - [Reweighted RIF Regression
      Decomposition](#reweighted-rif-regression-decomposition-1)
  - [Replication of Firpo, Fortin, and Lemieux
    (2018)](#replication-of-firpo-fortin-and-lemieux-2018)
    - [Loading Data](#loading-data)
    - [Reweighted RIF Regression Decomposition (Table
      4)](#reweighted-rif-regression-decomposition-table-4)
    - [Presentation of Results](#presentation-of-results)
  - [Credits](#credits)
  - [References](#references)

<!-- README.md is generated from README.Rmd. Please edit that file -->

# `ddecompose`: Detailed Distributional Decompositions of Between-Group Differences in R

## Overview

**ddecompose** implements the [Oaxaca
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

`dfl_decompose()` divides differences in distributional statistics into an
composition effect and a structure effect using inverse probability
weighting as introduced by [DiNardo, Fortin, and Lemieux
(1996)](https://www.jstor.org/stable/2171954). The function also allows
to sequentially decompose the composition effect into the contribution
of single covariates.

The package contains generic summary, print and plot functions for the
results and computes standard errors. This documentation provides a
brief overview of the functions implemented in the package. For a more
detailed discussion of decomposition methods including their respective
assumptions and limitations, refer to [Fortin, Lemieux, and Firpo
(2011)](https://economics.ubc.ca/wp-content/uploads/sites/38/2013/05/pdf_paper_nicole-fortin-decomposition-methods.pdf)
and Firpo et al. (2018).

## Installation

You can either install the CRAN version of `ddecompose`

``` r
install.packages("ddecompose")
```

or the latest development version from GitHub:

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
overcomes both shortcomings. Instead of modeling the conditional mean,
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
by industry) of group 0 but rather combine the marginal of group 0 with
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
16-9](https://docs.iza.org/dp6397.pdf) or [Rothe, 2015:
328](https://doi.org/10.1080/07350015.2014.948959)).

### Inference

Analytical standard errors are straightforward for the original OB
decomposition (see Jann,
[2005](https://boris.unibe.ch/69506/1/oaxaca_se_handout.pdf) and
[2008](https://journals.sagepub.com/doi/abs/10.1177/1536867X0800800401)).
[Firpo
(2007)](https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1468-0262.2007.00738.x)
and [Firpo and Pinto (2016)](https://www.jstor.org/stable/26609621)
develop asymptotic theory for the reweighting estimator. [Firpo, Fortin,
and Lemieux
(2009b)](https://www.econometricsociety.org/publications/econometrica/2009/05/01/unconditional-quantile-regressions/supp/6822_extensions_0.pdf)
and Firpo et al. (2018) do the same for the RIF estimator and for the
RIF decomposition terms, respectively. The authors propose to bootstrap
standard errors. `ddecompose` allows to bootstrap standard errors in both
`ob_decompose()` and `dfl_decompose()`. For the classical OB decomposition,
analytical standard errors are available, too.

## Examples (unfinished)

The following examples illustrate the workings of the main decomposition
functions in `ddecompose`. We use a sample of the National Longitudinal
Survey (NLSY) 79 containig wage data from 2000 of workers who were aged
35 to 43 in that year. The data is from O’Neill and O’Neill (2006) and
is used as an illustration of the Oxaca-Blinder mean decomposition in
Fortin, Lemieux, and Firpo (2011). The data contains 2655 male and 2654
female observations, respectively.

``` r
library(ddecompose)
data("nlys00")
```

### Oaxaca-Blinder Decomposition

We are interested in the change of mean log hourly wages from 1988-1990
and 2014-2016. The change is decomposed into a composition and wage
structure effect.

The parameter `statistic` specifies quantiles as our distributional
statistic of interest, while `probs` defines the probabilities of the
quantiles.

``` r
  

  model <- log(wage) ~ age + region + black + hispanic + education + years_worked_civilian + part_time + industry

  ob_decompose_estimate <- ob_decompose(formula = model,
                                      data = nlys00,
                                      group = female)
```

The `summary` method returns the estimated coefficients and the standard
errors. Per default, `summary` returns heteroscedasticity-consistent
standard errors estimated if the variance is not bootstrapped.

``` r
summary(ob_decompose_estimate)
```

`ddecompose` comes with a convenient plot function. It illustrates the
estimated UQPE across the distribution. The confidence bands are based
on the same standard errors as returned by `summary`.

-\> DO WE NEED PLOT FUNCTION FOR CLASSIC OB DECO?

``` r
plot(ob_decompose_estimate, varselect = c("education"))
```

### Reweighting Decomposition

Fortin, Lemieux, and Firpo (FLF, 2011, p. 79-88) decompose the increase
in US male wage inequality between the early 1980s and the early 2000s
using the CPS data. In this example, we replicate their results and
treat the observations from 1983 to 1985 as reference group that is
reweighted.

``` r
library(ddecompose)
data("men8305")
flf_model <- log(wage) ~ union*(education + experience) + education*experience
flf_male_inequality  <- dfl_decompose(flf_model,
                                  data = men8305,
                                  weights = weights,
                                  group = year)
```

We can summarize the results:

``` r
summary(flf_male_inequality)
```

Using `plot()`, we can illustrate the decomposition accross different
quantiles.

``` r
plot(flf_male_inequality)
```

### Reweighted RIF Regression Decomposition

*The influence function of quantiles, the mean, the variance, the Gini
coefficient, the interquantile range and the quantile ratio are
available in `ddecompose`. Additionally, the package enables calculating the
RIF for additional statistics with user-written functions (see example
below). [Cowell and Flachaire
(2007)](https://doi.org/10.1016/j.jeconom.2007.01.001), [Essama-Nssah &
Lambert (2012)](https://doi.org/10.1108/S1049-2585(2012)0000020009), and
[Rios-Avila (2020)](https://doi.org/10.1177/1536867X20909690) derive the
influence functions for an array of distributional statistics. More
information about RIF regressions can be found in the documentation of
the companion package [`rifreg`](https://github.com/samumei/rifreg).*

The same function is used to decompose RIF regressions. Setting
`rifreg = TRUE` and indicating `rifreg_statistic`, decomposes the gap at
the specified statistic. The default is `rifreg_statistic = "quantiles"`
and `rifreg_probs = c(1:9)/10`, meaning the gap is decomposed at each
decile.

``` r
  model <- log(wage) ~ age + region + black + hispanic + education + years_worked_civilian + part_time + industry

  ob_decompose_estimate <- ob_decompose(formula = model,
                              data = nlys00,
                              group = female, 
                              reweighting = TRUE, 
                              rifreg = TRUE,
                              bootstrap = TRUE,
                              bootstrap_iteration = 10)
```

To decompose the the difference in the Gini coefficient of two groups,
`rifreg_statistic = "gini"` needs to be indicated accordingly. The RIF
functions for the following functions are currently implemented:
`"quantiles"`, `"mean"`, `"variance"`, `"gini"`,
`"interquantile_range"`, and `"interquantile_ratio"`. However, `ob_decompose`
also allows to pass a custom RIF function for the decomposition, by
setting `rifreg_statistic = "custom"` and passing the custom function to
`custom_rif_function`.

#### Maybe add Reweighting, Weights, and Changing the Reference Group examples

To add reweighting in Setting `bootstrap=TRUE` bootstraps standard
errors by resampling from all observations and reestimating both the RIF
and the regression in every iteration. We can set number of
`bootstrap_iterations` and the number of `cores`.

``` r
  model <- log(wage) ~ age + region + black + hispanic + education + years_worked_civilian + part_time + industry

  ob_decompose_estimate <- ob_decompose(formula = model,
                              data = nlys00,
                              group = female, 
                              reweighting = TRUE)
```

## Replication of Firpo, Fortin, and Lemieux (2018)

To validate the functions and provide users with an additional example,
we replicate the reweighted RIF regression decomposition estimates in
Firpo et al. (2018, p. 30). In their empirical example, Firpo et
al. focus in changes in wage inequality in the US between 1988 and 2016.
Using a large sample of male log wages in 1988-1990 (268’494
observations) and in 2014-2016 (236’296 observation) based on the
Outgoing Rotation Group (ORG) supplement of the Current Population
Survey (CPS), they attribute changes in variance, gini, and
interquantile ranges between the two time periods to changes in
composition (e.g., union coverage, education, and experience) and
changes in wage structure.

This replication follows the Stata replication code, that the authors
published alongside the paper and that is available on their website.
(INSERT LINK). The highest wages in the public CPS data set are
unavailable due to privacy concerns. The Firpo et al. impute these wages
using a random uniform Pareto distribution (CHECK WORDING). Since random
numbers are generated differently in R and Stata, even with an
equivalent seed. Therefore, we perform the entire data preparation in
Stata up until the estimation of the decomposition, i.e., the “oaxaca”
commands in Stata. This ensures that changes in the estimation results
between our function and the orignal paper do not stem from different
input data.

To reproduce the code below, make sure you download the entire data set
and Stata code from the journal’s website, execute the code up to the
oaxaca command, save the data and load it into your R environment.

### Loading Data

``` r

set.seed(987421)


# Make sure you execute the Stata Code until the oaxaca commands 
# and save the data in the appropriate folder.
men8816_t4 <- readstata13::read.dta13("data-raw/FFL_2018/usmen8816_t4.dta")

# Removing redundant observations - we replicate the reweighting within the function
men8816_t4 <- men8816_t4[men8816_t4$time <= 1,] 
```

### Reweighted RIF Regression Decomposition (Table 4)

The model is specified as in the Stata files, containing weights,
computing bootstrapped standard errors with 100 iterations, and setting
a fixed bandwidth of 0.06 for the epanechnikov kernel density
estimation.

``` r

men8816_t4 <- men8816
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
                paste(grep("^marr", names(men8816_t4), value = TRUE), collapse = " + "), "+",
                paste(c("ed0", "ed1", "ed3", "ed4", "ed5"), collapse = " + "), "+",
                paste(grep("^ex[1-4]|^ex[6-9]", names(men8816_t4), value = TRUE), collapse = " + "), "+",
                paste(grep("^uned", names(men8816_t4), value = TRUE), collapse = " + "), "+",
                paste(grep("^unex", names(men8816_t4), value = TRUE), collapse = " + "), "+",
                paste(grep("^ex[1-9]ed", names(men8816_t4), value = TRUE), collapse = " + "), "+",
                paste(grep("^pub", names(men8816_t4), value = TRUE), collapse = " + "), "+",
                paste(grep("^indd(1|1e|[3-9]|10|11|13|14)(?!2)", 
                           names(men8816_t4), perl = TRUE, value = TRUE), collapse = " + "), "+",
                paste(grep("^occd", names(men8816_t4), value = TRUE), collapse = " + "))))


# Interquantile Ratio 90-10 
decompose_90_10  <- ddecompose::ob_decompose(formula = ffl_model_with_reweighting,
                                data = men8816_t4,
                                weights = eweight,
                                group = time,
                                reference_0 = TRUE,
                                rifreg = TRUE,
                                rifreg_statistic = "interquantile_range",
                                rifreg_probs = c(0.9, 0.1),
                                bw = 0.065,
                                kernel = "epanechnikov",
                                reweighting = TRUE,
                              reweighting_method = "fast_glm",
                                trimming = TRUE,
                                trimming_threshold = 100,
                                bootstrap = TRUE,
                                bootstrap_iterations = 100, 
                              cores = 4)

## Variance
decompose_variance  <- ddecompose::ob_decompose(formula = var_model,
                                 data = men8816_t4,
                                 weights = eweight,
                                 group = time,
                                 reference_0 = TRUE,
                                 rifreg = TRUE,
                                 rifreg_statistic = "variance",
                                 reweighting = TRUE,
                              reweighting_method = "fast_glm",
                                 trimming = TRUE,
                                 trimming_threshold = 100,
                                 bootstrap = TRUE,
                                 bootstrap_iterations = 100, 
                              cores = 4)

# Gini 
ffl_model_with_reweighting_gini <- update.formula(ffl_model_with_reweighting, exp(.) ~ .)

decompose_gini  <- ddecompose::ob_decompose(formula = ffl_model_with_reweighting_gini,
                               data = men8816_t4,
                               weights = eweight,
                               group = time,
                               rifreg = TRUE,
                               rifreg_statistic = "gini",
                               reweighting = TRUE,
                              reweighting_method = "fast_glm",
                               trimming = TRUE,
                               trimming_threshold = 100,
                               bootstrap = TRUE,
                               bootstrap_iterations = 100, 
                              cores = 4)
```

### Presentation of Results

Describe the results and to what extent they are similar to those of ffl
and where they differ and why.

Also note that the replication with the same imputed data and code to
replicate table 1-4 is available upon request.

MAKE SURE THAT YOU write down where the differences are arising:

- Data (reweighting)
- seed
- rounding
- density estimation
- SE (includes Reweighting)

MAYBE ALSO NOTE HOW MUCH THE DIFFERENCE IS, IN THE OTHER TABLES ETC. AND
THAT THESE FILES ARE AVAILABLE UPON REQUEST BZW: IN THE TEST_FILE XXX
-\> PUT IN A SEPARATE TEST FILE!

Looking at the plots, we see that our plots correspond to those
presented by Firpo, Fortin, and Lemieux (2009a: 965). The plot example
illustrates, how the plot from the generic rifreg::plot() function can
be further enhanced, for instance with a horizontal line indicating the
OLS coefficient for comparison.

This validation example illustrates that the rifreg package works as
intended in computing RIF regressions and reliably yields the expected
results.

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

Jann, Ben, 2005. “Standard errors for the Blinder-Oaxaca decomposition.”
*3rd German Stata Users’ Group Meeting 2005*. Available from
<https://boris.unibe.ch/69506/1/oaxaca_se_handout.pdf>.

Jann, Ben, 2008. “The Blinder–Oaxaca Decomposition for Linear Regression
Models”. *Stata Journal* 8: 435–479.

Kline, Pat, 2009. “Blinder-Oaxaca as a Reweighting Estimator”. *UC
Berkeley mimeo.*

Rios-Avila, Fernando (2020): “Recentered influence functions (RIFs) in
Stata: RIF regression and RIF decomposition.” *The Stata Journal* 20(1):
51-94.

Rothe, Christoph (2012): “Decomposing the Composition Effect. The Role
of Covariates in Determining Between-Group Differences in Economic
Outcomes.” *IZA Discussion Paper* No. 6397.

Rothe, Christoph (2015): “Decomposing the Composition Effect. The Role
of Covariates in Determining Between-Group Differences in Economic
Outcomes.” *Journal of Business & Economic Statistics* 33(3): 323-337.
