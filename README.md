
- [`ddeco`: Detailed decompositions of between-group differences in
  R](#ddeco-detailed-decompositions-of-between-group-differences-in-r)
  - [Overview](#overview)
  - [Installation](#installation)
  - [Background](#background)
    - [Oacaxa-Blinder Decomposition](#oacaxa-blinder-decomposition)
    - [Reweighting decomposition](#reweighting-decomposition)
    - [Sequential decomposition](#sequential-decomposition)
    - [‘Doubly Robust’ Oaxca-Blinder](#doubly-robust-oaxca-blinder)
    - [RIF decomposition](#rif-decomposition)
  - [Inference](#inference)
  - [Examples (unfinished)](#examples-unfinished)
    - [Oaxaca-Blinder Decomposition](#oaxaca-blinder-decomposition)
    - [Reweighting Decomposition](#reweighting-decomposition-1)
    - [Reweighted RIF Regression
      Decomposition](#reweighted-rif-regression-decomposition)
    - [Validation](#validation)
  - [Credits](#credits)
  - [References](#references)

<!-- README.md is generated from README.Rmd. Please edit that file -->

# `ddeco`: Detailed decompositions of between-group differences in R

## Overview

**ddeco** implements the [Oaxaca
(1973)](https://www.jstor.org/stable/2525981)-[Blinder
(1973)](https://www.jstor.org/stable/144855) decomposition method and
generalizations of it that decompose differences in distributional
statistics beyond the mean.

`ob_deco()` decomposes differences in the mean outcome between two
groups into one part explained by different covariates (composition
effect) and into another part due to differences in the way covariates
are linked to the outcome variable (wage structure effect). The function
further divides the two effects into the contribution of each covariate
and allows for weighted ‘doubly robust’ decompositions. For
distributional statistics beyond the mean, the function performs the RIF
decomposition proposed by [Firpo, Fortin, and Lemieux
(2018)](https://doi.org/10.3390/econometrics6020028).

`dfl_deco()` divides differences in distributional statistics into an
composition effect and a wage structure effect using inverse probability
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

You can either install the CRAN version of `ddeco`

``` r
install.packages("ddeco")
```

or the latest development version from GitHub:

``` r
devtools::install_github("samumei/ddeco")
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
i.e.  $E( \varepsilon_{g,i} | X_{1,i}, \ldots ,X_{k,i}) = 0$, and that
there is an overlap in observable characteristics across groups (‘common
support’).

The coefficients are estimated with OLS. Together with the sample means
of the covariates, one can derive a counterfactual mean  
$$\overline Y_C = \widehat \beta_{0,0} + \sum^K_{k=1}\overline X_{1,k} \widehat \beta_{0,k}$$
that would be observed if group $1$ had the same coefficients like group
$0$. By adding and subtracting the counterfactual, the observed
difference
$$\widehat\Delta^\mu_O = (\overline Y_1 - \overline Y_C) + (\overline Y_C - \overline Y_0) = \widehat\Delta^\mu_S + \widehat\Delta^\mu_C, $$
is divided into the aggregate wage structure effect
$$\widehat\Delta^\mu_S  = (\widehat \beta_{1,0} - \widehat \beta_{0,0}) + \sum^K_{k=1}\overline X_{1,k}(\widehat \beta_{1,k} - \widehat \beta_{0,k}),$$

that captures outcome differences due to different coefficients, and the
composition effect
$$\widehat\Delta^\mu_X = \sum^K_{k=1} (\overline X_{1,k} - \overline X_{0,k})\widehat \beta_{0,k}, $$

which accounts for mean differences of the covariates. Note that we
could also combine coefficients from group 1 with covariates of group 0
to define a counterfactual. Such a change in the “reference” group
generally leads to different results.

$\widehat\Delta^\mu_S$ and $\widehat\Delta^\mu_X$ denote the aggregate
composition terms. Since the method hinges on an additive linear model,
the terms can be further decomposed into the contribution of each
covariate $k$ in a *detailed decomposition*. For instance, the
contribution of covariate $k = 1$ to the structure effect is
$\overline X_{1,1}(\widehat \beta_{1,1} - \widehat \beta_{0,1})$ and
$(\overline X_{1,1} - \overline X_{0,1})\widehat \beta_{0,1}$ to the
composition effect.

### Reweighting decomposition

The Oaxaca-Blinder decomposition is limited to the mean. Moreover, the
decomposition terms are biased if the expectation of the outcome
conditional on the covariates is not linear. DiNardo, Fortin, and
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

can be rewritten where $P(g)$ and $P(g|x)$ correspond to the
(conditional) probabilities of belonging to group $g$. We can estimate
the reweighting factor using sample probabilities of each group in the
joint sample and conditional probability models (e.g., logit). The
fitted factors are then used to estimate weighted distributional
statistics of interest (e.g., mean, quantiles or Gini coefficient) in
the reference sample – group 0 in the present example. The resulting
counterfactual distributional statistic,
$\widehat\nu_C=\widehat\nu(F_{Y_C})$, is then contrasted with the
observed difference
$$\widehat\Delta_O^{\nu} = (\widehat\nu_1 - \widehat\nu_C) + (\widehat\nu_C - \widehat\nu_0) = \widehat\Delta_S^\nu + \widehat\Delta_X^\nu,$$

which yields again an aggregate wage structure effect and aggregate
composition effect.

The two decomposition terms account for the contribution of the
covariates and the conditional outcome distribution, assuming common
support and ignorability. The latter condition asserts that the
distribution of unobserved covariates, conditional on observed
covariates, is independent of the group.

### Sequential decomposition

In contrast to the Oaxaca-Blinder decomposition, where the contributions
of each covariate simply add up, detailed decompositions are not
straightforward in the reweighting framework. However, DFL show that we
can sequentially alter the covariates distributions to decompose the
composition effect into the contribution of single covariates. For
instance, assume we want to distinguish the effect of covariate $X_1$
(e.g., union status) from that of covariate $X_2$ (e.g., industry). We
begin again with the counterfactual distribution of group 0 with the
covariates distribution of group 1
$$F_{Y_{C}}(y) = \iint F_{Y_0}(y|x_1,x_2)dF_{X_{1,1}}(x_1|x_2)dF_{X_{1,2}}(x_2)$$

and introduce a second counterfactual where we combine the conditional
outcome distribution of group 0 as well as the conditional covariate
distribution of $X_1$ given $X_2$ (e.g., union coverage by industry) of
group 0 with the covariates distribution $X_2$ (e.g., industry
distribution) of group 1
$$F_{Y_{C,X_2}}(y) = \iint F_{Y_0}(y|x_1,x_2)dF_{X_{0,1}}(x_1|x_2)dF_{X_{1,2}}(x_2) $$

which can be expressed as the outcome distribution 0
$$F_{Y_{C,X_2}}(y) = \iint F_{Y_0}(y|x_1,x_2)dF_{X_{0,1}}(x_1|x_2)\Psi_{X_2}(x_2)dF_{X_{0,2}}(x_2),$$
reweighted by the factor
$$\Psi_{X_2}(x_2) = \frac{dF_{X_{1,2}}(x_2)}{dF_{X_{0,2}}(x_2)} =  \frac{P(g=0)P(g=1|x_2)}{P(g=1)P(g=0|x_2)}.$$

With the distributional statistics of the additional counterfactual, we
can divide the aggregate decomposition effect into the contribution of
each covariate
$$\widehat \Delta_X^{\nu} =  (\widehat \nu_C - \widehat \nu_{C,X_2}) + (\widehat \nu_{C,X_2} - \widehat \nu_0) = \widehat \Delta_{X_1}^\nu + \widehat \Delta_{X_2}^\nu,$$

However, sequential decompositions are path dependent because the
detailed composition effects attributed to single covariates depend on
the order of which we include the variables into the sequence. For
instance, it matters if we reweight union coverage by industry
($F_{X_{0,1}}(x_1 \mid x_2)$) or the industry given union coverage
($F_{X_{0,2}}(x_2 \mid x_1)$).

Moreover, we get different results if we do not manipulate the marginal
covariate distribution $X_2$ (e.g., industry distribution) but the
conditional distribution function of $X_1$ given $X_2$ (e.g., union
density by industry) to derive the counterfactual, e.g.,
$$F_{Y_{C,X_1}}(y) = \iint F_{Y_0}(y|x_1,x_2)dF_{X_{1,1}}(x_1|x_2)dF_{X_0,2}(x_2).$$
where we would reweight with a slightly different factor
$$\Psi_{X_1}(x_1,x_2) = \frac{dF_{X_{1,1}}(x_1|x_2)}{dF_{X_{0,1}}(x_1|x_2)} =  \frac{P(g=0|x_2)P(g=1|x_1,x_2)}{P(g=1|x_2)P(g=0|x_1,x_2)}.$$

### ‘Doubly Robust’ Oaxca-Blinder

A robust and path independent alternative for decompositions at the mean
has been introduced by [Barsky et
al. (2002)](https://www.jstor.org/stable/3085702). To produce a robust
estimate for the case that conditional mean function is non linear, they
propose combining the DFL reweighting with the Oaxaca-Blinder
decomposition. Additionally, this approach has the valuable side effect
of accounting for potential errors introduced by an incomplete inverse
probability weighting and the linear model specification, respectively.

The estimation proceeds in two steps. First, the reweighting function
$\widehat\Psi_X(x)$ that matches the characteristics of group $0$ ($1$)
to those of group $1$ ($0$) is derived. Second, the linear model and the
covariate means are estimated in the two observed samples as well as in
the reweighted sample $0$ ($1$), yielding $\widehat \beta_{C,0}$,
$\widehat \beta_{C,k}$, and $\overline X_{C,k}$.

The mean of the reweighted sample builds the main counterfactual to
derive the aggregate structure effect and composition effect,
respectively. The detailed decomposition terms are also evaluated with
respect to the statistics of the reweigthed sample, i.e. 

$$\widehat\Delta^\mu_{O,R} = (\widehat \beta_{1,0} - \widehat \beta_{C,0}) + \sum^K_{k=1} (\overline X_{1,k}\widehat \beta_{1,k} - \overline X_{C,k}\widehat \beta_{C,k}) + \sum^K_{k=1} (\overline X_{C,k}\beta_{C,k} - \overline X_{0,k}\widehat \beta_{0,k}) = \widehat\Delta^\mu_{S,R}  + \widehat\Delta^\mu_{X,R}.$$
These decomposition terms can be further decomposed into a structure and
into a composition effect, respectively. This separates the errors from
reweighting $\widehat\Delta^\mu_{S,e}$ and the linear specification
$\widehat\Delta^\mu_{X,e}$, respectively, and yields a “pure”
composition effect $\widehat\Delta^\mu_{X,p}$ and the “pure” structure
effect $\widehat\Delta^\mu_{S,p}$. Concretely, the additional
decomposition of the structure effect reads as

$$\widehat\Delta^\mu_{S,R} = (\widehat \beta_{1,0} - \widehat \beta_{C,0}) + \sum^K_{k=1}\overline X_{1,k}(\widehat \beta_{1,k} - \widehat \beta_{C,k}) + \sum^K_{k=1} (\overline X_{1,k} - \overline X_{C,k})\widehat \beta_{C,k}= \widehat\Delta^\mu_{S,p} + \widehat\Delta^\mu_{S,e}.$$

Similarly, the composition effect can be written as

$$\widehat\Delta^\mu_{X,R} = \sum^K_{k=1} (\overline X_{C,k} - \overline X_{0,k})\widehat \beta_{0,k} + (\widehat \beta_{C,0} - \widehat \beta_{0,0}) + \sum^K_{k=1}\overline X_{C,k}(\widehat \beta_{C,k} - \widehat \beta_{0,k}) = \widehat\Delta^\mu_{X,p} + \widehat\Delta^\mu_{X,e}.$$

**Text müsste hier noch konkreter werden, man könnte ihn auch kürzen.
Meine Intuition: Die “puren” Effekte sind ja v.a. relevant für die
detailierte Dekomposition. Das müsste man erwähnen. Der
Spezifikationsfehler misst inwiefern sich die geschätzen Koeffzienten
bei gleichbleinder Lohnstruktur alleine aufgrund der
Kovariatenverteilung verändern und lassen dadurch den “puren” Effekt für
jede Variable isolieren. Der “reweighting” Fehler wiederum zeigt an,
inwiefen die geschätzte Lohnstruktur anders auffällt, weil man die
Kovariateverteilung ungenau umgewichtet. Dadurch kann man den “puren”
Effekt der Struktur messen.**

**Weiter könnte man Kline (2010,
<https://eml.berkeley.edu/~pkline/papers/oaxaca6.pdf>) erwähnen, der
argumentiert weshalb es nun ‘doubly robust’. Grundsätzlich heisst es so,
weil a) man sich sowohl gegen die Misspezifikation des linearen Modelles
absichert und b) des Propensity Scores Schätzers. Wenn das lineare
Modell richtig ist, braucht es keine Umgewichtung bzw. die Umgewichtung
kann falsch spezifiziert sein. Umgekehrt hält der Schätzer bei richtiger
Umgewichtung, auch wenn das lineare Modell misspezifiziert ist. **

The reweighted regression decomposition improves the original OB method
in two aspects. Firstly, the estimation of the structure effect involves
comparing $\widehat \beta_{1,k}$ with the reweighted group $0$
estimates, $\widehat \beta_{C,k}$, instead of $\widehat \beta_{0,k}$.
This allows for a more accurate representation of the true underlying
structural differences between groups $0$ and $1$ regarding wage
structure effects (Fortin et al. 2011). The reweighting error indicates
how accurate the reweighting function $\widehat\Psi(X_i)$ is and will be
zero if the composition of group $1$ is equal to the reweighted group
$0$.

Second, the pure composition effect displays the differences between the
composition of group $0$ and the reweighted group $0$. If the model is
truly linear and correctly specified, the specification error is zero
and $\widehat\Delta^\mu_{X,p} = \widehat\Delta^\mu_{X,R}$. If the
specification error is not close to zero, it is advisable to revise the
model specification. While this decomposition method is doubly robust
and path independent, its main limitation is that it only allows for
decomposition at the mean.

### RIF decomposition

A path independent decomposition method that goes beyond the mean is the
reweighted RIF regression decomposition as proposed by Firpo et
al. (2018). Firpo et al. propose to approximate the conditional
expectation of the RIF given the explanatory variables with a linear
regression. The regression coefficients can be consistent estimates of
the average derivatives $\widehat{\alpha}(\nu)$ if the conditional
expectations of the RIF are linear in $X$ (see [Firpo et al.,
2009a](https://www.econometricsociety.org/publications/econometrica/2009/05/01/unconditional-quantile-regressions),
[Rothe 2015: 328](https://doi.org/10.1080/07350015.2014.948959)).

In the reweighted RIF regression decompositions, first the RIF of the
outcome variable $Y$ and a distributional statistic of interest $\nu$
are computed. Then, an OLS regression of the transformed outcome
variable $Y$ is run on the explanatory variables $X$. Thereafter, the
decomposition method is analogous to the original OB method:

$$\widehat\Delta^\nu_{O,R} = (\widehat \beta_{1,0} - \widehat \beta_{C,0}) + \sum^K_{k=1}\overline X_{1,k}(\widehat \beta_{1,k} - \widehat \beta_{C,k}) + \sum^K_{k=1} (\overline X_{1,k} - \overline X_{C,k})\widehat \beta_{C,k} + \sum^K_{k=1} (\overline X_{C,k} - \overline X_{0,k})\widehat \beta_{0,k} + (\widehat \beta_{C,0} - \widehat \beta_{0,0}) + \sum^K_{k=1}\overline X_{C,k}(\widehat \beta_{C,k} - \widehat \beta_{0,k}) = \widehat\Delta^\nu_{S,p} + \widehat\Delta^\nu_{S,e} + \widehat\Delta^\nu_{X,p} +  \widehat\Delta^\nu_{X,e}.  $$
with $\widehat \beta$ being the coefficients of the OLS regression of
the transformed outcome variable $Y$, the RIF variable, run on the
explanatory variables $X$ of group $0$ and $1$.

The RIF of quantiles, the mean, the variance, the Gini coefficient, the
interquantile range and the quantile ratio are available in `ddeco`.
Additionally, the package enables calculating the RIF for additional
statistics with user-written functions (see example below). [Cowell and
Flachaire (2007)](https://doi.org/10.1016/j.jeconom.2007.01.001),
[Essama-Nssah & Lambert
(2012)](https://doi.org/10.1108/S1049-2585(2012)0000020009), and
[Rios-Avila (2020)](https://doi.org/10.1177/1536867X20909690) derive the
influence functions for an array of distributional statistics. More
information about RIF regressions can be found in the documentation of
the [`rifreg`](https://github.com/samumei/rifreg) package written by the
same authors as `ddeco`.

It is possible to compute the RIF regression decomposition without
reweighting. However, it is strongly advisable to use reweighting and
assess the reweighting and specification error, as large errors
highlight inconsistent estimations.

## Inference

`ddeco` allows to bootstrap standard errors in both the `ob_deco()` and
`dfl_deco()` functions. Furthermore, we implemented analytical standard
errors for the original OB decomposition as introduced by Jann
([2005](https://boris.unibe.ch/69506/1/oaxaca_se_handout.pdf),
[2008](https://journals.sagepub.com/doi/abs/10.1177/1536867X0800800401).
Analytical standard errors can be nontrivial when the RIF introduces an
additional estimation step. In particular, this is the case for
quantiles where the density has to be estimated (see [Firpo, Fortin, and
Lemieux,
2009b](https://www.econometricsociety.org/publications/econometrica/2009/05/01/unconditional-quantile-regressions/supp/6822_extensions_0.pdf)).

## Examples (unfinished)

The following examples illustrate the workings of the main decomposition
functions in `ddeco`. We use a sample of the National Longitudinal
Survey (NLSY) 79 containig wage data from 2000 of workers who were aged
35 to 43 in that year. The data is from O’Neill and O’Neill (2006) and
is used as an illustration of the Oxaca-Blinder mean decomposition in
Fortin, Lemieux, and Firpo (2011). The data contains 2655 male and 2654
female observations, respectively.

``` r
library(ddeco)
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

  ob_deco_estimate <- ob_deco(formula = model,
                                      data = nlys00,
                                      group = female)
```

The `summary` method returns the estimated coefficients and the standard
errors. Per default, `summary` returns heteroscedasticity-consistent
standard errors estimated if the variance is not bootstrapped.

``` r
summary(ob_deco_estimate)
```

`ddeco` comes with a convenient plot function. It illustrates the
estimated UQPE across the distribution. The confidence bands are based
on the same standard errors as returned by `summary`.

-\> DO WE NEED PLOT FUNCTION FOR CLASSIC OB DECO?

``` r
plot(ob_deco_estimate, varselect = c("education"))
```

### Reweighting Decomposition

Fortin, Lemieux, and Firpo (FLF, 2011, p. 79-88) decompose the increase
in US male wage inequality between the early 1980s and the early 2000s
using the CPS data. In this example, we replicate their results and
treat the observations from 1983 to 1985 as reference group that is
reweighted.

``` r
library(ddeco)
data("men8305")
flf_model <- log(wage) ~ union*(education + experience) + education*experience
flf_male_inequality  <- dfl_deco(flf_model,
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

The same function is used to decompose RIF regressions. Setting
`rifreg = TRUE` and indicating `rifreg_statistic`, decomposes the gap at
the specified statistic. The default is `rifreg_statistic = "quantiles"`
and `rifreg_probs = c(1:9)/10`, meaning the gap is decomposed at each
decile.

``` r
  model <- log(wage) ~ age + region + black + hispanic + education + years_worked_civilian + part_time + industry

  ob_deco_estimate <- ob_deco(formula = model,
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
`"interquantile_range"`, and `"interquantile_ratio"`. However, `ob_deco`
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

  ob_deco_estimate <- ob_deco(formula = model,
                              data = nlys00,
                              group = female, 
                              reweighting = TRUE)
```

### Validation

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

Firpo, Sergio, Nicole M. Fortin, and Thomas Lemieux. 2007a.
“Unconditional Quantile Regressions.” *Technical Working Paper 339,
National Bureau of Economic Research*. Cambridge, MA.

Firpo, Sergio, Nicole M. Fortin, and Thomas Lemieux. 2009a.
“Unconditional Quantile Regressions.” *Econometrica* 77(3): 953-973.

Firpo, Sergio, Nicole M. Fortin, and Thomas Lemieux. 2009b. “Supplement
to ‘Unconditional Quantile Regressions’.” *Econometrica Supplemental
Material*, 77.

Fortin, Nicole M., Thomas Lemieux, and Sergio Firpo. 2011.
“Decomposition Methods in Economics.” *National Bureau of Economic
Research - Working Paper Series*, 16045.

Jann, Ben, 2005. “Standard errors for the Blinder-Oaxaca decomposition.”
*3rd German Stata Users’ Group Meeting 2005*. Available from
<https://boris.unibe.ch/69506/1/oaxaca_se_handout.pdf>.

Jann, Ben, 2008. “The Blinder–Oaxaca Decomposition for Linear Regression
Models”. *Stata Journal* 8: 435–479.

Rios-Avila, Fernando (2020): “Recentered influence functions (RIFs) in
Stata: RIF regression and RIF decomposition.” *The Stata Journal* 20(1):
51-94.

Rothe, Christoph (2015): “Decomposing the Composition Effect. The Role
of Covariates in Determining Between-Group Differences in Economic
Outcomes.” *Journal of Business & Economic Statistics* 33(3): 323-337.
