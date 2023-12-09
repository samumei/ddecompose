
- [`ddeco`: Detailed decompositions of between-group differences in
  R](#ddeco-detailed-decompositions-of-between-group-differences-in-r)
  - [Overview](#overview)
  - [Installation](#installation)
  - [Background](#background)
    - [Oaxaca-Blinder Decomposition (and Reweighted
      Variant)](#oaxaca-blinder-decomposition-and-reweighted-variant)
  - [Background (By DG)](#background-by-dg)
    - [Reweighted RIF Regression
      Decomposition](#reweighted-rif-regression-decomposition)
    - [Decomposition Through Reweighting (DiNardo, Fortin, and Lemieux
      1996)](#decomposition-through-reweighting-dinardo-fortin-and-lemieux-1996)
  - [Reweighting decompositions with
    `dfl_deco()`](#reweighting-decompositions-with-dfl_deco)
    - [Aggregate decomposition](#aggregate-decomposition)
    - [Sequential decomposition](#sequential-decomposition)
    - [Example](#example)
    - [Inference](#inference)
  - [Examples](#examples)
    - [Oaxaca-Blinder Decomposition](#oaxaca-blinder-decomposition)
    - [Reweighted RIF Regression
      Decomposition](#reweighted-rif-regression-decomposition-1)
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

`ob_deco()` divides differences in the mean outcome between two groups
into one part explained by different covariates (composition effect) and
into another part due to differences in the way covariates are linked to
the outcome variable (wage structure effect). The function further
divides the two effects into the contribution of each covariate and
allows for weighted ‘doubly robust’ decompositions. For distributional
statistics beyond the mean, the function performs the RIF decomposition
proposed by [Firpo, Fortin, and Lemieux
(2018)](https://doi.org/10.3390/econometrics6020028).

`dfl_deco()` divides differences in distributional statistics into an
composition effect and a wage structure effect using inverse probability
weighting as introduced by [DiNardo, Fortin, and Lemieux
(1996)](https://www.jstor.org/stable/2171954). The function also allows
to sequentially decompose the composition effect into the contribution
of single covariates.

The package contains generic summary, print and plot functions for the
results and computes bootstrapped standard errors. This documentation
provides a brief overview over the functions implemented in the package.
For a more detailed discussion of decomposition methods including their
respective assumptions and limitations, refer to [Fortin, Lemieux, and
Firpo
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

### Oaxaca-Blinder Decomposition (and Reweighted Variant)

The original decomposition method introduced by Oaxaca (1973) and
Blinder (1973) decomposes group differences in the mean of an outcome
variable (e.g. wage) <span style="color:red">between two groups
$g = A, B$.</span> The group differences in the mean
$\widehat\Delta^\mu_O$ are decomposed into two parts.One part, the
“composition effect” $\widehat\Delta^\mu_X$, is explained by differences
in the observable covariates (e.g. educational level or experience). The
other part, the “structure effect” $\widehat\Delta^\mu_S$, is due to
different returns to these covariates.

<span style="color:red"> The method linearly models the relationship
between outcome $Y$ and covariates $X$, i.e.
$$Y_{g,i} = X_{g,i}\beta_g + \varepsilon_{g,i}.$$ Moreover, it assumes
error term $\varepsilon$ is conditionally independent of $X$
($E( \varepsilon_{g,i} | X, g) = 0$). The coefficients are then
estimated with OLS. </span>

<span style="color:red">~~Given a group $A$ and $B$, the group
difference in the outcome variable $Y$ can be written as~~ Since
coefficients of linear regressions relate the mean of a covariate
directly to the mean of the outcome, the difference in the mean </span>

$$\widehat\Delta^\mu_O = \overline Y_B - \overline Y_A,$$

<span style="color:red"> can be expressed as  
$$\widehat\Delta^\mu_O = \widehat \beta_{B0} + \sum^K_{k=1}\overline X_{Bk} \widehat \beta_{Bk} - \beta_{A,0} + \sum^K_{k=1}\overline X_{A,k} \widehat \beta_{A,k}$$
where $\beta_{g}$ are the intercept and slope coefficents of covariates
$k = 1,..., K$. Adding and substracting the counterfactual mean
$$\overline Y_C = \beta_{A,0} + \sum^K_{k=1}\overline X_{B,k} \widehat \beta_{A,k}$$
the observed difference is diveded into the aggregate wage structure
effect \$^\_S \$ and and the aggregate composition effect
$\widehat\Delta^\mu_X$  
$$\widehat\Delta^\mu_O = \widehat\Delta^\mu_S + \widehat\Delta^\mu_C, $$
that can be again written as the product of estimated coefficients and
covariates means
$$\widehat\Delta^\mu_S  = (\widehat \beta_{B0} - \widehat \beta_{A0}) + \sum^K_{k=1}\overline X_{Bk}(\widehat \beta_{Bk} - \widehat \beta_{Ak})$$
$$\widehat\Delta^\mu_X = \sum^K_{k=1} (\overline X_{Bk} - \overline X_{Ak})\widehat \beta_{Ak}. $$
</span>

<span style="color:red">\~~with the composition and structure effect as:

$$\widehat\Delta^\mu_O = \underbrace{(\widehat \beta_{B0} - \widehat \beta_{A0}) + \sum^K_{k=1}\overline X_{Bk}(\widehat \beta_{Bk} - \widehat \beta_{Ak})}_{\widehat\Delta^\mu_S \text{ (Unexplained)}} + \underbrace{\sum^K_{k=1} (\overline X_{Bk} - \overline X_{Ak})\widehat \beta_{Ak}}_{\widehat\Delta^\mu_X \text{ (Explained)}} $$
where $\beta_{g}$ are the intercept and slope coefficents of covariates
$k = 1,..., K$ <span style="color:red">~~, in the OLS regression models
of group $g = A, B$. In this notation~~. \<\> The counterfactual
reference scenario is $\overline X_{Bk}\beta_{Ak}$, the covariates’ mean
of group $B$, mutliplied by the covariates’ coefficients of group $A$.
However, the reference group could also be the set as
$\overline X_{Ak}\beta_{Bk}$. \~~</span>

The terms \$^\_S \$ and $\widehat\Delta^\mu_X$ denote the , i.e. the
overall structure and composition effect. In the , the contribution of
each covariate $k$ to the structure and composition effect is estimated.
For instance, the contribution of covariate $k = 1$ to the structure
effect is $\overline X_{B1}(\widehat \beta_{B1} - \widehat \beta_{A1})$
and $(\overline X_{B1} - \overline X_{A1})\widehat \beta_{A1}$ to the
composition effect.

## Background (By DG)

The **ddeco** package allows to decompose observed differences in a
distributional statistics $\nu_t=\nu(F_{Y_t})$ of an outcome variable
$Y_t$ (e.g., hourly wages) between two groups $t = 0,1$.
$$\Delta_O^{\nu} = \nu_1 - \nu_0.$$ We are interested to what extent the
difference can explained by a different composition of observable
covariates and by differences in structure linking covariates to the
outcome variables. For this purpose, we define a conterfactual outcome
statistic that would be observed if group 0 had the same composition
like group 1 but otherwise the same structure linking covariates to the
outcome variable. By contrasting the observed statistics to the
counterfactual, we can indentify the composition effect $\Delta_S^\nu$
and the (wage) structure effect $\Delta_X^\nu$ that capture differencs
in the structure linking covariates to the outcome and in the difference
in the composition of covariates, respectively:
$$\Delta_O^{\nu} = (\nu_1 - \nu_C) + (\nu_C - \nu_0) = \Delta_S^\nu + \Delta_X^\nu.$$
Note, we could specifiy the counterfactual the other way round.

Under certain assumptions (see Fortin, Lemieux, and Firpo 2011)
(ignorability/CIA, overlapping support), we can derive…

Standard OB case. …

#### Reweighted Regression Decompostion

[Barsky et al. (2002)](https://www.jstor.org/stable/3085702) pointed out
that OB decompositions do not provide consistent estimates of the wage
structure and composition effect if the conditional mean function is non
linear. They propose using a reweighting approach as in DiNardo et
al. (1996) in the decomposition. In this reweighting approach, the
reweighting function $\widehat\Psi(X_i)$ makes the characteristics of
group $A$ similar to those of group $B$, defined as $X^C_A$. In
addition, the OLS regression coefficients are estimated for the
reweighted group $\overline X^C_A$, yielding $\widehat \beta^C_A$.
Analogous, group A can be used as reference group, thereby reweighting
the characteristics of group $B$ would be reweighted to reflect the
characteristics of group $A$.

To decompose the overall group difference in the reweighted regression
decomposition, $\overline X^C_A\widehat \beta^C_A$ is used as
counterfactual scenario. The overall gap is then decomposed as:

$$
\begin{aligned} \widehat\Delta^\mu_{O,R} &= (\widehat \beta_{B0} - \widehat \beta^C_{A0}) + \sum^K_{k=1} (\overline X_{Bk}\widehat \beta_{Bk} - \overline X^C_{Ak}\widehat \beta^C_{Ak})) + \sum^K_{k=1} (\overline X^C_{Ak}\beta^C_{Ak} - \overline X_{Ak}\widehat \beta_{Ak}) \\
 &= \widehat\Delta^\mu_{S,R}  + \widehat\Delta^\mu_{X,R}. 
\end{aligned}
$$ These decompostion parts are further divided into pure structure and
composition effects, $\widehat\Delta^\mu_{S,p}$ and
$\widehat\Delta^\mu_{X,p}$ respectively, and into reweighting error
$\widehat\Delta^\mu_{S,e}$ and specification error
$\widehat\Delta^\mu_{X,e}$. The structure effect
$\widehat\Delta^\mu_{S,R}$ can be written as

$$\widehat\Delta^\mu_{S,R} = \underbrace{(\widehat \beta_{B0} - \widehat \beta^C_{A0}) + \sum^K_{k=1}\overline X_{Bk}(\widehat \beta_{Bk} - \widehat \beta^C_{Ak})}_{\widehat\Delta^\mu_{S,p} \text{ (Pure structure effect)}} + \underbrace{\sum^K_{k=1} (\overline X_{Bk} - \overline X^C_{Ak})\widehat \beta^C_{Ak}.}_{\widehat\Delta^\mu_{S,e} \text{ (Reweighting error)}}$$

Similarly, the composition effect can be written as

$$\widehat\Delta^\mu_{X,R} = \underbrace{\sum^K_{k=1} (\overline X^C_{Ak} - \overline X_{Ak})\widehat \beta_{Ak}}_{\widehat\Delta^\mu_{X,p} \text{ (Pure composition effect)}} + \underbrace{(\widehat \beta^C_{A0} - \widehat \beta_{A0}) + \sum^K_{k=1}\overline X^C_{Ak}(\widehat \beta^C_{Ak} - \widehat \beta_{Ak}).}_{\widehat\Delta^\mu_{X,e} \text{ (Specification error)}}$$

The reweighted regression decomposition enhances the original OB
decomposition method in two aspects. First, the structure effect is
estimated by comparing $\widehat \beta_{B0}$ with the reweighted
estimates of group $A$, $\widehat \beta^C_{A0}$ instead of
$\widehat \beta_{A0}$. Therefore, the pure wage structure effect entails
the true underlying structural differences in group $A$ and $B$ (Fortin
et al. 2011). The reweighting error indicates how accurate the
reweighting function $\widehat\Psi(X_i)$ and will be zero if the
composition of group $B$ is equal to the reweighted group $A$.

Second, the pure composition effect displays the differences between the
composition of group $A$ and the reweighted group $A$. If the model is
truly linear and correctly specified, the specification error is zero
and $\widehat\Delta^\mu_{X,p} = \widehat\Delta^\mu_{X,R}$. If the
specification error is not close to zero, it is advisable to revise the
model specification. While the original OB decomposition only allows for
decomposition at the mean, the following methods provide means to
decompose other distributional metrics.

### Reweighted RIF Regression Decomposition

To estimate decompositions beyond the mean, we implemented a reweighted
RIF regression decomposition as proposed by Firpo et al. (2018). Firpo
et al. propose to approximate the conditional expectation of the RIF
given the explanatory variables with a linear regression. The regression
coefficients can be consistent estimates of the average derivatives
$\widehat{\alpha}(\nu)$ if the conditional expectations of the RIF are
linear in $X$ (see [Firpo et al.,
2009a](https://www.econometricsociety.org/publications/econometrica/2009/05/01/unconditional-quantile-regressions),
[Rothe 2015: 328](https://doi.org/10.1080/07350015.2014.948959)).

In the reweighted RIF regression decompositions, first the RIF of the
outcome variable $Y$ and a distributional statistic of interest $\nu$ is
computed. Then, an OLS regression of the transformed outcome variable
$Y$ is run on the explanatory variables $X$. Thereafter, the
decompostion method is analogous to the original OB method:

$$\widehat\Delta^\nu_{O,R} = \underbrace{(\widehat \beta_{B0} - \widehat \beta^C_{A0}) + \sum^K_{k=1}\overline X_{Bk}(\widehat \beta_{Bk} - \widehat \beta^C_{Ak})}_{\widehat\Delta^\nu_{S,p} \text{ (Pure structure effect)}} + \underbrace{\sum^K_{k=1} (\overline X_{Bk} - \overline X^C_{Ak})\widehat \beta^C_{Ak}}_{\widehat\Delta^\nu_{S,e} \text{ (Reweighting error)}} + \underbrace{\sum^K_{k=1} (\overline X^C_{Ak} - \overline X_{Ak})\widehat \beta_{Ak}}_{\widehat\Delta^\nu_{X,p} \text{ (Pure composition effect)}} + \underbrace{(\widehat \beta^C_{A0} - \widehat \beta_{A0}) + \sum^K_{k=1}\overline X^C_{Ak}(\widehat \beta^C_{Ak} - \widehat \beta_{Ak}).}_{\widehat\Delta^\nu_{X,e} \text{ (Specification error)}}$$
with $\widehat \beta$ being the coefficients of the OLS regression of
the transformed outcome variable $Y$ run on the explanatory variables
$X$ of group $A$ and $B$.

By default, the RIF of quantiles, the mean, the variance, the Gini
coefficient, the interquantile range and the quantile ratio are
available in `ddeco`. Moreover, the package allows to calculate the RIF
for additional statistics with user-written functions (see example
below). [Cowell and Flachaire
(2007)](https://doi.org/10.1016/j.jeconom.2007.01.001), [Essama-Nssah &
Lambert (2012)](https://doi.org/10.1108/S1049-2585(2012)0000020009), and
[Rios-Avila (2020)](https://doi.org/10.1177/1536867X20909690) derive the
influence functions for an array of distributional statistics. More
further information on RIF regressions can be found in the documentation
of the [`rifreg`](https://github.com/samumei/rifreg) package written by
the same authors as `ddeco`.

It is possible to compute the RIF regression decomposition without
reweighting. However, it is strongly advisable to use reweighting and
assess the reweighting and specification error, as large errors
highlight inconsistent estimations.

### Decomposition Through Reweighting (DiNardo, Fortin, and Lemieux 1996)

Another method to decompose group differences in an outcome variable
beyond the mean is using inverse propensity reweighting as proposed by
DiNardo, Fortin, and Lemieux (1996). We implement this decompostion
method in the function `dfl_deco()`. The procedure reweights the sample
distribution of a reference group such that the group’s covariates
distribution matches the covariates distribution of a counterfactual
group. Reweighting factors are derived by modelling the probability of
belonging to the one group instead of the other conditional on
covariates. The function allows detailed decompositions of the
composition effect by sequentially reweighting (conditional) covariate
distributions.

If group $A$ is set as reference, it is reweighted to match the
distribution of characteristics of group $B$. First, a logit model is
run to estimate the probability of belonging to group $B$ conditional on
covariates:

$$ Pr(D_B = 1 | X) = 1 - Pr(D_A = 1 | X ). $$

Next, the reweighting factor $\widehat\Psi(X)$ for observations in group
$A$ is estimated:

$$ \widehat\Psi(X_i) = \frac{\widehat Pr(D_B = 1 | X) / \widehat Pr(D_B = 1) }{\widehat Pr(D_B = 0 | X) / \widehat Pr(D_B = 0)}.$$
Note that the same estimation method for $\widehat\Psi(X)$ is applied
when computing the reweighting factor for the reweighted OB and RIF
regression decomposition discussed above. Using $\widehat\Psi(X)$ to
reweight the sample obersavtions in group $A$, the counterfactual
statistic of interest then be computed. For instance, to estimate the
overall wage gap at a specific quantile $Q_X$ and decompose it into
aggregate structure and composition effect following approach is
applied:

$$
\begin{aligned}
\widehat\Delta^{Q_X}_O &= Q_{B,X} - Q_{A,X} \\
&= \underbrace{Q^C_{A,X}- Q_{A,X} }_{\widehat\Delta^{Q_X}_S \text{ (Unexplained)}} + \underbrace{Q_{B,X} - Q^C_{A,X}}_{\widehat\Delta^{Q_X}_X \text{ (Explained)}}
\end{aligned}
$$

with $Q_{gX}$ being the estimated Quantile in groups $g = A, B$ and
$Q^C_AX$ being the quantile estimate in the reweighted group $A$. For
further details on this method, including the detailed decomposition,
refer to DiNardo et al. (1996) and Fortin et al. (2011, p. 63-69).

## Reweighting decompositions with `dfl_deco()`

### Aggregate decomposition

`dfl_deco()` uses inverse probability weighting to estimate
counterfactual distributional statistics. The counterfactual combines
the conditional outcome distribution of the reference group, here
defined as group 0, with the covariates distribution of the comparison
group, and can expressed as reweighted version of the reference
distribution
$$F_{Y_C}(y) = \int F_{Y_0}(y|x)dF_{X_1} (x)= \int F_{Y_0}(y|x)\Psi(x)dF_{X_0}(x),$$

where the reweighting factor
$\Psi(x) = \frac{P(t=0)P(t=1|X)}{P(t=1)P(t=0|X)}$ captures the inverse
probability of belonging to group 1. The reweighting factor can be
written as function of the probabilities of the binary group variable
$t$ and can, thus, readily estimated in the pooled sample of the two
group unsing conditional probability models. The counterfactual
statistics are then estimated with weighted estimator using the observed
data of group 0 and the fitted reweighting factors.

### Sequential decomposition

We can extend the the reweighting approach to sequentially decompose the
composition effect into the contribution of single covariates $X$ and
$Z$, respectively
$$\Delta_X^{\nu} =  (\nu_C - \nu_{C,X}) + (\nu_{C,X} - \nu_0) = \Delta_{X,X}^\nu + \Delta_{X,Z}^\nu,$$

where we define $\nu_C$ as the statistic of the counterfactual where we
combine the conditional outcome of group 0 with the covariates
distribution of group 1, i.e.
$$F_{Y_{C}}(y) = \iint F_{Y_0}(y|x,z)dF_{X_1}(x|z)dF_{Z_1}(x).$$

For the decomposition work, we require an additional counterfactual
$\nu_{C,X}$ that combines the conditional outcome distribution of group
0 and the distribution of covariate $X$ given $Z$ of the same group with
the marginal distribution of $Z$ of group 1, i.e.
$$F_{Y_{C,X}}(y) = \iint F_{Y_0}(y|x,z)dF_{X_0}(x|z)dF_{Z_1}(x).$$ Note,
sequential decompositions are path-dependent, i.e. the detailed
composition effects attributed to single covariates depend on the order
we include the variables into the sequence (i.e. if we integrate over
the conditional distribution of $X$ or that of $Z$, respecitvely).
Moreover, we get different results if we derive $\nu_{C,X}$ using the
conditional covariate distribution from the other group, e.g.
$$F_{Y_{C,X}}(y) = \iint F_{Y_0}(y|x,z)dF_{X_1}(x|z)dF_{Z_0}(x).$$

### Example

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

### Inference

`ddeco` allows to bootstrap standard errors in both the `ob_deco()` and
`dfl_deco()` function. Analytical standard errors can be nontrivial when
the RIF introduces an additional estimation step. In particular, this is
the case for quantiles where the density has to be estimated (see
[Firpo, Fortin, and Lemieux,
2009b](https://www.econometricsociety.org/publications/econometrica/2009/05/01/unconditional-quantile-regressions/supp/6822_extensions_0.pdf)).

## Examples

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

Rios-Avila, Fernando (2020): “Recentered influence functions (RIFs) in
Stata: RIF regression and RIF decomposition.” *The Stata Journal* 20(1):
51-94.

Rothe, Christoph (2015): “Decomposing the Composition Effect. The Role
of Covariates in Determining Between-Group Differences in Economic
Outcomes.” *Journal of Business & Economic Statistics* 33(3): 323-337.
