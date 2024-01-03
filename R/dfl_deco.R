#' DFL reweighting decomposition
#'
#' @description \code{dfl_deco} decomposes between-group differences in distributional
#' statistics of an outcome variable into a structure effect and a composition
#' effect. Following DiNardo, Fortin, and Lemieux (1996), the procedure reweights
#' the sample distribution of a reference group such that the group's covariates
#' distribution matches the covariates distribution of a comparison group.
#'
#' The function derives counterfactual distributions with inverse probability
#' weigthing. Reweighting factors are estimate by modelling the probability of
#' belonging to the comparison group conditional on covariates.
#'
#' The function allows detailed decompositions of the composition effect by
#' sequentially reweighting (conditional) covariate distributions. Standard
#' errors can be bootstrapped.
#'
#' @param formula a \code{formula} object with an outcome variable Y on the left-hand side
#' and the covariates X on the right-hand side. For sequential decompositions,
#' the sequence of covariates X are distinguished by the \code{|} operator. Covariates
#' are used to estimate the conditional probabilities for the reweighting factors.
#' @param data a \code{data.frame} containing all variables and observations of
#' both groups.
#' @param weights name of the observation weights variable or vector of
#' observation weights.
#' @param na.action a function to filter missing data (default \code{na.exclude}).
#' @param group name of a binary variable (numeric or factor) identifying the
#' two groups for which the differences are to be decomposed. The group
#' identified by the lower ranked value in \code{group} (i.e., 0 in the
#' case of a dummy variable or the first level of factor variable) is defined
#' as group 0. Per default, group 0 is the reference group (see \code{reference_0}).
#' @param reference_0 boolean: if \code{TRUE} (default), then the group 0 -- i.e.,
#' the group identified by the lower ranked value in \code{group} -- will be defined
#' as reference group. The reference group will be reweighted to match the
#' covariates distribution of the sample of the comparison group.
#' @param subtract_1_from_0 boolean: By default (`FALSE`), the distributional statistic
#' of group 0 is subtracted from the one of group 1 to compute the overall difference.
#' Setting `subtract_1_from_0` to `TRUE` merely changes the sign of the decomposition
#' results.
#' @param right_to_left determines the direction of a sequential decomposition.
#' If \code{TRUE} (default), the sequential decomposition starts right and reweights
#' first the reference group using only the variables entered last into the
#' \code{formula} sequence. Sequentially, the other variables are added. Otherwise,
#' the decomposition start left and using all variables entered into
#' \code{formula} object from the start, sequentially removing variables.
#' @param method specifies the method fit and predict conditional probabilities
#' used to derive the reweighting factor. At the moment, \code{"logit"}, \code{"fastglm"},
#' and \code{"random forests"} are available.
#' @param estimate_statistics boolean: if \code{TRUE} (default), then distributional
#' statistics are estimated and the decomposition is performed. If \code{FALSE},
#' the function only returns the fitted inverse propensity weights.
#' @param statistics a character vector that defines the distributional statistics
#' for which the decomposition is perforemd. Per default,
#' \code{c("quantiles", "mean", "variance", "gini", "iq_range_p90_p10", "iq_range_p90_p50", "iq_range_p50_p10")}
#' are estimated and decomposed. Also implemented are `c("iq_ratio_p90_p10", "iq_ratio_p90_p50", "iq_ratio_p50_p10")`.
#' Note: The function calculates the Gini coefficient for the untransformed variable
#' (i.e., \code{exp(log(Y))}), if the logarithm of a variable Y is set as outcome variable
#' in \code{formula}).
#' @param probs a vector of length 1 or more with the probabilities of the quantiles
#' to be estimated with default \code{c(1:9)/10}.
#' @param custom_statistic_function a custom statistic function to pass as argument.
#' @param trimming boolean: If \code{TRUE}, observations with dominant reweighting factor
#' values are trimmend according to rule of Huber, Lechner, and Wunsch (2013). Per
#' default, trimming is set to \code{FALSE}.
#' @param trimming_threshold numeric: threshold defining the maximal accepted
#' relative weight of the reweighting factor value (i.e., inverse probability weight)
#' of a single observation. If \code{NULL}, the threshold is set to \eqn{sqrt(N)/N},
#' where \eqn{N} is the number of observations in the reference group.
#' @param return_model boolean: If \code{TRUE} (default), the object(s) of the model
#' fit(s) used to predict the conditional probabilities for the reweighting factor(s)
#' are returned.
#' @param bootstrap boolean: If \code{FALSE} (default), then the estimation is not boostrapped
#' and no standard errors are calculated.
#' @param bootstrap_iterations positive integer with default \code{100} indicating the
#' number of bootstrap iterations to be executed.
#' @param bootstrap_robust boolean: if \code{FALSE} (default), then bootstrapped standard
#' errores are estimated as the standard deviations of the bootstrapp estimates.
#' Otherwise, the function uses the bootstrap interquartile range rescaled by the
#' interquantile range of the standard distribution to estimate standard errors.
#' @param cores positive integer with default \code{1} indicating the number of cores
#' to use when computing bootstrap standard errors.
#' @param ... other parameters passed to the function estimating the conditional probabilities.
#'
#' @details
#' The observed difference to be decomposed equals the difference between the values
#' of the distributional statistic of \code{group} 1 and \code{group} 0, respectively:
#'
#' \deqn{\DeltaO = \nu1 - \nu0,}
#'
#' where \eqn{\nut = \nu(F_g)} denotes the statistics of the outcome distribution
#' \eqn{F_g} of group \eqn{g}. Group 0 is identified by the lower ranked value
#' of the \code{group} variable.
#'
#' If \code{reference_0=TRUE}, then group 0 is the reference group and its observations
#' are reweighted such that they match the covariates distribution of group 1, the
#' comparison group. The counterfactual combines the covariates distribution
#' \eqn{F_1(x)} of group 1 with the conditional outcome distribution \eqn{F_0(y|x)}
#' of group 0 and is derived by reweighting group 0
#'
#' \deqn{F_C(y) = \int F_0(y|x)dF_1(x) = \int F_0(y|x)\Psi(x)dF_0(x),}
#'
#' where \eqn{\Psi(x)} is the reweighting factor, i.e. the inverse probabilities
#' of belonging to the comparison group conditional on covariates x.
#'
#' The distributional statistic of the counterfactual distribution,
#' \eqn{\nuC = \nu(F_C)}, allows to decompose the observed difference into
#' a (wage) structure effect (\eqn{\DeltaS = \nu1 - \nuC}) and a
#' composition effect (\eqn{\DeltaC = \nuC - \nu0}).
#'
#' If \code{reference_0=FALSE}, then the counterfactual is derived by combining
#' the covariates distribution of group 0 with the conditional outcome
#' distribution of group 1 and, thus, reweighting group 1
#'
#' \deqn{F_C(y) = \int F_1(y|x)dF_0(x) = \int F_1(y|x)\Psi(x)dF_1(x).}
#'
#' The composition effect becomes \eqn{\DeltaC = \nu1 - \nuC} and the
#' structure effect \eqn{\DeltaS = \nuC - \nu0}, respectively.
#'
#' The covariates are defined in \code{formula}. The reweighting factor is
#' estimated in the pooled sample with observations from both groups. \code{method = "logit"}
#' uses a logit model to fit the conditional probabilities. \code{method = "fastglm"}
#' also fits a logit model but with a faster algorithm from \strong{fastglm}.
#' \code{method = "random forests"} uses the \strong{Ranger} implementation of
#' the random forests classifier.
#'
#' The counterfactual statistics are then estimated with the observed data of
#' the reference group and the fitted reweighting factors.
#'
#' \code{formula} allows to specify interaction terms in the conditional
#' probability models. If you are interested in an aggregate decomposition,
#' then all covariates have to be entered at once, e.g. \code{Y ~ X + Z}.
#'
#' The procedure allows for sequential decomposition of the composition effect.
#' In this case, more than one reweighting factor based on different sets of
#' covariates are estimated.
#'
#' If you are interested in a sequential decomposition, the decomposition
#' sequence has to be distinguished by the \code{|} operator in the \code{formula}
#' object. For instance, \code{Y ~ X | Z} would decompose the aggregate composition
#' effect into the contribution of covariate(s) X and the one of covariate(s) Z,
#' respectively.
#'
#' In this two-fold sequential decomposition, we have the detailed composition
#' effects
#'
#' \deqn{\DeltaC_X = \nu1 - \nuCX}
#' and
#' \deqn{\DeltaC_Z = \nuCX - \nu_C}
#'
#' which sum up to the aggregate composition effect \eqn{\DeltaC}.
#'  \eqn{\nuC} is definied as above. It captures the contribution of all
#' covariates (i.e., X and Z). In contrast, \eqn{\nuCX} corresponds
#' to the statistic of the counterfactual distribution isolating the contribution
#' of covariate(s) X in contrast to the one of covariate(s) Z.
#'
#' If \code{right_to_left=TRUE}, then the counterfactual is defined as
#' \deqn{F_CX(y) = \iint F_0(y|x,z)dF_0(x|z)dF_1(z),}
#'
#' where \eqn{F_1(x|z)} is the conditional distribution of X given Z of
#' group 1 and \eqn{F_0(z)} the distribution of Z. If \code{right_to_left=FALSE},
#' we have
#' \deqn{F_CX(y) = \iint F_{0}(y|x,z)dF_1(x|z)dF_0(z).}
#'
#' Note that it is possible to specify the detailed models in every part of \code{formula}.
#' This is useful if you want to estimate in every step a fully saturated model,
#' e.g., \code{Y ~ X * Z | Z}. If not further specified, the variables are
#' additively included in the model used to derived the aggregate reweighting
#' factor.
#'
#' The detailed decomposition terms are path-dependent. The results depend on the sequence
#' the covariates enter the decomposition (e.g, \code{Y ~ X | Z} yields different
#' detailed decomposition terms than \code{Y ~ Z | X}) . Even for the same sequence,
#' the results differ depending on the 'direction' of the decomposition. In
#' the example above using \code{right_to_left=TRUE}, the contribution of Z is evaluated
#' using the conditional distribution of X given Z from group 0. If we use
#' \code{right_to_left=FALSE} instead, the same contribution is evaluated using
#' the conditional distribution from group 1.
#'
#' Per default, the distributional statistics for which the between group ifferences
#' are decomposed are quantiles, the mean, the variance, the Gini coefficient
#' and the interquantile range between the 9th and the 1st decile, the 9th decile
#' and the median, and the median and the first decile, respectively. The interquantile
#' ratios between the same quantiles are implemented, as well.
#'
#' The quantiles can be specified by \code{probs} that sets the corresponding
#' probabilities of the quantiles of interest. For other distributional statistics,
#' please use \code{custom_statistic_function}
#'
#' @return an object of class \code{dfl_deco} containing a data.frame with the
#' decomposition results for the quantiles and for the other distributional
#' statistics, respectively, a data.frame with the estimated reweighting factor
#' for every observation, a data.frame with sample quantiles of the reweighting
#' factors and a list with standard errors for the decomposition terms, the
#' quantiles of the reweighting factor as well as the bootstrapped
#' Kolmogorov-Smirnov distribution to construct uniform confidence bands for
#' quantiles.
#'
#' @references
#' DiNardo, John, Nicole M. Fortin, and Thomas Lemieux. 1996. "Labor Market
#' Institutions and the Distribution of Wages, 1973-1992: A Semiparametric Approach."
#' \emph{Econometrica}, 64(5), 1001-1044.
#'
#' Firpo, Sergio P., Nicole M. Fortin, and Thomas Lemieux. 2018. "Decomposing Wage
#' Distributions Using Recentered Influence Function Regressions."
#' \emph{Econometrics} 6(2), 28.
#'
#' Fortin, Nicole M., Thomas Lemieux, and Sergio Firpo. 2011. "Decomposition methods in economics."
#' In Orley Ashenfelter and David Card, eds., \emph{Handbook of Labor Economics}. Vol. 4. Elsevier, 1-102.
#'
#' Firpo, Sergio P., and Cristine Pinto. 2016. "Identification and Estimation of
#' Distributional Impacts of Interventions Using Changes in Inequality Measures."
#' \emph{Journal of Applied Econometrics}, 31(3), 457-486.
#'
#' Huber, Martin, Michael Lechner, and Conny Wunsch. 2013. "The performance of
#' estimators based on the propensity score." \emph{Journal of Econometrics},
#' 175(1), 1-21.
#'
#' @export
#'
#' @examples
#' ## Example from handbook chapter of Fortin, Lemieux, and Firpo (2011: 67)
#' ## with a sample of the original data
#'
#' data("men8305")
#'
#' flf_model <- log(wage) ~ union*(education + experience) + education*experience
#'
#' # Reweighting sample from 1983-85
#' flf_male_inequality  <- dfl_deco(flf_model,
#'                                  data = men8305,
#'                                  weights = weights,
#'                                  group = year)
#'
#' # Summarize results
#' summary(flf_male_inequality)
#'
#' # Plot decomposition of quantile differences
#' plot(flf_male_inequality)
#'
#' # Use alternative reference group (i.e., reweight sample from 2003-05)
#' flf_male_inequality_reference_0305  <- dfl_deco(flf_model,
#'                                                 data = men8305,
#'                                                 weights = weights,
#'                                                 group = year,
#'                                                 reference_0 = FALSE)
#' summary(flf_male_inequality_reference_0305)
#'
#' # Bootstrap standard errors (using smaller sample for the sake of illustration)
#' \dontrun{
#' set.seed(123)
#' flf_male_inequality_boot  <- dfl_deco(flf_model,
#'                                  data = men8305[1:1000, ],
#'                                  weights = weights,
#'                                  group = year,
#'                                  bootstrap = TRUE,
#'                                  bootstrap_iterations = 100,
#'                                  cores = 1)
#'
#' # Get standard errors and confidence intervals
#' summary(flf_male_inequality_boot)
#'
#' # Plot quantile differences with pointwise confidence intervals
#' plot(flf_male_inequality_boot)
#'
#' # Plot quantile differences with uniform confidence intervals
#' plot(flf_male_inequality_boot, uniform_bands=TRUE)
#'
#' }
#'
#'
#' ## Sequential decomposition
#'
#' # Here we distinguish the contribution of education and experience
#' # from the contribution of unionization conditional on education and experience.
#'
#' model_sequential <- log(wage) ~ union*(education + experience) +
#'                                   education*experience |
#'                                   education*experience
#'
#' # First variant:
#' # Contribution of union is evaluated using composition of
#' # education and experience from 2003-2005 (group 1)
#'
#' male_inequality_sequential  <- dfl_deco(model_sequential,
#'                                         data = men8305,
#'                                         weights = weights,
#'                                         group = year)
#'
#' # Summarize results
#' summary(male_inequality_sequential)
#'
#' # Second variant:
#' # Contribution of union is evaluated using composition of
#' # education and experience from 1983-1985 (group 0)
#'
#' male_inequality_sequential_2  <- dfl_deco(model_sequential,
#'                                         data = men8305,
#'                                         weights = weights,
#'                                         group = year,
#'                                         right_to_left = FALSE)
#'
#' # Summarize results
#' summary(male_inequality_sequential_2)
#'
#' # The domposition effects associated with (conditional) unionization for deciles
#' cbind(male_inequality_sequential$decomposition_quantiles$prob,
#'       male_inequality_sequential$decomposition_quantiles$`Comp. eff. X1|X2`,
#'       male_inequality_sequential_2$decomposition_quantiles$`Comp. eff. X1|X2`)
#'
#'
#' ## Trim observations with weak common support
#' ## (i.e. observations with relative factor weights > \sqrt(N)/N)
#'
#' set.seed(123)
#' data_weak_common_support <- data.frame(d = factor(c(c("A", "A", rep("B", 98)),
#'                                                   c(rep("A", 90), rep("B", 10)))),
#'                                        group = rep(c(0, 1), each = 100))
#' data_weak_common_support$y <- ifelse(data_weak_common_support$d == "A", 1, 2) +
#'                               data_weak_common_support$group +
#'                               rnorm(200, 0, 0.5)
#'
#' deco_results_trimmed <- dfl_deco(y ~ d,
#'                                  data_weak_common_support,
#'                                  group = group,
#'                                  trimming = TRUE)
#'
#' identical(deco_results_trimmed$trimmed_observations,
#'           which(data_weak_common_support$d == "A"))
#'
dfl_deco <-  function(formula,
                      data,
                      weights,
                      group,
                      na.action = na.exclude,
                      reference_0 = TRUE,
                      subtract_1_from_0 = FALSE,
                      right_to_left = TRUE,
                      method = "logit",
                      estimate_statistics = TRUE,
                      statistics = c("quantiles", "mean", "variance", "gini", "iq_range_p90_p10", "iq_range_p90_p50", "iq_range_p50_p10"),
                      probs = c(1:9)/10,
                      custom_statistic_function = NULL,
                      trimming = FALSE,
                      trimming_threshold = NULL,
                      return_model = TRUE,
                      bootstrap = FALSE,
                      bootstrap_iterations = 100,
                      bootstrap_robust = FALSE,
                      cores = 1,
                      ...){

  ## Get model.frame
  function_call <- match.call()
  data_arguments_index <- match(c("formula", "data", "weights", "group"), names(function_call), 0)
  data_arguments <- function_call[c(1, data_arguments_index)]
  data_arguments$drop.unused.levels <- TRUE

  data_arguments[[1]] <- as.name("model.frame")
  formula <- Formula::as.Formula(formula)
  data_arguments$formula <- formula
  data_used <- eval.parent(data_arguments)
  data_used <- na.action(data_used)
  function_terms <- attr(data_used, "terms")
  dep_var <- model.response(data_used, "numeric")

  ## Get weights
  weights <- model.weights(data_used)
  if (!is.null(weights) && !is.numeric(weights)) {
    stop("'weights' must be a numeric vector")
  }
  if (is.null(weights)) {
    weights <- rep(1, length(dep_var))
  }
  if (sum(weights) <= 1) {
    weights <- weights * length(weights)
  }

  ## Get group variable and reference level
  group_variable_name <- as.character(data_arguments[["group"]])
  remove(group)
  group_variable <- data_used[, ncol(data_used)]
  names(data_used)[ncol(data_used)] <- "group_variable"
  check_group_variable <- is.numeric(group_variable)&length(unique(group_variable))==2|is.factor(group_variable)&length(unique(group_variable))==2
  if(check_group_variable==FALSE){
    stop("Group variable must either be a binary numeric variable or a binary factor variable.")
  }

  if(is.numeric(group_variable)){
    data_used$group_variable <- group_variable <- as.factor(group_variable)
  }
  reference_group <- ifelse(reference_0, 0, 1)
  reference_group_print <- levels(group_variable)[reference_group + 1]

  ## Check if statistics are implemented
  statistics_implemented <-  c("quantiles", "mean", "variance", "gini",
                               "iq_range_p90_p10", "iq_range_p90_p50", "iq_range_p50_p10",
                               "iq_ratio_p90_p10", "iq_ratio_p90_p50", "iq_ratio_p50_p10")
  statistics_not_implemented <- setdiff(statistics, statistics_implemented)
  statistics <- setdiff(statistics,statistics_not_implemented)
  if(length(statistics_not_implemented)>0){
    cat("Warning:", paste0("Selected statistics (", paste0(statistics_not_implemented, collapse=", "),")"), "not implemented! \n")
    cat("Implemented statistics:", paste0(statistics_implemented, collapse=", "), "\n \n")
  }
  if("quantiles" %in% statistics & length(probs)==0){
    probs <- seq(5,95,5)/100
  }
  if(length(statistics)==0){
    estimate_statistics <- FALSE
  }

  if(method %in% c("logit", "random forests") == FALSE){
    stop("Only 'logit' and 'random forests' are available methods to estimate reweighting factors")
  }

  results <- dfl_deco_estimate(formula = formula,
                               dep_var = dep_var,
                               data_used = data_used ,
                               weights = weights,
                               group_variable = group_variable,
                               reference_group = reference_group,
                               method = method,
                               estimate_statistics = estimate_statistics,
                               statistics = statistics,
                               probs = probs,
                               right_to_left = right_to_left,
                               trimming = trimming,
                               trimming_threshold = trimming_threshold,
                               return_model = return_model,
                               ...)


  if(bootstrap){
    cat("Bootstrapping standard errors...\n")
    if(cores == 1) {
      bootstrap_estimates <- pbapply::pblapply(1:bootstrap_iterations,
                                               function(x) dfl_deco_bootstrap(formula = formula,
                                                                              dep_var = dep_var,
                                                                              data_used = data_used,
                                                                              weights = weights,
                                                                              group_variable = group_variable,
                                                                              reference_group = reference_group,
                                                                              method=method,
                                                                              estimate_statistics = estimate_statistics,
                                                                              statistics = statistics,
                                                                              probs = probs,
                                                                              right_to_left = right_to_left,
                                                                              trimming = trimming,
                                                                              trimming_threshold = trimming_threshold,
                                                                              return_model = FALSE,
                                                                              ...))
    } else {
      cores <- min(cores, parallel::detectCores() - 1)
      cluster <- parallel::makeCluster(cores)
      parallel::clusterSetRNGStream(cluster, round(runif(1, 0, 100000)))
      parallel::clusterExport(cl = cluster,
                              varlist = ls(),
                              envir = environment())
      parallel::clusterExport(cl = cluster,
                              c("dfl_deco_bootstrap",
                                "dfl_deco_estimate",
                                "fit_and_predict_probabilities",
                                "select_observations_to_be_trimmed",
                                "get_distributional_statistics",
                                "estimate_iq_range",
                                "estimate_iq_ratio"
                              ))
      bootstrap_estimates <- pbapply::pblapply(1:bootstrap_iterations,
                                               function(x) dfl_deco_bootstrap(formula = formula,
                                                                              dep_var = dep_var,
                                                                              data_used = data_used,
                                                                              weights = weights,
                                                                              group_variable = group_variable,
                                                                              reference_group = reference_group,
                                                                              method = method,
                                                                              estimate_statistics = estimate_statistics,
                                                                              statistics = statistics,
                                                                              probs = probs,
                                                                              right_to_left = right_to_left,
                                                                              trimming = trimming,
                                                                              trimming_threshold = trimming_threshold,
                                                                              return_model = FALSE,
                                                                              ...),
                                               cl = cluster)
      parallel::stopCluster(cluster)
    }

    if("quantiles" %in% statistics){
      bootstrapped_quantiles <- as.data.frame(do.call("rbind", lapply(bootstrap_estimates, function(x) x[["decomposition_quantiles"]])))
      bootstrapped_quantiles$iteration <- rep(1:bootstrap_iterations, each=length(probs))

      bs_se_deco_quantiles <-  stats::reshape(bootstrapped_quantiles,
                                              idvar = c("probs","iteration"),
                                              times = setdiff(names( bootstrapped_quantiles),c("probs","iteration")),
                                              timevar="effect",
                                              varying = list(setdiff(names( bootstrapped_quantiles),c("probs","iteration"))),
                                              direction = "long",
                                              v.names = "value")

      bs_se_deco_quantiles <- lapply(split(bs_se_deco_quantiles, bs_se_deco_quantiles[,c("probs","effect")]),
                                     function(x) data.frame(probs=x$probs[1],
                                                            effect=x$effect[1],
                                                            se=ifelse(bootstrap_robust,
                                                                      (quantile(x$value, 0.75) - quantile(x$value, 0.25))/(qnorm(0.75) - qnorm(0.25)),
                                                                      sqrt(var(x$value))))
      )

      bs_se_deco_quantiles <- do.call("rbind",  bs_se_deco_quantiles)
      bs_se_deco_quantiles <- stats::reshape(bs_se_deco_quantiles,
                                             idvar = c("probs"),
                                             timevar="effect",
                                             direction = "wide")
      names(bs_se_deco_quantiles) <- gsub("se[.]","",names( bs_se_deco_quantiles))


      bs_se_deco_quantiles <- bs_se_deco_quantiles[, names(results$decomposition_quantiles)]

      bs_kolmogorov_smirnov_stat  <-  stats::reshape(bootstrapped_quantiles,
                                                     idvar = c("probs","iteration"),
                                                     times = setdiff(names(bootstrapped_quantiles),c("probs","iteration")),
                                                     timevar="effect",
                                                     varying = list(setdiff(names(bootstrapped_quantiles),c("probs","iteration"))),
                                                     direction = "long",
                                                     v.names = "value")

      bs_kolmogorov_smirnov_stat <- lapply(split(bs_kolmogorov_smirnov_stat, bs_kolmogorov_smirnov_stat[, c("probs","effect")]),
                                           function(x) data.frame(probs=x$probs[1],
                                                                  effect=x$effect[1],
                                                                  iteration=x$iteration,
                                                                  value=x$value,
                                                                  se=ifelse(bootstrap_robust,
                                                                            (quantile(x$value, 0.75) - quantile(x$value, 0.25))/(qnorm(0.75) - qnorm(0.25)),
                                                                            sqrt(var(x$value)))
                                           )
      )
      bs_kolmogorov_smirnov_stat <- do.call("rbind",  bs_kolmogorov_smirnov_stat)
      bs_kolmogorov_smirnov_stat$value_over_se <- bs_kolmogorov_smirnov_stat$value/bs_kolmogorov_smirnov_stat$se
      bs_kolmogorov_smirnov_stat <- lapply(split(bs_kolmogorov_smirnov_stat, bs_kolmogorov_smirnov_stat[, c("effect","iteration")]),
                                           function(x) data.frame(effect=x$effect[1],
                                                                  iteration=x$iteration[1],
                                                                  kms_t_value = max(x$value_over_se)
                                           )
      )


      bs_kolmogorov_smirnov_stat <- do.call("rbind", bs_kolmogorov_smirnov_stat)

    } else {
      bs_se_deco_quantiles <- NULL
    }

    if(length(setdiff(statistics,"quantiles"))>1){
      bootstrapped_other_statistics <- as.data.frame(do.call("rbind", lapply(bootstrap_estimates, function(x) x[["decomposition_other_statistics"]])))
      bootstrapped_other_statistics$iteration <- rep(1:bootstrap_iterations, each=length(unique(bootstrapped_other_statistics$statistic)))

      bs_se_deco_other_statistics <- stats::reshape(bootstrapped_other_statistics,
                                                    idvar = c("statistic","iteration"),
                                                    ids=unique(bootstrapped_other_statistics$statistic),
                                                    times = setdiff(names(bootstrapped_other_statistics),c("statistic","iteration")),
                                                    timevar="effect",
                                                    varying = list(setdiff(names(bootstrapped_other_statistics),c("statistic","iteration"))),
                                                    direction = "long",
                                                    v.names = "value")

      bs_se_deco_other_statistics <- lapply(split(bs_se_deco_other_statistics, bs_se_deco_other_statistics[,c("statistic","effect")]),
                                            function(x) data.frame(statistic=x$statistic[1],
                                                                   effect=x$effect[1],
                                                                   se=ifelse(bootstrap_robust,
                                                                             (quantile(x$value, 0.75) - quantile(x$value, 0.25))/(qnorm(0.75) - qnorm(0.25)),
                                                                             sqrt(var(x$value))))
      )

      bs_se_deco_other_statistics <- do.call("rbind",  bs_se_deco_other_statistics)

      bs_se_deco_other_statistics <- stats::reshape(bs_se_deco_other_statistics,
                                                    idvar = c("statistic"),
                                                    timevar="effect",
                                                    direction = "wide")
      names(bs_se_deco_other_statistics) <- gsub("se[.]","",names( bs_se_deco_other_statistics ))
      bs_se_deco_other_statistics <- as.data.frame(bs_se_deco_other_statistics[match(results$decomposition_other_statistics$statistic,
                                                                                     bs_se_deco_other_statistics$statistic),
                                                                               names(results$decomposition_other_statistics)])
    }else{
      bs_se_deco_other_statistics <- NULL
    }

    bootstrapped_quantiles_reweighting_factor <- as.data.frame(do.call("rbind", lapply(bootstrap_estimates, function(x) x[["quantiles_reweighting_factor"]])))
    bs_se_quantiles_reweighting_factor <- tidyr::pivot_longer(bootstrapped_quantiles_reweighting_factor,
                                                              col=names(bootstrapped_quantiles_reweighting_factor)[-1],
                                                              names_to="effect")
    bs_se_quantiles_reweighting_factor <- dplyr::summarise(dplyr::group_by(bs_se_quantiles_reweighting_factor,
                                                                           probs,
                                                                           effect),
                                                           se=ifelse(bootstrap_robust,
                                                                     (quantile(value, 0.75) - quantile(value, 0.25))/(qnorm(0.75) - qnorm(0.25)),
                                                                     sqrt(var(value))),
                                                           .groups = "keep")
    bs_se_quantiles_reweighting_factor <- tidyr::pivot_wider(bs_se_quantiles_reweighting_factor,
                                                             id_cols = "probs",
                                                             names_from = "effect",
                                                             values_from = "se")
    bs_se_quantiles_reweighting_factor <- as.data.frame(bs_se_quantiles_reweighting_factor[, names(results$quantiles_reweighting_factor)])
    bs_se_quantiles_reweighting_factor[which(bs_se_quantiles_reweighting_factor$probs %in% c(0,1)),
                                       2:ncol(bs_se_quantiles_reweighting_factor)] <- NA

    bootstrap_se <- list(decomposition_quantiles=bs_se_deco_quantiles,
                         decomposition_other_statistics=bs_se_deco_other_statistics,
                         quantiles_reweighting_factor=bs_se_quantiles_reweighting_factor,
                         decomposition_quantiles_kms_distribution=bs_kolmogorov_smirnov_stat)
  }
  else {
    bootstrap_se <- NULL

  }

  if(subtract_1_from_0 & estimate_statistics){
    if(!is.null(results$decomposition_quantiles)) results$decomposition_quantiles[,-1] <- results$decomposition_quantiles[,-1] * -1
    if(!is.null(results$decomposition_other_statistics)) results$decomposition_other_statistics[,-1] <- results$decomposition_other_statistics[,-1] * -1
  }

  add_to_results <- list(bootstrapped_standard_errors = bootstrap_se,
                         group_variable_name = group_variable_name,
                         group_variable_levels = levels(group_variable),
                         reference_group = reference_group_print,
                         subtract_1_from_0 = subtract_1_from_0)
  results <- c(results, add_to_results)

  class(results) <- "dfl_deco"
  return(results)

}


#' Estimate DFL reweighting decomposition
#'
#' This function performs the actual DFL decomposition. It derives the
#' reweighting factors, estimates the distributional statistics and
#' calculates the decomposition terms.
#'
dfl_deco_estimate <- function(formula,
                              dep_var,
                              data_used ,
                              weights,
                              group_variable,
                              reference_group,
                              method,
                              estimate_statistics,
                              statistics,
                              probs,
                              right_to_left,
                              trimming,
                              trimming_threshold,
                              return_model,
                              ...){

  # Estimate probabilities -------------------------------------------------------
  mod <- group_variable ~ 1
  p1 <- mean(fit_and_predict_probabilities(mod, data_used, weights, method = "logit", return_model = FALSE)[[1]])
  p0 <- 1-p1
  estimated_probabilities <- rep(p0/p1, nrow(data_used))

  formula <- Formula::as.Formula(formula)
  nvar <- length(formula)[2] # Number of detailed decomposition effects
  covariates_labels <- fitted_models <- vector("list", nvar)

  for(i in nvar:1){
    mod <- update(stats::formula(formula, rhs=nvar:i, collapse=TRUE), group_variable ~ .)
    fitted_model <- fit_and_predict_probabilities(mod,
                                                  data_used,
                                                  weights,
                                                  method = method,
                                                  return_model = return_model,
                                                  ...)

    p1 <- fitted_model[[1]]
    p0 <- 1 - p1
    estimated_probabilities <- cbind(estimated_probabilities, p0/p1)

    covariates_labels[[i]] <- attr(terms(mod), "term.labels")[which(attr(terms(mod), "order") == 1)]
    fitted_models[i] <- fitted_model[2]
  }

  names(fitted_models) <- sapply(nvar:1, function(i) paste0("P(g=1|", paste0(paste0("X", nvar:i), collapse = ","), ")"))


  # Collect covariates' labels ---------------------------------------------------

  covariates_labels <- rev(covariates_labels)
  all_covariates <- paste0(unique(do.call("c", covariates_labels)), collapse = ", ")
  if(nvar > 1){
    covariate_index <- nvar:1
    for(i in nvar:2){
      add_vars <- setdiff(covariates_labels[[i]], covariates_labels[[i-1]])
      add_index <- paste0(paste0("X", covariate_index[i]), "|", paste0(paste0("X",covariate_index[(i-1):1]), collapse = ","))
      covariates_labels[[i]] <- paste0("Detailed effect ",
                                       add_index,
                                       ": ",
                                       paste0(add_vars, collapse = ", "),
                                       " | ",
                                       paste0(covariates_labels[[i-1]], collapse = ", "))
    }
    covariates_labels[[1]] <- paste0("Detailed effect X",
                                     nvar,
                                     ": ",
                                     paste0(covariates_labels[[1]], collapse = ", "))

    covariates_labels <- c(paste0(c("Aggregate effect: ", all_covariates), collapse =""), rev(covariates_labels))
  }else{
    covariates_labels[[1]] <- all_covariates
  }


  # Derive reweighting factors ---------------------------------------------------

  # e.g., in the case of Y ~ X1 + X2 + X3 | X2 + X3 | X3
  # Matrix probs contains nvar+1 columns:
  # first column  [P(g=0)/P(g=1)]
  # second column [P(g=0|X3)/P(g=1|X3)]
  # third column  [P(g=0|X2,X3)/P(g=1|X2,X3)]
  # fourth column [P(g=0|X1,X2,X3)/P(g=1|X1,X2,X3)]

  psi <- NULL

  if(right_to_left){
    # e.g., in the case of Y ~ X1 + X2 + X3 | X2 + X3 | X3
    # if right_to_left==TRUE,
    # then the matrix psi has the following columns:
    # first column:  [P(g=1)/P(g=0)]*[P(g=0|X3)/P(g=1|X3)]
    # second column: [P(g=1)/P(g=0)]*[P(g=0|X2,X3)/P(g=1|X2,X3)]
    # third column:  [P(g=1)/P(g=0)]*[P(g=0|X1,X2,X3)/P(g=1|X1,X2,X3)]
    for(i in 1:nvar){
      psi <- cbind(psi, (estimated_probabilities[, 1]^-1) * estimated_probabilities[, i+1])
    }
    psi <- as.data.frame(psi)
    names(psi) <- paste0("Psi_",
                         sapply(nvar:1, function(i) paste0(paste0("X", i:nvar), collapse=",")))

    #names_decomposition_terms <- nvar:1
    names_decomposition_terms <- paste0("X", nvar)
    if(nvar>1){
      names_decomposition_terms <- c(names_decomposition_terms,
                                     sapply((nvar-1):1,
                                            function(i) paste0(paste0(paste0("X", i), collapse=","), "|", paste0(paste0("X", (i+1):nvar), collapse=","))))

    }
  }else{
    # if right_to_left==FALSE,
    # then the matrix psi has the following columns:
    # first column  [P(g=1|X2,X3)/P(g=0|X2,X3)]*[P(g=0|X1,X2,X3)/P(g=1|X1,X2,X3)]
    # second column [P(g=1|X3)/P(g=0|X3)]*[P(g=0|X1,X2,X3)/P(g=1|X1,X2,X3)]
    # third column  [P(g=1)/P(g=0)]*[P(g=0|X1,X2,X3)/P(g=1|X1,X2,X3)]
    for(i in nvar:1){
      psi <- cbind(psi, (estimated_probabilities[, i]^-1) * estimated_probabilities[, nvar+1])
    }
    psi <- as.data.frame(psi)
    names(psi)[nvar] <- paste0("Psi_", paste0(paste0("X", 1:nvar), collapse=","))

    names_decomposition_terms <- paste0("X", nvar)
    if(nvar>1){
      names(psi)[1:(nvar-1)] <- paste0("Psi_", sapply(1:(nvar-1),
                                                      function(i) paste0(paste0(paste0("X", 1:i), collapse=","), "|", paste0(paste0("X", (i+1):nvar), collapse=","))))
      names_decomposition_terms <- c(sapply(1:(nvar-1),
                                            function(i) paste0(paste0(paste0("X", i), collapse=","), "|", paste0(paste0("X", (i+1):nvar), collapse=","))),
                                     names_decomposition_terms)
    }
    #names_decomposition_terms <- 1:nvar


  }

  # Take inverse of reweighting factor if reference group is 0 -----------------

  psi_power <- ifelse(reference_group == 1, 1, -1)
  psi <- psi^psi_power

  # Trimming: Set weights of reweighting_factors above trimming threshold to zero

  if(trimming){

    observations_to_be_trimmed <- list()

    for(j in 1:ncol(psi)){
      observations_to_be_trimmed[[j]] <- select_observations_to_be_trimmed(reweighting_factor = psi[, j],
                                                                           group_variable = group_variable,
                                                                           group = reference_group,
                                                                           trimming_threshold = trimming_threshold)
    }

    observations_to_be_trimmed <- unique(do.call("c", observations_to_be_trimmed))
    weights[observations_to_be_trimmed] <- 0

  }else{

    observations_to_be_trimmed <- NULL

  }


  # Estimate distributional statistics and perform decomposition -----------------
  if(estimate_statistics){

    log_transformed <- grepl(pattern = "log[(]", strsplit(as.character(formula), split="~")[[2]])
    nu1 <- get_distributional_statistics(dep_var,
                                         weights,
                                         group_variable,
                                         group = 1,
                                         statistics = statistics,
                                         probs = probs,
                                         log_transformed = log_transformed)

    nu0 <- get_distributional_statistics(dep_var,
                                         weights,
                                         group_variable,
                                         group = 0,
                                         statistics = statistics,
                                         probs = probs,
                                         log_transformed = log_transformed)

    #if reference group==0, take inverse of rw factors

    nuC <- NULL

    for(i in 1:nvar){
      nuC <- cbind(nuC,
                   get_distributional_statistics(dep_var,
                                                 weights*psi[,i],
                                                 group_variable,
                                                 group=reference_group,
                                                 statistics=statistics,
                                                 probs=probs,
                                                 log_transformed=log_transformed))
    }
    nuC <- as.matrix(nuC)



    # Aggregate decomposition ----------------------------------------------------
    Delta <- cbind(nu1 - nu0,
                   nu1 - nuC[, nvar],
                   nuC[, nvar] - nu0)
    if(reference_group==1){
      colnames(Delta) <- c("Observed difference", "Composition effect", "Structure effect")
    } else {
      colnames(Delta) <- c("Observed difference", "Structure effect", "Composition effect")
    }

    Delta <- Delta[, match(c("Observed difference", "Composition effect", "Structure effect"), colnames(Delta))]

    # Detailed decomposition -----------------------------------------------------
    if(nvar>1){
      if(reference_group==1){
        Delta <- cbind(Delta,
                       nu1-nuC[, 1])
        colnames(Delta)[length(colnames(Delta))] <- paste("Comp. eff. ", names_decomposition_terms[1], sep="")
        for(i in 2:nvar){
          Delta <- cbind(Delta, nuC[,i-1]-nuC[,i])
          colnames(Delta)[length(colnames(Delta))] <- paste("Comp. eff. ", names_decomposition_terms[i], sep="")
        }
      }else{
        for(i in nvar:2){
          Delta <- cbind(Delta, nuC[,i]-nuC[,i-1])
          colnames(Delta)[length(colnames(Delta))] <- paste("Comp. eff. ", names_decomposition_terms[i], sep="")
        }
        Delta <- cbind(Delta, nuC[,1]-nu0)
        colnames(Delta)[length(colnames(Delta))] <- paste("Comp. eff. ", names_decomposition_terms[1], sep="")
      }

      Delta <- Delta[, c(1:3,
                         order(colnames(Delta)[-c(1:3)])+3)]
    }

    Delta <- as.data.frame(Delta)

    # Save quantiles and other statistics in different objects -------------------

    if("quantiles" %in% statistics){
      decomposition_quantiles <- Delta[1:length(probs), ]
      cn <- names(decomposition_quantiles)
      decomposition_quantiles$probs <- probs
      decomposition_quantiles <- decomposition_quantiles[,c("probs", cn)]
    }else{
      decomposition_quantiles <- NULL
    }
    if("quantiles" %in% statistics & length(statistics) > 1){
      decomposition_other_statistics <- Delta[(length(probs)+1):nrow(Delta), ]
    }else if("quantiles" %in% statistics == FALSE){
      decomposition_other_statistics <- Delta
    }else{
      decomposition_other_statistics <- NULL
    }

    if(is.null(decomposition_other_statistics)==FALSE){
      cn <- names(decomposition_other_statistics)
      decomposition_other_statistics$statistic <- rownames(decomposition_other_statistics)
      decomposition_other_statistics <- decomposition_other_statistics[,c("statistic", cn)]
    }


  }else{
    decomposition_quantiles <- NULL
    decomposition_other_statistics <- NULL
    Delta <- NULL
  }

  # Compute sample quantiles of reweighting factors ----------------------------
  quantiles_reweighting_factor <- data.frame(probs=c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1))
  rownames(quantiles_reweighting_factor) <- c("Min.", "10%-quantile", "25%-quantile", "50%-quantile", "75%-quantile", "90%-quantile", "Max.")
  factors_to_be_considered <- setdiff(1:nrow(psi), observations_to_be_trimmed)
  factors_to_be_considered <- intersect(factors_to_be_considered,
                                        which(group_variable == levels(group_variable)[reference_group + 1]))

  for(i in 1:ncol(psi)){
    quantiles_reweighting_factor <- cbind(quantiles_reweighting_factor,
                                          quantile(psi[factors_to_be_considered, i], probs=quantiles_reweighting_factor$probs, na.rm=TRUE))
    names(quantiles_reweighting_factor)[i+1] <- names(psi)[i]
  }

  # Export results
  results <- list(decomposition_quantiles = decomposition_quantiles,
                  decomposition_other_statistics = decomposition_other_statistics,
                  reweighting_factor = psi,
                  quantiles_reweighting_factor = quantiles_reweighting_factor,
                  trimmed_observations = observations_to_be_trimmed,
                  covariates_labels = covariates_labels,
                  fitted_models = fitted_models)
  return(results)
}


#' Bootstrap function
dfl_deco_bootstrap <- function(formula,
                               dep_var,
                               data_used ,
                               weights,
                               group_variable,
                               reference_group,
                               estimate_statistics,
                               statistics,
                               probs,
                               right_to_left,
                               trimming,
                               trimming_threshold,
                               ...){

  sampled_observations <- sample(1:nrow(data_used),
                                 nrow(data_used),
                                 replace=TRUE,
                                 prob=weights/sum(weights,na.rm=TRUE))

  deco_estimates <- dfl_deco_estimate(formula = formula,
                                      dep_var = dep_var[sampled_observations],
                                      data_used = data_used[sampled_observations,],
                                      weights = (weights[sampled_observations]/sum(weights[sampled_observations],na.rm=TRUE))*sum(weights,na.rm=TRUE),
                                      group_variable=group_variable[sampled_observations],
                                      reference_group=reference_group,
                                      estimate_statistics=estimate_statistics,
                                      statistics = statistics,
                                      probs = probs,
                                      right_to_left = right_to_left,
                                      trimming = trimming,
                                      trimming_threshold = trimming_threshold,
                                      ...)

  deco_estimates$reweighting_factor <- NULL

  return(deco_estimates)
}


#' Predict conditional probabilities
#'
#' This function fits a binary choice model and predicts probabilities for every
#' observations.
#'
#' @param formula \code{formula} object specifying the conditional probability model.
#' @param data_used \code{data.frame} with data.
#' @param method method to estimate conditional probabilities
#' @param  boolean: If \code{TRUE} (default), the object(s) of the model
#' fit(s) used to predict the conditional probabilities for the reweighting factor(s)
#' are returned.
#' @param newdata \code{data.frame} with data to be used for predictions.
#' @param ... other parameters passed to estimation function.
#'
fit_and_predict_probabilities <- function(formula,
                                          data_used,
                                          weights,
                                          method = "logit",
                                          return_model = FALSE,
                                          newdata = NULL,
                                          ...){

  data_used$`(weights)` <- weights

  if(method=="logit"){

    model_fit <- glm(formula,
                     data = data_used,
                     weights = `(weights)`,
                     family = quasibinomial(link = "logit"),
                     ...)

    fitted_probabilities  <- predict.glm(model_fit,
                                         newdata = newdata,
                                         type="response",
                                         na.action = na.exclude)
  }

  if(method=="fastglm"){

    x <- model.matrix(formula, data_used)
    y <- model.response(model.frame(formula, data_used))
    if(is.factor(y)){
      y <- as.numeric(y == max(levels(y)))
    }else if(is.numeric(y)){
      y <- as.numeric(y == max(y))
    }
    model_fit <- fastglm::fastglm(x = x,
                                  y = y,
                                  quasibinomial(link = "logit"),
                                  weights = weights,
                                  ...
    )

    if(is.null(newdata)){
      newdata <- data_used
    }
    newdata_matrix <- model.matrix(formula, newdata)

    fitted_probabilities  <- predict(model_fit,
                                      newdata = newdata_matrix,
                                      type="response")


  }


  ### Random forestes predictions with ranger
  if(method=="random forests"){

    model_fit <- ranger::ranger(formula,
                                data = data_used,
                                case.weights = data_used$`(weights)`,
                                classification = TRUE,
                                probability = TRUE,
                                ...)

    if(is.null(newdata)){
      newdata <- data_used
    }

    fitted_probabilities  <-  predict(model_fit,
                                      data = newdata)$predictions[,2]
  }


  if(as.logical(return_model) != TRUE){
    model_fit <- NULL
  }

  results <- list(fitted_probabilities,
                  model_fit)

  return(results)
}


#' Select observations with little common support to be trimmed
#'
#' This function implements the trimming rule proposed by Huber, Lechner,
#' and Wunsch (2014). Observations above the trimming threshold are trimmed in
#' the reference group and in the comparison group. Per default, the timming
#' is set to sqrt(N)/N, where N is the number of observation in the  reweighted
#' reference group. The function returns vector index of observation to be trimmed.
#'
#' @param reweighting_factor Estimated reweigting factor
#' @param group_variable Variable identifying the reference and comparison group, respectively.
#' @param group Identifier of reference group
#' @param trimming_threshold threshold defining the maximal accepted relative weight of a reweighting factor/observation. If `NULL`, the threshold is set to `sqrt(N)/N`, where N is the number of observations in the reference group.
#'

select_observations_to_be_trimmed <- function(reweighting_factor,
                                              group_variable,
                                              group,
                                              trimming_threshold = NULL){

  group <- levels(group_variable)[group + 1]

  reweighting_factor_control <- reweighting_factor[which(group_variable==group)]
  if(is.null(trimming_threshold)){
    nobs <- length(reweighting_factor_control)
    trimming_threshold <- sqrt(nobs) / nobs
  }
  reweighting_factor_in_trimming_range <- reweighting_factor_control[which(reweighting_factor_control / sum(reweighting_factor_control) > trimming_threshold)]
  trim_value <- ifelse(length(reweighting_factor_in_trimming_range) == 0,
                       sum(reweighting_factor),
                       min(reweighting_factor_in_trimming_range))

  observations_to_be_trimmed <- which(reweighting_factor >= trim_value)

  return(observations_to_be_trimmed)
}


