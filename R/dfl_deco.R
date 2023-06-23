#' DFL reweighting decomposition
#'
#' @description `dfl_deco` decomposes between-group differences in distributional
#' statistics of an outcome variable into a structure effect and a composition effect.
#' As proposed by DiNardo, Fortin, and Lemieux (1996), the procedure
#' reweights the sample distribution of a reference group such that the group's
#' covariates distribution matches the covariates distribution of a counterfactual
#' group.
#'
#' Reweighting factors are derived by modelling the probability of belonging to
#' the one group instead of the other conditional on covariates. The function
#' allows detailed decompositions of the composition effect by sequentially
#' reweighting (conditional) covariate distributions. Standard errors can be
#' bootstrapped.
#'
#' @param formula a `formula` object with an outcome variable Y on the left-hand side
#' and the covariates X on the right-hand side. For sequential decompositions,
#' the sequence of covariates X are distinguished by the `|` operator. Covariates
#' are used to estimate the reweighting factors.
#' @param data a `data.frame` containing all variables and observations for
#' both groups.
#' @param weights name of the observation weights variable or vector of
#' observation weights.
#' @param na.action a function to filter missing data (default `na.exclude`).
#' @param group name of the a binary variable (numeric or factor)
#' identifying the two groups for which the differences are to be decomposed.
#' The group identified by the lower ranked value in `group` (i.e., 0 in the
#' case of a dummy variable or the first level of factor variable) is defined
#' as group 0.
#' @param reference_0 boolean: if `TRUE` (default), then the group 0 -- i.e.,
#' the group identified by the lower ranked value in `group` -- will be defined
#' as reference group. The reference group will be reweighted to match the
#' covariates distribution of the sample of the comparison group.
#' @param reweight_marginals If `TRUE` (default), then the sequential decomposition
#' reweights first the marginal (joint) distribution of the last covariate (covariates)
#' entered into `formula` sequence. Otherwise, the conditional distribution of
#' first covariate(s) entered are reweigthed.
#' @param method specifies the method fit and predict conditional probabilities
#' used to derive the reweighting factor. At the moment, `logit` is the only method
#' available.
#' @param estimate_statistics boolean: if `TRUE` (default), then distributional
#' statistics are estimated and the decomposition is performed. If `FALSE`,
#' the function only returns the propensity weights.
#' @param statistics a character vector that defines the distributional statistics
#' for which the decomposition is perforemd. Per default, `c("quantiles", "mean", "variance", "gini", "iq_range_p90_p10", "iq_range_p90_p50", "iq_range_p50_p10")`
#' are estimated and decomposed. Also implemented are `c("iq_ratio_p90_p10", "iq_ratio_p90_p50", "iq_ratio_p50_p10")`.
#' Note: The function calculates the Gini coefficient for untransformed variable
#' (i.e., exp(log(Y))), if the logarithm of a variable Y is set as outcome variable
#' in `formula`.
#' @param probs a vector of length 1 or more with the probabilities of the quantiles
#' to be estimated with default `c(1:9)/10`.
#' @param custom_statistic_function a custom statistic function to pass as argument.
#' @param trimming boolean: If `TRUE`, observations with dominant reweighting factor values are trimmend according to rule of Huber, Lechner, and Wunsch (2014). Per default, trimming is set to `FALSE`.
#' @param trimming_threshold numeric: threshold defining the maximal accepted relative weight of the reweighting factor value of a single observation. If `NULL`, the threshold is set to `sqrt(N)/N`, where N is the number of observations in the reference group.
#' @param bootstrap boolean: If `FALSE`, then the estimation is not boostrapped and no
#' standard errors are calculated.
#' @param bootstrap_iterations positive integer indicating the number of bootstrap
#'  iterations to execute. Only required if \code{bootstrap = TRUE}.
#' @param bootstrap_robust boolean: if `FALSE` (default), then bootstrapped standard
#' errores are estimated as the standard deviations of the bootstrapp estimates.
#' Otherwise, the function uses the bootstrap interquartile range rescaled by the
#' interquantile range of the standard distribution to estimate standard errors.
#' @param cores positive integer indicating the number of cores to use when
#' computing bootstrap standard errors. Only required if \code{bootstrap = TRUE}.
#' @param ... other parameters passed to the function estimating the conditional probabilities.
#'
#' @details
#' The covariates entered into `formula` will be used to estimate the reweighting
#' factors, i.e. the propensities of belonging to one of the groups. `formula` also
#' allows to specify interaction terms in the conditional probability models.
#'
#' If you are interest in an aggregate decomposition and no sequential decomposition
#' should be performend, the all covariates have to be entered at once, e.g.
#' `Y ~ X1 + X2 + X3`.
#'
#' If you are interested in a sequential decomposition, the the decomposition sequence
#' has to be distinguished by the `|` operator. For instance, `Y ~ X1 | X2 + X3`
#' would decomposes the aggregate composition effect into the contribution of X1
#' and of the remaining two covariates, respectively. The function sequentially
#' collapses the multiple parts of the `formula` object. For instance, if
#' `reweight_marginals=TRUE` (see below), then it estimtates in a first step
#' the reweighting factor based on `X2 + X3` and then in a second step the one
#' based on `X1 + X2 + X3`.
#'
#' You can also specify the detailed models in every part of `formula`. This is
#' useful if you want to estimate in every step a fully saturated model, e.g.
#' `Y ~ X1*X2*X3 | X2*X3`.
#'
#' The observed difference to be decomposed is equals the value of the
#' distributioal statistic of `group` 0 substracted from the value of the same
#' statistic of `group` 1.
#'
#' If `reference_0=TRUE`, then group 0 is the reference group and its observations
#' are reweighted such that they match the covariates distribution of group 1, the
#' counterfacutal group. Group 0 is identified by the lower ranked value of the
#' `group` variable. The composition effect is evaluated using the structure of
#' the reference group. The wage structure effect will be evaluated using the
#' covariates' distribution of the counterfactual group.
#'
#' If `reweight_marginals=TRUE`, the sequential decomposition is performed by
#' reweighting the marginal distribution of the last covariate or the joint
#' distribution of the last covariates in the sequence. For instance, if we have
#' `Y ~ X1 | X2 + X3`, then the reference group would be first reweighted such
#' that its distribution of covariates X2 and X3 matches the one of the
#' counterfactual group. Thus, we use the marginal distribution of X2 and X3 of
#' the counterfactual group but the conditional distribution of X1 given X2
#' and X3 as well as the (wage) structure of the reference group to evaluated
#' the first composition effect.
#'
#' If `reweight_marginals=FALSE`, the sequential decomposition is performed by
#' reweighting the conditional distribution of the first covariate or the joint
#' distribution of the first covariates in the sequence. For instance, if we have
#' `Y ~ X1 | X2 + X3`, then the reference group would be first reweighted such
#' that its conditional distribution of X1 given covariates X2 and X3 matches
#' the one of the counterfacutal group. Thus, we use the conditional distribution
#' of X1 given X2 and X3 of the counterfactual group but the joint distribution
#' of X2 and X3 as well as the (wage) structure of the reference group to evaluate
#' the first composition effect.
#'
#' `method="logit"` uses a logit model to estimate the conditional probabilities
#' used to derive the reweighting factors. You can specify interaction terms
#' within the `formula` object.
#'
#' Per default, the function decomposes the between-group differences for
#' quantiles, the mean, the variance, the Gini coefficient and the interquantile
#' range between the 9th and the 1st decile, the 9th decile and the median, and
#' the median and the first decile, respectively. Implemented are also the
#' interquantile ratios for between the same quantiles.
#'
#' The quantiles can be specified by `probs` that sets the corresponding
#' probabilities of the quantiles of interest. For other distributional statistics,
#' please use `custom_statistic_function`.
#'
#' @return an object of class `dfl_deco` containing a data.frame with the
#' decomposition results for the quantiles and for the other distributional
#' statistics, respectively, a data.frame with the estimated reweighting factor
#' for every observation, a data.frame with sample quantiles of the reweighting
#' factors and a list with standard errors for the decomposition terms, the quantiles
#' of the reweighting factor as well as the bootstrapped Kolmogorov-Smirnov
#' distribution to construct uniform confidence bands for quantiles.
#'
#' @references
#' DiNardo, John, Nicole M. Fortin, and Thomas Lemieux. 1996. "Labor Market
#' Institutions and the Distribution of Wages, 1973-1992: A Semiparametric Approach."
#' \emph{Econometrica}, 64(5), 1001-1044.
#'
#' Firpo, Sergio P., Nicole M. Fortin, and Thomas Lemieux. 2018. “Decomposing Wage
#' Distributions Using Recentered Influence Function Regressions.”
#' \emph{Econometrics} 6(2), 28.
#'
#' Fortin, Nicole M., Thomas Lemieux, and Sergio Firpo. 2011. "Decomposition methods in economics."
#' In Orley Ashenfelter and David Card, eds., \emph{Handbook of labor economics}. Vol. 4. Elsevier, 1-102.
#'
#' Firpo, Sergio P., and Cristine Pinto. 2016. "Identification and Estimation of
#' Distributional Impacts of Interventions Using Changes in Inequality Measures."
#' \emph{Journal of Applied Econometrics}, 31(3), 457– 486.
#'
#' @export
#'
#' @examples
#' ## Replicate example from handbook chapter of Fortin, Lemieux, and Firpo (FLF, 2011)
#' data("men8305")
#' men8305$weights <- men8305$weights/sum(men8305$weights) * length(men8305$weights)
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
#' # Plot decomposition of qunatile differences
#' plot(flf_male_inequality)
#'
#' # Use alternative reference group (i.e. reweight sample from 2003-05)
#' flf_male_inequality_reference_0305  <- dfl_deco(flf_model,
#'                                                 data = men8305,
#'                                                 weights = weights,
#'                                                 group = year,
#'                                                 reference_0 = FALSE)
#'
#' # Summarize results
#' summary(flf_male_inequality_reference_0305)
#'
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
#' # Trim observations
#' flf_male_inequality_trimming  <- dfl_deco(flf_model,
#'                               data = men8305,
#'                               weights = weights,
#'                               group = year,
#'                               reference_0 = TRUE,
#'                               trimming = TRUE,
#'                               trimming_threshold = 0.00005)
#' summary(flf_male_inequality_trimming)
#'
#'
#' ## Sequential decomposition
#' # Distinguishing the contribution of the marginal distribution of education
#' # and experience from the contribution of unionization conditional on
#' # education and experience
#'
#' model_sequential <- log(wage) ~ union*(education + experience) + education*experience | education*experience
#'
#' male_inequality_sequential  <- dfl_deco(model_sequential,
#'                                         data = men8305,
#'                                         weights = weights,
#'                                         group = year)
#'
#' # Summarize results
#' summary(male_inequality_sequential)
#'
dfl_deco <-  function(formula,
                      data,
                      weights,
                      group,
                      na.action = na.exclude,
                      reference_0 = TRUE,
                      reweight_marginals = TRUE,
                      method = "logit",
                      estimate_statistics = TRUE,
                      statistics = c("quantiles", "mean", "variance", "gini", "iq_range_p90_p10", "iq_range_p90_p50", "iq_range_p50_p10"),
                      probs = c(1:9)/10,
                      custom_statistic_function = NULL,
                      trimming = FALSE,
                      trimming_threshold = NULL,
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
  cat("Reweighted reference group:",  paste0(group_variable_name, " == '", reference_group_print ,"'"), "\n \n")

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

  results <- dfl_deco_estimate(formula = formula,
                               dep_var = dep_var,
                               data_used = data_used ,
                               weights = weights,
                               group_variable = group_variable,
                               reference_group = reference_group,
                               method=method,
                               estimate_statistics = estimate_statistics,
                               statistics = statistics,
                               probs = probs,
                               reweight_marginals = reweight_marginals,
                               trimming = trimming,
                               trimming_threshold = trimming_threshold,
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
                                                                              reweight_marginals = reweight_marginals,
                                                                              trimming = trimming,
                                                                              trimming_threshold = trimming_threshold,
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
                                "fit_and_predict_probabilities",
                                "dfl_deco_estimate",
                                "select_observations_to_be_trimmed",
                                "get_distributional_statistics",
                                "estimate_iq_range",
                                "estimate_iq_ratio"
                                ))
      #parallel::clusterEvalQ(cl = cluster, library("ddeco"))
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
                                                                              reweight_marginals = reweight_marginals,
                                                                              trimming = trimming,
                                                                              trimming_threshold = trimming_threshold,
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

  add_to_results <- list(bootstrapped_standard_errors = bootstrap_se,
                         group_variable_name = group_variable_name,
                         group_variable_levels = levels(group_variable),
                         reference_group = reference_group_print)
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
                              reweight_marginals,
                              trimming,
                              trimming_threshold,
                              ...){

# Estimate probabilities -------------------------------------------------------

  mod <- group_variable ~ 1
  p1 <- mean(fit_and_predict_probabilities(mod, data_used, weights, method = "logit"))
  p0 <- 1-p1
  estimated_probabilities <- rep(p0/p1, nrow(data_used))

  nvar <- length(formula)[2] # Number of detailed decomposition effects
  covariates_labels <- vector("list", nvar)
  for(i in nvar:1){
    mod <- update(stats::formula(formula, rhs=nvar:i, collapse=TRUE), group_variable ~ .)
    p1 <- fit_and_predict_probabilities(mod, data_used, weights, method = method, ...)
    p0 <- 1 - p1
    estimated_probabilities <- cbind(estimated_probabilities, p0/p1)
    covariates_labels[[i]] <- attr(terms(mod), "term.labels")[which(attr(terms(mod), "order") == 1)]
  }


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

  if(reweight_marginals){
    # e.g., in the case of Y ~ X1 + X2 + X3 | X2 + X3 | X3
    # if reweight_marginals==TRUE,
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
    # if reweight_marginals==FALSE,
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

# Trimming: Set weights of reweighting_factors above trimming threshold to zero

  if(trimming){

   observations_to_be_trimmed <- list()

   for(j in 1:ncol(psi)){
     observations_to_be_trimmed[[j]] <- select_observations_to_be_trimmed(psi[, j],
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
                      probs=probs,
                      log_transformed = log_transformed)

    nu0 <- get_distributional_statistics(dep_var,
                      weights,
                      group_variable,
                      group = 0,
                      statistics = statistics,
                      probs = probs,
                      log_transformed = log_transformed)

    #if reference group==0, take inverse of rw factors
    psi_power <- ifelse(reference_group == 1, 1, -1)
    nuC <- NULL

    for(i in 1:nvar){
      psi[,i] <- psi[,i]^psi_power
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
                  covariates_labels = covariates_labels)
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
                               reweight_marginals,
                               trimming,
                               trimming_threshold,
                               ...){
  sampled_observations <- sample(1:nrow(data_used),
                                 nrow(data_used),
                                 replace=TRUE,
                                 prob=weights/sum(weights,na.rm=TRUE))
  deco_estimates <- dfl_deco_estimate(formula=formula,
                                          dep_var = dep_var[sampled_observations],
                                          data_used = data_used[sampled_observations,],
                                          weights = (weights[sampled_observations]/sum(weights[sampled_observations],na.rm=TRUE))*sum(weights,na.rm=TRUE),
                                          group_variable=group_variable[sampled_observations],
                                          reference_group=reference_group,
                                          estimate_statistics=estimate_statistics,
                                          statistics = statistics,
                                          probs = probs,
                                          reweight_marginals = reweight_marginals,
                                          trimming,
                                          trimming_threshold,
                                          ...)

  deco_estimates$reweighting_factor <- NULL
  return(deco_estimates)
}


#' Predict conditional probabilities
#'
#' This function fits a binary choice model and predicts probabilities for every
#' observations.
fit_and_predict_probabilities <- function(formula,
                                          data_used,
                                          weights,
                                          method="logit",
                                          newdata=NULL,
                                          ...){

  if(method=="logit"){

  data_used$`(weights)` <- weights
  model_fit <- glm(formula,
                   data = data_used,
                   weights = `(weights)`,
                   family = quasibinomial(link = "logit"),
                   ...)

  p_X_1  <- predict.glm(model_fit,
                        newdata=newdata,
                        type="response",
                        na.action = na.exclude)
  }

  ### Include here prediction with ranger::ranger!

  return(p_X_1)
}


#' Select observations with little common support to be trimmed
#'
#' This function implements the trimming rule proposed by Huber, Lechner,
#' and Wunsch (2014). Observations above the trimming threshold are trimmed in
#' the reference group and in the comparison group. Per default, the timming
#' is set to sqrt(N)/N, where N is the number of observation in the  reweighted
#' reference group. The function returns vector index of observation to be trimmed.
#'
#' @param reweigting_factor Estimated reweigting factor
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
    trimming_threshold <- sqrt(nobs) / obs
  }
  reweighting_factor_in_trimming_range <- reweighting_factor_control[which(reweighting_factor_control / sum(reweighting_factor_control) > trimming_threshold)]
  trim_value <- ifelse(length(reweighting_factor_in_trimming_range) == 0,
                       sum(reweighting_factor),
                       min(reweighting_factor_in_trimming_range))

  observations_to_be_trimmed <- which(reweighting_factor >= trim_value)

  return(observations_to_be_trimmed)
}


