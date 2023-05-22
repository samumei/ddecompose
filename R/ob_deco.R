# oaxaca package ---------------------------------------------------------------
# library("oaxaca")
# ?oaxaca
# data("chicago")
# oaxaca.results.1 <- oaxaca(ln.real.wage ~ age + female + LTHS + some.college +
#                              college + advanced.degree | foreign.born,
#                            data = chicago, R = 30)
# print(oaxaca.results.1)
# plot(oaxaca.results.1)

# deco_results <- ob_deco(formula = ln.real.wage ~ age + female + LTHS + some.college +
#                             college + advanced.degree, data = chicago, group = foreign.born)
# deco_results


# Start actual function -------------------------------------------------------

# Data
# set.seed(125)
# library("AER")
# data("CPS1985")
# formula <- log(wage) ~ education + experience + union + ethnicity
# data_used <- CPS1985
# data_used$weights <- runif(nrow(CPS1985), 0.5, 1.5)
# # data_used[1,1] <- NA
# data_used <- get_all_vars(formula, data_used, weights=weights, group=gender)
# # cluster <- data_used$cluster <- rbinom(nrow(CPS1985), 10, 0.5)
# cluster <- NULL
#
# na.action = na.exclude
# reference_0 <- TRUE
# vcov=stats::vcov
# normalize_factors=FALSE
# bootstrap = FALSE
# bootstrap_iterations = 100


#' Oaxaca-Blinder decomposition
#'
#' Following Oaxaca (1973) and Blinder (1973), this function decomposes between-
#' group differences in the mean of an outcome variable into a part explained
#' by mean-differences in the observable covariates ("composition effect") and
#' another part due to the returns to these covariates ("structure effect"). The
#' function estimates the returns to covariates using OLS separately for the
#' two groups.
#'
#' @param formula an object of class "formula". See \link[stats]{lm} for further details.
#' @param data a data frame containing the variables in the model.
#' @param weights numeric vector of non-negative observation weights, hence of same length as \code{dep_var}.
#'                The default (\code{NULL}) is equivalent to \code{weights = rep(1, length(dep_var))}.
#' @param na.action generic function that defines how NAs in the data should be handled.
#'                  Default is \code{na.omit}, leading to exclusion of observations that contain one or more missings.
#'                  See \link[stats]{na.action} for further details.
#' @param group name of the a binary variable (numeric or factor)
#' identifying the two groups that will be compared. The group identified by the
#' lower ranked value in `group` (i.e., 0 in the case of a dummy variable or the
#' first level of factor variable) is defined as group 0.
#' @param reference_0 boolean: if `TRUE` (default), then the group 0 -- i.e.,
#' the group identified by the lower ranked value in `group` -- will be defined
#' as reference group. The reference group will be reweighted to match the
#' covariates distribution of the counterfactual sample.
#' @param normalize_factors boolean: If `TRUE`, then factor variables are normalized as
#' proposed by Gardeazabal/Ugidos (2004) and results are not dependent on the factor's
#' reference group. Per default (\code{normalize_factors  = FALSE}), factors are
#' normalized.
#' @param bootstrap boolean: If `FALSE`, then the estimation is not boostrapped and no
#' standard errors are calculated.
#' @param bootstrap_iterations positive integer indicating the number of bootstrap
#'  iterations to execute. Only required if \code{bootstrap = TRUE}.
#' @param bootstrap_robust boolean: if `FALSE` (default), then bootstrapped standard
#' errors are estimated as the standard deviations of the bootstrapp estimates.
#' Otherwise, the function uses the bootstrap interquartile range rescaled by the
#' interquantile range of the standard distribution to estimate standard errors.
#' @param cluster numeric vector of same length as \code{dep_var} indicating the
#' clustering of observations. If \code{cluster = NULL} (default), no clustering
#' is a assumend and bootstrap procedure resamples individual observations. Other
#' wise bootstrap procedure resamples clusters.
#' @param cores positive integer indicating the number of cores to use when
#' computing bootstrap standard errors. Only required if \code{bootstrap = TRUE}.
#' @param vcov function estimating covariance matrix of regression coefficients if
#' standard errors are not bootstrapped (i.e., \code{bootstrap = FALSE}). By default,
#' \link[stats]{vcov} is used assuming homoskedastic errors.
#'
#' @references
#' Fortin, Nicole, Thomas Lemieux, and Sergio Firpo. 2011. "Decomposition methods in economics."
#' In Orley Ashenfelter and David Card, eds., \emph{Handbook of labor economics}. Vol. 4. Elsevier, 1-102.
#'
#' Gardeazabal, Javier, and Arantza Ugidos. 2004. "More on identification in detailed wage decompositions."
#' \emph{Review of Economics and Statistics} 86(4): 1034-1036.
#'
#' @examples
#'
#' ## Decompose gender wage gap
#' ## with NLYS79 data like in Fortin, Lemieux, & Firpo (2011: 41)
#'
#' load("data/nlys00.rda")
#'
#' mod1 <- log(wage) ~ age + central_city + msa + region + black +
#' hispanic + education + afqt + family_responsibility + years_worked_civilian +
#' years_worked_military + part_time + industry
#'
#' # Using female coefficients (reference_0 = TRUE) to estimate counterfactual mean
#' deco_female_as_reference <- ob_deco(formula = mod1, data = nlys00, group = female, reference_0 = TRUE)
#' deco_female_as_reference
#'
#' # Using male coefficients (reference_0 = FALSE)
#' deco_male_as_reference <- ob_deco(formula = mod1, data = nlys00, group = female, reference_0 = FALSE)
#' deco_male_as_reference
#'
#' # Replicate first and third column in Table 3 in Fortin, Lemieux, & Firpo (2011: 41)
#' # Define aggregation of decomposition terms
#' custom_aggregation <- list(`Age, race, region, etc.` = c("age", "blackyes", "hispanicyes", "regionNorth-central", "regionSouth", "regionWest", "central_cityyes", "msayes"),
#'                    `Education` = c("education<10 yrs", "educationHS grad (diploma)", "educationHS grad (GED)", "educationSome college", "educationBA or equiv. degree", "educationMA or equiv. degree", "educationPh.D or prof. degree"),
#'                    `AFTQ` = "afqt",
#'                    `L.T. withdrawal due to family` =  "family_responsibility",
#'                    `Life-time work experience` = c("years_worked_civilian", "years_worked_military", "part_time"),
#'                    `Industrial sectors` = c("industryManufacturing", "industryEducation, Health, Public Admin.", "industryOther services"))
#'
#' # First column
#' summary(deco_male_as_reference, custom_aggregation = custom_aggregation)
#'
#' # Third column
#' summary(deco_female_as_reference, custom_aggregation = custom_aggregation)
#'
#' # Compare bootstrapped standard errors...
#' deco_female_as_reference_bs <- ob_deco(formula = mod1,
#'                                         data = nlys00,
#'                                         group = female,
#'                                         bootstrap = TRUE,
#'                                         bootstrap_iterations = 100)
#' summary(deco_female_as_reference_bs, custom_aggregation = custom_aggregation)
#'
#' # ... to analytical standard errors
#' summary(deco_female_as_reference, custom_aggregation = custom_aggregation)
#'
ob_deco <- function(formula,
                    data,
                    weights,
                    na.action = na.exclude,
                    group,
                    reference_0=TRUE,
                    normalize_factors = FALSE,
                    bootstrap = FALSE,
                    bootstrap_iterations = 100,
                    bootstrap_robust = FALSE,
                    cluster = NULL,
                    cores = 1,
                    vcov=stats::vcov){

  ## Get model.frame
  function_call <- match.call()
  data_arguments_index <- match(c("formula", "data", "weights", "group"), names(function_call), 0)
  data_arguments <- function_call[c(1, data_arguments_index)]
  #data_arguments$drop.unused.levels <- TRUE

  data_arguments[[1]] <- as.name("get_all_vars") #as.name("model.frame")
  data_used <- eval.parent(data_arguments)
  data_used <- na.action(data_used)
  data_used <- lapply(list(data_used), na.action)[[1]]

  ## Get weights
  weights <- model.weights(data_used)
  if (!is.null(weights) && !is.numeric(weights)) {
    stop("'weights' must be a numeric vector")
  }
  if (is.null(weights)) {
    data_used$weights <- rep(1, nrow(data_used))
  }

  ## Check group variable
  group_variable_name <- data_arguments[["group"]]
  group_variable <- data_used[, "group"]
  check_group_variable <- is.numeric(group_variable)&length(unique(group_variable))==2|is.factor(group_variable)&length(unique(group_variable))==2
  if(check_group_variable==FALSE){
    stop("Group variable must either be a binary numeric variable or a binary factor variable.")
  }
  if(is.numeric(group_variable)){
    data_used[, "group"] <- as.factor(group_variable)
  }

  reference_group <- ifelse(reference_0, 0, 1)
  reference_group_print <- levels(data_used[, "group"])[reference_group + 1]

  compute_analytical_se <- ifelse(bootstrap, FALSE, TRUE)
  estimated_decomposition <- estimate_ob_deco(formula = formula,
                                              data_used = data_used,
                                              reference_0 = reference_0,
                                              normalize_factors = normalize_factors,
                                              compute_analytical_se = compute_analytical_se,
                                              return_model_fit = TRUE,
                                              vcov = vcov)

  if(bootstrap){
    if(is.null(cluster)==FALSE){
      if(length(cluster) != nrow(data_used)){
        stop("Vector `cluster` must have the same length as number of observations in `data`.")
      }
      cluster_weights <- do.call("rbind",lapply(split(data_used, cluster), function(x) data.frame(cluster_weights=sum(x$weight))))
      data_used$cluster <- cluster
      data_used$cluster_weights <-  cluster_weights[match(as.character(cluster),rownames(cluster_weights)),]
    }

    cat("Bootstrapping standard errors...\n")

    if(cores == 1) {
      bootstrap_estimates <- pbapply::pblapply(1:bootstrap_iterations,
                                               function(x) bootstrap_estimate_ob_deco(formula = formula,
                                                                                      data_used = data_used,
                                                                                      reference_0 = reference_0,
                                                                                      normalize_factors = normalize_factors,
                                                                                      cluster = cluster))
    } else {
      cores <- min(cores, parallel::detectCores() - 1)
      core_cluster <- parallel::makeCluster(cores)
      parallel::clusterSetRNGStream(core_cluster, round(runif(1, 0, 100000)))
      parallel::clusterExport(cl = core_cluster,
                              varlist = ls(),
                              envir = environment())
      bootstrap_estimates <- pbapply::pblapply(1:bootstrap_iterations,
                                               function(x) bootstrap_estimate_ob_deco(formula = formula,
                                                                                      data_used = data_used,
                                                                                      reference_0 = reference_0,
                                                                                      normalize_factors = normalize_factors,
                                                                                      cluster = cluster),
                                               cl = core_cluster)
      parallel::stopCluster(core_cluster)
    }

    bootstrap_estimates <- do.call("rbind", bootstrap_estimates)
    bootstrap_estimates$iteration <- rep(1:bootstrap_iterations, each=length(unique(bootstrap_estimates$Variable)))
    bootstrap_estimates <- stats::reshape(bootstrap_estimates,
                                                  idvar = c("Variable","iteration"),
                                                  ids=unique(bootstrap_estimates$variable),
                                                  times = setdiff(names(bootstrap_estimates),c("Variable","iteration")),
                                                  timevar="effect",
                                                  varying = list(setdiff(names(bootstrap_estimates),c("Variable","iteration"))),
                                                  direction = "long",
                                                  v.names = "value")
    bootstrap_estimates <- stats::reshape(bootstrap_estimates,
                                                  idvar = c("iteration","effect"),
                                                  timevar="Variable",
                                                  direction = "wide")
    names(bootstrap_estimates) <- gsub("value[.]", "", names(bootstrap_estimates))
    bootstrap_estimates <- lapply(split(bootstrap_estimates, bootstrap_estimates$effect), function(x) stats::cov(x[, -c(1:2)]))

    decomposition_terms_se <- as.data.frame(do.call("cbind", lapply(bootstrap_estimates, function(x) sqrt(diag(x)))))
    decomposition_terms_se$Variable <- rownames(decomposition_terms_se)
    decomposition_terms_se <- decomposition_terms_se[, c("Variable", "Observed_difference", "Composition_effect", "Structure_effect")]

    decomposition_terms_vcov <- lapply(bootstrap_estimates, function(x) x[-1,-1])

    estimated_decomposition$decomposition_vcov$decomposition_terms_se <- decomposition_terms_se
    estimated_decomposition$decomposition_vcov$decomposition_terms_vcov <- decomposition_terms_vcov
  }

  add_to_results <- list(group_variable_name=group_variable_name,
                         group_variable_levels=levels(group_variable),
                         reference_group=reference_group_print,
                         normalize_factors=normalize_factors,
                         bootstrap=bootstrap)

  estimated_decomposition <- c(estimated_decomposition,
                               add_to_results)

  class(estimated_decomposition) <- "ob_deco"
  return(estimated_decomposition)

}

#' Estimate OB decomposition
#'
estimate_ob_deco <- function(formula,
                             data_used,
                             reference_0 = TRUE,
                             normalize_factors = FALSE,
                             compute_analytical_se = TRUE,
                             return_model_fits = TRUE,
                             vcov=stats::vcov){

  group0 <- levels(data_used[, "group"])[1]

  obs_0 <- which(data_used[, "group"] == group0)
  obs_1 <- which(data_used[, "group"] != group0)

  weights0 <- data_used[obs_0, "weights"]
  weights1 <- data_used[obs_1, "weights"]

  if(normalize_factors){
    normalized_data <- GU_normalization(formula=formula,
                                        data=data_used,
                                        weights=weights,
                                        group=group)
    formula <- normalized_data$formula
    data_used <- normalized_data$data
    adjusted_coefficient_names <- normalized_data$adjusted_coefficient_names
    X0 <- normalized_data$regressors_for_prediction[obs_0, ]
    X1 <- normalized_data$regressors_for_prediction[obs_1, ]

    #Insert here different X0 and X1 for predictions!
  }else{
    X0 <- model.matrix(formula, data_used[obs_0, ])
    X1 <- model.matrix(formula, data_used[obs_1, ])

    adjusted_coefficient_names <- NULL
  }

  #browser()
  fit0 <- lm(formula, data = subset(data_used, group == group0), weights = weights)
  fit1 <- lm(formula, data = subset(data_used, group != group0), weights = weights)

  beta0 <- coef(fit0)
  beta1 <- coef(fit1)

  if(normalize_factors){
    beta0 <- GU_normalization_get_coefficients(coef_names = adjusted_coefficient_names,
                                               est_coef = beta0)
    beta1 <- GU_normalization_get_coefficients(coef_names = adjusted_coefficient_names,
                                               est_coef = beta1)
  }

  estimated_deco_terms <- ob_deco_calculate_terms(beta0 = beta0,
                                                  beta1 = beta1,
                                                  X0 = X0,
                                                  X1 = X1,
                                                  weights0 = weights0,
                                                  weights1 = weights1,
                                                  reference_0 = reference_0)

  if(compute_analytical_se){
    Cov_beta0 <- vcov(fit0) #lapply(list(fit0), vcov)[[1]]
    Cov_beta1 <- vcov(fit1) #lapply(list(fit1), vcov)[[1]]
    if(normalize_factors){
      Cov_beta0 <- GU_normalization_get_vcov(coef_names = adjusted_coefficient_names,
                                             Cov_beta = Cov_beta0)
      Cov_beta1 <- GU_normalization_get_vcov(coef_names = adjusted_coefficient_names,
                                             Cov_beta = Cov_beta1)
    }

    estimated_deco_vcov <-  ob_deco_calculate_vcov(beta0 = beta0,
                                                   beta1 = beta1,
                                                   X0 = X0,
                                                   X1 = X1,
                                                   weights0 = weights0,
                                                   weights1 = weights1,
                                                   Cov_beta0  =  Cov_beta0,
                                                   Cov_beta1  =  Cov_beta1,
                                                   reference_0 = reference_0)
  }else{
    estimated_deco_vcov <- NULL
  }

  if(return_model_fits){
    model_fits <- list(fit_group_0 = fit0,
                       fit_group_1 = fit1)
  }else{
    model_fits <- NULL
  }

  results <- list(decomposition_terms = estimated_deco_terms,
                  decomposition_vcov = estimated_deco_vcov,
                  model_fits = model_fits,
                  GU_normalized_coefficient_names = adjusted_coefficient_names)
  return(results)
}

#' Estimate OB decomposition in bootstrap replications
#'
bootstrap_estimate_ob_deco <- function(formula,
                                       data_used,
                                       reference_0 = TRUE,
                                       normalize_factors = FALSE,
                                       cluster = NULL){
  if(is.null(cluster)){
  sampled_observations <- sample(1:nrow(data_used),
                                 size = nrow(data_used),
                                 replace = TRUE,
                                 prob = data_used$weights/sum(data_used$weights, na.rm=TRUE))
  } else {
  unique_cluster <- unique(data_used$cluster)
  cluster_weights <- data_used[match(unique_cluster, data_used$cluster), "cluster_weights"]
  sampled_cluster <- sample(unique_cluster,
                            size = length(unique_cluster),
                            replace = TRUE,
                            prob = cluster_weights/sum(data_used$weights, na.rm=TRUE))
  sampled_observations <- do.call("c",sapply(sampled_cluster, function(x) which(data_used$cluster %in% x)))
  data_used$weights <- data_used$weights * sum(data_used[ ,"weights"], na.rm=TRUE) / sum(data_used[sampled_observations,"weights"], na.rm=TRUE)
  }

  deco_estimates <- estimate_ob_deco(formula = formula,
                                     data_used = data_used[sampled_observations, ],
                                     reference_0 = reference_0,
                                     normalize_factors = normalize_factors,
                                     compute_analytical_se = FALSE,
                                     return_model_fit = FALSE)

  deco_estimates <- deco_estimates[["decomposition_terms"]]

  return(deco_estimates)
}



#' Calculate decomposition terms based on model.matrix and estimated OLS coefficients
#'
ob_deco_calculate_terms <- function(beta0,
                                    beta1,
                                    X0,
                                    X1,
                                    weights0,
                                    weights1,
                                    reference_0){

  X0 <- apply(X0, 2, weighted.mean, w=weights0)
  X1 <- apply(X1, 2, weighted.mean, w=weights1)

  Xb0 <- X0*beta0
  Xb1 <- X1*beta1

  observed_diff <- Xb1 - Xb0
  if(reference_0){
    XbC <- X1*beta0
    structure_effect <- Xb1 - XbC
    composition_effect <- XbC - Xb0
  }else{
    XbC <- X0*beta1
    composition_effect <- Xb1 - XbC
    structure_effect <- XbC - Xb0
  }

  agg_observed_diff <- sum(observed_diff)
  agg_composition_effect <- sum(composition_effect)
  agg_structure_effect <- sum(structure_effect)

  decomposition_terms <- data.frame(Variable = c("Total", names(observed_diff)),
                                    Observed_difference = c(agg_observed_diff, observed_diff),
                                    Composition_effect = c(agg_composition_effect, composition_effect),
                                    Structure_effect = c(agg_structure_effect, structure_effect))
  return(decomposition_terms)

}

#' Estimate covariance matrix for OB decomposition terms
#'
#' assuming indenpendence between groups
#' see: https://www.stata.com/meeting/3german/jann.pdf
ob_deco_calculate_vcov  <- function(beta0,
                                    beta1,
                                    X0,
                                    X1,
                                    weights0,
                                    weights1,
                                    Cov_beta0,
                                    Cov_beta1,
                                    reference_0){


  Cov_X0 <- stats::cov.wt(X0, wt=weights0)$cov
  Cov_X1 <- stats::cov.wt(X1, wt=weights1)$cov
  Cov_X_diff <- Cov_X1 + Cov_X0

  Cov_beta_diff <- Cov_beta1 + Cov_beta0

  X0 <- apply(X0, 2, weighted.mean, w=weights0)
  X1 <- apply(X1, 2, weighted.mean, w=weights1)
  X_diff <- X1 - X0

  Xb0 <- X0*beta0
  Xb1 <- X1*beta1

  beta_diff <- beta1 - beta0

  if(reference_0){
    XbC <- X1*beta0
    structure_effect <- Xb1 - XbC
    composition_effect <- XbC - Xb0
    beta_reference <- beta0
    X_counterfactual <- X1
    Cov_beta_reference <-   Cov_beta0
    Cov_X_counterfactual <- Cov_X1
  }else{
    XbC <- X0*beta1
    composition_effect <- Xb1 - XbC
    structure_effect <- XbC - Xb0
    beta_reference <- beta1
    X_counterfactual <- X0
    Cov_beta_reference <-   Cov_beta1
    Cov_X_counterfactual <- Cov_X0
  }

  Var_agg_observed_diff <- t(X1) %*% Cov_beta1 %*% X1 +
    t(beta1) %*% Cov_X1 %*% beta1 +
    sum(diag(Cov_X1 %*% Cov_beta1)) +
    t(X0) %*% Cov_beta0 %*% X0 +
    t(beta0) %*% Cov_X0 %*% beta0 +
    sum(diag(Cov_X0 %*% Cov_beta0))

  Cov_observed_diff <- Cov_beta1 * X1 %*% t(X1) +
    Cov_X1 * beta1 %*% t(beta1) +
    Cov_X1 %*% Cov_beta1 +
    Cov_beta0 * X0 %*% t(X0) +
    Cov_X0 * beta0 %*% t(beta0) +
    Cov_X0 %*% Cov_beta0

  Var_agg_composition_effect <- t(X_diff) %*% Cov_beta_reference %*% X_diff +
    t(beta_reference) %*% Cov_X_diff %*% beta_reference +
    sum(diag(Cov_X_diff %*% Cov_beta_reference))

  Cov_composition_effect <- Cov_beta_reference * X_diff %*% t(X_diff) +
    Cov_X_diff * beta_reference %*% t(beta_reference) +
    Cov_X_diff %*% Cov_beta_reference

  Var_agg_structure_effect <- t(X_counterfactual) %*% Cov_beta_diff %*% X_counterfactual +
    t(beta_diff) %*% Cov_X_counterfactual %*% beta_diff +
    sum(diag(Cov_X_counterfactual %*% Cov_beta_diff))

  Cov_structure_effect <- Cov_beta_diff * X_counterfactual %*% t(X_counterfactual) +
    Cov_X_counterfactual * beta_diff %*% t(beta_diff) +
    Cov_X_counterfactual %*% Cov_beta_diff

  decomposition_terms_se <- data.frame(Variable = c("Total", names(X0)),
                                    Observed_difference = sqrt(c(Var_agg_observed_diff, diag(Cov_observed_diff))),
                                    Composition_effect = sqrt(c( Var_agg_composition_effect , diag(Cov_composition_effect))),
                                    Structure_effect = sqrt(c(Var_agg_structure_effect, diag(Cov_structure_effect))))

  vcov_list <- list(Observed_difference=Cov_observed_diff,
                    Composition_effect= Cov_composition_effect,
                    Structure_effect= Cov_structure_effect)

  results <- list(decomposition_terms_se=decomposition_terms_se,
                  vcov=vcov_list)
  return(results)
}

