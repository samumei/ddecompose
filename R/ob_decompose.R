#' Oaxaca-Blinder decomposition
#'
#' @description \code{ob_deco} implements the Oaxaca-Blinder decomposition that
#' divides differences in the mean outcome between two groups into one part explained
#' by different covariate means (composition effect) and into another part due to
#' differences in linear regression coefficients linking covariates to the outcome
#' variable (structure effect).
#'
#' The function allows for 'doubly robust' decompositions where the sample of one
#' group is reweighted such that it matches the covariates distribution of the
#' other group before the regression coefficients are estimated.
#'
#' For distributional statistics beyond the mean, the function performs the RIF
#' regression decomposition proposed by Firpo, Fortin, and Lemieux (2018).
#'
#' @param formula a \code{formula} object with an outcome variable Y on the left-hand side
#' and the covariates X on the right-hand side. If \code{reweighting = TRUE}, the same
#' covariates are used to estimate the conditional probabilities for the reweighting factor.
#' A different model for estimating the conditional probabilities can be defined
#' after a \code{|} operator on the right-hand side.
#' @param data a data frame containing the variables in the model.
#' @param weights numeric vector of non-negative observation weights, hence of same length as \code{dep_var}.
#'                The default (\code{NULL}) is equivalent to \code{weights = rep(1, length(dep_var))}.
#'                If no weights are used, make sure you do not define this parameter (e.g. with \code{weights = NULL}).
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
#' By default, the composition effect is computed as \code{(X1 - X0) * b0} and
#' the structure effect as \code{X1 * (b1 - b0)}. Putting \code{reference_0 = FALSE} changes
#' the reference structure. Hence, the composition effect is computed as \code{(X1 - X0) * b1} and
#' the structure effect as \code{X0 * (b1 - b0)}.
#' @param subtract_1_from_0 boolean: By default (`FALSE`), X0 is subtracted from X1 and beta0 from beta1 (X1b1 - X0b0)
#' to compute the overall difference. Setting `subtract_1_from_0` to `TRUE` merely changes the sign of the decomposition results.
#' This means the composition effect is computed as \code{(X0 - X1) * b1} and
#' the structure effect as \code{X0 * (b0 - b1)}.
#' @param rifreg_statistic string containing the distributional statistic for which to compute the RIF.
#'                         If `NULL` (default), no RIF regression decomposition is computed.
#'                         If an available statistic is selected, `ob_deco` estimates a RIF regression decomposition.
#'                         The `rifreg_statistic` can be one of
#'                  "quantiles", "mean", "variance", "gini", "interquantile_range", "interquantile_ratio", or "custom".
#'                   If "custom" is selected, a \code{custom_rif_function} needs to be provided.
#' @param rifreg_probs a vector of length 1 or more with probabilities of quantiles. Each quantile is indicated with a value between 0 and 1.
#'              Default is \code{c(1:9)/10}. If \code{statistic = "quantiles"}, a single RIF regression for every quantile in \code{probs}
#'              is estimated. An interquantile ratio (range) is defined by the ratio (difference) between the \code{max(probs)}-quantile and
#'              the \code{min(probs)}-quantile.
#' @param custom_rif_function the RIF function to compute the RIF of the custom distributional statistic.
#'                            Default is NULL. Only needs to be provided if \code{statistic = "custom"}.
#'                            Every custom_rif_function needs the parameters \code{dep_var}, \code{weights} and \code{probs}.
#'                            If they are not needed, they must be set to NULL in the function definition (e.g. \code{probs = NULL}).
#'                            A custom function must return a data frame containing at least a "rif" and "weights" column.
#'                            See \code{examples} for further details.
#' @param reweighting boolean: if `TRUE`, then the decomposition is performed with
#' with respect to reweighted reference group yielding either a 'doubly robust'
#' Oaxaca-Blinder decomposition or a reweighted RIF decomposition.
#' @param reweighting_method  specifies the method fit and predict conditional probabilities
#' used to derive the reweighting factor. Currently, \code{"logit"}, \code{"fastglm"},
#' and \code{"random forests"} are available.
#' @param trimming boolean: If \code{TRUE}, observations with dominant reweighting factor
#' values are trimmend according to rule of Huber, Lechner, and Wunsch (2013). Per
#' default, trimming is set to \code{FALSE}.
#' @param trimming_threshold numeric: threshold defining the maximal accepted
#' relative weight of the reweighting factor value (i.e., inverse probability weight)
#' of a single observation. If \code{NULL}, the threshold is set to \eqn{sqrt(N)/N},
#' where \eqn{N} is the number of observations in the reference group.
#' @param normalize_factors boolean: If `TRUE`, then factor variables are normalized as
#' proposed by Gardeazabal/Ugidos (2004) and results are not dependent on the factor's
#' reference group. Per default (\code{normalize_factors  = FALSE}) and factors are not
#' normalized.
#' @param bootstrap boolean: If `FALSE` (default), then no bootstrapped standard
#' errors are calculated and, in the case of a standard Oaxaca-Blinder decomposition,
#' analytical standard errors are estimated (assuming independence between groups).
#' @param bootstrap_iterations positive integer indicating the number of bootstrap
#' iterations to execute. Only required if \code{bootstrap = TRUE}.
#' @param bootstrap_robust boolean: if `FALSE` (default), then bootstrapped standard
#' errors are estimated as the standard deviations of the bootstrapp estimates.
#' Otherwise, the function uses the bootstrap interquartile range rescaled by the
#' interquantile range of the standard distribution to estimate standard errors.
#' @param cluster numeric vector of same length as \code{dep_var} indicating the
#' clustering of observations. If \code{cluster = NULL} (default), no clustering
#' is a assumend and bootstrap procedure resamples individual observations. Otherwise
#' bootstrap procedure resamples clusters.
#' @param cores positive integer indicating the number of cores to use when
#' computing bootstrap standard errors. Only required if \code{bootstrap = TRUE}.
#' @param vcov function estimating covariance matrix of regression coefficients if
#' standard errors are not bootstrapped (i.e., \code{bootstrap = FALSE}). By default,
#' \link[stats]{vcov} is used assuming homoscedastic errors.
#' @param ... additional parameters passed to the custom_rif_function.
#' Apart from dep_var, weights and probs they must have a different name than the the ones in rifreg.
#' For instance, if you want to pass a parameter statistic to the custom_rif_function, name it custom_statistic.
#' Additional parameters can also be passed to the \link[stats]{density} function used
#' to estimate the RIF of quantiles.
#'
#' @references
#' Fortin, Nicole, Thomas Lemieux, and Sergio Firpo. 2011. "Decomposition methods in economics."
#' In Orley Ashenfelter and David Card, eds., \emph{Handbook of labor economics}. Vol. 4. Elsevier, 1-102.
#'
#' Gardeazabal, Javier, and Arantza Ugidos. 2004. "More on identification in detailed wage decompositions."
#' \emph{Review of Economics and Statistics} 86(4): 1034-1036.
#'
#' @export
#'
#' @examples
#'
#' ## Oxaca-Blinder decomposition of gender wage gap
#' ## with NLYS79 data like in Fortin, Lemieux, & Firpo (2011: 41)
#'
#' data("nlys00")
#'
#' mod1 <- log(wage) ~ age + central_city + msa + region + black +
#'   hispanic + education + afqt + family_responsibility + years_worked_civilian +
#'   years_worked_military + part_time + industry
#'
#' # Using female coefficients (reference_0 = TRUE) to estimate counterfactual mean
#' deco_female_as_reference <- ob_decompose(formula = mod1,
#'                                     data = nlys00,
#'                                     group = female,
#'                                     reference_0 = TRUE)
#' deco_female_as_reference
#'
#' # Using male coefficients (reference_0 = FALSE)
#' deco_male_as_reference <- ob_decompose(formula = mod1,
#'                                   data = nlys00,
#'                                   group = female,
#'                                   reference_0 = FALSE)
#' deco_male_as_reference
#'
#' # Replicate first and third column in Table 3 in Fortin, Lemieux, & Firpo (2011: 41)
#' # Define aggregation of decomposition terms
#' custom_aggregation <- list(`Age, race, region, etc.` = c("age",
#'                                                          "blackyes",
#'                                                          "hispanicyes",
#'                                                          "regionNorth-central",
#'                                                          "regionSouth",
#'                                                          "regionWest",
#'                                                          "central_cityyes",
#'                                                          "msayes"),
#'                    `Education` = c("education<10 yrs",
#'                                    "educationHS grad (diploma)",
#'                                    "educationHS grad (GED)",
#'                                    "educationSome college",
#'                                    "educationBA or equiv. degree",
#'                                    "educationMA or equiv. degree",
#'                                    "educationPh.D or prof. degree"),
#'                    `AFTQ` = "afqt",
#'                    `L.T. withdrawal due to family` =  "family_responsibility",
#'                    `Life-time work experience` = c("years_worked_civilian",
#'                                                    "years_worked_military",
#'                                                    "part_time"),
#'                    `Industrial sectors` = c("industryManufacturing",
#'                                             "industryEducation, Health, Public Admin.",
#'                                             "industryOther services"))
#'
#' # First column
#' summary(deco_male_as_reference, custom_aggregation = custom_aggregation)
#'
#' # Third column
#' summary(deco_female_as_reference, custom_aggregation = custom_aggregation)
#'
#' ## Compare bootstrapped standard errors...
#' deco_female_as_reference_bs <- ob_decompose(formula = mod1,
#'                                        data = nlys00,
#'                                        group = female,
#'                                        bootstrap = TRUE,
#'                                        bootstrap_iterations = 100)
#' summary(deco_female_as_reference_bs, custom_aggregation = custom_aggregation)
#'
#' # ... to analytical standard errors (assuming independence between groups and
#' # homoskedasticity)
#' deco_female_as_reference <- ob_decompose(formula = mod1,
#'                                     data = nlys00,
#'                                     group = female,
#'                                     reference_0 = TRUE)
#' summary(deco_female_as_reference, custom_aggregation = custom_aggregation)
#'
#' # Return standard errors for all detailed terms
#' summary(deco_female_as_reference, aggregate_factors = FALSE)
#'
#'
#' ## 'Doubly robust' Oxaca-Blinder decomposition of gender wage gap
#'  mod2 <- log(wage) ~ age + central_city + msa + region + black +
#'   hispanic + education + afqt + family_responsibility + years_worked_civilian +
#'   years_worked_military + part_time + industry | age  + (central_city + msa) * region + (black +
#'   hispanic) * (education + afqt) + family_responsibility * (years_worked_civilian +
#'   years_worked_military) + part_time * industry
#' deco_male_as_reference_robust <- ob_decompose(formula = mod2,
#'                                          data = nlys00,
#'                                          group = female,
#'                                          reference_0 = FALSE,
#'                                          reweighting = TRUE)
#'
#'
ob_decompose <- function(formula,
                    data,
                    group,
                    subtract_1_from_0 = FALSE,
                    weights = NULL,
                    rifreg_statistic = NULL,
                    rifreg_probs = c(1:9)/10,
                    custom_rif_function = NULL,
                    reweighting = FALSE,
                    reweighting_method = "logit",
                    reference_0 = TRUE,
                    na.action = na.omit,
                    normalize_factors = FALSE,
                    trimming = FALSE,
                    trimming_threshold = NULL,
                    bootstrap = FALSE,
                    bootstrap_iterations = 100,
                    bootstrap_robust = FALSE,
                    cluster = NULL,
                    cores = 1,
                    vcov=stats::vcov,
                    ...){

  rifreg <- ifelse(is.null(rifreg_statistic), FALSE, TRUE)

  if(rifreg & !reweighting) warning("If you want to decompose rifregression it is highly recommended to apply reweighting!
                                    See references for further information.")
  if(!typeof(data) == "list") stop("The provided \"data\" is not a data frame.")


  ## Get model.frame
  function_call <- match.call()
  data_arguments_index <- match(c("formula", "data", "weights", "group"), names(function_call), 0)
  data_arguments <- function_call[c(1, data_arguments_index)]
  #data_arguments$drop.unused.levels <- TRUE

  data_arguments[[1]] <- as.name("get_all_vars") #as.name("model.frame")
  data_used <- eval.parent(data_arguments)
  data_used <- na.action(data_used)
  #data_used <- lapply(list(data_used), na.action)[[1]]

  ## Check group variable
  group_variable_name <- data_arguments[["group"]]
  group_variable <- data_used[, "group"]
  check_group_variable <- is.numeric(group_variable)&length(unique(group_variable))==2|is.factor(group_variable)&length(unique(group_variable))==2
  if(!check_group_variable){
    stop("Group variable must either be a binary numeric variable or a binary factor variable.")
  }
  if(is.numeric(group_variable)){
    data_used[, "group"] <- as.factor(group_variable)
  }

  reference_group <- ifelse(reference_0, 0, 1)
  reference_group_print <- levels(data_used[, "group"])[reference_group + 1]

  # Get formula(s)
  formula <- Formula::as.Formula(formula)
  #data_arguments$formula <- formula

  nvar <- length(formula)[2] # Number of detailed decomposition effects
  if(nvar == 1) {
    if(reweighting) {
      print("The same model specification is used for decomposition and to compute reweighting weights.")
    }
    formula_decomposition <- formula
    formula_reweighting <- formula
  }
  else {
    if(!nvar == 2) stop("Cannot parse formula. See documentation and examples for further details.")
    if(!reweighting) warning("Parameter \"reweighting\" is set to FALSE. No reweighting is applied and given reweighting formula is ignored.")
    formula_decomposition <- stats::formula(formula, rhs=1, collapse=TRUE)
    formula_reweighting <- stats::formula(formula, rhs=2, collapse=TRUE)
  }


  ## Get weights
  if (!is.null(data_used$weights) && !is.numeric(data_used$weights)) {
    stop("'weights' must be a numeric vector")
  }
  if (is.null(data_used$weights)) {
    data_used$weights <- rep(1, nrow(data_used))
  }



  if(!bootstrap & !rifreg & !reweighting) compute_analytical_se <- TRUE
  else compute_analytical_se <- FALSE

  if(reweighting) {
    dfl_deco_results <- dfl_decompose(formula = formula_reweighting,
                                 data = data_used,
                                 weights = weights,
                                 group = group,
                                 reference_0 = reference_0,
                                 method = reweighting_method,
                                 estimate_statistics = FALSE,
                                 trimming = trimming,
                                 trimming_threshold = trimming_threshold)

    reweighting_factor <- dfl_deco_results$reweighting_factor$Psi_X1
    data_used$weights_and_reweighting_factors <- data_used[, "weights"] * reweighting_factor
  }

  if(rifreg && rifreg_statistic == "quantiles" & length(rifreg_probs) > 1) {
      estimated_decomposition <- lapply(rifreg_probs, estimate_ob_decompose,
                                        formula = formula_decomposition, data_used = data_used,
                                        reference_0 = reference_0, normalize_factors = normalize_factors,
                                        compute_analytical_se = compute_analytical_se,
                                        return_model_fit = TRUE,
                                        reweighting = reweighting,
                                        rifreg = rifreg,
                                        rifreg_statistic = rifreg_statistic,
                                        custom_rif_function = custom_rif_function,
                                        na.action = na.action,
                                        vcov = vcov,
                                        ... = ...)

      for (i in seq_along(estimated_decomposition)) {
        estimated_decomposition[[i]][["decomposition_terms"]][, 2:6] <- estimated_decomposition[[i]][["decomposition_terms"]][, 2:6] * -1
      }

      names(estimated_decomposition) <- paste0("quantile_", as.character(rifreg_probs))

  }
  else {
    estimated_decomposition <- estimate_ob_decompose(formula = formula_decomposition,
                                                data_used = data_used,
                                                reference_0 = reference_0,
                                                normalize_factors = normalize_factors,
                                                compute_analytical_se = compute_analytical_se,
                                                return_model_fit = TRUE,
                                                reweighting = reweighting,
                                                rifreg = rifreg,
                                                rifreg_statistic = rifreg_statistic,
                                                rifreg_probs = rifreg_probs,
                                                custom_rif_function = custom_rif_function,
                                                na.action = na.action,
                                                vcov = vcov,
                                                ... = ...)

    if(subtract_1_from_0) estimated_decomposition$decomposition_terms[2:6] <- estimated_decomposition$decomposition_terms[2:6] *-1

    estimated_decomposition <- list(estimated_decomposition)

    # set name
    if(rifreg) {
      if(rifreg_statistic == "quantiles") {
        deco_name <- paste0("quantile_", as.character(rifreg_probs))
      }
      else {
        deco_name <- rifreg_statistic
      }
    }
    else {
      if(reweighting) {
        deco_name <- "reweighted_ob_decompose"
      }
      else{
        deco_name <- "ob_decompose"
      }
    }
    names(estimated_decomposition) <- deco_name
  }


  if(bootstrap){
    if(!is.null(cluster)){
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
                                               function(x) bootstrap_estimate_ob_decompose(formula_decomposition = formula_decomposition,
                                                                                      formula_reweighting = formula_reweighting,
                                                                                      data_used = data_used,
                                                                                      group = group,
                                                                                      reference_0 = reference_0,
                                                                                      normalize_factors = normalize_factors,
                                                                                      reweighting = reweighting,
                                                                                      reweighting_method = reweighting_method,
                                                                                      trimming = trimming,
                                                                                      trimming_threshold = trimming_threshold,
                                                                                      rifreg = rifreg,
                                                                                      rifreg_statistic = rifreg_statistic,
                                                                                      rifreg_probs = rifreg_probs,
                                                                                      custom_rif_function = custom_rif_function,
                                                                                      na.action = na.action,
                                                                                      cluster = cluster,
                                                                                      ... = ...))
    } else {
      cores <- min(cores, parallel::detectCores() - 1)
      core_cluster <- parallel::makeCluster(cores)
      parallel::clusterSetRNGStream(core_cluster, round(runif(1, 0, 100000)))
      parallel::clusterExport(cl = core_cluster,
                              varlist = ls(),
                              envir = environment())
      bootstrap_estimates <- pbapply::pblapply(1:bootstrap_iterations,
                                               function(x) bootstrap_estimate_ob_decompose(formula_decomposition = formula_decomposition,
                                                                                      formula_reweighting = formula_reweighting,
                                                                                      data_used = data_used,
                                                                                      group = group,
                                                                                      reference_0 = reference_0,
                                                                                      normalize_factors = normalize_factors,
                                                                                      reweighting = reweighting,
                                                                                      reweighting_method = reweighting_method,
                                                                                      trimming,
                                                                                      trimming_threshold,
                                                                                      rifreg = rifreg,
                                                                                      rifreg_statistic = rifreg_statistic,
                                                                                      rifreg_probs = rifreg_probs,
                                                                                      custom_rif_function = custom_rif_function,
                                                                                      na.action = na.action,
                                                                                      cluster = cluster,
                                                                                      ... = ...),
                                               cl = core_cluster)
      parallel::stopCluster(core_cluster)
    }

    if(rifreg && rifreg_statistic == "quantiles" & length(rifreg_probs) > 1) {
      for(i in 1:length(rifreg_probs)) {
        current_bootstrap_estimates <- lapply(bootstrap_estimates, function(x) x[[i]])
        bootstrap_vcov <- retrieve_bootstrap_vcov(current_bootstrap_estimates, bootstrap_iterations)
        estimated_decomposition[[i]][["decomposition_vcov"]][["decomposition_terms_se"]] <- bootstrap_vcov$decomposition_terms_se
        estimated_decomposition[[i]][["decomposition_vcov"]][["decomposition_terms_vcov"]] <- bootstrap_vcov$decomposition_terms_vcov
      }
    }
    else {
        current_bootstrap_estimates <- lapply(bootstrap_estimates, function(x) x[[1]])
        bootstrap_vcov <- retrieve_bootstrap_vcov(current_bootstrap_estimates, bootstrap_iterations)
        estimated_decomposition[[1]][["decomposition_vcov"]][["decomposition_terms_se"]] <- bootstrap_vcov$decomposition_terms_se
        estimated_decomposition[[1]][["decomposition_vcov"]][["decomposition_terms_vcov"]] <- bootstrap_vcov$decomposition_terms_vcov
    }


# estimated_decomposition$decomposition_vcov$decomposition_terms_se <- bootstrap_vcov$decomposition_terms_se
# estimated_decomposition$decomposition_vcov$decomposition_terms_vcov <- bootstrap_vcov$decomposition_terms_vcov

    # bootstrap_vcov <- lapply(bootstrap_estimates, retrieve_bootstrap_vcov,
    #                          bootstrap_iterations = bootstrap_iterations)


  }


  if(reweighting) {
    reweighting_estimates <- dfl_deco_results[c(3, 5:12)]
  }
  else {
    reweighting_estimates <- NA
  }

  add_to_results <- list(group_variable_name=group_variable_name,
                         group_variable_levels=levels(group_variable),
                         reference_group=reference_group_print,
                         reweighting_estimates=reweighting_estimates,
                         input_parameters = list(rifreg_statistic = rifreg_statistic,
                                                 rifreg_probs = rifreg_probs,
                                                 reweighting = reweighting,
                                                 reweighting_method = reweighting_method,
                                                 reference_0 = reference_0,
                                                 subtract_1_from_0 = subtract_1_from_0,
                                                 normalize_factors=normalize_factors,
                                                 bootstrap=bootstrap,
                                                 bootstrap_iterations = bootstrap_iterations))

  estimated_decomposition <- c(estimated_decomposition,
                               add_to_results)

  class(estimated_decomposition) <- "ob_decompose"
  return(estimated_decomposition)

}


# Estimate OB decomposition
#
estimate_ob_decompose <- function(formula,
                             data_used,
                             reference_0,
                             normalize_factors,
                             compute_analytical_se,
                             return_model_fit,
                             reweighting,
                             rifreg,
                             rifreg_statistic,
                             rifreg_probs,
                             custom_rif_function,
                             na.action,
                             vcov,
                             ...){

  group0 <- levels(data_used[, "group"])[1]

  obs_0 <- which(data_used[, "group"] == group0)
  obs_1 <- which(data_used[, "group"] != group0)

  weights0 <- data_used[obs_0, "weights"]
  weights1 <- data_used[obs_1, "weights"]

  if(normalize_factors){
    if(reweighting) {
      # store reweighting weights
      weights_and_reweighting_factors <- data_used$weights_and_reweighting_factors
    }

    normalized_data <- GU_normalization(formula=formula,
                                        data=data_used,
                                        weights=weights,
                                        group=group)
    formula <- normalized_data$formula
    data_used <- normalized_data$data

    if(reweighting) {
      # attach reweighting weights again
      data_used$weights_and_reweighting_factors <- weights_and_reweighting_factors
    }

    adjusted_coefficient_names <- normalized_data$adjusted_coefficient_names
    X0 <- normalized_data$regressors_for_prediction[obs_0, ]
    X1 <- normalized_data$regressors_for_prediction[obs_1, ]

    #Insert here different X0 and X1 for predictions!
  }else{
    X0 <- model.matrix(formula, data_used[obs_0, ])
    X1 <- model.matrix(formula, data_used[obs_1, ])

    adjusted_coefficient_names <- NULL
  }

  if(rifreg) {
    fit0 <- rifreg::rifreg(formula = formula,
                                  data = subset(data_used, group == group0),
                                  statistic = rifreg_statistic,
                                  weights = weights,
                                  probs = rifreg_probs,
                                  custom_rif_function = custom_rif_function,
                                  na.action = na.action,
                           bootstrap = FALSE,
                           ... = ...)

    fit1 <- rifreg::rifreg(formula = formula,
                           data = subset(data_used, group != group0),
                           statistic = rifreg_statistic,
                           weights = weights,
                           probs = rifreg_probs,
                           custom_rif_function = custom_rif_function,
                           na.action = na.action,
                           bootstrap = FALSE,
                           ... = ...)

    beta0 <- fit0$estimates[,1]
    beta1 <- fit1$estimates[,1]
  }
  else {
    fit0 <- lm(formula, data = subset(data_used, group == group0), weights = weights)
    fit1 <- lm(formula, data = subset(data_used, group!= group0), weights = weights)

    beta0 <- coef(fit0)
    beta1 <- coef(fit1)
  }

  if(normalize_factors){
    beta0 <- GU_normalization_get_coefficients(coef_names = adjusted_coefficient_names,
                                               est_coef = beta0)
    beta1 <- GU_normalization_get_coefficients(coef_names = adjusted_coefficient_names,
                                               est_coef = beta1)
  }

  if(reweighting) {
    if(rifreg) {
      if(reference_0) {
        fit_reweighted <- rifreg::rifreg(formula = formula,
                                                data = subset(data_used, group == group0),
                                                statistic = rifreg_statistic,
                                                weights = weights_and_reweighting_factors,
                                                probs = rifreg_probs,
                                                custom_rif_function = custom_rif_function,
                                                na.action = na.action,
                                         bootstrap = FALSE,
                                         ... = ...)
      }
      else {
        fit_reweighted <- rifreg::rifreg(formula = formula,
                                                data = subset(data_used, group != group0),
                                                statistic = rifreg_statistic,
                                                weights = weights_and_reweighting_factors,
                                                probs = rifreg_probs,
                                                custom_rif_function = custom_rif_function,
                                                na.action = na.action,
                                         bootstrap = FALSE,
                                         ... = ...)
        }


      beta_reweighted <- fit_reweighted$estimates[,1]

    }
    else {
      if(reference_0) {
        fit_reweighted <- lm(formula, data = subset(data_used, group == group0), weights = weights_and_reweighting_factors)
      }
      else {
        fit_reweighted <- lm(formula, data = subset(data_used, group != group0), weights = weights_and_reweighting_factors)
      }
      beta_reweighted <- coef(fit_reweighted)
    }


    if(normalize_factors){
      beta_reweighted <- GU_normalization_get_coefficients(coef_names = adjusted_coefficient_names,
                                                           est_coef = beta_reweighted)
    }

    if(reference_0) {
      deco_results_group0_group_reweighted <-
        ob_deco_calculate_terms(beta0 = beta0,
                                beta1 = beta_reweighted,
                                X0 = X0,
                                X1 = X0,
                                weights0 = weights0,
                                weights1 = data_used[obs_0, "weights_and_reweighting_factors"],
                                reference_0 = TRUE)

      deco_results_group_reweighted_group_1 <-
        ob_deco_calculate_terms(beta0 = beta_reweighted,
                                beta1 = beta1,
                                X0 = X0,
                                X1 = X1,
                                weights0 = data_used[obs_0, "weights_and_reweighting_factors"],
                                weights1 = weights1,
                                reference_0 = TRUE)

      deco_results <- deco_results_group0_group_reweighted
      deco_results$Composition_effect <- deco_results_group0_group_reweighted$Composition_effect
      deco_results$Specification_error <- deco_results_group0_group_reweighted$Structure_effect

      deco_results$Structure_effect <- deco_results_group_reweighted_group_1$Structure_effect
      deco_results$Reweighting_error <- deco_results_group_reweighted_group_1$Composition_effect

      deco_results$Observed_difference <- deco_results_group0_group_reweighted$Observed_difference +
        deco_results_group_reweighted_group_1$Observed_difference

    }
    else {
      deco_results_group0_group_reweighted <-
        ob_deco_calculate_terms(beta0 = beta0,
                                beta1 = beta_reweighted,
                                X0 = X0,
                                X1 = X1,
                                weights0 = weights0,
                                weights1 = data_used[obs_1, "weights_and_reweighting_factors"],
                                reference_0 = FALSE)

      deco_results_group_reweighted_group_1 <-
        ob_deco_calculate_terms(beta0 = beta_reweighted,
                                beta1 = beta1,
                                X0 = X1,
                                X1 = X1,
                                weights0 = data_used[obs_1, "weights_and_reweighting_factors"],
                                weights1 = weights1,
                                reference_0 = FALSE)

      deco_results <- deco_results_group0_group_reweighted

      deco_results$Composition_effect <- deco_results_group_reweighted_group_1$Composition_effect
      deco_results$Specification_error <- deco_results_group_reweighted_group_1$Structure_effect

      deco_results$Structure_effect <- deco_results_group0_group_reweighted$Structure_effect
      deco_results$Reweighting_error <- deco_results_group0_group_reweighted$Composition_effect

      deco_results$Observed_difference <- deco_results_group0_group_reweighted$Observed_difference +
        deco_results_group_reweighted_group_1$Observed_difference
    }

    estimated_deco_vcov <- NULL

  }
  else {
    deco_results <- ob_deco_calculate_terms(beta0 = beta0,
                                            beta1 = beta1,
                                            X0 = X0,
                                            X1 = X1,
                                            weights0 = weights0,
                                            weights1 = weights1,
                                            reference_0 = reference_0)
    deco_results$Specification_error <- NA
    deco_results$Reweighting_error <- NA
    fit_reweighted <- NA


    if(compute_analytical_se) {
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
  }


  if(return_model_fit){
    model_fits <- list(fit_group_0 = fit0,
                       fit_group_1 = fit1,
                       fit_group_reweighted = fit_reweighted)
  }else{
    model_fits <- NULL
  }

  results <- list(decomposition_terms = deco_results,
                  decomposition_vcov = estimated_deco_vcov,
                  model_fits = model_fits,
                  GU_normalized_coefficient_names = adjusted_coefficient_names)
  return(results)
}

# Estimate OB decomposition in bootstrap replications
#
bootstrap_estimate_ob_decompose <- function(formula_decomposition,
                                       formula_reweighting,
                                       data_used,
                                       group,
                                       reference_0,
                                       normalize_factors,
                                       reweighting,
                                       reweighting_method,
                                       trimming,
                                       trimming_threshold,
                                       rifreg,
                                       rifreg_statistic,
                                       rifreg_probs,
                                       custom_rif_function,
                                       na.action,
                                       cluster = NULL,
                                       ...){
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

  tryCatch({
    sink(nullfile())  # Start suppressing output

    if(reweighting) {
    dfl_deco_results <- suppressWarnings(dfl_decompose(formula = formula_reweighting,
                                 data = data_used[sampled_observations, ],
                                 weights = weights,
                                 group = group,
                                 reference_0 = reference_0,
                                 method = reweighting_method,
                                 estimate_statistics = FALSE,
                                 trimming = trimming,
                                 trimming_threshold = trimming_threshold))

    reweighting_factor <- dfl_deco_results$reweighting_factor$Psi_X1
    data_used[sampled_observations, "weights_and_reweighting_factors"] <- data_used[sampled_observations, "weights"] * reweighting_factor
  }

    if(rifreg && rifreg_statistic == "quantiles" & length(rifreg_probs) > 1) {
      deco_estimates <- suppressWarnings(lapply(rifreg_probs, estimate_ob_decompose,
                                        formula = formula_decomposition,
                               data_used = data_used[sampled_observations, ],
                                        reference_0 = reference_0,
                               normalize_factors = normalize_factors,
                                        compute_analytical_se = FALSE,
                                        return_model_fit = FALSE,
                                        reweighting = reweighting,
                                        rifreg = rifreg,
                                        rifreg_statistic = rifreg_statistic,
                                        custom_rif_function = custom_rif_function,
                                        na.action = na.action,
                                        vcov = NULL,
                               ... = ...))

      deco_estimates <- lapply(deco_estimates, function(x) x$decomposition_terms)

  }
  else {
    deco_estimates <- suppressWarnings(estimate_ob_decompose(formula = formula_decomposition,
                                       data_used = data_used[sampled_observations, ],
                                       reference_0 = reference_0,
                                       normalize_factors = normalize_factors,
                                       reweighting = reweighting,
                                       rifreg = rifreg,
                                       rifreg_statistic = rifreg_statistic,
                                       rifreg_probs = rifreg_probs,
                                       custom_rif_function = custom_rif_function,
                                       na.action = na.action,
                                       compute_analytical_se = FALSE,
                                       vcov = NULL,
                                       return_model_fit = FALSE,
                                       ... = ...))

    deco_estimates <- list(deco_estimates[["decomposition_terms"]])
  }
  }, error = function(e) {
    print(e)
    stop(e)
  }, finally = {
    sink()
  })




  return(deco_estimates)
}

retrieve_bootstrap_vcov <- function(bootstrap_estimates, bootstrap_iterations) {
  bootstrap_estimates_as_dataframe <- do.call("rbind", bootstrap_estimates)
  bootstrap_estimates_as_dataframe$iteration <- rep(1:bootstrap_iterations, each=length(unique(bootstrap_estimates_as_dataframe$Variable)))
  bootstrap_estimates_long <- stats::reshape(bootstrap_estimates_as_dataframe,
                                        idvar = c("Variable", "iteration"),
                                        ids=unique(bootstrap_estimates_as_dataframe$variable),
                                        times = setdiff(names(bootstrap_estimates_as_dataframe), c("Variable", "iteration")),
                                        timevar="effect",
                                        varying = list(setdiff(names(bootstrap_estimates_as_dataframe), c("Variable", "iteration"))),
                                        direction = "long",
                                        v.names = "value")
  bootstrap_estimates_wide <- stats::reshape(bootstrap_estimates_long,
                                        idvar = c("iteration", "effect"),
                                        timevar= "Variable",
                                        direction = "wide")

  names(bootstrap_estimates_wide) <- gsub("value[.]", "", names(bootstrap_estimates_wide))
  bootstrap_estimates <- lapply(split(bootstrap_estimates_wide, bootstrap_estimates_wide$effect), function(x) stats::cov(x[, -c(1:2)]))

  decomposition_terms_se <- as.data.frame(do.call("cbind", lapply(bootstrap_estimates, function(x) sqrt(diag(x)))))
  decomposition_terms_se$Variable <- rownames(decomposition_terms_se)
  decomposition_terms_se <- decomposition_terms_se[, c("Variable",
                                                       "Observed_difference",
                                                       "Composition_effect",
                                                       "Structure_effect",
                                                       "Specification_error",
                                                       "Reweighting_error")]

  decomposition_terms_vcov <- lapply(bootstrap_estimates, function(x) x[-1,-1])

  return(list(decomposition_terms_se = decomposition_terms_se,
              decomposition_terms_vcov = decomposition_terms_vcov))
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

  Xb0 <- X0 * beta0
  Xb1 <- X1 * beta1

  if(reference_0){
    observed_diff <- Xb1 - Xb0
    XbC <- X1*beta0
    composition_effect <- XbC - Xb0
    structure_effect <- Xb1 - XbC
  }else{
    observed_diff <- Xb1 - Xb0
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
#' assuming independence between groups
#' see also: https://www.stata.com/meeting/3german/jann.pdf
#'
#' Var(xb) = Var(x)Var(b) + Var(x)E(b)^2 + Var(b)E(x)^2
#' Cov(x1b1, x2b2) = Cov(x1,x2)Cov(b1,b2) + Cov(x1,x2)E(b1)E(b2) + Cov(b1,b2)E(x1)E(x2)
#' Matrix notation Sigma_b * Sigma_X + Sigma_b * mu %*% mu' + Sigma_X * b %*% b'
ob_deco_calculate_vcov  <- function(beta0,
                                    beta1,
                                    X0,
                                    X1,
                                    weights0,
                                    weights1,
                                    Cov_beta0,
                                    Cov_beta1,
                                    reference_0){

  ### New

  Cov_X0 <- stats::cov.wt(X0, wt=weights0)$cov / nrow(X0)
  Cov_X1 <- stats::cov.wt(X1, wt=weights1)$cov / nrow(X1)

  X0 <- apply(X0, 2, weighted.mean, w=weights0)
  X1 <- apply(X1, 2, weighted.mean, w=weights1)

  Cov_Xb0 <-  Cov_X0 * Cov_beta0 + Cov_beta0 * (X0 %*% t(X0)) + Cov_X0 * (beta0 %*% t(beta0))
  Cov_Xb1 <-  Cov_X1 * Cov_beta1 + Cov_beta1 * (X1 %*% t(X1)) + Cov_X1 * (beta1 %*% t(beta1))

  Cov_observed_diff <-  Cov_Xb1 + Cov_Xb0

  if(reference_0){
    Cov_structure_effect <- Cov_X1 * (Cov_beta1 + Cov_beta0) +
      (Cov_beta1 + Cov_beta0) * (X1 %*% t(X1)) +
      Cov_X1 * ((beta1 - beta0) %*% t((beta1 - beta0)))
    Cov_composition_effect <-  (Cov_X1 + Cov_X0) * Cov_beta0 +
      Cov_beta0 * ((X1 - X0)%*% t((X1 - X0))) +
      (Cov_X1 + Cov_X0) * (beta0 %*% t(beta0))
  }else{
    Cov_structure_effect <- Cov_X0 * (Cov_beta1 + Cov_beta0) +
      (Cov_beta1 + Cov_beta0) * (X0 %*% t(X0)) +
      Cov_X0 * ((beta1 - beta0) %*% t((beta1 - beta0)))
    Cov_composition_effect <-  (Cov_X1 + Cov_X0) * Cov_beta1 +
      Cov_beta1 * ((X1 - X0)%*% t((X1 - X0))) +
      (Cov_X1 + Cov_X0) * (beta1 %*% t(beta1))
  }

  Var_agg_observed_diff <-  sum(Cov_observed_diff)
  Var_agg_composition_effect <- sum(Cov_composition_effect)
  Var_agg_structure_effect <- sum(Cov_structure_effect)

  ### Old

  # Cov_X0 <- stats::cov.wt(X0, wt=weights0)$cov
  # Cov_X1 <- stats::cov.wt(X1, wt=weights1)$cov
  # Cov_X_diff <- Cov_X1 + Cov_X0
  #
  # Cov_beta_diff <- Cov_beta1 + Cov_beta0
  #
  # X0 <- apply(X0, 2, weighted.mean, w=weights0)
  # X1 <- apply(X1, 2, weighted.mean, w=weights1)
  # X_diff <- X1 - X0
  #
  # Xb0 <- X0*beta0
  # Xb1 <- X1*beta1
  #
  # beta_diff <- beta1 - beta0
  #
  # if(reference_0){
  #   XbC <- X1*beta0
  #   structure_effect <- Xb1 - XbC
  #   composition_effect <- XbC - Xb0
  #   beta_reference <- beta0
  #   X_counterfactual <- X1
  #   Cov_beta_reference <-   Cov_beta0
  #   Cov_X_counterfactual <- Cov_X1
  # }else{
  #   XbC <- X0*beta1
  #   composition_effect <- Xb1 - XbC
  #   structure_effect <- XbC - Xb0
  #   beta_reference <- beta1
  #   X_counterfactual <- X0
  #   Cov_beta_reference <-   Cov_beta1
  #   Cov_X_counterfactual <- Cov_X0
  # }
  #
  # Var_agg_observed_diff <- t(X1) %*% Cov_beta1 %*% X1 +
  #   t(beta1) %*% Cov_X1 %*% beta1 +
  #   sum(diag(Cov_X1 %*% Cov_beta1)) +
  #   t(X0) %*% Cov_beta0 %*% X0 +
  #   t(beta0) %*% Cov_X0 %*% beta0 +
  #   sum(diag(Cov_X0 %*% Cov_beta0))
  #
  # Cov_observed_diff <- Cov_beta1 * X1 %*% t(X1) +
  #   Cov_X1 * beta1 %*% t(beta1) +
  #   Cov_X1 %*% Cov_beta1 +
  #   Cov_beta0 * X0 %*% t(X0) +
  #   Cov_X0 * beta0 %*% t(beta0) +
  #   Cov_X0 %*% Cov_beta0
  #
  # Var_agg_composition_effect <- t(X_diff) %*% Cov_beta_reference %*% X_diff +
  #   t(beta_reference) %*% Cov_X_diff %*% beta_reference +
  #   sum(diag(Cov_X_diff %*% Cov_beta_reference))
  #
  # Cov_composition_effect <- Cov_beta_reference * X_diff %*% t(X_diff) +
  #   Cov_X_diff * beta_reference %*% t(beta_reference) +
  #   Cov_X_diff %*% Cov_beta_reference
  #
  # Var_agg_structure_effect <- t(X_counterfactual) %*% Cov_beta_diff %*% X_counterfactual +
  #   t(beta_diff) %*% Cov_X_counterfactual %*% beta_diff +
  #   sum(diag(Cov_X_counterfactual %*% Cov_beta_diff))
  #
  # Cov_structure_effect <- Cov_beta_diff * X_counterfactual %*% t(X_counterfactual) +
  #   Cov_X_counterfactual * beta_diff %*% t(beta_diff) +
  #   Cov_X_counterfactual %*% Cov_beta_diff

  decomposition_terms_se <- data.frame(Variable = c("Total", names(X0)),
                                       Observed_difference = sqrt(c(Var_agg_observed_diff,
                                                                    diag(Cov_observed_diff))),
                                       Composition_effect = sqrt(c(Var_agg_composition_effect,
                                                                   diag(Cov_composition_effect))),
                                       Structure_effect = sqrt(c(Var_agg_structure_effect,
                                                                 diag(Cov_structure_effect))))
  rownames(decomposition_terms_se)[1] <- "Total"

  vcov_list <- list(Observed_difference = Cov_observed_diff,
                    Composition_effect = Cov_composition_effect,
                    Structure_effect = Cov_structure_effect)

  results <- list(decomposition_terms_se = decomposition_terms_se,
                  decomposition_terms_vcov = vcov_list)
  return(results)
}
