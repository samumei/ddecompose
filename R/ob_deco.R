


library("AER")
data("CPS1985")
mod <- log(wage) ~ education + experience + union + ethnicity


# oaxaca package ---------------------------------------------------------------
library("oaxaca")
?oaxaca
data("chicago")
oaxaca.results.1 <- oaxaca(ln.real.wage ~ age + female + LTHS + some.college + 
                             college + advanced.degree | foreign.born, 
                           data = chicago, R = 30)
print(oaxaca.results.1)
plot(oaxaca.results.1)

# GUnormalization

# Aproach I: Normalize coefficient after estimation 
#- downside: estimating vcov is mess
# 1, Estimate coefficients
# 2, Create empty model.matrix()
# 3. Check for every variable if factor
#    if factor: 
#    create model.matrix without intercept, i.e. model.matrix(dep_var ~ factor - 1, data)
#    normalize coefficient (incl. adjusting intercept), extend coefficient vector

# Approach II: Normalize factors/dummies before estimation
#- downside: you need `newdata` in est_ob_decomposition function.
#            i.e. normalized model.frame for estimation, plus de-normalized normalized model.frame
#- positive: function exists already, 
#- positive: vcov is drictly estimated

adjust.dummies <- function(beta, dummies) {
  out <- beta
  c <- sum(beta[dummies])/(length(dummies) + 1)
  out[dummies] <- beta[dummies] - c
  out[1] <- beta[1] + c
  return(out)
}

mf <- model.frame(mod, CPS1985)
mt <- attr(mf, "terms")

names_of_regressors <- attr(mt, "term.labels")


set.seed(125)
CPS1985$weights <- runif(nrow(CPS1985), 0.5,1.5)

data0 <- subset(CPS1985, gender=="female")
fit0 <- lm(mod, data=data0, weights=weights)

data1 <- subset(CPS1985, gender=="male")
fit1 <- lm(mod, data=data1, weights=weights)

mod2 <- update(mod, . ~ (.)*gender )
fit2 <- lm(mod2, CPS1985 )

X0 <- apply(model.matrix(mod, fit0$model), 2, weighted.mean, w=fit0$model[,"(weights)"])
X1 <- apply(model.matrix(mod, fit1$model), 2, weighted.mean, w=fit1$model[,"(weights)"])

beta0 <- coef(fit0)
beta1 <- coef(fit1)

Xb0 <- X0*beta0
Xb1 <- X1*beta1
XbC <- X1*beta0

observed_diff <- Xb1 - Xb0
composition_effect <- XbC - Xb0
structure_effect <- Xb1 - XbC

agg_observed_diff <- sum(observed_diff)
agg_composition_effect <- sum(composition_effect)
agg_structure_effect <- sum(structure_effect)


data.frame(variable = c("Totel", names(observed_diff)), 
           observed_diff = c(agg_observed_diff, observed_diff), 
           composition_effect = c(agg_composition_effect, composition_effect), 
           structure_effect = c(agg_structure_effect, structure_effect))



weighted.mean(model.matrix(mod, fit0$model)  %*% coef(fit0), w=fit0$model[,"(weights)"])  
weighted.mean(x=fit0$model[,1], w=fit0$model[,"(weights)"])


# Standardfehler
est_ob_decomposition(fit0, fit1)
est_ob_vcov(fit0, fit1)$decomposition_terms_se



# ob_deco function -------------------------------------------------------------

data0 <- subset(CPS1985, gender=="female")
fit0 <- lm(mod, data=data0, weights=weights)

data1 <- subset(CPS1985, gender=="male")
fit1 <- lm(mod, data=data1, weights=weights)



# Start actual function -------------------------------------------------------


# Data 
set.seed(125)
library("AER")
data("CPS1985")
formula <- log(wage) ~ education + experience + union + ethnicity 
data_used <- CPS1985
data_used$weights <- runif(nrow(CPS1985), 0.5, 1.5)
# data_used[1,1] <- NA
data_used <- get_all_vars(formula, data_used, weights=weights, group=gender)
# cluster <- data_used$cluster <- rbinom(nrow(CPS1985), 10, 0.5) 
cluster <- NULL

#data_used <- model.frame(formula, data_used, weights=weights, group=gender)


na.action = na.exclude
reference_0 <- TRUE
vcov=stats::vcov
normalize_factors=FALSE
bootstrap_iterations = 100
# function


#' @example 
#' set.seed(125)
#' library("AER")
#' data("CPS1985")
#' formula <- log(wage) ~ education + experience + union + ethnicity 
#' data_used <- CPS1985
#' data_used$weights <- runif(nrow(CPS1985), 0.5,1.5)
#' ob_deco(formula, data_used, group=gender, weights=weights)
#' 
#' ob_deco(formula, data_used, group=gender, weights=weights, bootstrap = TRUE )
#' 
#' 
ob_deco <- function(formula, 
                    data, 
                    weights, 
                    na.action = na.exclude,
                    group,
                    reference_0=TRUE,
                    normalize_factors=FALSE,
                    bootstrap=FALSE,
                    bootstrap_iterations=100, 
                    bootstrap_robust=FALSE,
                    cluster=NULL,
                    cores=1,
                    vcov=stats::vcov){
  
  ## Get model.frame
  function_call <- match.call()
  data_arguments_index <- match(c("formula", "data", "weights", "na.action", "group"), names(function_call), 0)
  data_arguments <- function_call[c(1, data_arguments_index)]
  data_arguments$drop.unused.levels <- TRUE
  
  data_arguments[[1]] <- as.name("get_all_vars") #as.name("model.frame")
  data_used <- eval.parent(data_arguments)
  data_used <- lapply(list(data_used), na.action)[[1]]
  
  ## Get weights
  if (!is.null(data_used[, "weights"]) && !is.numeric(data_used[, "weights"])) {
    stop("'weights' must be a numeric vector")
  }
  if (is.null(data_used[, "weights"])) {
    data_used[, "weights"] <- rep(1, nrow(data_used))
  }
  
  ## Check group variable
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
  
  estimated_decomposition <- estimate_ob_deco(formula = formula,
                                              data_used = data_used,
                                              reference_0 = reference_0,
                                              normalize_factors = normalize_factors,
                                              compute_analytical_se=TRUE,
                                              return_model_fit = TRUE,
                                              vcov=vcov)
  
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
    names(bootstrap_estimates) <- gsub("value[.]","", names(bootstrap_estimates))
    bootstrap_estimates <- lapply(split(bootstrap_estimates, bootstrap_estimates$effect), function(x) stats::cov(x[, -c(1:2)]))
  
    decomposition_terms_se <- as.data.frame(do.call("cbind", lapply(bootstrap_estimates, function(x) sqrt(diag(x)))))
    decomposition_terms_se$Variable <- rownames(decomposition_terms_se)
    decomposition_terms_se <- decomposition_terms_se[, c("Variable", "Observed_difference", "Composition_effect", "Structure_effect")]
    
    decomposition_terms_vcov <- lapply(bootstrap_estimates, function(x) x[-1,-1])
    
    estimated_decomposition$decomposition_vcov$decomposition_terms_se <- decomposition_terms_se
    estimated_decomposition$decomposition_vcov$decomposition_terms_vcov <- decomposition_terms_vcov 
  }
  
  class(estimated_decomposition) <- "ob_deco"
  return(estimated_decomposition)
  
}



#' Estimate OB decomposition
#'  
estimate_ob_deco <- function(formula, 
                             data_used,
                             reference_0=TRUE,
                             normalize_factors=FALSE,
                             compute_analytical_se=TRUE,
                             return_model_fits = TRUE,
                             vcov=stats::vcov){
  
  group0 <- levels(data_used[, "group"])[1]
  
  obs_0 <- which(data_used[, "group"] == group0)
  obs_1 <- which(data_used[, "group"] != group0)
  
  weights0 <- data_used[obs_0, "weights"]
  weights1 <- data_used[obs_1, "weights"]
  
  if(normalize_factors){
    normalized_data <- GU_normalization(formula=formula,
                                        data=data,
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
  }
    
  fit0 <- lm(formula, data = data_used, subset = group == group0, weights = weights)
  fit1 <- lm(formula, data = data_used, subset = group != group0, weights = weights)
  
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
                                                  reference_0 = TRUE)
 
  if(compute_analytical_se){ 
    Cov_beta0 <- lapply(list(fit0), vcov)[[1]]
    Cov_beta1 <- lapply(list(fit1), vcov)[[1]]
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
                                                  reference_0 = TRUE)
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
                 model_fits = model_fits)
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
  data_used$weights <- data_used$weights * sum(data_used[ ,"weights"], na.rm=TRUE)/sum(data_used[sampled_observations,"weights"], na.rm=TRUE)
  }
  
  deco_estimates <- estimate_ob_deco(formula = formula,
                                     data_used = data_used[sampled_observations, ],
                                     reference_0 = reference_0,
                                     normalize_factors = normalize_factors,
                                     compute_analytical_se=FALSE,
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
                                    reference_0=TRUE){
  
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
                                    reference_0=TRUE){
  
  
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

