#' Gardeazabal and Ugidos normalization of factor variables
#'


# load("data/men8305.rda")
# formula <- log(wage) ~ union + married + nonwhite + education + experience
# data <- get_all_vars(formula, men8305, weights=weights, group=year)
# data_used <- model.frame(formula, data, weights=weights, group=group)
#

# set.seed(125)
# library("AER")
# data("CPS1985")
# formula <- log(wage) ~ education + experience + union + ethnicity
# data_used <- CPS1985
# data_used$weights <- runif(nrow(CPS1985), 0.5, 1.5)
# # data_used[1,1] <- NA
# data <- get_all_vars(formula, data_used, weights=weights, group=gender)
# data_used <- model.frame(formula, data, weights=weights, group=group)



#' Perform normalization of factor variables as proposed by Gardeazabal and Ugidos (2004)
#' @param formula an object of class "formula". See \link[stats]{lm} for further details.
#' @param data a data frame containing the variables in the model.
#' @param weights numeric vector of non-negative observation weights, hence of same length as \code{dep_var}.
#' @param group name of the a binary variable (numeric or factor) identifying the two groups that will be compared.
#' @export
#' @examples
#' mod1 <- log(wage) ~ union + married + nonwhite + education + experience
#' normalized_data <- GU_normalization(formula = mod1,
#'                                    data = men8305,
#'                                    weights = weights,
#'                                    group = year)
#'
#'
GU_normalization <- function(formula, data, weights, group){

  function_call <- match.call()
  data_arguments_index = match(c("formula", "data", "weights", "group"), names(function_call), 0)
  data_arguments <- function_call[c(1, data_arguments_index)]

  # Assign formula to model.frame attributes
  data_arguments[[1]] <- as.name("model.frame")
  data_arguments$formula <- formula
  data_used <- eval.parent(data_arguments)
  function_terms <- attr(data_used, "terms")

  # Stop if there are interactions terms:
  if(sum(attr(function_terms, "order"))>length(attr(function_terms, "order"))){
    stop("GU normalization does not allow interactions.")
  }
  cat("\nFactor variables are normalized as proposed by Gardeazabal & Ugidos (2004)\n")

  term_labels <- attr(function_terms, "term.labels")


  unadjusted_regressors <- data_used[, -c(1, which(names(data_used) %in% c("(weights)","(group)")))]
  adjusted_regressors <- data.frame(matrix(nrow=nrow(unadjusted_regressors), ncol=0))


  regressors_for_prediction <- matrix(nrow=nrow(unadjusted_regressors), ncol=0)
  if(attr(function_terms, "intercept")==1){
    # regressors_for_prediction <- model.matrix(object = stats::update(formula, . ~ . - 1), data =  data)
    regressors_for_prediction <- cbind(rep(1, nrow(regressors_for_prediction)), regressors_for_prediction)
    colnames(regressors_for_prediction)[1] <- "(Intercept)"
  }
  # }else{
  #   regressors_for_prediction <- model.matrix(object = formula, data = data)
  # }
  #colnames(regressors_for_prediction) <- make.names(colnames(regressors_for_prediction))


  adjusted_coefficient_names <- list()

  for(i in 1:ncol(unadjusted_regressors)){

    regressor_i <- unadjusted_regressors[,i]
    regressor_name_i <- colnames(unadjusted_regressors)[i]

    # # if it is a factor variable with more than two levels: normalize
    # if(is.factor(regressor_i)&&length(levels(regressor_i))>2){
    if(is.factor(regressor_i)){

      formula_i <- as.formula(paste("~", regressor_name_i, "+ 0", sep=""))
      mod_matrix_i <- model.matrix(formula_i, unadjusted_regressors)
      regressors_for_prediction <- cbind(regressors_for_prediction, mod_matrix_i)

      #adjusted_coefficient_names_i <- colnames(mod_matrix_i) <- make.names(colnames(mod_matrix_i))
      adjusted_coefficient_names_i <- colnames(mod_matrix_i)
      mod_matrix_i <- mod_matrix_i[, 2:ncol(mod_matrix_i)] - mod_matrix_i[,1]
      mod_matrix_i <- as.data.frame(mod_matrix_i)
      names(mod_matrix_i) <- adjusted_coefficient_names_i[-1]
      adjusted_regressors <- cbind(adjusted_regressors, mod_matrix_i)

      term_labels <- c(term_labels[-pmatch(grep(regressor_name_i,
                                                term_labels,
                                                value=TRUE),
                                           term_labels)],
                       adjusted_coefficient_names_i[-1])


      select_to_add_parenthesis <- which(make.names(adjusted_coefficient_names_i) != adjusted_coefficient_names_i)
      adjusted_coefficient_names_i[select_to_add_parenthesis] <- paste0("`", adjusted_coefficient_names_i[select_to_add_parenthesis], "`")
      adjusted_coefficient_names[[i]] <- adjusted_coefficient_names_i
      names(adjusted_coefficient_names)[i] <- regressor_name_i

    # }else if(is.factor(regressor_i)){
    #
    #   # if it is a two level factor adjust normally with dummies
    #   formula_i <- as.formula(paste("~", regressor_name_i, sep=""))
    #   mod_matrix_i <- model.matrix(formula_i, unadjusted_regressors)
    #   adjusted_regressors <- cbind(adjusted_regressors, mod_matrix_i[,-1])
    #   colnames(adjusted_regressors)[ncol(adjusted_regressors)] <- colnames(mod_matrix_i)[2]
    #
    #   term_labels <- c(colnames(mod_matrix_i)[2],
    #                    term_labels[-pmatch(grep(regressor_name_i,
    #                                             term_labels,
    #                                             value=TRUE),
    #                                        term_labels)])

    }else{

      # if not: pass unadjusted variable thru
      adjusted_regressors <- cbind(adjusted_regressors, regressor_i)
      regressors_for_prediction <- cbind(regressors_for_prediction, regressor_i)

      colnames(adjusted_regressors)[ncol(adjusted_regressors)] <- regressor_name_i
      colnames(regressors_for_prediction)[ncol(regressors_for_prediction)] <- regressor_name_i
      names(adjusted_coefficient_names)[i] <- adjusted_coefficient_names[[i]] <- regressor_name_i

    }
  } #end loop

  adjusted_data <- as.data.frame(cbind(data[,1],
                                       adjusted_regressors,
                                       data[, which(names(data) %in% c("weights","group"))]))
  names(adjusted_data)[1] <- names(data)[1]


  adjusted_formula <-  update(formula, as.formula(paste0(". ~ ", paste0(paste0("`", term_labels, "`"), collapse=" + "))))

  if(attr(function_terms, "intercept")==1){
    adjusted_coefficient_names <- c(list(`(Intercept)`="(Intercept)"), adjusted_coefficient_names)
  }

  return(list(formula=adjusted_formula,
              data=adjusted_data,
              regressors_for_prediction=regressors_for_prediction,
              adjusted_coefficient_names=adjusted_coefficient_names))

} #end normalize GU



#' Sum coefficients for GU normalization
#'
#' After adjusting model.matrix and estimating the regression, this function sums
#' the coefficients for reference group of GU normalized factor variables.
#'
GU_normalization_sum_coefficients <- function(coef_names, est_coef){

  est_coef <- est_coef[which(names(est_coef) %in% coef_names)]
  if(length(coef_names)>1){
    coefficient_reference_group <- -sum(est_coef)
    est_coef <- c(coefficient_reference_group, est_coef)
    names(est_coef)[1] <- coef_names[1]
  }
  return(est_coef)
}


#' Get coefficients for GU normalization
#'
#' After adjusting model.matrix and estimating the regression, this function computes
#' the coefficients for reference group of factors variables for GU normalization.
#'
GU_normalization_get_coefficients <- function(coef_names, est_coef){
  est_coef <- do.call("c",lapply(coef_names,
                                 GU_normalization_sum_coefficients,
                                 est_coef = est_coef))
  names(est_coef) <- do.call("c", lapply(strsplit(names(est_coef), split="[.]"),
                                         function(x) x[-1]))
  return(est_coef)
}


#' Sum covariance matrix for GU normalization
#'
#' After adjusting model.matrix and estimating the regression, this function computes
#' the covariance matrix containing aggregate coefficients for the reference group of
#' GU normalized factor variables.
#'
GU_normalization_sum_vcov <- function(coef_names, Cov_beta){

  if(length(coef_names)>1){
    Cov_beta_adjusted <- matrix(NA, nrow=nrow(Cov_beta)+1, ncol=ncol(Cov_beta)+1)
    index_factors <- which(colnames(Cov_beta) %in% coef_names)
    ri <- range(index_factors)
    if(ri[1]!=1){
        index_lower <- 1:(ri[1]-1)
    }else{
        index_lower <- NULL
    }
    index_greater <- ri[1]:ncol(Cov_beta)

    cov_reference_group_i <- Cov_beta[, index_factors]
    cov_reference_group_i <- -rowSums(cov_reference_group_i)
    var_reference_group_i <- -sum(Cov_beta_i[index_factors])
    cov_reference_group_i <- c(cov_reference_group_i[index_lower],
                               var_reference_group_i,
                               cov_reference_group_i[index_greater])
    names(cov_reference_group_i)[ri[1]] <- coef_names[1]

    Cov_beta_adjusted[ri[1],] <- cov_reference_group_i
    Cov_beta_adjusted[,ri[1]] <- cov_reference_group_i

    Cov_beta_adjusted[index_lower, index_lower] <- Cov_beta[index_lower, index_lower]
    Cov_beta_adjusted[index_lower, index_greater+1] <- Cov_beta[index_lower, index_greater]
    Cov_beta_adjusted[index_greater+1, index_lower] <- Cov_beta[index_greater, index_lower]
    Cov_beta_adjusted[index_greater+1, index_greater+1] <- Cov_beta[index_greater, index_greater]

    rownames(Cov_beta_adjusted) <- colnames(Cov_beta_adjusted) <- names(cov_reference_group_i)
  }else{
    Cov_beta_adjusted <- Cov_beta
  }
  return(Cov_beta_adjusted)
}

#' Sum covariance matrix for GU normalization
#'
#' After adjusting model.matrix and estimating the regression, this function computes
#' the covariance matrix containing aggregate coefficients for the reference group of
#' GU normalized factor variables.
GU_normalization_get_vcov <- function(coef_names, Cov_beta){
  for(i in 1:length(coef_names)){
    Cov_beta <- GU_normalization_sum_vcov(coef_names[[i]], Cov_beta)
  }
  return(Cov_beta)
}
