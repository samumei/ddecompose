#' summary method for class "dfl_deco"
#'
#' @param x an object of class "dfl_deco", a result of a call to [dfl_deco()].
#' @param confidence_level numeric value between 0 and 1 (default = 0.95) that defines the
#'              confidence level of the printed confidence intervals. Pointwise confidences bands
#'              are defined as \code{qnorm((1-confidence_level)/2)} * standard error. Uniform bands
#'              are constructed by multiplying the standard error with \code{confidence_level}-quantile
#'              of the bootstrapped Kolmogorov-Smirnov statistic.
#' @param digits number of digits to be printed.
#' @param ... other parameters to be passed through to printing functions.
#' @return The function \code{summary.dfl_deco()} displays the decompositions
#' terms save in \code{x}. If standard errors have been bootstrapped, standard
#' errors and confidence bands are given. Uniform confidence bands for quantiles
#' are constructed based on the bootstrappend Kolmogorov-Smirnov distribution.
#'
#' @export
#'
summary.dfl_deco <- function(x, confidence_level=0.95, digits=4, ...){

  if(x$subtract_1_from_0 == FALSE){
  cat("Decomposition of difference between",
      paste0(x$group_variable_name, " == '",x$group_variable_levels[2],"'"),
      "(group 1) and\n",
      paste0(x$group_variable_name, " == '",x$group_variable_levels[1],"'"),
      "(group 0)\n\n")
  }else{
    cat("Decomposition of difference between",
        paste0(x$group_variable_name, " == '",x$group_variable_levels[1],"'"),
        "(group 0) and\n",
        paste0(x$group_variable_name, " == '",x$group_variable_levels[2],"'"),
        "(group 1)\n\n")
  }

  cat("Reweighted reference group:",  paste0(x$group_variable_name, " == '", x$reference_group,"' (group ", ifelse(x$group_variable_levels[1]==x$reference_group,0,1), ")"), "\n \n")

  if(length(x$covariates_labels) == 1){
    cat("Composition effect accounts for between-group differences\nin the distribution of the following covariates:\n\n")
  }else{
    cat("Composition effects of the sequential decomposition account \nfor between-group differences in the distribution of the\nfollowing covariates:\n\n")
  }
  for(i in 1:length(x$covariates_labels)){
    cat(x$covariates_labels[[i]], "\n")
  }
  cat("\n")
  cat("---------------------------------------------------------------------------------\n")

  if(is.null(x$decomposition_quantiles)==FALSE){
    cat("Decomposition of difference at conditional quantiles:\n\n")
    decomposition_quantiles <- x$decomposition_quantiles

    if(is.null(x$bootstrapped_standard_errors)==FALSE){
      decomposition_quantiles_se <- x$bootstrapped_standard_errors$decomposition_quantiles
      kolmogorov_smirnov_stat <- x$bootstrapped_standard_errors$decomposition_quantiles_kms_distribution
      kolmogorov_smirnov_stat <- lapply(split(kolmogorov_smirnov_stat, kolmogorov_smirnov_stat$effect),
                                        function(x) data.frame(effect = x$effect[1],
                                                               t_value = quantile(x$kms_t_value, confidence_level)))
      kolmogorov_smirnov_stat <- do.call("rbind", kolmogorov_smirnov_stat)
      # kolmogorov_smirnov_stat <- dplyr::summarise(dplyr::group_by(x$bootstrapped_standard_errors$decomposition_quantiles_kms_distribution,
      #                                                             effect),
      #                                             t_value = quantile(kms_t_value, confidence_level))
      kolmogorov_smirnov_stat <- as.data.frame(kolmogorov_smirnov_stat[match(names(decomposition_quantiles), kolmogorov_smirnov_stat$effect), ])

      for(i in 2:ncol(decomposition_quantiles)){
        cat(paste0(names(decomposition_quantiles)[i],":"), "\n")
        cat("---------------------------------------------------------------------------------\n")
        results_to_print <- as.data.frame(cbind(decomposition_quantiles[,c(1, i)],
                                                decomposition_quantiles_se[,i]))
        results_to_print$ci_p_low <- results_to_print[,2] - qnorm(1-(1-confidence_level)/2)*results_to_print[,3]
        results_to_print$ci_p_high <- results_to_print[,2] + qnorm(1-(1-confidence_level)/2)*results_to_print[,3]
        results_to_print$ci_u_low <- results_to_print[,2] - kolmogorov_smirnov_stat[i,2]*results_to_print[,3]
        results_to_print$ci_u_high <- results_to_print[,2] + kolmogorov_smirnov_stat[i,2]*results_to_print[,3]

        names(results_to_print) <- c("Quantile", "Estimate", "Std. Error","Pointwise CI: [low", "high]", "Uniform CI: [low", "high]")
        print(results_to_print, digits = digits, row.names = FALSE)
        cat("\n")
      }

    }else{
      print(decomposition_quantiles[,-1])
      cat("\n")
      cat("---------------------------------------------------------------------------------\n")
    }

  }
  if(is.null(x$decomposition_other_statistics)==FALSE){
    cat("Decomposition of difference for other distributional statistics\n\n")
    decomposition_other_statistics  <- x$decomposition_other_statistics
    if(is.null(x$bootstrapped_standard_errors)==FALSE){
      decomposition_other_statistics_se <- x$bootstrapped_standard_errors$decomposition_other_statistics

      for(i in 2:ncol(decomposition_other_statistics)){
        cat(paste0(names(decomposition_other_statistics)[i],":"), "\n")
        cat("---------------------------------------------------------------------------------\n")
        results_to_print <- as.data.frame(cbind(decomposition_other_statistics[,c(1, i)],
                                                decomposition_other_statistics_se[,i]))
        results_to_print$ci_p_low <- results_to_print[,2] - qnorm(1-(1-confidence_level)/2)*results_to_print[,3]
        results_to_print$ci_p_high <- results_to_print[,2] + qnorm(1-(1-confidence_level)/2)*results_to_print[,3]

        names(results_to_print) <- c("Statistic", "Estimate", "Std. Error","CI [low", "high]")
        print(results_to_print, digits = digits, row.names = FALSE)
        cat("\n")
      }

    }else{
      print(x$decomposition_other_statistics[, -1])
      cat("\n")
      cat("---------------------------------------------------------------------------------\n")
    }
  }

  cat("Summary statistics of reweighting factors\n\n")

  cat(paste0("Number of trimmed observations (not included in statistics): ",
             length(x$trimmed_observations),
             " (",
             round(length(x$trimmed_observations)/nrow(x$reweighting_factor)*100, 1),
             "%)\n\n"))

  if(is.null(x$bootstrapped_standard_errors)==FALSE){
    quantiles_reweighting_factor <- x$quantiles_reweighting_factor
    quantiles_reweighting_factor_se <- x$bootstrapped_standard_errors$quantiles_reweighting_factor
    for(i in 2:ncol(quantiles_reweighting_factor_se)){
      cat(paste0(names(quantiles_reweighting_factor)[i],":"), "\n")
      cat("---------------------------------------------------------------------------------\n")
      results_to_print <- cbind(quantiles_reweighting_factor[,c(1,i)], quantiles_reweighting_factor_se[,c(1,i)])[,c(2,4)]
      names(results_to_print) <- c("Estimate","Std. Error")
      print(results_to_print, digits = digits)
      cat("\n")
    }

  }else{
    quantiles_reweighting_factor <- as.data.frame(x$quantiles_reweighting_factor[,-1])
    names(quantiles_reweighting_factor) <- names(x$quantiles_reweighting_factor)[-1]
    rownames(quantiles_reweighting_factor ) <- rownames(x$quantiles_reweighting_factor)
    print(quantiles_reweighting_factor)
    cat("\n")
  }
}


#' summary method for class "ob_deco"
#'
#' Apart from displaying the (detailed) decomposition results with standard
#' errors, \code{summary.ob_deco()} allows to customize the aggregation of the
#' detailed decomposition terms.
#'
#' @param x an object of class "ob_deco", usually , a result of a call to [ob_deco()].
#' @param aggregate_factors boolean, if `TRUE` (default) terms associated with detailed factor
#' levels are aggregated to a single term for every factor variable.
#' @param custom_aggregation list specifying the aggregation of detailed decomposition
#' terms. The parameter `custom_aggregation` overrides the parameter `aggregate_factors`.
#' If `NULL` (default), then either all detailed terms or all terms associated with
#' a single variable are returned.
#' @param confidence_level numeric value between 0 and 1 (default = 0.95) that defines the printed confidence interval.
#' @param ... other parameters to be passed through to print functions.
#'
#' @return The function \code{summary.ob_deco()} summarizes the decompositions terms saved in \code{x}.
#'
#' @export
#'
#' @examples
#' data("nlys00")
#' mod1 <- log(wage) ~ age + education + years_worked_civilian +
#'   years_worked_military + part_time + industry
#'
#' deco_results <- ob_deco(formula = mod1,
#'                         data = nlys00,
#'                         group = female,
#'                         reference_0 = TRUE)
#'
#' # Print standard errors
#' summary(deco_results)
#'
#' # Aggregate decomposition terms associated with factor levels
#' summary(deco_results, aggregate_factors = TRUE)
#'
#' # custom aggregation of decomposition terms
#' custom_aggregation <-
#'   list(`Age` = c("age"),
#'        `Education` = c("education<10 yrs",
#'                        "educationHS grad (diploma)",
#'                        "educationHS grad (GED)",
#'                        "educationSome college",
#'                        "educationBA or equiv. degree",
#'                        "educationMA or equiv. degree",
#'                        "educationPh.D or prof. degree"),
#'        `Life-time work experience` = c("years_worked_civilian",
#'                                        "years_worked_military",
#'                                        "part_time"),
#'        `Industrial sectors` = c("industryManufacturing",
#'                                 "industryEducation, Health, Public Admin.",
#'                                 "industryOther services"))
#' summary(deco_results, custom_aggregation = custom_aggregation)
#'
summary.ob_deco <- function(x,
                            aggregate_factors = TRUE,
                            custom_aggregation = NULL,
                            confidence_level = 0.95,
                            ...){

  reweighting <- ifelse(x$input_parameters$reweighting, TRUE, FALSE)

  if(!x$input_parameters$rifreg) {
    if(!reweighting) {
      decomposition_type <- "\n\nOaxaca-Blinder decomposition of mean difference\nbetween"
    }
    else{
      decomposition_type <- "\n\n'Doubly robust' Oaxaca-Blinder decomposition of mean difference\nbetween"
    }
  }
  else {
    if(!reweighting) {
      decomposition_type <- paste0("\n\nRIF regression decomposition of difference in ",
                                   x$input_parameters$rifreg_statistic  ,
                                   "\nbetween")
    }
    else{
      decomposition_type <- paste0("\n\nReweighted RIF regression decomposition of difference in ",
                                   x$input_parameters$rifreg_statistic  ,
                                   "\nbetween")
    }
  }
  cat(decomposition_type,
      paste0(x$group_variable_name, " == '", x$group_variable_levels[2], "'"),
      "(group 1) and",
      paste0(x$group_variable_name, " == '", x$group_variable_levels[1], "'"),
      "(group 0). \nThe reference group is", paste0("'",x$reference_group,"'."), "\n\n")


  if(!reweighting) {


    cat("Group 0:", paste0(x$group_variable_name, " == '", x$group_variable_levels[1], "'"),
        paste0("(",length(x[[1]]$model_fits$fit_group_0$residuals)," observations)"),
        "\nGroup 1:", paste0(x$group_variable_name, " == '", x$group_variable_levels[2], "'"),
        paste0("(",length(x[[1]]$model_fits$fit_group_1$residuals)," observations)"),"\n\n")

    if(x$input_parameters$reference_0) {
      x_reference <- "X1"
      b_reference <- "b0"
    }
    else {
      x_reference <- "X0"
      b_reference <- "b1"
    }
    if(!x$input_parameters$subtract_1_from_0) {
      x_subtraction <- "(X1 - X0)"
      b_subtraction <- "(b1 - b0)"
    }
    else {
      x_subtraction <- "(X0 - X1)"
      b_subtraction <- "(b0 - b1)"
    }

    composition_formula <- paste("Composition Effect:", x_subtraction, "*", b_reference)
    structure_formula <- paste(" Structure Effect:", x_reference, "*", b_subtraction)

    cat(composition_formula, "\n",
        structure_formula, "\n\n")
  }
  else {
    # With Reweighting
    cat("Group 0:", paste0(x$group_variable_name, " == '", x$group_variable_levels[1], "'"),
        paste0("(",length(x[[1]]$model_fits$fit_group_0$residuals)," observations)"),
        "\nGroup 1:", paste0(x$group_variable_name, " == '", x$group_variable_levels[2], "'"),
        paste0("(",length(x[[1]]$model_fits$fit_group_1$residuals)," observations)"),
        "\nGroup C: The reweighted group C is group", paste0("'",x$reference_group,"'"), "(reference group) reweighted
         to match the characteristics of the other group", paste0("(",length(x[[1]]$model_fits$fit_group_reweighted$residuals)," observations).\n\n"))

    if(x$input_parameters$reference_0) {
      x_p_reference <- "b0"
      x_e_reference <- "XC"
      s_p_reference <- "X1"
      s_e_reference <- "bC"
    }
    else {
      x_p_reference <- "bC"
      x_e_reference <- "X0"
      s_p_reference <- "XC"
      s_e_reference <- "b1"
    }
    if(!x$input_parameters$subtract_1_from_0) {
      x_p_subtraction <- "(XC - X0)"
      x_e_subtraction <- "(bC - b0)"
      s_p_subtraction <- "(b1 - bC)"
      s_e_subtraction <- "(X1 - XC)"
    }
    else {
      x_p_subtraction <- "(X0 - XC)"
      x_e_subtraction <- "(b0 - bC)"
      s_p_subtraction <- "(bC - b1)"
      s_e_subtraction <- "(XC - X1)"
    }


    composition_formula <- paste("Pure Composition Effect:", x_p_subtraction, "*", x_p_reference)
    structure_formula <- paste(" Pure Structure Effect:", x_e_reference, "*", x_e_subtraction)
    spec_err_formula <- paste("   Specification Error:", s_e_subtraction, "*", s_e_reference)
    rw_err_formula <- paste("     Reweighting Error:", s_p_reference, "*", s_p_subtraction)

    cat(composition_formula, "\n",
        structure_formula, "\n",
        spec_err_formula, "\n",
        rw_err_formula, "\n\n")
  }

  n_decompositions <- length(x) - 5

  for(i in 1:n_decompositions) {

    if(x$input_parameters$rifreg &
       x$input_parameters$rifreg_statistic == "quantiles") {
      cat("\n*** Quantile:",  x$input_parameters$rifreg_probs[i], "***")
      cat("\n\n")
    }

    if(aggregate_factors | !is.null(custom_aggregation)){
      results <- aggregate_terms(x[[i]],
                           aggregate_factors = aggregate_factors,
                           custom_aggregation = custom_aggregation,
                           reweighting = reweighting)
    }
    else {
      results <- x[[i]][1:4]

      no_se <- ifelse(is.null(results$decomposition_vcov), TRUE, FALSE)

      if(no_se) {
        results$decomposition_vcov$decomposition_terms_se <- results$decomposition_terms
        results$decomposition_vcov$decomposition_terms_vcov <- results$decomposition_terms[, -1]
        results$decomposition_vcov$decomposition_terms_se[] <- NA
        results$decomposition_vcov$decomposition_terms_vcov[] <- NA
      }

      if(length(results[["decomposition_vcov"]][["decomposition_terms_vcov"]]) == 3){ # no reweighting
        results$decomposition_terms <- results$decomposition_terms[, 1:4]
      }
    }

    decomposition_terms <- results$decomposition_terms[,-1]
    decomposition_terms_se <- results$decomposition_vcov$decomposition_terms_se[,-1]
    names(decomposition_terms) <- gsub("_", " ", names(decomposition_terms))
    aggregate_decomposition <- decomposition_terms[1, ]
    aggregate_decomposition_se <- decomposition_terms_se[1, ]
    detailed_decomposition <-  decomposition_terms[-1, ]
    detailed_decomposition_se <-  decomposition_terms_se[-1, ]

    aggregate_decomposition <- data.frame(Effect = names(aggregate_decomposition),
                                          Estimate = as.numeric(aggregate_decomposition[1, ]),
                                          se = as.numeric(aggregate_decomposition_se[1, ]))
    aggregate_decomposition$low <-  aggregate_decomposition$Estimate - aggregate_decomposition$se * qnorm(1 - (1 - confidence_level)/2)
    aggregate_decomposition$high <-  aggregate_decomposition$Estimate + aggregate_decomposition$se * qnorm(1 - (1 - confidence_level)/2)
    names(aggregate_decomposition) <- c("Effect", "Estimate", "Std. Error", "CI [Low", "High]")
    rownames(aggregate_decomposition) <- aggregate_decomposition$Effect
    aggregate_decomposition$Effect <- NULL

    detailed_decomposition_observed <-  data.frame(Estimate = detailed_decomposition[, c(which(names(detailed_decomposition) == "Observed difference"))],
                                                   se = detailed_decomposition_se[, c(which(names(detailed_decomposition) == "Observed difference"))])
    detailed_decomposition_composition <- data.frame(Estimate = detailed_decomposition[, c(which(names(detailed_decomposition) == "Composition effect"))],
                                                     se = detailed_decomposition_se[, c(which(names(detailed_decomposition) == "Composition effect"))])
    detailed_decomposition_structure <-  data.frame(Estimate = detailed_decomposition[, c(which(names(detailed_decomposition) == "Structure effect"))],
                                                    se = detailed_decomposition_se[, c(which(names(detailed_decomposition) == "Structure effect"))])

    detailed_decomposition_observed$low <-  detailed_decomposition_observed$Estimate - detailed_decomposition_observed$se * qnorm(1 - (1 - confidence_level)/2)
    detailed_decomposition_observed$high <-  detailed_decomposition_observed$Estimate + detailed_decomposition_observed$se * qnorm(1 - (1 - confidence_level)/2)

    detailed_decomposition_composition$low <-  detailed_decomposition_composition$Estimate - detailed_decomposition_composition$se * qnorm(1 - (1 - confidence_level)/2)
    detailed_decomposition_composition$high <-  detailed_decomposition_composition$Estimate + detailed_decomposition_composition$se * qnorm(1 - (1 - confidence_level)/2)

    detailed_decomposition_structure$low <-  detailed_decomposition_structure$Estimate - detailed_decomposition_structure$se * qnorm(1 - (1 - confidence_level)/2)
    detailed_decomposition_structure$high <-  detailed_decomposition_structure$Estimate + detailed_decomposition_structure$se * qnorm(1 - (1 - confidence_level)/2)

    names(detailed_decomposition_observed) <- names(detailed_decomposition_structure) <- names(detailed_decomposition_composition) <- c("Estimate", "Std. Error", "CI [Low", "High]")
    rownames(detailed_decomposition_observed) <-  rownames(detailed_decomposition_structure) <-  rownames(detailed_decomposition_composition) <- rownames(detailed_decomposition)


    if(reweighting) {
      detailed_decomposition_spec_error <- data.frame(Estimate = detailed_decomposition[, c(which(names(detailed_decomposition) == "Specification error"))],
                                                      se = detailed_decomposition_se[, c(which(names(detailed_decomposition) == "Specification error"))])
      detailed_decomposition_rw_error <-  data.frame(Estimate = detailed_decomposition[, c(which(names(detailed_decomposition) == "Reweighting error"))],
                                                     se = detailed_decomposition_se[, c(which(names(detailed_decomposition) == "Reweighting error"))])

      detailed_decomposition_spec_error$low <-  detailed_decomposition_spec_error$Estimate - detailed_decomposition_spec_error$se * qnorm(1 - (1 - confidence_level)/2)
      detailed_decomposition_spec_error$high <-  detailed_decomposition_spec_error$Estimate + detailed_decomposition_spec_error$se * qnorm(1 - (1 - confidence_level)/2)

      detailed_decomposition_rw_error$low <-  detailed_decomposition_rw_error$Estimate - detailed_decomposition_rw_error$se * qnorm(1 - (1 - confidence_level)/2)
      detailed_decomposition_rw_error$high <-  detailed_decomposition_rw_error$Estimate + detailed_decomposition_rw_error$se * qnorm(1 - (1 - confidence_level)/2)

      names(detailed_decomposition_spec_error) <- names(detailed_decomposition_rw_error) <- names(detailed_decomposition_observed)
      rownames(detailed_decomposition_spec_error) <- rownames(detailed_decomposition_rw_error) <- rownames(detailed_decomposition)

    }


    #rownames(aggregate_decomposition) <- paste0("Total difference ", paste0(rep(" ",  max(nchar(rownames(detailed_decomposition)))-nchar("Total difference ")), collapse=""))
    cat("Aggregate decomposition:\n\n")
    print(aggregate_decomposition, ...)
    cat("\n")
    cat("\n")
    cat("Observed difference:\n\n")
    print(detailed_decomposition_observed, ...)
    cat("\n")
    cat("\n")

    if(!reweighting) {
      cat("Structure effect:\n\n")
      print(detailed_decomposition_structure, ...)
      cat("\n")
      cat("\n")
      cat("Composition effect:\n\n")
      print(detailed_decomposition_composition, ...)
      cat("\n")
    }
    else {
      cat("Pure structure effect:\n\n")
      print(detailed_decomposition_structure, ...)
      cat("\n")
      cat("\n")
      cat("Pure composition effect:\n\n")
      print(detailed_decomposition_composition, ...)
      cat("\n")
      cat("\n")
      cat("Specification error:\n\n")
      print(detailed_decomposition_spec_error, ...)
      cat("\n")
      cat("\n")
      cat("Reweighting error:\n\n")
      print(detailed_decomposition_rw_error, ...)
      cat("\n")
    }
  }
}



#' Aggregate decomposition terms
#'
#' The function aggregates decomposition terms and calculates
#' their covariance matrix based on detailed decomposition results.
#'
#' @param x an object of class "ob_deco", usually , a result of a call to [ob_deco()].
#' @param aggregate_factors boolean, if `TRUE` (default) terms associated with detailed factor
#' levels are aggregated to a single term for every factor variable.
#' @param custom_aggregation list specifying the aggregation of detailed decomposition
#' terms. The parameter `custom_aggregation` overrides the parameter `aggregate_factors`.
#' If `NULL` (default), then either all detailed terms or all terms associated with
#' a single variable are returned.
#' @param reweighting boolean, if `TRUE` the decompostion in `x` contains reweighting
#' (i.e. specification and reweighting error)
#'
#' @return The function returns an updated object of class "ob_deco" containing
#' the aggregated decomposition terms.
#'
aggregate_terms <- function(x,
                            aggregate_factors = TRUE,
                            custom_aggregation = NULL,
                            reweighting){



  no_se <- ifelse(is.null(x$decomposition_vcov), TRUE, FALSE)

  if(no_se) {
    x$decomposition_vcov$decomposition_terms_se <- x$decomposition_terms
    x$decomposition_vcov$decomposition_terms_vcov <- x$decomposition_terms[, -1]
    x$decomposition_vcov$decomposition_terms_se[] <- NA
    x$decomposition_vcov$decomposition_terms_vcov[] <- NA

    decomposition_vcov <- x$decomposition_vcov
  }
  else{
    decomposition_vcov <- x$decomposition_vcov
  }

  decomposition_terms <- x$decomposition_terms


    if(aggregate_factors & is.null(custom_aggregation)){

      if(is.null(x$GU_normalized_coefficient_names)){

        model_variables <- all.vars(x$model_fits[[1]]$terms)[-1]
        if(names(x$model_fits[[1]]$coefficients[1])=="(Intercept)") {
          model_variables <- c("(Intercept)", model_variables)
        }

        factor_levels <- x$model_fits[[1]]$xlevels
        number_of_factors <- length(factor_levels)
        factor_variables <- names(factor_levels)

        custom_aggregation <- list()
        for(i in 1:length(model_variables)){
          sel_factor_variable <- match(model_variables[i], factor_variables)
          if(is.na(sel_factor_variable) == FALSE){
            custom_aggregation[[i]] <- paste0(model_variables[i], factor_levels[[sel_factor_variable]])[-1]
          }else{
            custom_aggregation[[i]] <-  model_variables[i]
          }
          names(custom_aggregation)[i] <-  model_variables[i]
        }
      }
      else{
        custom_aggregation <- x$GU_normalized_coefficient_names
      }

    }else if(!is.null(custom_aggregation)){

      missing_variables <- setdiff(do.call("c", custom_aggregation), decomposition_terms$Variable[-1])
      if(length(missing_variables) == 1){
        stop(paste0("Cannot aggregate terms. A variable (", missing_variables, ") is not defined."))
      }else if(length(missing_variables) > 1){
        stop(paste0("Cannot aggregate terms. Some variables (", paste0(missing_variables, collapse=", ") ,") are not defined."))
      }
      other_variables <- setdiff(decomposition_terms$Variable[-1], do.call("c", custom_aggregation))
      if(length(other_variables)>0){
        custom_aggregation <- c(custom_aggregation, list(other_variables))
        names(custom_aggregation)[length(custom_aggregation)] <- "(Other variables)"
      }

    }

    # Aggregate terms and vcov
    aggregated_terms <- decomposition_terms[1, ]
    aggregated_vcov_Observed_difference <- decomposition_vcov$decomposition_terms_vcov$Observed_difference
    aggregated_vcov_Composition_effect <- decomposition_vcov$decomposition_terms_vcov$Composition_effect
    aggregated_vcov_Structure_effect <- decomposition_vcov$decomposition_terms_vcov$Structure_effect

    if(reweighting) {
      aggregated_vcov_Spec_error <- decomposition_vcov$decomposition_terms_vcov$Specification_error
      aggregated_vcov_RW_error <- decomposition_vcov$decomposition_terms_vcov$Reweighting_error
    }

    for(i in 1:length(custom_aggregation)){

      # Decomposition terms
      add <- subset(decomposition_terms, Variable %in% custom_aggregation[[i]])
      add[1, -1] <- colSums(add[,-1])
      add[1, "Variable"] <- names(custom_aggregation)[i]
      add <- add[1, ]
      rownames(add) <- names(custom_aggregation)[i]
      aggregated_terms <- rbind(aggregated_terms, add)

      # vcov
      if(!no_se) {
        sel_terms <- which(colnames(aggregated_vcov_Observed_difference) %in% custom_aggregation[[i]])
        if(length(sel_terms) > 1){
          aggregated_vcov_Observed_difference[sel_terms[1], ] <- colSums(aggregated_vcov_Observed_difference[sel_terms, ])
          aggregated_vcov_Observed_difference[, sel_terms[1]] <- rowSums(aggregated_vcov_Observed_difference[, sel_terms])
          aggregated_vcov_Observed_difference <- aggregated_vcov_Observed_difference[-setdiff(sel_terms, sel_terms[1]),-setdiff(sel_terms, sel_terms[1])]

          aggregated_vcov_Composition_effect[sel_terms[1], ] <- colSums(aggregated_vcov_Composition_effect[sel_terms, ])
          aggregated_vcov_Composition_effect[, sel_terms[1]] <- rowSums(aggregated_vcov_Composition_effect[, sel_terms])
          aggregated_vcov_Composition_effect <- aggregated_vcov_Composition_effect [-setdiff(sel_terms, sel_terms[1]),-setdiff(sel_terms, sel_terms[1])]

          aggregated_vcov_Structure_effect[sel_terms[1], ] <- colSums(aggregated_vcov_Structure_effect[sel_terms, ])
          aggregated_vcov_Structure_effect[, sel_terms[1]] <- rowSums(aggregated_vcov_Structure_effect[, sel_terms])
          aggregated_vcov_Structure_effect <- aggregated_vcov_Structure_effect [-setdiff(sel_terms, sel_terms[1]),-setdiff(sel_terms, sel_terms[1])]

          if(reweighting) {
            aggregated_vcov_Spec_error[sel_terms[1], ] <- colSums(aggregated_vcov_Spec_error[sel_terms, ])
            aggregated_vcov_Spec_error[, sel_terms[1]] <- rowSums(aggregated_vcov_Spec_error[, sel_terms])
            aggregated_vcov_Spec_error <- aggregated_vcov_Spec_error [-setdiff(sel_terms, sel_terms[1]),-setdiff(sel_terms, sel_terms[1])]

            aggregated_vcov_RW_error[sel_terms[1], ] <- colSums(aggregated_vcov_RW_error[sel_terms, ])
            aggregated_vcov_RW_error[, sel_terms[1]] <- rowSums(aggregated_vcov_RW_error[, sel_terms])
            aggregated_vcov_RW_error <- aggregated_vcov_RW_error [-setdiff(sel_terms, sel_terms[1]),-setdiff(sel_terms, sel_terms[1])]
          }


        }

        colnames(aggregated_vcov_Observed_difference)[sel_terms[1]] <- names(custom_aggregation)[i]
        rownames(aggregated_vcov_Observed_difference)[sel_terms[1]] <- names(custom_aggregation)[i]

        colnames(aggregated_vcov_Composition_effect)[sel_terms[1]] <- names(custom_aggregation)[i]
        rownames(aggregated_vcov_Composition_effect)[sel_terms[1]] <- names(custom_aggregation)[i]

        colnames(aggregated_vcov_Structure_effect)[sel_terms[1]] <- names(custom_aggregation)[i]
        rownames(aggregated_vcov_Structure_effect)[sel_terms[1]] <- names(custom_aggregation)[i]

        if(reweighting) {
          colnames(aggregated_vcov_Spec_error)[sel_terms[1]] <- names(custom_aggregation)[i]
          rownames(aggregated_vcov_Spec_error)[sel_terms[1]] <- names(custom_aggregation)[i]

          colnames(aggregated_vcov_RW_error)[sel_terms[1]] <- names(custom_aggregation)[i]
          rownames(aggregated_vcov_RW_error)[sel_terms[1]] <- names(custom_aggregation)[i]
        }
      }
    }


    if(all(is.na(aggregated_terms$Specification_error)) &
       all(is.na(aggregated_terms$Reweighting_error))) {
      aggregated_terms <- aggregated_terms[, 1:4]
      decomposition_vcov$decomposition_terms_se <- decomposition_vcov$decomposition_terms_se[, 1:4]
    }

    aggregated_terms_se <- aggregated_terms

    if(no_se) {
      aggregated_terms_se[] <- NA
      x$decomposition_vcov$decomposition_terms_se <- aggregated_terms_se
      x$decomposition_vcov$decomposition_terms_vcov <- aggregated_terms_se[, -1]
    }
    else {
      aggregated_vcov_Observed_difference <- aggregated_vcov_Observed_difference[names(custom_aggregation), names(custom_aggregation)]
      aggregated_vcov_Composition_effect <- aggregated_vcov_Composition_effect[names(custom_aggregation), names(custom_aggregation)]
      aggregated_vcov_Structure_effect <- aggregated_vcov_Structure_effect[names(custom_aggregation), names(custom_aggregation)]

      aggregated_terms_se[1, ] <- decomposition_vcov$decomposition_terms_se[1, ]
      aggregated_terms_se[-1, "Observed_difference"] <- sqrt(diag(aggregated_vcov_Observed_difference))
      aggregated_terms_se[-1, "Composition_effect"] <- sqrt(diag(aggregated_vcov_Composition_effect))
      aggregated_terms_se[-1, "Structure_effect"] <- sqrt(diag(aggregated_vcov_Structure_effect))

      x$decomposition_vcov$decomposition_terms_se <- aggregated_terms_se

      x$decomposition_vcov$decomposition_terms_vcov$Observed_difference <- aggregated_vcov_Observed_difference
      x$decomposition_vcov$decomposition_terms_vcov$Composition_effect <- aggregated_vcov_Composition_effect
      x$decomposition_vcov$decomposition_terms_vcov$Structure_effect <- aggregated_vcov_Structure_effect


      if(reweighting) {
        aggregated_vcov_Spec_error <- aggregated_vcov_Spec_error[names(custom_aggregation), names(custom_aggregation)]
        aggregated_vcov_RW_error <- aggregated_vcov_RW_error[names(custom_aggregation), names(custom_aggregation)]

        aggregated_terms_se[-1, "Specification_error"] <- sqrt(diag(aggregated_vcov_Spec_error))
        aggregated_terms_se[-1, "Reweighting_error"] <- sqrt(diag(aggregated_vcov_RW_error))

        x$decomposition_vcov$decomposition_terms_vcov$Specifiation_error <- aggregated_vcov_Spec_error
        x$decomposition_vcov$decomposition_terms_vcov$RW_error <- aggregated_vcov_Structure_effect
      }
    }

    x$decomposition_terms <- aggregated_terms


  return(x)
}

