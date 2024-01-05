#' Estimate distributional statistics
#'
#' Estimate weighted distributional statistics for the reference or
#' the counterfactual group.
#'
#' @param dep_var vector of outcome variable
#' @param weights vector of observations weights
#' @param group_variable vector of group assignment
#' @param group identifier of group for which distributional statistics are calculated
#' @param statistics vector of statistics to be calculated
#' @param custom_statistic_function a custom statistic function to be evaluated
#' @param probs probabilities of quantiles to be calculated
#' @param log_transformed indicator if outcome variable is log transformed
#'
get_distributional_statistics <- function(dep_var,
                                          weights,
                                          group_variable,
                                          group,
                                          statistics,
                                          custom_statistic_function = NULL,
                                          probs=1:9/10,
                                          log_transformed){

  group <- levels(group_variable)[group + 1]
  dep_var <- dep_var[which(group_variable == group & weights != 0)]
  weights <- weights[which(group_variable == group & weights != 0)]

  results <- NULL
  if("quantiles" %in% statistics){
    results <- c(results, Hmisc::wtd.quantile(x=dep_var, weights=weights, probs=probs))
    names(results) <- paste0(names(results),"-quantile")
  }

  if("mean" %in% statistics){
    results  <- c(results, weighted.mean(x=dep_var, w=weights))
    names(results)[length(results)] <- "Mean"
  }

  if("variance" %in% statistics){
    results  <- c(results, Hmisc::wtd.var(x=dep_var, weight=weights))
    names(results)[length(results)] <- "Variance"
  }

  if("gini" %in% statistics){
    if(log_transformed){
      results  <- c(results, rifreg::compute_gini(dep_var=exp(dep_var), weights = weights))
      names(results)[length(results)] <- "Gini of untransformed Y (=exp(log(Y)))"
    }else{
      results  <- c(results, rifreg::compute_gini(dep_var=dep_var, weights = weights))
      names(results)[length(results)] <- "Gini"
    }
  }

  if("iq_range_p90_p10" %in% statistics){
    results  <- c(results, estimate_iq_range(dep_var=dep_var, weights=weights, probs=c(0.9,0.1)))
    names(results)[length(results)] <- "Interquantile range p90-p10"
  }

  if("iq_range_p90_p50" %in% statistics){
    results  <- c(results, estimate_iq_range(dep_var=dep_var, weights=weights, probs=c(0.9,0.5)))
    names(results)[length(results)] <- "Interquantile range p90-p50"
  }

  if("iq_range_p50_p10" %in% statistics){
    results  <- c(results, estimate_iq_range(dep_var=dep_var, weights=weights, probs=c(0.5,0.1)))
    names(results)[length(results)] <- "Interquantile range p50-p10"
  }

  if("iq_ratio_p90_p10" %in% statistics){
    results  <- c(results, estimate_iq_ratio(dep_var=dep_var, weights=weights, probs=c(0.9,0.1)))
    names(results)[length(results)] <- "Interquantile ratio p90-p10"
  }

  if("iq_ratio_p90_p50" %in% statistics){
    results  <- c(results, estimate_iq_ratio(dep_var=dep_var, weights=weights, probs=c(0.9,0.5)))
    names(results)[length(results)] <- "Interquantile ratio p90/p50"
  }

  if("iq_ratio_p90_p10" %in% statistics){
    results  <- c(results, estimate_iq_ratio(dep_var=dep_var, weights=weights, probs=c(0.5,0.1)))
    names(results)[length(results)] <- "Interquantile ratio p50/p10"
  }

  if(is.null(custom_statistic_function) == FALSE){
    results  <- c(results, custom_statistic_function(dep_var = dep_var, weights = weights))
    names(results)[length(results)] <- "Custom statistic"
  }

  return(results)
}

#' Interquantile ratio
#'
#' @param dep_var numeric vector of outcome variable
#' @param weights numeric vector of weights
#' @param probs a vector with probabilities whose range defines the interquantile range
#' @export
#'
estimate_iq_ratio <- function(dep_var,
                              weights,
                              probs=c(0.1,0.9)){
  probs <- range(probs)
  iqr <- Hmisc::wtd.quantile(x = dep_var, weights = weights, probs = probs)
  iqr <- iqr[2]/iqr[1]
  return(iqr)
}

#' Interquantile range
#'
#' @param dep_var numeric vector of outcome variable
#' @param weights numeric vector of weights
#' @param probs a vector with probabilities whose range defines the interquantile range
#'
#' @export
#'
estimate_iq_range <- function(dep_var,
                              weights,
                              probs=c(0.1,0.9)){
  probs <- range(probs)
  iqr <- Hmisc::wtd.quantile(x = dep_var, weights = weights, probs = probs)
  iqr <- iqr[2] - iqr[1]
  return(iqr)
}

