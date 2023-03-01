#' Estimate distributional statistics
#'
#' Estimate weighted distributional statistics for the reference or
#' the counterfactual group.
#' 
get_distributional_statistics <- function(dep_var,
                                          weights,
                                          group_variable,
                                          group,
                                          statistics,
                                          probs=1:9/10,
                                          log_transformed){
  
  group <- levels(group_variable)[group+1]
  dep_var <- dep_var[which(group_variable==group)]
  weights <- weights[which(group_variable==group)]
  
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
      results  <- c(results, estimate_gini(dep_var=exp(dep_var), weights=weights))
      names(results)[length(results)] <- "Gini of untransformed Y (=exp(log(Y)))"
    }else{
      results  <- c(results, estimate_gini(dep_var=dep_var, weights=weights))
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
  
  return(results)
}


#' Estimate the Gini coefficient
#'
#' #' @author Rothe (2015)
#' estimate_gini <- function (dep_var, weights) {
#'   n <- length(dep_var)
#'   weights <- weights/sum(weights)
#'   G <- sum(dep_var[order(dep_var)] * 1:n * weights[order(dep_var)])
#'   G <- 2 * G/(n*sum(dep_var[order(dep_var)]  * weights[order(dep_var)]))
#'   G <- G - 1 - (1/n)
#'   return(G)
#' }
#' 
#' @export
#' 
estimate_gini <- function (dep_var, weights) {
  weights <- weights/sum(weights)
  weighted_mean <- weighted.mean(x = dep_var, w = weights)
  integrated_generalized_lorenz_curve <- integrate_generalized_lorenz_curve(dep_var, weights)
  gini_coef <- 1 - (2/weighted_mean)*integrated_generalized_lorenz_curve
  return(gini_coef)
}


#' Integrate generalize Lorenz curve
integrate_generalized_lorenz_curve <- function (dep_var, weights) {
  weights <- weights/sum(weights)
  weighted_ecdf <- cumsum(weights[order(dep_var)])
  generalized_lorenz_ordinates <- cumsum(dep_var[order(dep_var)]*weights[order(dep_var)])
  lorenz_curve <- approxfun(c(0,weighted_ecdf),  c(0,generalized_lorenz_ordinates))
  integrated_lorenz_curve <- integrate(lorenz_curve, 0,1)$value
  return(integrated_lorenz_curve )
} 

#' Interquantile ratio
#' 
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
#' #' @export
#' 
estimate_iq_range <- function(dep_var, 
                              weights,
                              probs=c(0.1,0.9)){
  probs <- range(probs)
  iqr <- Hmisc::wtd.quantile(x = dep_var, weights = weights, probs = probs)
  iqr <- iqr[2] - iqr[1]
  return(iqr)
}

