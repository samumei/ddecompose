#' summary method for class "dfl_deco"
#'
#' @param x an object of class "dfl_deco", a result of a call to [dfl_deco()].
#' @param confidence_level numeric value between 0 and 1 (default = 0.95) that defines the
#'              confidence level of the printed confidence intervals. Pointwise confidences bands
#'              are defined as \code{qnorm((1-confidence_level)/2)} * standard error. Uniform bands
#'              are constructed by multiplying the standard error with \code{confidence_level}-quantile 
#'              of the bootstrapped Kolmogorov-Smirnov statistic.
#' @param ... other parameters to be passed through to printing functions.
#' @return The function \code{summary.dfl_deco()} displays the decompositions 
#' terms save in \code{x}. If standard errors have been bootstrapped, standard 
#' errors and confidence bands are given. Uniform confidence bands for quantiles
#' are constructed based on the bootstrappend Kolmogorov-Smirnov distribution.
#' 
#' @export
#' 
summary.dfl_deco <- function(x, confidence_level=0.95, digits=4, ...){
  cat("Decomposition of difference between",
      paste0(x$group_variable_name, " == '",x$group_variable_levels[2],"'"), 
      "(group 1) and",
      paste0(x$group_variable_name, " == '",x$group_variable_levels[1],"'"),
      "(group 0)\n\n")
  
  
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
      print(decomposition_quantiles)
      cat("\n")
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
                                                decomposition_other_statistics[,i]))
        results_to_print$ci_p_low <- results_to_print[,2] - qnorm(1-(1-confidence_level)/2)*results_to_print[,3]
        results_to_print$ci_p_high <- results_to_print[,2] + qnorm(1-(1-confidence_level)/2)*results_to_print[,3]
        
        names(results_to_print) <- c("Statistic", "Estimate", "Std. Error","CI [low", "high]")
        print(results_to_print, digits = digits, row.names = FALSE)
        cat("\n")
      }
      
    }else{
      cat("Decomposition of differences in other statistics\n\n")
      print(x$decomposition_other_statistics[, -1])
      cat("\n")
    }
  }
  
  cat("Summary statistics of reweighting factors\n\n")  
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