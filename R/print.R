#' print method for class "dfl_deco"
#'
#' @param x an object of class "dfl_deco", usually , a result of a call to [dfl_deco()].
#' @param ... other parameters to be passed through to printing functions.
#'
#' @return The function \code{print.dfl_deco()} displays the decompositions terms saved in \code{x}.
#' 
#' @export
#' 
print.dfl_deco <- function(x, ...){
  cat("Decomposition of difference between",
      paste0(x$group_variable_name, " == '",x$group_variable_levels[2],"'"), 
      "(group 1) and",
      paste0(x$group_variable_name, " == '",x$group_variable_levels[1],"'"),
      "(group 0)\n\n")
  if(is.null(x$decomposition_quantiles)==FALSE){
    cat("Decomposition of difference at conditional quantiles:\n\n")
    print(x$decomposition_quantiles[, -1])
    cat("\n")
  }
  if(is.null(x$decomposition_other_statistics)==FALSE){
    cat("Decomposition of differences in other statistics\n\n")
    print(x$decomposition_other_statistics[, -1])
    cat("\n")
  }
  
  # cat("Summary statistics of reweighting factors\n\n")
  # quantiles_reweighting_factor <- as.data.frame(x$quantiles_reweighting_factor[,-1])
  # names(quantiles_reweighting_factor) <- names(x$quantiles_reweighting_factor)[-1]
  # rownames(quantiles_reweighting_factor ) <- rownames(x$quantiles_reweighting_factor)
  # print(quantiles_reweighting_factor )
  # cat("\n")
  
}