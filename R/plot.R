#' Plot decomposition terms for quantiles
#'
#' The function plots decomposition terms for quantiles estimtated
#' with \code{dfl_deco} over the  unit interval.
#'
#' @param x an object of class "dfl_deco", usually, a result of a call to [dfl_deco()] with code{statistics = "quantiles"}.
#' @param confidence_bands If `TRUE` (default) and if standard errors have been bootstrapped, confidence bands are plotted.
#' @param confidence_level numeric value between 0 and 1 (default = 0.95) that defines the confidence interval
#'              plotted as a ribbon and defined as \code{qnorm((1-confidence_level)/2)} * standard error.
#' @param uniform_bands If `FALSE` (default), pointsise confidence bands are computed. Otherwise, uniform bands are constructed
#'              based on the bootstrapped Kolmogrov-Smirnov statistic.
#' @param ... other parameters to be passed through to plotting functions.
#'
#' @return a ggplot illustrating the decomposition terms for quantiles.
#' @export
#'
#' @examples
#' # NOT YET PROVIDED
#'
plot.dfl_deco <- function(x, confidence_bands=TRUE, confidence_level = 0.95, uniform_bands=FALSE, ...){

  # decomposition_quantiles <- tidyr::pivot_longer(x$decomposition_quantiles,
  #                                                cols=-c("probs"),
  #                                                names_to="effect",
  #                                                values_to="value")
  decomposition_quantiles <-  stats::reshape(x$decomposition_quantiles,
                                                idvar = c("probs"),
                                                times = setdiff(names(x$decomposition_quantiles),c("probs")),
                                                timevar="effect",
                                                varying = list(setdiff(names(x$decomposition_quantiles),c("probs"))),
                                                direction = "long",
                                                v.names = "value")
  confidence_bands <- ifelse(confidence_bands == TRUE
                             & is.null(x$bootstrapped_standard_errors) ==FALSE,
                             TRUE,
                             FALSE)
  if(confidence_bands){
    decomposition_quantiles_se <- x$bootstrapped_standard_errors$decomposition_quantiles
    decomposition_quantiles_se <-  stats::reshape(decomposition_quantiles_se,
                                                  idvar = c("probs"),
                                                  times = setdiff(names(decomposition_quantiles_se),c("probs")),
                                                  timevar="effect",
                                                  varying = list(setdiff(names(decomposition_quantiles_se),c("probs"))),
                                                  direction = "long",
                                                  v.names = "se")
    # decomposition_quantiles_se <- tidyr::pivot_longer(x$bootstrapped_standard_errors$decomposition_quantiles,
    #                                                   cols=-c("probs"),
    #                                                   names_to="effect",
    #                                                   values_to="se")
    kolmogorov_smirnov_stat <- x$bootstrapped_standard_errors$decomposition_quantiles_kms_distribution
    kolmogorov_smirnov_stat <- lapply(split(kolmogorov_smirnov_stat, kolmogorov_smirnov_stat$effect),
                                      function(x) data.frame(effect = x$effect[1],
                                                             t_value = quantile(x$kms_t_value, confidence_level)))
    kolmogorov_smirnov_stat <- do.call("rbind", kolmogorov_smirnov_stat)
    rn <- names(decomposition_quantiles_se)
    decomposition_quantiles_se <- cbind(decomposition_quantiles_se,
                                        kolmogorov_smirnov_stat[match(decomposition_quantiles_se$effect,  kolmogorov_smirnov_stat$effect), "t_value"])
    names(decomposition_quantiles_se) <- c(rn, "t_value")
    decomposition_quantiles_se$effect_probs <- paste0(decomposition_quantiles_se$effect,decomposition_quantiles_se$probs)
    decomposition_quantiles$effect_probs <- paste0(decomposition_quantiles$effect,decomposition_quantiles$probs)
    decomposition_quantiles <- cbind(decomposition_quantiles,
                                     decomposition_quantiles_se[match(decomposition_quantiles$effect_probs,  decomposition_quantiles_se$effect_probs), c("se","t_value")])
    decomposition_quantiles$effect_probs <- NULL

    # kolmogorov_smirnov_stat <- dplyr::summarise(dplyr::group_by(x$bootstrapped_standard_errors$decomposition_quantiles_kms_distribution,
    #                                                             effect),
    #                                             t_value = quantile(kms_t_value, confidence_level))
    # decomposition_quantiles_se <- dplyr::left_join(decomposition_quantiles_se,
    #                                                kolmogorov_smirnov_stat,
    #                                                by="effect")
    # decomposition_quantiles <- dplyr::left_join(decomposition_quantiles,
    #                                             decomposition_quantiles_se,
    #                                             by=c("probs","effect"))

    if(uniform_bands==FALSE){
      decomposition_quantiles$t_value <- stats::qnorm(1-(1-confidence_level)/2)
    }
    decomposition_quantiles$effect <- relevel(as.factor(decomposition_quantiles$effect), ref="Observed difference")
    plot <-  ggplot(data=decomposition_quantiles, aes(x=probs,
                                                      y=value,
                                                      #col=effect,
                                                      #fill=effect,
                                                      ymin=value-t_value*se,
                                                      ymax=value+t_value*se
                                                      )) +
      geom_hline(yintercept=0, col ="darkgrey", linewidth=.75) +
      geom_ribbon(alpha=0.2, col=NA, fill="red") +
      geom_line(col="red") +
      geom_point(col="red") +
      facet_wrap(~ effect) +
      labs(y="Difference", x="Quantile rank")
  }else{
    decomposition_quantiles$effect <- relevel(as.factor(decomposition_quantiles$effect), ref="Observed difference")
    plot <- ggplot(decomposition_quantiles, aes(probs, value, col=effect, shape=effect)) +
      geom_hline(yintercept=0, col ="darkgrey", linewidth=.75) +
      geom_line() +
      geom_point() +
      labs(y="Difference", x="Quantile rank")
  }

  return(plot)
}
