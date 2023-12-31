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
#' data("men8305")
#' flf_model <- log(wage) ~ union*(education + experience) + education*experience
#' flf_male_inequality  <- dfl_deco(flf_model,
#'                                  data = men8305,
#'                                  weights = weights,
#'                                  group = year)
#' plot(flf_male_inequality)
#'
plot.dfl_deco <- function(x, confidence_bands=TRUE, confidence_level = 0.95, uniform_bands=FALSE, ...){

  decomposition_quantiles <-  stats::reshape(x$decomposition_quantiles,
                                                idvar = c("probs"),
                                                times = setdiff(names(x$decomposition_quantiles),c("probs")),
                                                timevar="effect",
                                                varying = list(setdiff(names(x$decomposition_quantiles),c("probs"))),
                                                direction = "long",
                                                v.names = "value")
  confidence_bands <- ifelse(confidence_bands
                             & !is.null(x$bootstrapped_standard_errors),
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




#' Plot decomposition terms for quantiles
#'
#' The function plots decomposition terms for quantiles estimtated
#' with \code{ob_deco} over the  unit interval.
#'
#' @param x an object of class "ob_deco", usually, a result of a call to [ob_deco()] with code{statistics = "quantiles"}.
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
plot.ob_deco <- function(x, confidence_bands=TRUE, confidence_level = 0.95, uniform_bands=FALSE, ...){


  n_quantiles <- length(x) - 5
  col_names <- c("probs", names(x[[1]][["decomposition_terms"]][-1]))

  na_matrix <- matrix(NA, nrow = n_quantiles, ncol = length(col_names))
  deco_results <- as.data.frame(na_matrix)
  colnames(deco_results) <- col_names

  deco_results$probs <- x$input_parameters$rifreg_probs

  for(i in 1:n_quantiles) {
    deco_results[i, 2:length(col_names)] <- x[[i]]$decomposition_terms[1, 2:length(col_names)]
  }


    decomposition_quantiles <-  stats::reshape(deco_results,
                                               idvar = c("probs"),
                                               times = setdiff(names(deco_results),c("probs")),
                                               timevar="effect",
                                               varying = list(setdiff(names(deco_results),c("probs"))),
                                               direction = "long",
                                               v.names = "value")

  confidence_bands <- ifelse(confidence_bands
                             & x$input_parameters$bootstrap,
                             TRUE,
                             FALSE)
  if(confidence_bands){
#
#     deco_se <- as.data.frame(na_matrix)
#     colnames(deco_se) <- col_names
#
#     deco_se$probs <- x$input_parameters$rifreg_probs
#
#     for(i in 1:n_quantiles) {
#       deco_se[i, 2:length(col_names)] <- x[[i]]$decomposition_vcov$decomposition_terms_se[2:length(col_names)]
#     }
#
#     decomposition_quantiles_se <-  stats::reshape(deco_se,
#                                                   idvar = c("probs"),
#                                                   times = setdiff(names(deco_se),c("probs")),
#                                                   timevar="effect",
#                                                   varying = list(setdiff(names(deco_se),c("probs"))),
#                                                   direction = "long",
#                                                   v.names = "se")
#
#     kolmogorov_smirnov_stat <- x$bootstrapped_standard_errors$decomposition_quantiles_kms_distribution
#     kolmogorov_smirnov_stat <- lapply(split(kolmogorov_smirnov_stat, kolmogorov_smirnov_stat$effect),
#                                       function(x) data.frame(effect = x$effect[1],
#                                                              t_value = quantile(x$kms_t_value, confidence_level)))
#     kolmogorov_smirnov_stat <- do.call("rbind", kolmogorov_smirnov_stat)
#     rn <- names(decomposition_quantiles_se)
#     decomposition_quantiles_se <- cbind(decomposition_quantiles_se,
#                                         kolmogorov_smirnov_stat[match(decomposition_quantiles_se$effect,  kolmogorov_smirnov_stat$effect), "t_value"])
#     names(decomposition_quantiles_se) <- c(rn, "t_value")
#     decomposition_quantiles_se$effect_probs <- paste0(decomposition_quantiles_se$effect,decomposition_quantiles_se$probs)
#     decomposition_quantiles$effect_probs <- paste0(decomposition_quantiles$effect,decomposition_quantiles$probs)
#     decomposition_quantiles <- cbind(decomposition_quantiles,
#                                      decomposition_quantiles_se[match(decomposition_quantiles$effect_probs,  decomposition_quantiles_se$effect_probs), c("se","t_value")])
#     decomposition_quantiles$effect_probs <- NULL
#
#     if(uniform_bands==FALSE){
#       decomposition_quantiles$t_value <- stats::qnorm(1-(1-confidence_level)/2)
#     }
#     decomposition_quantiles$effect <- relevel(as.factor(decomposition_quantiles$effect), ref="Observed difference")
#     plot <-  ggplot(data=decomposition_quantiles, aes(x=probs,
#                                                       y=value,
#                                                       #col=effect,
#                                                       #fill=effect,
#                                                       ymin=value-t_value*se,
#                                                       ymax=value+t_value*se
#     )) +
#       geom_hline(yintercept=0, col ="darkgrey", linewidth=.75) +
#       geom_ribbon(alpha=0.2, col=NA, fill="red") +
#       geom_line(col="red") +
#       geom_point(col="red") +
#       facet_wrap(~ effect) +
#       labs(y="Difference", x="Quantile rank")
  }
  else{

    decomposition_quantiles$effect <- relevel(as.factor(decomposition_quantiles$effect), ref="Observed_difference")
    plot <- ggplot(decomposition_quantiles, aes(probs, value, col=effect, shape=effect)) +
      geom_hline(yintercept=0, col ="darkgrey", linewidth=.75) +
      geom_line() +
      geom_point() +
      labs(y="Difference", x="Quantile rank")
  }

  return(plot)
}


#############################################################
### Plot function for composition effect results
get_all_names <- function(formula, df){
  mf <- get_all_vars(formula,df)[,-1]
  mm <- colnames(model.matrix(formula,df))[1]
  vn <- list()
  k <- 1
  if( colnames(model.matrix(formula,df))[1]=="(Intercept)"){
    vn[[k]] <- "(Intercept)"
    names(vn)[k] <- "(Intercept)"
    k <- k+1
  }
  for(i in names(mf)){
    if(is.factor(mf[,i])){
      vn[[k]] <- paste0(i,levels(mf[,i])[-1])
    }else{
      vn[[k]] <- i
    }
    names(vn)[k] <- i
    k <- k+1
  }
  return(vn)
}

#############################################################
### Sum coefficients
aggregate_terms <- function(result,vn){
  x <- result[,1:2]
  for(i in 1:length(vn)){
    if(length(vn[[i]])==1){
      x <- cbind(x,result[,vn[[i]]])
    }else{
      x <- cbind(x,rowSums(result[,vn[[i]]]))
    }
    names(x)[i+2] <- names(vn)[i]
  }
  return(x)
}


#############################################################
### Plot function for composition effect results
rifreg_deco_plot <- function(estimations, type=c(1,2,3,4),
                             effect=c("all","X","S"),
                             varselect=NULL, ylim=c(-1, 1)){

  # Title for plot types 2 and 3
  title <- "Detailed rifreg decomposition (absolute)"

  ###########################################################
  ## Aggregate plot (type 4)
  if(type==4){
    result <- NULL
    ntau <- nrow(estimations[[2]])
    for(i in 2:6){
      effect <- rep(names(estimations[[2]])[i],ntau)
      resultn <- cbind(estimations[[2]][,c(1,i)],effect)
      names(resultn) <- c("tau","delta","effect")
      result <- rbind(result,resultn)
    }

    title <- "Rifreg decomposition: Aggregate decomposition terms"
    plot <- ggplot(result, aes(tau,delta, colour=effect))  +
      geom_line() + geom_point() +
      geom_hline(yintercept = 0, colour="grey") +
      ggtitle(title)

  } else {

    tl <- estimations[[3]]
    result <- estimations[[1]]

    ###########################################################
    # Relativ share of every variable

    if(type==2){
      result <- result[1:(nrow(result)/2),]
      result[,3:ncol(result)] <- result[,3:ncol(result)]/rowSums(result[,3:ncol(result)])
      title <- "Detailed rifreg decomposition (relative)"
    }

    ###########################################################
    # Plot aggregating factors

    if(type==3){
      # Aggregating categorical variables
      resultVars <- result[,3:ncol(result)]
      result <- result[,1:3]

      for(i in 1:length(tl)){
        if(length(pmatch(grep(tl[i], names(resultVars), value=TRUE),names(resultVars)))==1){
          result <- cbind(result,resultVars[,pmatch(grep(tl[i], names(resultVars), value=TRUE),names(resultVars))])
        } else {
          result <- cbind(result,rowSums(resultVars[,pmatch(grep(tl[i], names(resultVars), value=TRUE),names(resultVars))]))
        }
        names(result)[ncol(result)] <- tl[i]
      }

    }


    ###########################################################
    # Create df for plot

    newname <- c("delta","variable")
    vars <- names(result)[-(1:2)]
    if(is.null(varselect)==FALSE){
      vars <- vars[varselect]
    }

    dep <- c("tau","effect")
    df2 <- cbind(result[,c(dep,vars[1])],vars[1])
    select <- (length(df2)-1):length(df2)
    names(df2)[select] <- newname

    # Add more variable if more than one variable is selected
    if(length(vars)>1){
      for(i in 2:length(vars)){
        df_add <- cbind(result[,c(dep,vars[i])],vars[i])
        names(df_add)[select] <- newname
        df2 <- rbind(df2,df_add)
      }
    }

    if(effect=="X"){
      df2 <- df2[which(is.element(df2$effect,c("Delta_X","Spec_error"))),]
    } else if(effect=="S"){
      df2 <- df2[which(is.element(df2$effect,c("Delta_S","RW_error"))),]
    }

    ###########################################################
    # Actual plot

    plot <- ggplot(df2, aes(tau,delta, colour=effect)) +
      geom_line() + geom_point() +
      geom_hline(yintercept = 0, colour="grey") +
      ggtitle(title) + facet_wrap( ~ variable)

    if(type==2){
      plot <- plot + coord_cartesian(ylim=ylim) +  geom_hline(yintercept = 0, color="lightgray")
    }

  }
  return(plot)

}

#############################################################
### Rifreg coefficient plot
rifreg_coef_plot <- function(result, varselect=NULL){

  if(is.null(varselect)){
    varselect=1:ncol(result$estimates$b0)
  }

  tau1 <- result$estimates$tau
  coefb0 <- as.matrix(result$estimates$b0[,varselect])
  coefb1 <- as.matrix(result$estimates$b1[,varselect])
  coefbc <- as.matrix(result$estimates$bC[,varselect])
  group_levels <- result$group_levels
  ref <- result$reference+1

  varnames <-  colnames(result$estimates$b0)[varselect]
  coefs <- NULL
  variable <- NULL
  tau <- NULL
  group <- NULL
  for(i in 1:ncol(coefb0)){
    coefs <- c(coefs,coefb0[,i],coefb1[,i],coefbc[,i])
    tau <- c(tau,rep(tau1,3))
    variable <- c(variable,rep(varnames[i],3*length(tau1)))
    group <- c(group, rep(c(group_levels[1],group_levels[2],paste(group_levels[ref],"reweighted",sep=" ")),each=length(tau1)))
  }
  df <- data.frame(coefs,tau,variable,group)

  #Actual plot
  if(length(varselect)==1){
    plot <- ggplot(df, aes(tau,coefs, shape=group, colour=group)) +
      geom_point() + geom_line() +
      geom_hline(yintercept = 0, colour="grey") +
      ggtitle(paste("Rifreg coefficients",varnames,sep=": "))
  } else {
    plot <- ggplot(df, aes(tau,coefs, shape=group, colour=group)) +
      geom_point() + geom_line() +
      geom_hline(yintercept = 0, colour="grey") +
      ggtitle("Rifreg coefficients") + facet_wrap( ~ variable)
  }
  print(plot)
}
