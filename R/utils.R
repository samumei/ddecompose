# Packages to import into NAMESPACE

#' @import ggplot2
#' @import pbapply
#' @import stats
#' @import Formula
#' @import parallel
#' @import sandwich
#' @import fastglm
NULL

utils::globalVariables(c("Variable",
                         "effect",
                         "group",
                         "(weights)",
                         "probs",
                         "value",
                         "t_value",
                         "se"))
