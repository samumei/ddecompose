# test_that("ob_deco() does not throw an error", {
#   set.seed(43825081)
#   library("AER")
#   data("CPS1985")
#   formula <- log(wage) ~ education + experience + union + ethnicity
#   data_used <- CPS1985
#   data_used$weights <- runif(nrow(CPS1985), 0.5, 1.5)
#   # data_used[1,1] <- NA
#   data_used <- get_all_vars(formula, data_used, weights=weights, group=gender)
#
#   deco_results <- ob_deco(formula = formula,
#                           data = data_used,
#                           group = group)
#
#   testthat::expect_error(deco_results, NA)
# })
#
# test_that("ob_deco() returns the same results as R-pacakge oaxaca", {
#   set.seed(43825081)
#   library("oaxaca")
#   data("chicago")
#
#   oaxaca_results <- oaxaca(ln.real.wage ~ age + female + LTHS + some.college +
#                                college + advanced.degree | foreign.born,
#                              data = chicago, R = 30)
#
#   ddeco_results <- ob_deco(formula = ln.real.wage ~ age + female +
#                              LTHS + some.college + college + advanced.degree,
#                            data = chicago,
#                            group = foreign.born,
#                            bootstrap = TRUE,
#                            bootstrap_iterations = 30)
#
#   # no errors
#   testthat::expect_error(oaxaca_results, NA)
#   testthat::expect_error(ddeco_results, NA)
#
#   # same regression results
#   testthat::expect_equal(ddeco_results$model_fits$fit_group_0$coefficients,
#                          oaxaca_results$reg$reg.A$coefficients)
#   testthat::expect_equal(ddeco_results$model_fits$fit_group_0$residuals,
#                          oaxaca_results$reg$reg.A$residuals)
#   testthat::expect_equal(ddeco_results$model_fits$fit_group_0$fitted.values,
#                          oaxaca_results$reg$reg.A$fitted.values)
#   testthat::expect_equal(ddeco_results$model_fits$fit_group_0$effects,
#                          oaxaca_results$reg$reg.A$effects)
#   testthat::expect_equal(ddeco_results$model_fits$fit_group_0$rank,
#                          oaxaca_results$reg$reg.A$rank)
#   testthat::expect_equal(ddeco_results$model_fits$fit_group_0$assign,
#                          oaxaca_results$reg$reg.A$assign)
#   testthat::expect_equal(ddeco_results$model_fits$fit_group_0$df.residual,
#                          oaxaca_results$reg$reg.A$df.residual)
#
#   testthat::expect_equal(ddeco_results$model_fits$fit_group_1$coefficients,
#                          oaxaca_results$reg$reg.B$coefficients)
#   testthat::expect_equal(ddeco_results$model_fits$fit_group_1$residuals,
#                          oaxaca_results$reg$reg.B$residuals)
#   testthat::expect_equal(ddeco_results$model_fits$fit_group_1$fitted.values,
#                          oaxaca_results$reg$reg.B$fitted.values)
#   testthat::expect_equal(ddeco_results$model_fits$fit_group_1$effects,
#                          oaxaca_results$reg$reg.B$effects)
#   testthat::expect_equal(ddeco_results$model_fits$fit_group_1$rank,
#                          oaxaca_results$reg$reg.B$rank)
#   testthat::expect_equal(ddeco_results$model_fits$fit_group_1$assign,
#                          oaxaca_results$reg$reg.B$assign)
#   testthat::expect_equal(ddeco_results$model_fits$fit_group_1$df.residual,
#                          oaxaca_results$reg$reg.B$df.residual)
#
#
#   # same decomposition results
#
#   # DIFFERENT VALUES ARE RETURNED -> OAXACA AUTOMATICALLY RETURNS THREEFOLD
#   # HENCE, CHECK/COMPARE RESULTS LATER
#   # MAYBE ADD MORE INFORMATION TO OB_DECO RESULTS TOO
#
#
# })



test_that("ob_deco() provides expected output", {

  # make manual ob-decomposition by hand and check if the results are the same
  # -> if it works, you can do that for additional parameters to test, them.
  # only look at the actual function, if error fails (or of course after it passes)
  # once, you finished everything here, finish the oaxaca comparison test.

  set.seed(43825081)
  library("AER")
  data("CPS1985")
  formula <- wage ~ education + experience + union
  data_used <- CPS1985
  data_used <- get_all_vars(formula, data_used, group=gender)
  data_used$union <- ifelse(data_used$union == "no", 0, 1)

  deco_results <- ob_deco(formula = formula,
                          data = data_used,
                          group = group)


  ## manual calculation
  data_female <- data_used[data_used$group == "female", ]
  data_male <- data_used[data_used$group == "male", ]

  # regression coefs
  reg_female <- lm(formula, data_female)
  reg_male <- lm(formula, data_male)

  # group means
  mean_female <- colMeans(data_female[, 1:4])
  mean_male <- colMeans(data_male[, 1:4])


  # check that regression outputs are the same
  testthat::expect_equal(deco_results$model_fits$fit_group_0$coefficients,
                         reg_male$coefficients)
  testthat::expect_equal(deco_results$model_fits$fit_group_0$residuals,
                         reg_male$residuals)
  testthat::expect_equal(deco_results$model_fits$fit_group_0$fitted.values,
                         reg_male$fitted.values)

  testthat::expect_equal(deco_results$model_fits$fit_group_1$coefficients,
                         reg_female$coefficients)
  testthat::expect_equal(deco_results$model_fits$fit_group_1$residuals,
                         reg_female$residuals)
  testthat::expect_equal(deco_results$model_fits$fit_group_1$fitted.values,
                         reg_female$fitted.values)


  # check that means are the same

  # YET TO IMPLEMENT in ob_deco

  # check that overall diff is as expected
  overall_diff <- mean_female[1] - mean_male[1]

  female_observed <- sum(mean_female[2:4] * reg_female$coefficients[2:4],
                         reg_female$coefficients[1])
  testthat::expect_equal(female_observed, mean_female[[1]])

  male_observed <- sum(mean_male[2:4] * reg_male$coefficients[2:4],
                         reg_male$coefficients[1])
  testthat::expect_equal(male_observed, mean_male[[1]])
  testthat::expect_equal(overall_diff[[1]], female_observed - male_observed)
  testthat::expect_equal(overall_diff[[1]], deco_results$decomposition_terms$Observed_difference[1])


  # check that decomposition results are as expected
  female_counterfactual <- sum(mean_female[2:4] * reg_male$coefficients[2:4],
                             reg_male$coefficients[1])

  unexplained <- female_observed - female_counterfactual
  explained <- female_counterfactual - male_observed

  testthat::expect_equal(deco_results$decomposition_terms$Composition_effect[1],
                         explained)
  testthat::expect_equal(deco_results$decomposition_terms$Structure_effect[1],
                         unexplained)
  testthat::expect_equal(deco_results$decomposition_terms$Structure_effect[2:5],
                         unname(c(reg_female$coefficients[1] - reg_male$coefficients[1],
                             (mean_female[2:4] * reg_female$coefficients[2:4]) -
                               (mean_female[2:4] * reg_male$coefficients[2:4]))))

  testthat::expect_equal(deco_results$decomposition_terms$Composition_effect[2:5],
                         unname(c(0, (mean_female[2:4] * reg_male$coefficients[2:4]) -
                                    (mean_male[2:4] * reg_male$coefficients[2:4]))))


})

#
# test_that("ob_deco() provides correct output", {
#
#   # make manual ob-decomposition by hand and check if the results are the same
#   # then test, if the results are the same as in oaxaca (maybe test that first)
#   # -> if it works, you can do that for additional parameters to test, them.
#   # only look at the actual function, if error fails (or of course after it passes)
#
#   set.seed(43825081)
#   library("AER")
#   data("CPS1985")
#   formula <- log(wage) ~ education + experience + union + ethnicity
#   data_used <- CPS1985
#   data_used$weights <- runif(nrow(CPS1985), 0.5, 1.5)
#   # data_used[1,1] <- NA
#   data_used <- get_all_vars(formula, data_used, weights=weights, group=gender)
#
#   # cluster <- data_used$cluster <- rbinom(nrow(CPS1985), 10, 0.5)
#   # cluster <- NULL
#   # na.action = na.exclude
#   # reference_0 <- TRUE
#   # vcov=stats::vcov
#   # normalize_factors=FALSE
#   # bootstrap_iterations = 100
#
#   browser()
#   deco_results <- ob_deco(formula = formula,
#                           data = data_used,
#                           group = group)
#
#   testthat::expect_error(deco_results, NA)
# })



