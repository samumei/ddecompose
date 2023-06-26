# test_that("ob_deco() does not throw an error", {
#   set.seed(43825081)
#   data("nlys00")
#
#   mod1 <- log(wage) ~ age + central_city + msa + region + black +
#   hispanic + education + afqt + family_responsibility + years_worked_civilian +
#   years_worked_military + part_time + industry
#
#   # Using female coefficients (reference_0 = TRUE) to estimate counterfactual mean
#   deco_female_as_reference <- ob_deco(formula = mod1,
#                                       data = nlys00,
#                                       group = female,
#                                       reference_0 = TRUE)
#   testthat::expect_error(deco_female_as_reference, NA)
# })
#
# test_that("ob_deco() returns the same results as R-pacakge oaxaca", {
#   set.seed(43825081)
#   library("oaxaca")
#   data("chicago")
#
#   oaxaca_results <- oaxaca(ln.real.wage ~ age + female | foreign.born,
#                              data = chicago, R = 100)
#
#   ddeco_results <- ob_deco(formula = ln.real.wage ~ age + female,
#                            data = chicago,
#                            group = foreign.born,
#                            bootstrap = TRUE,
#                            bootstrap_iterations = 100)
#
#   # no errors
#   testthat::expect_error(oaxaca_results, NA)
#   testthat::expect_error(ddeco_results, NA)
#
#   # same regression results
#   testthat::expect_equal(ddeco_results$ob_deco$model_fits$fit_group_0$coefficients,
#                          oaxaca_results$reg$reg.A$coefficients)
#   testthat::expect_equal(ddeco_results$ob_deco$model_fits$fit_group_0$residuals,
#                          oaxaca_results$reg$reg.A$residuals)
#   testthat::expect_equal(ddeco_results$ob_deco$model_fits$fit_group_0$fitted.values,
#                          oaxaca_results$reg$reg.A$fitted.values)
#   testthat::expect_equal(ddeco_results$ob_deco$model_fits$fit_group_0$effects,
#                          oaxaca_results$reg$reg.A$effects)
#   testthat::expect_equal(ddeco_results$ob_deco$model_fits$fit_group_0$rank,
#                          oaxaca_results$reg$reg.A$rank)
#   testthat::expect_equal(ddeco_results$ob_deco$model_fits$fit_group_0$assign,
#                          oaxaca_results$reg$reg.A$assign)
#   testthat::expect_equal(ddeco_results$ob_deco$model_fits$fit_group_0$df.residual,
#                          oaxaca_results$reg$reg.A$df.residual)
#
#   testthat::expect_equal(ddeco_results$ob_deco$model_fits$fit_group_1$coefficients,
#                          oaxaca_results$reg$reg.B$coefficients)
#   testthat::expect_equal(ddeco_results$ob_deco$model_fits$fit_group_1$residuals,
#                          oaxaca_results$reg$reg.B$residuals)
#   testthat::expect_equal(ddeco_results$ob_deco$model_fits$fit_group_1$fitted.values,
#                          oaxaca_results$reg$reg.B$fitted.values)
#   testthat::expect_equal(ddeco_results$ob_deco$model_fits$fit_group_1$effects,
#                          oaxaca_results$reg$reg.B$effects)
#   testthat::expect_equal(ddeco_results$ob_deco$model_fits$fit_group_1$rank,
#                          oaxaca_results$reg$reg.B$rank)
#   testthat::expect_equal(ddeco_results$ob_deco$model_fits$fit_group_1$assign,
#                          oaxaca_results$reg$reg.B$assign)
#   testthat::expect_equal(ddeco_results$ob_deco$model_fits$fit_group_1$df.residual,
#                          oaxaca_results$reg$reg.B$df.residual)
#
#
#   ## same decomposition results
#   # overall
#   testthat::expect_equal(ddeco_results$ob_deco$decomposition_terms$Composition_effect[1],
#                          unname(oaxaca_results$twofold$overall[2,2]))
#   testthat::expect_equal(ddeco_results$ob_deco$decomposition_terms$Structure_effect[1],
#                          unname(oaxaca_results$twofold$overall[2,4]))
#
#   # detailed
#   testthat::expect_equal(ddeco_results$ob_deco$decomposition_terms$Composition_effect[2:4],
#                          unname(oaxaca_results$twofold$variables[[2]][,2]))
#   testthat::expect_equal(ddeco_results$ob_deco$decomposition_terms$Structure_effect[2:4],
#                          unname(oaxaca_results$twofold$variables[[2]][,4]))
#
#   # SE similar
#   testthat::expect_equal(ddeco_results$ob_deco$decomposition_vcov$decomposition_terms_se$Composition_effect[1],
#                          unname(oaxaca_results$twofold$overall[2,3]),
#                          tolerance = 0.03)
#   testthat::expect_equal(ddeco_results$ob_deco$decomposition_vcov$decomposition_terms_se$Structure_effect[1],
#                          unname(oaxaca_results$twofold$overall[2,5]),
#                          tolerance = 0.1)
#
#   testthat::expect_equal(ddeco_results$ob_deco$decomposition_vcov$decomposition_terms_se$Composition_effect[2:4],
#                          unname(oaxaca_results$twofold$variables[[2]][,3]),
#                          tolerance = 0.03)
#   testthat::expect_equal(ddeco_results$ob_deco$decomposition_vcov$decomposition_terms_se$Structure_effect[2:4],
#                          unname(oaxaca_results$twofold$variables[[2]][,5]),
#                          tolerance = 0.1)
#
#   # compare analytical se
#   ddeco_results_analytical_se <- ob_deco(formula = ln.real.wage ~ age + female,
#                                          data = chicago,
#                                          group = foreign.born)
#
#   testthat::expect_equal(ddeco_results_analytical_se$ob_deco$decomposition_vcov$decomposition_terms_se$Composition_effect[1],
#                          unname(oaxaca_results$twofold$overall[2,3]),
#                          tolerance = 0.03)
#   testthat::expect_equal(ddeco_results_analytical_se$ob_deco$decomposition_vcov$decomposition_terms_se$Structure_effect[1],
#                          unname(oaxaca_results$twofold$overall[2,5]),
#                          tolerance = 0.1)
#
#   testthat::expect_equal(ddeco_results_analytical_se$ob_deco$decomposition_vcov$decomposition_terms_se$Composition_effect[2:4],
#                          unname(oaxaca_results$twofold$variables[[2]][,3]),
#                          tolerance = 0.03)
#   testthat::expect_equal(ddeco_results_analytical_se$ob_deco$decomposition_vcov$decomposition_terms_se$Structure_effect[2:4],
#                          unname(oaxaca_results$twofold$variables[[2]][,5]),
#                          tolerance = 0.15)
# })
#
#
# test_that("ob_deco() provides expected output", {
#
#   set.seed(43825081)
#   library("AER")
#   data("CPS1985")
#   formula <- wage ~ education + experience + union
#   data_used <- CPS1985
#   data_used <- get_all_vars(formula, data_used, group=gender)
#   data_used$union <- ifelse(data_used$union == "no", 0, 1)
#
#   deco_results <- ob_deco(formula = formula,
#                           data = data_used,
#                           group = group)
#
#
#   ## manual calculation
#   data_female <- data_used[data_used$group == "female", ]
#   data_male <- data_used[data_used$group == "male", ]
#
#   # regression coefs
#   reg_female <- lm(formula, data_female)
#   reg_male <- lm(formula, data_male)
#
#   # group means
#   mean_female <- colMeans(data_female[, 1:4])
#   mean_male <- colMeans(data_male[, 1:4])
#
#   # check that regression outputs are the same
#   testthat::expect_equal(deco_results$ob_deco$model_fits$fit_group_0$coefficients,
#                          reg_male$coefficients)
#   testthat::expect_equal(deco_results$ob_deco$model_fits$fit_group_0$residuals,
#                          reg_male$residuals)
#   testthat::expect_equal(deco_results$ob_deco$model_fits$fit_group_0$fitted.values,
#                          reg_male$fitted.values)
#
#   testthat::expect_equal(deco_results$ob_deco$model_fits$fit_group_1$coefficients,
#                          reg_female$coefficients)
#   testthat::expect_equal(deco_results$ob_deco$model_fits$fit_group_1$residuals,
#                          reg_female$residuals)
#   testthat::expect_equal(deco_results$ob_deco$model_fits$fit_group_1$fitted.values,
#                          reg_female$fitted.values)
#
#
#   # check that means are the same
#
#   # check that overall diff is as expected
#   overall_diff <- mean_male[1] - mean_female[1]
#
#   female_observed <- sum(mean_female[2:4] * reg_female$coefficients[2:4],
#                          reg_female$coefficients[1])
#   testthat::expect_equal(female_observed, mean_female[[1]])
#
#   male_observed <- sum(mean_male[2:4] * reg_male$coefficients[2:4],
#                          reg_male$coefficients[1])
#   testthat::expect_equal(male_observed, mean_male[[1]])
#   testthat::expect_equal(overall_diff[[1]],male_observed -  female_observed)
#   testthat::expect_equal(overall_diff[[1]], deco_results$ob_deco$decomposition_terms$Observed_difference[1])
#
#
#   # check that decomposition results are as expected
#   female_counterfactual <- sum(mean_female[2:4] * reg_male$coefficients[2:4],
#                              reg_male$coefficients[1])
#
#   unexplained <- female_counterfactual - female_observed
#   explained <- male_observed -female_counterfactual
#
#   testthat::expect_equal(deco_results$ob_deco$decomposition_terms$Composition_effect[1],
#                          explained)
#   testthat::expect_equal(deco_results$ob_deco$decomposition_terms$Structure_effect[1],
#                          unexplained)
# })
#
#
# test_that("ob_deco() analytical and bootstrapped se are the same", {
#
#   set.seed(43825081)
#   library("AER")
#   data("CPS1985")
#   formula <- log(wage) ~ education + experience + union + ethnicity
#   data_used <- CPS1985
#   data_used$weights <- runif(nrow(CPS1985), 0.5, 1.5)
#   data_used <- get_all_vars(formula, data_used, weights=weights, group=gender)
#
#   deco_analytical_se <- ob_deco(formula = formula,
#                                 data = data_used,
#                                 group = group)
#   deco_bootstrapped_se <- ob_deco(formula = formula,
#                                   data = data_used,
#                                   group = group,
#                                   bootstrap = TRUE,
#                                   bootstrap_iterations = 400)
#
#   testthat::expect_equal(deco_analytical_se$ob_deco$decomposition_vcov$decomposition_terms_se,
#                          deco_bootstrapped_se$ob_deco$decomposition_vcov$decomposition_terms_se[1:4],
#                          tolerance = 0.04)
# })
#
#
# test_that("reweighted ob decomposition does not throw an error", {
#
#   set.seed(43825081)
#   data("nlys00")
#
#   mod1 <- log(wage) ~ age + central_city + msa + region + black +
#     hispanic + education + afqt + family_responsibility + years_worked_civilian +
#     years_worked_military + part_time + industry
#
#   reweighted_deco_results <- ob_deco(formula = mod1,
#                                       data = nlys00,
#                                       group = female,
#                                       reference_0 = TRUE,
#                                       reweighting = TRUE)
#
#   testthat::expect_error(reweighted_deco_results, NA)
# })
#
# test_that("reweighted ob decomposition aggregate results are as expected", {
#
#   set.seed(43825081)
#   data("nlys00")
#
#   mod1 <- log(wage) ~ age + central_city + msa + region + black +
#     hispanic + education + afqt + family_responsibility + years_worked_civilian +
#     years_worked_military + part_time + industry
#
#   reweighted_deco_results <- ob_deco(formula = mod1,
#                                      data = nlys00,
#                                      group = female,
#                                      reference_0 = TRUE,
#                                      reweighting = TRUE)
#
#   deco_results <- ob_deco(formula = mod1,
#                           data = nlys00,
#                           group = female,
#                           reference_0 = TRUE,
#                           reweighting = FALSE)
#
#   testthat::expect_equal(reweighted_deco_results$reweighted_ob_deco$decomposition_terms$Observed_difference,
#                          deco_results$ob_deco$decomposition_terms$Observed_difference)
#   testthat::expect_equal(reweighted_deco_results$reweighted_ob_deco$decomposition_terms$Composition_effect +
#                            reweighted_deco_results$reweighted_ob_deco$decomposition_terms$Specification_error +
#                            reweighted_deco_results$reweighted_ob_deco$decomposition_terms$Structure_effect +
#                            reweighted_deco_results$reweighted_ob_deco$decomposition_terms$Reweighting_error,
#                          deco_results$ob_deco$decomposition_terms$Observed_difference)
#
# })
#
# test_that("reweighted ob decomposition aggregate results are as expected with male reference", {
#
#   set.seed(43825081)
#   data("nlys00")
#
#   mod1 <- log(wage) ~ age + central_city + msa + region + black +
#     hispanic + education + afqt + family_responsibility + years_worked_civilian +
#     years_worked_military + part_time + industry
#
#   reweighted_deco_results <- ob_deco(formula = mod1,
#                                      data = nlys00,
#                                      group = female,
#                                      reference_0 = FALSE,
#                                      reweighting = TRUE)
#
#   deco_results <- ob_deco(formula = mod1,
#                           data = nlys00,
#                           group = female,
#                           reference_0 = FALSE,
#                           reweighting = FALSE)
#
#   testthat::expect_equal(reweighted_deco_results$reweighted_ob_deco$decomposition_terms$Observed_difference,
#                          deco_results$ob_deco$decomposition_terms$Observed_difference)
#   testthat::expect_equal(reweighted_deco_results$reweighted_ob_deco$decomposition_terms$Composition_effect +
#                            reweighted_deco_results$reweighted_ob_deco$decomposition_terms$Specification_error +
#                            reweighted_deco_results$reweighted_ob_deco$decomposition_terms$Structure_effect +
#                            reweighted_deco_results$reweighted_ob_deco$decomposition_terms$Reweighting_error,
#                          deco_results$ob_deco$decomposition_terms$Observed_difference)
#
# })
#
#
# test_that("reweighted ob decomposition aggregate results are as expected with normalization", {
#
#   set.seed(43825081)
#   data("nlys00")
#
#   mod1 <- log(wage) ~ age + central_city + msa + region + black +
#     hispanic + education + afqt + family_responsibility + years_worked_civilian +
#     years_worked_military + part_time + industry
#
#   reweighted_deco_results <- ob_deco(formula = mod1,
#                                      data = nlys00,
#                                      group = female,
#                                      reference_0 = TRUE,
#                                      reweighting = TRUE,
#                                      normalize_factors = TRUE)
#
#   deco_results <- ob_deco(formula = mod1,
#                           data = nlys00,
#                           group = female,
#                           reference_0 = TRUE,
#                           reweighting = FALSE,
#                           normalize_factors = TRUE)
#
#   testthat::expect_equal(reweighted_deco_results$reweighted_ob_deco$decomposition_terms$Observed_difference,
#                          deco_results$ob_deco$decomposition_terms$Observed_difference)
#   testthat::expect_equal(reweighted_deco_results$reweighted_ob_deco$decomposition_terms$Composition_effect +
#                            reweighted_deco_results$reweighted_ob_deco$decomposition_terms$Specification_error +
#                            reweighted_deco_results$reweighted_ob_deco$decomposition_terms$Structure_effect +
#                            reweighted_deco_results$reweighted_ob_deco$decomposition_terms$Reweighting_error,
#                          deco_results$ob_deco$decomposition_terms$Observed_difference)
#
# })
#
#
# test_that("reweighted ob decomposition computes bootstrap SE", {
#
#   set.seed(43825081)
#   data("nlys00")
#
#   mod1 <- log(wage) ~ age + central_city + msa + region + black +
#     hispanic + education + afqt + family_responsibility + years_worked_civilian +
#     years_worked_military + part_time + industry
#
#   reweighted_deco_results <- ob_deco(formula = mod1,
#                                      data = nlys00,
#                                      group = female,
#                                      reference_0 = TRUE,
#                                      reweighting = TRUE,
#                                      normalize_factors = TRUE,
#                                      bootstrap = TRUE)
#
#   testthat::expect_error(reweighted_deco_results, NA)
# })
#
# test_that("rifreg decomposition does not throw an error", {
#
#   set.seed(43825081)
#   data("nlys00")
#
#   mod1 <- log(wage) ~ age + central_city + msa + region + black +
#     hispanic + education + afqt + family_responsibility + years_worked_civilian +
#     years_worked_military + part_time + industry
#
#   rifreg_deco_results <- ob_deco(formula = mod1,
#                                      data = nlys00,
#                                      group = female,
#                                      reference_0 = TRUE,
#                                      reweighting = TRUE,
#                                      rifreg = TRUE,
#                                  rifreg_probs = 0.5,
#                                      normalize_factors = TRUE,
#                                      bootstrap = FALSE)
#
#   testthat::expect_error(rifreg_deco_results, NA)
#
# })
#
# test_that("rifreg decomposition does not throw an error with bootstrap", {
#
#   set.seed(43825081)
#   data("nlys00")
#
#   mod1 <- log(wage) ~ age + central_city + msa + region + black +
#     hispanic + education + afqt + family_responsibility + years_worked_civilian +
#     years_worked_military + part_time + industry
#
#   rifreg_deco_results <- ob_deco(formula = mod1,
#                                  data = nlys00,
#                                  group = female,
#                                  reference_0 = TRUE,
#                                  reweighting = TRUE,
#                                  rifreg = TRUE,
#                                  rifreg_probs = 0.5,
#                                  normalize_factors = TRUE,
#                                  bootstrap = TRUE)
#
#   testthat::expect_error(rifreg_deco_results, NA)
#
# })
#
# test_that("rifreg decomposition does not throw an error with multiple quantiles", {
#
#   set.seed(43825081)
#   data("nlys00")
#
#   mod1 <- log(wage) ~ age + central_city + msa + region + black +
#     hispanic + education + afqt + family_responsibility + years_worked_civilian +
#     years_worked_military + part_time + industry
#
#   rifreg_deco_results <- ob_deco(formula = mod1,
#                                  data = nlys00,
#                                  group = female,
#                                  reference_0 = TRUE,
#                                  reweighting = TRUE,
#                                  rifreg = TRUE,
#                                  rifreg_probs = c(1:9)/10,
#                                  normalize_factors = TRUE,
#                                  bootstrap = FALSE)
#
#   testthat::expect_error(rifreg_deco_results, NA)
#
# })
#
# test_that("rifreg decomposition does not throw an error with multiple quantiles and bootstrap", {
#
#   set.seed(43825081)
#   data("nlys00")
#
#   mod1 <- log(wage) ~ age + central_city + msa + region + black +
#     hispanic + education + afqt + family_responsibility + years_worked_civilian +
#     years_worked_military + part_time + industry
#
#   rifreg_deco_results <- ob_deco(formula = mod1,
#                                  data = nlys00,
#                                  group = female,
#                                  reference_0 = TRUE,
#                                  reweighting = TRUE,
#                                  rifreg = TRUE,
#                                  rifreg_probs = c(1:9)/10,
#                                  normalize_factors = TRUE,
#                                  bootstrap = TRUE,
#                                  bootstrap_iterations = 5)
#
#   testthat::expect_error(rifreg_deco_results, NA)
#
# })
#
#
# test_that("dfl_deco() replicates Table 5-D, p. 67, in FLF 2011 Handbook Chapter", {
# browser()
#   #load("data-raw/men8305_full.rda")
#   men8305_full$weights <- men8305_full$weights/sum(men8305_full$weights) * length(men8305_full$weights)
#   flf_model <- log(wage) ~ union*(education + experience) + education*experience
#
#   # Replicate statistics in table 5-D, p.67, in in FLF (2011)
#   flf_iq_range_p90_p10  <- ob_deco(flf_model,
#                                            data = men8305_full,
#                                            weights = weights,
#                                            group = year,
#                                            reference_0 = TRUE,
#                                            rifreg_statistic = "iq_range_p90_p10")
#   results_rifreg_deco <- flf_male_inequality_table_5$decomposition_terms
#
#   published_results_FLF_table_5 <- data.frame(statistic = results_dfl_deco$statistic,
#                                               `Observed difference` = c(0.0617, 0.1091, 0.1827, -0.0736),
#                                               `Structure effect` = c(0.0408, 0.0336, 0.1637, -0.1301),
#                                               `Composition effect` = c(0.0208, 0.0756, 0.0191, 0.0565))
#   rownames(results_dfl_deco) <- rownames(published_results_FLF_table_5) <- 1:4
#   names(published_results_FLF_table_5) <- names(results_dfl_deco)
#   testthat::expect_equal(results_dfl_deco,
#                          published_results_FLF_table_5,
#                          tolerance = 0.0075)
# })

# test_that("ob_deco() replicates Table 3, p. 41, in FLF 2011 Handbook Chapter", {
#
#   data("nlys00")
#   mod1 <- log(wage) ~ age + central_city + msa + region + black +
#       hispanic + education + afqt + family_responsibility + years_worked_civilian +
#       years_worked_military + part_time + industry
#
#     # Using female coefficients (reference_0 = TRUE) to estimate counterfactual mean
#     deco_female_as_reference <- ob_deco(formula = mod1,
#                                         data = nlys00,
#                                         group = female,
#                                         reference_0 = TRUE)
#   testthat::expect_equal(deco_female_as_reference$ob_deco$decomposition_terms$Observed_difference[1], 0.233, tolerance = 0.01)
#    testthat::expect_equal(deco_female_as_reference$ob_deco$decomposition_terms$Composition_effect[1], 0.136, tolerance = 0.01)
#    testthat::expect_equal(deco_female_as_reference$ob_deco$decomposition_terms$Structure_effect[1], 0.097, tolerance = 0.01)
#
#    deco_male_as_reference <- ob_deco(formula = mod1,
#                                        data = nlys00,
#                                        group = female,
#                                        reference_0 = FALSE)
#    testthat::expect_equal(deco_male_as_reference$ob_deco$decomposition_terms$Observed_difference[1], 0.233, tolerance = 0.01)
#    testthat::expect_equal(deco_male_as_reference$ob_deco$decomposition_terms$Composition_effect[1], 0.197, tolerance = 0.01)
#    testthat::expect_equal(deco_male_as_reference$ob_deco$decomposition_terms$Structure_effect[1], 0.036, tolerance = 0.01)
#
#
# })
# test_that("ob_deco() replicates Table 4-C, p. 62, in FLF 2011 Handbook Chapter", {
#   browser()
#
# #   Abweichungen: sind sie wirklich grösser als beim OLS?? -> ja, schon..
# #   NORMALizarione anschauen --> erst beim detaillierten Zeugs.
# #   Gewichte anschauen
# #   Ihre Berechnungen in Stata-File anschauen
# #   Umgewichtung haben sie starke umgewichtungen.
# #   RIFREG umgewichtung eigene Funktion -> wie bei DFL-Deco mit mehreren bedingungen.-> schauen, was du noch machen mussen
#
#
#   data("nlys00")
#   mod1 <- log(wage) ~ age + central_city + msa + region + black +
#     hispanic + education + afqt + family_responsibility + years_worked_civilian +
#     years_worked_military + part_time + industry
#
#   # Using male coefficients to estimate counterfactual mean
#   rifreg_deco <- ob_deco(formula = mod1,
#                                       data = nlys00,
#                                       group = female,
#                                       rifreg = TRUE,
#                                       rifreg_statistic = "quantiles",
#                                       rifreg_probs = 0.5, # c(0.1, 0.5, 0.9),
#                                       reference_0 = FALSE)
#
#   rifreg_deco_female <- ob_deco(formula = mod1,
#                          data = nlys00,
#                          group = female,
#                          rifreg = TRUE,
#                          rifreg_statistic = "quantiles",
#                          rifreg_probs = 0.5, #c(0.1, 0.5, 0.9),
#                          reference_0 = TRUE)
#   rifreg_deco$quantile_0.5$decomposition_terms$Observed_difference[1]
#
#   testthat::expect_equal(rifreg_deco$quantile_0.1$decomposition_terms$Observed_difference[1], 0.180, tolerance = 0.01)
#   testthat::expect_equal(rifreg_deco$quantile_0.1$decomposition_terms$Composition_effect[1], 0.274, tolerance = 0.01)
#   testthat::expect_equal(rifreg_deco$quantile_0.1$decomposition_terms$Structure_effect[1], -0.094, tolerance = 0.01)
#
#   testthat::expect_equal(rifreg_deco$quantile_0.5$decomposition_terms$Observed_difference[1], 0.241, tolerance = 0.01)
#   testthat::expect_equal(rifreg_deco$quantile_0.5$decomposition_terms$Composition_effect[1], 0.208, tolerance = 0.01)
#   testthat::expect_equal(rifreg_deco$quantile_0.5$decomposition_terms$Structure_effect[1], 0.033, tolerance = 0.01)
#
#   testthat::expect_equal(rifreg_deco$quantile_0.9$decomposition_terms$Observed_difference[1], 0.260, tolerance = 0.01)
#   testthat::expect_equal(rifreg_deco$quantile_0.9$decomposition_terms$Composition_effect[1], 0.136, tolerance = 0.01)
#   testthat::expect_equal(rifreg_deco$quantile_0.9$decomposition_terms$Structure_effect[1], 0.124, tolerance = 0.01)
#
#
#   browser()
#   # Replicate statistics in table 5-D, p.67, in in FLF (2011)
#   flf_variance  <- ob_deco(flf_model,
#                            data = men8305_full,
#                            weights = weights,
#                            group = year,
#                            reference_0 = TRUE,
#                            rifreg = TRUE,
#                            rifreg_statistic = "variance")
#   View(flf_variance[["variance"]][["decomposition_terms"]][1,])
#
#   flf_gini  <- ob_deco(flf_model,
#                        data = men8305_full,
#                        weights = weights,
#                        group = year,
#                        reference_0 = TRUE,
#                        rifreg = TRUE,
#                        rifreg_statistic = "gini")
#   View(flf_gini[["gini"]][["decomposition_terms"]][1,])
#
#
#   flf_iq_range_p90_p10  <- ob_deco(flf_model,
#                                    data = men8305_full,
#                                    weights = weights,
#                                    group = year,
#                                    reference_0 = TRUE,
#                                    rifreg = TRUE,
#                                    rifreg_statistic = "interquantile_range",
#                                    rifreg_probs = c(0.9, 0.1))
#   View(flf_iq_range_p90_p10[["interquantile_range"]][["decomposition_terms"]][1,])
#
#   rifreg_results <- rbind(flf_variance[["ob_deco"]][["decomposition_terms"]][1,],
#                           flf_gini[["ob_deco"]][["decomposition_terms"]][1,],
#                           flf_iq_range_p90_p10[["ob_deco"]][["decomposition_terms"]][1,])
#   results_rifreg_deco <- flf_male_inequality_table_5$decomposition_terms
#
#   published_results_FLF_table_5 <- data.frame(statistic = results_dfl_deco$statistic,
#                                               `Observed difference` = c(0.0617, 0.1091, 0.1827, -0.0736),
#                                               `Structure effect` = c(0.0408, 0.0336, 0.1637, -0.1301),
#                                               `Composition effect` = c(0.0208, 0.0756, 0.0191, 0.0565))
#   rownames(results_dfl_deco) <- rownames(published_results_FLF_table_5) <- 1:4
#   names(published_results_FLF_table_5) <- names(results_dfl_deco)
#   testthat::expect_equal(results_dfl_deco,
#                          published_results_FLF_table_5,
#                          tolerance = 0.0075)
# })
# test_that("dfl_deco() replicates Table 5-D, p. 67, in FLF 2011 Handbook Chapter", {
#
#   #load("data-raw/men8305_full.rda")
#   men8305_full$weights <- men8305_full$weights/sum(men8305_full$weights) * length(men8305_full$weights)
#   flf_model <- log(wage) ~ union*(education + experience) + education*experience
#
#   browser()
#   # Replicate statistics in table 5-D, p.67, in in FLF (2011)
#   flf_variance  <- ob_deco(flf_model,
#                            data = men8305_full,
#                            weights = weights,
#                            group = year,
#                            reference_0 = FALSE,
#                            rifreg = TRUE,
#                            rifreg_statistic = "variance")
#   View(flf_variance[["variance"]][["decomposition_terms"]][1,])
#
#   flf_gini  <- ob_deco(flf_model,
#                        data = men8305_full,
#                        weights = weights,
#                        group = year,
#                        reference_0 = TRUE,
#                        rifreg = TRUE,
#                        rifreg_statistic = "gini")
#   View(flf_gini[["gini"]][["decomposition_terms"]][1,])
#
#
#   flf_iq_range_p90_p10  <- ob_deco(flf_model,
#                                    data = men8305_full,
#                                    weights = weights,
#                                    group = year,
#                                    reference_0 = TRUE,
#                                    rifreg = TRUE,
#                                    rifreg_statistic = "interquantile_range",
#                                    rifreg_probs = c(0.9, 0.1))
#   View(flf_iq_range_p90_p10[["interquantile_range"]][["decomposition_terms"]][1,])
#
#   rifreg_results <- rbind(flf_variance[["ob_deco"]][["decomposition_terms"]][1,],
#                           flf_gini[["ob_deco"]][["decomposition_terms"]][1,],
#                           flf_iq_range_p90_p10[["ob_deco"]][["decomposition_terms"]][1,])
#   results_rifreg_deco <- flf_male_inequality_table_5$decomposition_terms
#
#   published_results_FLF_table_5 <- data.frame(statistic = results_dfl_deco$statistic,
#                                               `Observed difference` = c(0.0617, 0.1091, 0.1827, -0.0736),
#                                               `Structure effect` = c(0.0408, 0.0336, 0.1637, -0.1301),
#                                               `Composition effect` = c(0.0208, 0.0756, 0.0191, 0.0565))
#   rownames(results_dfl_deco) <- rownames(published_results_FLF_table_5) <- 1:4
#   names(published_results_FLF_table_5) <- names(results_dfl_deco)
#   testthat::expect_equal(results_dfl_deco,
#                          published_results_FLF_table_5,
#                          tolerance = 0.0075)
# })
#
#
# test_that("same results as in lecture exercises", {
#
#   # get and prepare data
#   #gsoep29 <- readstata13::read.dta13("data-raw/gsoep29.dta")
#   total_n <- nrow(gsoep29)
#   gsoep29$age <- 2012 - gsoep29$bcgeburt
#   gsoep29 <- subset(gsoep29, age >=25 & age <= 55)
#   total_n - nrow(gsoep29)
#
#   gsoep29$wage <- with(gsoep29, ifelse(labgro12 > 0 & bctatzeit > 0, labgro12 / (bctatzeit * 4.3), NA))
#   gsoep29$lnwage <- log(gsoep29$wage)
#
#   gsoep29$schooling <- with(gsoep29, ifelse(bcbilzeit > 0, bcbilzeit, NA))
#   gsoep29$ft_experience <- with(gsoep29, ifelse(expft12 >= 0, expft12, NA))
#   gsoep29$ft_experience2 <- with(gsoep29, ifelse(expft12 >= 0, expft12^2, NA))
#
#   summary(gsoep29[c("bcsex", "wage", "lnwage", "schooling", "ft_experience", "ft_experience2")])
#
#   gsoep29$tenure <- with(gsoep29, ifelse(bcerwzeit >= 0, bcerwzeit, NA))
#   gsoep29$ISEI <- with(gsoep29, ifelse(isei12 > 0, isei12, NA))
#
#   gsoep29 <- na.omit(gsoep29[, c("bcsex", "lnwage", "schooling", "ft_experience", "ft_experience2", "tenure", "ISEI", "bcphrf")])
#   browser()
#
#   # without weights
#   ob_deco_results <- ob_deco(formula = lnwage ~ schooling + ft_experience + ft_experience2,
#                              data = gsoep29,
#                              group = bcsex,
#                              weights = bcphrf,
#                              reference_0 = TRUE)
#
#   browser()
#
#   testthat::expect_equal(ob_deco_results$ob_deco$decomposition_terms$Observed_difference[1],
#                          0.250,
#                          tolerance = 0.0075)
#
#   testthat::expect_equal(ob_deco_results$ob_deco$decomposition_terms$Composition_effect[c(1, 3)],
#                          c(0.149, -0.00630),
#                          tolerance = 0.0075)
#
#   testthat::expect_equal(sum(ob_deco_results$ob_deco$decomposition_terms$Composition_effect[4:5]),
#                          0.155,
#                          tolerance = 0.0075)
#
#   testthat::expect_equal(ob_deco_results$ob_deco$decomposition_terms$Structure_effect[1:3],
#                          c(0.101, -0.121, 0.0856),
#                          tolerance = 0.0075)
#
#   testthat::expect_equal(sum(ob_deco_results$ob_deco$decomposition_terms$Structure_effect[4:5]),
#                          0.136,
#                          tolerance = 0.0075)
#   browser()
#   # with weights
#   ob_deco_results_weights <- ob_deco(formula = lnwage ~ schooling + ft_experience + ft_experience2,
#                                      data = gsoep29,
#                                      group = bcsex,
#                                      weights = bcphrf,
#                                      reference_0 = TRUE)
#
#   testthat::expect_equal(ob_deco_results_weights$ob_deco$decomposition_terms$Observed_difference[1],
#                          0.223,
#                          tolerance = 0.0075)
#
#   testthat::expect_equal(ob_deco_results_weights$ob_deco$decomposition_terms$Composition_effect[c(1, 3)],
#                          c(0.114, -0.0162),
#                          tolerance = 0.0075)
#
#   testthat::expect_equal(sum(ob_deco_results_weights$ob_deco$decomposition_terms$Composition_effect[4:5]),
#                          0.130,
#                          tolerance = 0.0075)
#
#   testthat::expect_equal(ob_deco_results_weights$ob_deco$decomposition_terms$Structure_effect[1:3],
#                          c(0.109, -0.126, 0.117),
#                          tolerance = 0.0075)
#
#   testthat::expect_equal(sum(ob_deco_results_weights$ob_deco$decomposition_terms$Structure_effect[4:5]),
#                          0.119,
#                          tolerance = 0.0075)
#
# })
