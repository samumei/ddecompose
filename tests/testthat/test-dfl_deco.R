
test_that("dfl_deco() does not throw an error", {
  set.seed(89342395)
  data("nlys00")
  formula <- log(wage) ~ age + central_city + msa + region + black +
    hispanic + education + afqt + family_responsibility + years_worked_civilian +
    years_worked_military + part_time + industry

  deco_results <- dfl_deco(formula = formula,
                       data = nlys00[1:300,],
                       weights = runif(nrow(nlys00[1:300,]), 0.5, 1.5),
                       group = female)

  testthat::expect_error(deco_results, NA)
})


test_that("dfl_deco_estimate() does not throw an error", {
  set.seed(89342395)
  data("men8385")
  men8305$weights <- men8305$weights/sum(men8305$weights) * length(men8305$weights)
  data_sample <- men8305[sample(1:nrow(men8305), size = 10000), ]
  formula <- Formula::as.Formula(log(wage) ~ union*(education + experience) + education*experience)
  data_used <- model.frame(formula, data_sample, weights = weights, group = year)
  dep_var <- model.response(data_used, "numeric")

  weights <- model.weights(data_used)
  group_variable_name <- "year"
  group_variable <- data_used[, ncol(data_used)]
  names(data_used)[ncol(data_used)] <- "group_variable"

  reference_0 <- TRUE
  reference_group <- ifelse(reference_0, 0, 1)

  statistics = c("quantiles", "mean", "variance", "gini", "iq_range_p90_p10", "iq_range_p90_p50", "iq_range_p50_p10")
  probs = c(1:9)/10
  trimming = TRUE
  trimming_threshold = 0.00005
  reweight_marginals = TRUE
  method = "logit"

  deco_results <- dfl_deco_estimate(formula = formula,
                    dep_var = dep_var,
                    data_used = data_used ,
                    weights = weights,
                    group_variable = group_variable,
                    reference_group = reference_group,
                    method=method,
                    estimate_statistics = estimate_statistics,
                    statistics = statistics,
                    probs = probs,
                    reweight_marginals = reweight_marginals,
                    trimming = trimming,
                    trimming_threshold = trimming_threshold)

  testthat::expect_error(deco_results, NA)

})



test_that("dfl_deco() does not throw an error without estimating statistics", {
  set.seed(89342395)
  data("nlys00")
  formula <- log(wage) ~ age + central_city + msa + region + black +
    hispanic + education + afqt + family_responsibility + years_worked_civilian +
    years_worked_military + part_time + industry

  deco_results <- dfl_deco(formula = formula,
                           data = nlys00[1:300,],
                           weights = runif(nrow(nlys00[1:300,]), 0.5, 1.5),
                           group = female,
                           estimate_statistics = FALSE)

  testthat::expect_error(deco_results, NA)
})


test_that("glm and surveyglm return the same coefficients", {
  set.seed(123)
  data("nlys00")
  formula <- female ~ education + afqt + industry + family_responsibility + part_time
  data_used <- get_all_vars(formula, nlys00)
  data_used$weights <- runif(nrow(data_used), 0.5, 1.5)

  design <- survey::svydesign(~0,
                              data=data_used,
                              weights=~weights)
  model_fit_surveyglm <- survey::svyglm(formula,
                                        data=data_used,
                                        design=design,family=quasibinomial(link="logit"))

  model_fit_glm <- glm(formula, data = data_used, weights = weights, family = binomial(link = "logit"))
  model_fit_quasiglm <- glm(formula, data = data_used, weights = weights, family = quasibinomial(link = "logit"))

  sglm <- summary(model_fit_glm)
  squasiglm <- summary(model_fit_quasiglm)
  ssurveryglm <- summary(model_fit_surveyglm)

  cbind(sglm$coefficients[,2],
        squasiglm$coefficients[,2],
        ssurveryglm$coefficients[,2])

  cbind(coef(model_fit_glm), coef(model_fit_surveyglm))
  testthat::expect_equal(coef(model_fit_glm),
                         coef(model_fit_surveyglm),
                         tolerance = 0.000001)

  cbind(coef(model_fit_quasiglm), coef(model_fit_surveyglm))
  testthat::expect_equal(coef(model_fit_quasiglm),
                         coef(model_fit_surveyglm),
                         tolerance = 0.000001)

})
