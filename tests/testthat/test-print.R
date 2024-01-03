test_that("Print function does not throw an error", {


  mod1 <- log(wage) ~ age + central_city + msa + region + black +
    hispanic + education + afqt + family_responsibility + years_worked_civilian +
    years_worked_military + part_time + industry

  deco_male_as_reference <- ob_deco(formula = mod1,
                                    data = nlys00,
                                    group = female,
                                    reference_0 = FALSE)

  testthat::expect_error(print(deco_male_as_reference), NA)
})




test_that("Print function does not throw an error with reweighting (no SE)", {

  data("nlys00")

  mod1 <- log(wage) ~ age + central_city + msa + region + black +
    hispanic + education + afqt + family_responsibility + years_worked_civilian +
    years_worked_military + part_time + industry

  deco_female_as_reference <- ob_deco(formula = mod1,
                                      data = nlys00,
                                      group = female,
                                      reference_0 = TRUE,
                                      reweighting = TRUE)

  print <- print(deco_female_as_reference)

  testthat::expect_error(print, NA)

  summary_detailed <- print(deco_female_as_reference)
  testthat::expect_error(summary_detailed, NA)

})



test_that("Print function does not throw an error with reweighting RIFREG and SE", {

  data("nlys00")

  mod1 <- log(wage) ~ age + central_city + msa + region + black +
    hispanic + education + afqt + family_responsibility + years_worked_civilian +
    years_worked_military + part_time + industry

  deco_female_as_reference <- ob_deco(formula = mod1,
                                      data = nlys00,
                                      group = female,
                                      reference_0 = TRUE,
                                      rifreg = TRUE,
                                      rifreg_statistic = "variance",
                                      reweighting = TRUE,
                                      bootstrap = TRUE,
                                      bootstrap_iterations = 50)

  print <- print(deco_female_as_reference)

  testthat::expect_error(print, NA)

})


test_that("Print function does not throw an error with multiple quantiles", {

  data("nlys00")

  mod1 <- log(wage) ~ age + central_city + msa + region + black +
    hispanic + education + afqt + family_responsibility + years_worked_civilian +
    years_worked_military + part_time + industry

  deco_female_as_reference <- ob_deco(formula = mod1,
                                      data = nlys00,
                                      group = female,
                                      reference_0 = TRUE,
                                      rifreg_statistic = "quantiles",
                                      rifreg_probs = c(0.1, 0.5, 0.9),
                                      reweighting = TRUE)

  print <- print(deco_female_as_reference)

  testthat::expect_error(print, NA)

})

