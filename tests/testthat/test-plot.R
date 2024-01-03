test_that("Plot function does not throw an error with dfl_deco", {
  set.seed(89342395)

  model_sequential <- log(wage) ~ union + experience + education | experience + education | education

  deco_results  <- dfl_deco(model_sequential,
                            data = men8305[1:1000, ],
                            weights = weights,
                            group = year)

  plot(deco_results)
  testthat::expect_error(plot(deco_results), NA)
})

test_that("Plot function does not throw an error with ob_deco", {

  data("nlys00")

  mod1 <- log(wage) ~ age + central_city + msa + region + black +
    hispanic + education + afqt + family_responsibility + years_worked_civilian +
    years_worked_military + part_time + industry

  deco_male_as_reference <- ob_deco(formula = mod1,
                                    data = nlys00,
                                    group = female,
                                    reweighting = TRUE,
                                    rifreg_statistic = "quantiles",
                                    rifreg_probs = 0.1,
                                    bootstrap = TRUE,
                                    bootstrap_iterations = 10,
                                    reference_0 = FALSE)

  testthat::expect_error(plot(deco_male_as_reference,
                              confidence_bands = FALSE), NA)
  testthat::expect_error(plot(deco_male_as_reference), NA)
})


