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

