test_that("Summary function does not throw an error", {

  data("nlys00")

  mod1 <- log(wage) ~ age + central_city + msa + region + black +
    hispanic + education + afqt + family_responsibility + years_worked_civilian +
    years_worked_military + part_time + industry

  deco_male_as_reference <- ob_deco(formula = mod1,
                                    data = nlys00,
                                    group = female,
                                    reference_0 = FALSE)

  testthat::expect_error(summary(deco_male_as_reference), NA)
})


test_that("Summary function does not throw an error with aggregation", {

  data("nlys00")

  mod1 <- log(wage) ~ age + central_city + msa + region + black +
    hispanic + education + afqt + family_responsibility + years_worked_civilian +
    years_worked_military + part_time + industry

  deco_male_as_reference <- ob_deco(formula = mod1,
                                    data = nlys00,
                                    group = female,
                                    reference_0 = FALSE)

  custom_aggregation <- list(`Age, race, region, etc.` = c("age",
                                                           "blackyes",
                                                           "hispanicyes",
                                                           "regionNorth-central",
                                                           "regionSouth",
                                                           "regionWest",
                                                           "central_cityyes",
                                                           "msayes"),
                             `Education` = c("education<10 yrs",
                                             "educationHS grad (diploma)",
                                             "educationHS grad (GED)",
                                             "educationSome college",
                                             "educationBA or equiv. degree",
                                             "educationMA or equiv. degree",
                                             "educationPh.D or prof. degree"),
                             `AFTQ` = "afqt",
                             `L.T. withdrawal due to family` = "family_responsibility",
                             `Life-time work experience` = c("years_worked_civilian",
                                                             "years_worked_military",
                                                             "part_time"),
                             `Industrial sectors` = c("industryManufacturing",
                                                        "industryEducation, Health, Public Admin.",
                                                        "industryOther services"))
  testthat::expect_error(summary(deco_male_as_reference,
                                 custom_aggregation = custom_aggregation),
                         NA)
})

# test_that("Summary function does not throw an error with aggregation and bootstrapped SE", {
#
#   data("nlys00")
#
#   mod1 <- log(wage) ~ age + central_city + msa + region + black +
#     hispanic + education + afqt + family_responsibility + years_worked_civilian +
#     years_worked_military + part_time + industry
#
#   deco_female_as_reference_bs <- ob_deco(formula = mod1,
#                                          data = nlys00,
#                                          group = female,
#                                          bootstrap = TRUE,
#                                          bootstrap_iterations = 100)
#
#   custom_aggregation <- list(`Age, race, region, etc.` = c("age",
#                                                            "blackyes",
#                                                            "hispanicyes",
#                                                            "regionNorth-central",
#                                                            "regionSouth",
#                                                            "regionWest",
#                                                            "central_cityyes",
#                                                            "msayes"),
#                              `Education` = c("education<10 yrs",
#                                              "educationHS grad (diploma)",
#                                              "educationHS grad (GED)",
#                                              "educationSome college",
#                                              "educationBA or equiv. degree",
#                                              "educationMA or equiv. degree",
#                                              "educationPh.D or prof. degree"),
#                              `AFTQ` = "afqt",
#                              `L.T. withdrawal due to family` = "family_responsibility",
#                              `Life-time work experience` = c("years_worked_civilian",
#                                                              "years_worked_military",
#                                                              "part_time"),
#                              `Industrial sectors` = c("industryManufacturing",
#                                                       "industryEducation, Health, Public Admin.",
#                                                       "industryOther services"))
#   testthat::expect_error(summary(deco_female_as_reference_bs,
#                                  custom_aggregation = custom_aggregation),
#                          NA)
# })


