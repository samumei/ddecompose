testthat::test_that("Test GU normalization", {

  data("men8305")
  data0 <- subset(men8305, year == "1983-1985")
  data1 <- subset(men8305, year != "1983-1985")

  m0 <- mean(data0$wage)
  m1 <- mean(data1$wage)

  m0_union <- mean(subset(data0, union == "yes")$wage)
  m0_non_union <- mean(subset(data0, union != "yes")$wage)

  m1_union <- mean(subset(data1, union == "yes")$wage)
  m1_non_union <- mean(subset(data1, union != "yes")$wage)

  union_share0 <- mean(as.numeric(data0$union)) - 1
  union_share1 <- mean(as.numeric(data1$union)) - 1

  mC <- m1_union * (union_share0) + m1_non_union * (1- union_share0)

  intercept_0 <- (m0_union + m0_non_union)/2
  intercept_1 <- (m1_union + m1_non_union)/2
  coef_union_0 <- m0_union - intercept_0
  coef_union_1 <- m1_union - intercept_1
  coef_non_union_0 <- m0_non_union - intercept_0
  coef_non_union_1 <- m1_non_union - intercept_1

  detailed_wage_structure_effect <- c((intercept_1 - intercept_0),
                                      (coef_non_union_1 - coef_non_union_0) * (1-union_share0),
                                      (coef_union_1 - coef_union_0) * union_share0)

  detailed_composition_effect <- c(0,
                                   coef_non_union_1 * ((1-union_share1) - (1-union_share0)),
                                   coef_union_1 * (union_share1 - union_share0))

  expected_decompose <- data.frame(`Composition_effect` = detailed_composition_effect,
                              `Structure_effect` = detailed_wage_structure_effect)
  rownames(expected_decompose) <- c("(Intercept)", "unionno", "unionyes")

  model_decompose <- wage ~ union
  deco_union <- ob_decompose(formula = model_decompose,
                        data = men8305,
                        group = year,
                        normalize_factors = TRUE,
                        reference_0 = FALSE)

  estimated_decompose <-
    deco_union$ob_decompose$decomposition_terms[
      which(rownames(deco_union$ob_decompose$decomposition_terms) %in%
              rownames(expected_decompose)), names(expected_decompose)]

  testthat::expect_equal(estimated_decompose,
                         expected_decompose,
                         tolerance = 0.0000000001)
})
