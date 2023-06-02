
# Compare glm and surveyglm
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
model_fit_glm <- glm(formula, data = data_used, weights = weights, family = quasibinomial(link = "logit"))

cbind(coef(model_fit_glm), coef(model_fit_surveyglm))
testthat::expect_equal(coef(model_fit_glm),
                       coef(model_fit_surveyglm),
                       tolerance = 0.000001)
