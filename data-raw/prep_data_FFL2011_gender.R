################################################
# Prepare data NLSY79 wage data
# of workers ages 35-43 in 2000
# used in  O'Neill and O'Neill (2006) and
# as an illustration in FFL (2011)

library(dplyr)
library("readstata13")

nlys00 <- read.dta13("data-raw/nlsy00_ind.dta")

# names(nlys00)

# Divide armed force qualification test score by 10 #AFQT percentile score
nlys00$afqt <- nlys00$afqtp89/10

# Industry factor
nlys00$industry <- NA
select_obs <- which(nlys00$indd1 == 1 | nlys00$indd2 == 1 | nlys00$indd7 == 1)
nlys00[select_obs,  "industry"] <- "Primary, Constr., Utilities"
select_obs <- which(nlys00$indd3 == 1 | nlys00$indd4 == 1)
nlys00[select_obs,  "industry"] <- "Manufacturing"
select_obs <- which(nlys00$indd11 == 1 | nlys00$indd13 == 1)
nlys00[select_obs,  "industry"] <- "Education, Health, Public Admin."
select_obs <- which(nlys00$indd5==1 | nlys00$indd6==1 | nlys00$indd8==1 | nlys00$indd9==1 | nlys00$indd10==1 | nlys00$indd12==1)
nlys00[select_obs,  "industry"] <- "Other services"

nlys00$industry <- factor(nlys00$industry, levels=c("Primary, Constr., Utilities",
                                                    "Manufacturing",
                                                    "Education, Health, Public Admin.",
                                                    "Other services"))


dim(nlys00)
sum(is.na(nlys00$industry))
nlys00 <- subset(nlys00,  ind3>0 & ind3<990)
sum(is.na(nlys00$industry))
dim(nlys00)
2655+2654


# Education variable
# Warning: There is an encoding error in the original data set. 4 observation
# are coded as having as HS diploma and a GED.
nlys00$education <- NA
select_obs <- which(nlys00$sch10_12 == 1)
nlys00[select_obs,  "education"] <- "10-12 yrs (no diploma or GED)"
select_obs <- which(nlys00$sch_10 == 1)
nlys00[select_obs,  "education"] <- "<10 yrs"
select_obs <- which(nlys00$diploma_hs == 1 & nlys00$ged_hs == 0)
nlys00[select_obs,  "education"] <- "HS grad (diploma)"
select_obs <- which(nlys00$diploma_hs == 0 & nlys00$ged_hs == 1)
nlys00[select_obs,  "education"] <- "HS grad (GED)"
select_obs <- which(nlys00$diploma_hs == 1 & nlys00$ged_hs == 1)
nlys00[select_obs,  "education"] <- "HS grad (GED)"
select_obs <- which(nlys00$smcol == 1)
nlys00[select_obs,  "education"] <- "Some college"
select_obs <- which(nlys00$bachelor_col == 1)
nlys00[select_obs,  "education"] <- "BA or equiv. degree"
select_obs <- which(nlys00$master_col == 1)
nlys00[select_obs,  "education"] <- "MA or equiv. degree"
select_obs <- which(nlys00$doctor_col == 1)
nlys00[select_obs,  "education"] <- "Ph.D or prof. degree"

nlys00$education <- factor(nlys00$education, levels=c("10-12 yrs (no diploma or GED)",
                                                      "<10 yrs",
                                                      "HS grad (diploma)",
                                                      "HS grad (GED)",
                                                      "Some college",
                                                      "BA or equiv. degree",
                                                      "MA or equiv. degree",
                                                      "Ph.D or prof. degree"))

nlys00 %>% group_by(education) %>%
  dplyr::summarise(school=length(which(sch10_12==1)),
                   less_school=length(which(sch_10==1)),
                   hs=length(which(diploma_hs==1)),
                   hs_ged=length(which(ged_hs==1)),
                   sc=length(which(smcol==1)),
                   ba=length(which(bachelor_col==1)),
                   ma=length(which(master_col==1)),
                   phd=length(which(doctor_col==1)))
subset(nlys00, diploma_hs==1 & ged_hs==1)

# Region
nlys00$region <- "East"
select_obs <- which(nlys00$north_central == 1)
nlys00[select_obs,  "region"] <- "North-central"
select_obs <- which(nlys00$south00 == 1)
nlys00[select_obs,  "region"] <- "South"
select_obs <- which(nlys00$west == 1)
nlys00[select_obs,  "region"] <- "West"

nlys00$region <- factor(nlys00$region, levels=c("East", "North-central", "South", "West"))

nlys00 %>% group_by(region) %>% dplyr::summarise(n=length(which(north_central==1)), s=length(which(south00==1)), w=length(which(west==1)))

# Lifetime withdrawal due to family
nlys00$family_responsibility <-  nlys00$famrspb #family responsibility?

# Life-time work experience
nlys00$years_worked_civilian <- nlys00$wkswk_18
nlys00$years_worked_military <- nlys00$yrsmil78_00
nlys00$part_time <- nlys00$pcntpt_22

# Wage and age
nlys00$wage <- exp(nlys00$lropc00)
nlys00$lwage <- nlys00$lropc00
nlys00$age <- nlys00$age00

# Gender variable
nlys00$female <- as.factor(nlys00$female)
levels(nlys00$female) <- c("no","yes")

# Other dummy variables
nlys00$central_city  <- as.factor(nlys00$ctrlcity)
levels(nlys00$central_city) <- c("no","yes")
nlys00$msa <- as.factor(nlys00$msa)
levels(nlys00$msa) <- c("no","yes")
nlys00$black <- as.factor(nlys00$black)
levels(nlys00$black) <- c("no","yes")
nlys00$hispanic <- as.factor(nlys00$hispanic)
levels(nlys00$hispanic) <- c("no","yes")


# Variable selection
sel_var_org <- c("lropc00", "female", "age00", "msa", "ctrlcity", "north_central",
                 "south00", "west", "hispanic", "black",
                 "sch_10",  "diploma_hs", "ged_hs", "smcol", "bachelor_col",
                 "master_col", "doctor_col", "afqtp89",
                 "famrspb", "wkswk_18", "yrsmil78_00", "pcntpt_22", "industry")

mod1 <- lropc00 ~ female + age00 + msa + ctrlcity + north_central + south00 + west+  hispanic+  black +
  sch_10 +  diploma_hs+  ged_hs + smcol +bachelor_col+ master_col+ doctor_col + afqtp89 +
  famrspb +wkswk_18+ yrsmil78_00+ pcntpt_22 + industry

fit1 <- lm(mod1, nlys00)
summary(fit1)


sel_var <- c("female", "lwage", "wage", "age", "central_city", "msa", "region", "black",
             "hispanic",  "education", "afqt", "family_responsibility",
             "years_worked_civilian", "years_worked_military", "part_time", "industry")
sel_var <- c(sel_var, sel_var_org)
nlys00 <- nlys00[, sel_var]


mod1 <- log(wage) ~ education
mod1 <- lwage ~ female + age + central_city + msa + region + black +
  hispanic + education + afqt + family_responsibility + years_worked_civilian +
  years_worked_military + part_time + industry
fit1 <- lm(mod1, nlys00)
summary(fit1)

mod2 <- lwage ~ age + central_city + msa + region + black +
  hispanic + education + afqt + family_responsibility + years_worked_civilian +
  years_worked_military + part_time + industry
fit2 <- lm(mod2, subset(nlys00, female=="no"))
summary(fit2)

fit3 <- lm(mod2, subset(nlys00, female!="no"))
summary(fit3)

round(cbind(apply(model.matrix(mod1, subset(nlys00, female != "yes")), 2, mean),
      apply(model.matrix(mod1, subset(nlys00, female == "yes")), 2, mean)),3)

nlys00 %>%
  group_by(female) %>%
  dplyr::summarise(mean(lwage)) %>% as.data.frame

nrow(subset(nlys00, female!="yes"))
2655

nrow(subset(nlys00, female=="yes"))
2654




# Save data
#save(nlys00, file="data/nlys00.rda")
usethis::use_data(nlys00, overwrite = TRUE)
