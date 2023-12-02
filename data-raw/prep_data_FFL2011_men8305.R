################################################
# Prepare data

# Setwd
#setwd("C:/Users/gallusse/switchdrive/Uni Basel/R Code/LSE/functions/dfl-rif-deco/")
# setwd("C:/Users/david/switchdrive/Uni Basel/Dissprojekte/7 Reweighting decomposition in R/rif-dfl-deco")

#load("data/men8305.rda")

#load data
library(dplyr)
library("readstata13")
df1 <- read.dta13("data-raw/usmen8385_occ.dta")
df2 <- read.dta13("data-raw/usmen0305_occ.dta")
df1$year <- "1983-1985"
df2$year <- "2003-2005"
men8305 <- rbind(df1, df2)
men8305$year <- as.factor(men8305$year)

men8305$experience <- 1
men8305[which(men8305$ex2==1),"experience"] <- 2
men8305[which(men8305$ex3==1),"experience"] <- 3
men8305[which(men8305$ex4==1),"experience"] <- 4
men8305[which(men8305$ex5==1),"experience"] <- 5
men8305[which(men8305$ex6==1),"experience"] <- 6
men8305[which(men8305$ex7==1),"experience"] <- 7
men8305[which(men8305$ex8==1),"experience"] <- 8
men8305[which(men8305$ex9==1),"experience"] <- 9


men8305$education <- 0
men8305[which(men8305$ed1==1),"education"] <- 1
men8305[which(men8305$ed2==1),"education"] <- 2
men8305[which(men8305$ed3==1),"education"] <- 3
men8305[which(men8305$ed4==1),"education"] <- 4
men8305[which(men8305$ed5==1),"education"] <- 5


men8305$education <- as.factor(men8305$education)
men8305$experience <- as.factor(men8305$experience)
men8305$experience <- dplyr::recode_factor(men8305$experience,
                                              "5"="20-24",
                                              "1"="0-4",
                                              "2"="5-9",
                                              "3"="10-14",
                                              "4"="15-19",
                                              "6"="25-29",
                                              "7"="30-34",
                                              "8"="35-39",
                                              "9"=">=40")

men8305$education <- dplyr::recode_factor(men8305$education,
                                             "2"="High School",
                                             "0"="Elementary",
                                             "1"="HS dropout",
                                             "3"="Some College",
                                             "4"="College",
                                             "5"="Post-graduate")


men8305$union <- as.factor(men8305$covered)
men8305$covered <- NULL
men8305$wage <- exp(men8305$lwage1)

men8305$married <- as.factor(men8305$marr)
men8305$marr <- NULL
men8305$nonwhite <- as.factor(men8305$nonwhite)
men8305$public <- as.factor(men8305$pub)
men8305$pub <- NULL
levels(men8305$married) <- levels(men8305$nonwhite) <- levels(men8305$public) <-  levels(men8305$union) <- c("no","yes")

# Rename weight variable
names(men8305)[names(men8305)=="eweight"] <- "weights"
names(men8305)[names(men8305)=="educ"] <- "education_in_years"


## Check variable and ranges
# do.call("rbind",lapply(split(men8305,men8305[,c("year","experience")]), function(x) range(x$exper)))
# do.call("rbind",lapply(split(men8305,men8305[,c("year","education")]), function(x) range(x$educ)))
#do.call("c",lapply(split(men8305,men8305[,c("year")]), function(x) wtd.mean(x$union,x$eweight)))

men8305 %>%
  group_by(year,experience) %>%
  dplyr::summarise(min=min(age - education_in_years - 6, na.rm=TRUE),
                   max=max(age - education_in_years - 6, na.rm=TRUE))

men8305 %>%
  group_by(year, education) %>%
  dplyr::summarise(min=min(education_in_years, na.rm=TRUE),
                   max=max(education_in_years, na.rm=TRUE))

# Select only variables used in FFL examples
sel <- c("wage", "union", "education", "experience", "married", "nonwhite", "year", "weights")
men8305 <- men8305[, sel]

men8305_full <- men8305
save(men8305_full, file="data-raw/men8305_full.rda")

# Save 10% sample of original data set
set.seed(123)
sel_obs <- sample(1:nrow(men8305), floor(nrow(men8305) / 10))
men8305 <- men8305[sel_obs, ]

#Gewichte normieren
men8305$weights <- men8305$weights/sum(men8305$weights) * length(men8305$weights)

# Save data
#save(men8305,file="data/men8305.rda")
usethis::use_data(men8305, overwrite = TRUE)
