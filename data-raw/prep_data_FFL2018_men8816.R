#### Prepare data ----


# Reading the data
library("dplyr")
library("readstata13")


men8816 <- read.dta13("data-raw/morgm_all8816.dta")

# Subsetting the data by year and non-missing values of lwage1 & uhrswk
#men8816 <- subset(men8816, year >= 88 & year <= 90)
men8816 <- subset(men8816, !is.na(lwage1) & !is.na(uhrswk))

# Generating covariates

## Unmarried
men8816$nmarr <- 1 - men8816$marr

## Public sector
men8816 <- men8816 %>%
  mutate(
    pub = case_when(
      class >= 1 & class <= 3 ~ 1,
      class >= 4 & class <= 5 ~ 0,
      class == 0 ~ 0,
      TRUE ~ NA_integer_
    )
  ) %>%
  mutate(
    pub = case_when(
      year == 88 & class == 1 ~ 0,
      year == 88 & class == 2 ~ 1,
      year >= 89 & year <= 90 & class == 1 ~ 0,
      year >= 89 & year <= 90 & class >= 2 & class <= 4 ~ 1,
      TRUE ~ pub
    )
  )



# Education
men8816$ed0 <- ifelse(men8816$educ < 9, 1, 0)
men8816$ed1 <- ifelse(men8816$educ >= 9 & men8816$educ < 12, 1, 0)
men8816$ed2 <- ifelse(men8816$educ >= 12 & men8816$educ < 13, 1, 0)
men8816$ed3 <- ifelse(men8816$educ >= 13 & men8816$educ <= 15, 1, 0)
men8816$ed4 <- ifelse(men8816$educ == 16, 1, 0)
men8816$ed5 <- ifelse(men8816$educ > 16, 1, 0)

# Experience
men8816$exper <- men8816$age - men8816$educ - 6
men8816$ex1 <- as.numeric(men8816$exper < 5)
men8816$ex2 <- as.numeric(men8816$exper >= 5 & men8816$exper < 10)
men8816$ex3 <- as.numeric(men8816$exper >= 10 & men8816$exper < 15)
men8816$ex4 <- as.numeric(men8816$exper >= 15 & men8816$exper < 20)
men8816$ex5 <- as.numeric(men8816$exper >= 20 & men8816$exper < 25)
men8816$ex6 <- as.numeric(men8816$exper >= 25 & men8816$exper < 30)
men8816$ex7 <- as.numeric(men8816$exper >= 30 & men8816$exper < 35)
men8816$ex8 <- as.numeric(men8816$exper >= 35 & men8816$exper < 40)
men8816$ex9 <- as.numeric(men8816$exper >= 40)



## Separate assignment for years 88-90
df_88_90 <- dplyr::filter(men8816, year >= 88, year <= 90)


# Public sector
df_88_90$pub <- NA

df_88_90$pub[df_88_90$class == 1 & df_88_90$year == 88] <- 0
df_88_90$pub[df_88_90$class == 2 & df_88_90$year == 88] <- 1
df_88_90$pub[df_88_90$class == 1 & df_88_90$year >= 89 & df_88_90$year <= 90] <- 0
df_88_90$pub[df_88_90$class >= 2 & df_88_90$class <= 4 & df_88_90$year >= 89 & df_88_90$year <= 90] <- 1


# Occupation
df_88_90$occd11 <- as.numeric(with(df_88_90, (occ3 >= 1 & occ3 <= 13) | (occ3 == 19)))
df_88_90$occd12 <- as.numeric(with(df_88_90, (occ3 >= 14 & occ3 <= 18) | (occ3 >= 20 & occ3 <= 37) | (occ3 >= 473 & occ3 <= 476)))
df_88_90$occd21 <- as.numeric(with(df_88_90, (occ3 >= 43 & occ3 <= 68) | (occ3 >= 213 & occ3 <= 218) | (occ3 == 229)))
df_88_90$occd22 <- as.numeric(with(df_88_90, (occ3 >= 69 & occ3 <= 83) | (occ3 >= 166 & occ3 <= 173) | (occ3 >= 223 & occ3 <= 225) | (occ3 == 235)))
df_88_90$occd23 <- as.numeric(with(df_88_90, (occ3 >= 113 & occ3 <= 165) | (occ3 >= 174 & occ3 <= 177) | (occ3 >= 183 & occ3 <= 199) | (occ3 == 234) | (occ3 == 228)))
df_88_90$occd24 <- as.numeric(with(df_88_90, (occ3 >= 84 & occ3 <= 85) | (occ3 >= 178 & occ3 <= 179)))
df_88_90$occd25 <- as.numeric(with(df_88_90, (occ3 >= 86 & occ3 <= 106) | (occ3 >= 203 & occ3 <= 208)))
df_88_90$occd30 <- as.numeric(with(df_88_90, occ3 >= 303 & occ3 <= 389))
df_88_90$occd40 <- as.numeric(with(df_88_90, (occ3 >= 243 & occ3 <= 252) | (occ3 >= 256 & occ3 <= 285)))
df_88_90$occd41 <- as.numeric(with(df_88_90, occ3 == 253 | occ3 == 254))
df_88_90$occd42 <- as.numeric(with(df_88_90, occ3 == 255))
df_88_90$occd50 <- as.numeric(with(df_88_90, occ3 >= 403 & occ3 <= 470))
df_88_90$occd60 <- as.numeric(with(df_88_90, occ3 >= 477 & occ3 <= 499))
df_88_90$occd70 <- as.numeric(with(df_88_90, (occ3 >= 503 & occ3 <= 617) | (occ3 >= 863 & occ3 <= 869)))
df_88_90$occd80 <- as.numeric(with(df_88_90, (occ3 >= 633 & occ3 <= 799) | occ3 == 873 | occ3 == 233))
df_88_90$occd90 <- as.numeric(with(df_88_90, occ3 == 803 | (occ3 >= 808 & occ3 <= 859) | (occ3 >= 875 & occ3 <= 889) | (occ3 >= 226 & occ3 <= 227)))
df_88_90$occd91 <- as.numeric(with(df_88_90, occ3 >= 804 & occ3 <= 806))


# Sum of occupation categories
df_88_90$occsum <- rowSums(df_88_90[, grep('occd', names(df_88_90))], na.rm = TRUE)

# Flag for any occupation
df_88_90$occun <- as.numeric(rowSums(df_88_90[, grep('occd', names(df_88_90))] == 1, na.rm = TRUE) >= 1)


# Industry
df_88_90$indd1 <- as.numeric(with(df_88_90, ind3 >= 10 & ind3 <= 50))
df_88_90$indd2 <- as.numeric(with(df_88_90, ind3 == 60))
df_88_90$indd3 <- as.numeric(with(df_88_90, ind3 == 310 | (ind3 >= 321 & ind3 <= 322) | (ind3 >= 340 & ind3 <= 372) | (ind3 >= 180 & ind3 <= 192) | (ind3 >= 210 & ind3 <= 212)))
df_88_90$indd4 <- as.numeric(with(df_88_90, (ind3 >= 100 & ind3 <= 162) | (ind3 >= 200 & ind3 <= 201) | (ind3 >= 220 & ind3 <= 301) | (ind3 >= 311 & ind3 <= 320) | (ind3 >= 331 & ind3 <= 332) | (ind3 >= 380 & ind3 <= 392)))
df_88_90$indd5 <- as.numeric(with(df_88_90, ind3 >= 500 & ind3 <= 571))
df_88_90$indd6 <- as.numeric(with(df_88_90, ind3 >= 580 & ind3 <= 691 & ind3 != 641))
df_88_90$indd7 <- as.numeric(with(df_88_90, (ind3 >= 400 & ind3 <= 432) | (ind3 >= 460 & ind3 <= 472)))
df_88_90$indd8 <- as.numeric(with(df_88_90, (ind3 >= 171 & ind3 <= 172) | ind3 == 852))
df_88_90$indd9 <- as.numeric(with(df_88_90, ind3 >= 700 & ind3 <= 712))
df_88_90$indd10 <- as.numeric(with(df_88_90, (ind3 >= 440 & ind3 <= 442) | (ind3 >= 732 & ind3 <= 740) | ind3 == 882))
df_88_90$indd11 <- as.numeric(with(df_88_90, (ind3 >= 721 & ind3 <= 731) | (ind3 >= 741 & ind3 <= 760) | ind3 == 890 | ind3 == 892))
df_88_90$indd12 <- as.numeric(with(df_88_90, (ind3 >= 812 & ind3 <= 872 & ind3 != 852) | ind3 == 891))
df_88_90$indd13 <- as.numeric(with(df_88_90, (ind3 >= 761 & ind3 <= 802) | (ind3 >= 880 & ind3 <= 881) | ind3 == 641))
df_88_90$indd14 <- as.numeric(with(df_88_90, ind3 >= 900 & ind3 <= 932))

# Sum of industry categories
df_88_90$indsum <- rowSums(df_88_90[, grep('indd', names(df_88_90))], na.rm = TRUE)

# Flag for any industry
df_88_90$indun <- as.numeric(rowSums(df_88_90[, grep('indd', names(df_88_90))] == 1, na.rm = TRUE) >= 1)

# Base case
df_88_90$base <- as.numeric(with(df_88_90, (covered == 0 & nonwhite == 0 & marr == 1 & ed2 == 1 & ex5 == 1 & occd70 == 1 & indd2 == 1)))



# Stochastically impute from Pareto for top-coded observations
set.seed(5783495)

# Year 88-90: Piketty and Saez (top 1%)
if (nrow(df_88_90) > 0) {
  df_88_90 <- df_88_90 %>%
    mutate(
      ranuni = ifelse(topcode == 1, runif(n()), NA),
      pareto = ifelse(!is.na(ranuni), 1 / (ranuni^(1 / 2.05)), NA),
      lwage2 = ifelse(topcode == 0, lwage, lwage + log(pareto))
    )
}


## Assignment for year 14-16
df_14_16 <- men8816 %>%
  dplyr::filter(year >= 114, year <= 116)

# Public sector
df_14_16$pub <- as.numeric(df_14_16$class >= 1 & df_14_16$class <= 3)
df_14_16$pub[df_14_16$class >= 4 & df_14_16$class <= 5] <- 0


# Occupation
df_14_16$occd11 <- as.numeric(with(df_14_16, (occ3 >= 10 & occ3 < 200) | (occ3 == 430)))
df_14_16$occd12 <- as.numeric(with(df_14_16, (occ3 >= 200 & occ3 < 1000) & occ3 != 430))
df_14_16$occd21 <- as.numeric(with(df_14_16, (occ3 >= 1000 & occ3 <= 1560)))
df_14_16$occd22 <- as.numeric(with(df_14_16, (occ3 >= 1600 & occ3 < 2000)))
df_14_16$occd23 <- as.numeric(with(df_14_16, (occ3 >= 2000 & occ3 < 2100) | (occ3 >= 2140 & occ3 < 3000)))
df_14_16$occd24 <- as.numeric(with(df_14_16, (occ3 >= 2100 & occ3 <= 2110) | (occ3 == 3010 | occ3 == 3060)))
df_14_16$occd25 <- as.numeric(with(df_14_16, (occ3 == 3000) | (occ3 >= 3030 & occ3 <= 3050) | (occ3 >= 3110 & occ3 <= 3540)))
df_14_16$occd30 <- as.numeric(with(df_14_16, (occ3 >= 5000 & occ3 <= 5930)))
df_14_16$occd40 <- as.numeric(with(df_14_16, (occ3 >= 4700 & occ3 <= 4960 & occ3 != 4810 & occ3 != 4820 & occ3 != 4920)))
df_14_16$occd41 <- as.numeric(with(df_14_16, (occ3 == 4810 | occ3 == 4920)))
df_14_16$occd42 <- as.numeric(with(df_14_16, (occ3 == 4820)))
df_14_16$occd50 <- as.numeric(with(df_14_16, (occ3 >= 3600 & occ3 < 4700)))
df_14_16$occd60 <- as.numeric(with(df_14_16, (occ3 >= 6000 & occ3 <= 6130)))
df_14_16$occd70 <- as.numeric(with(df_14_16, (occ3 >= 6200 & occ3 <= 7630)))
df_14_16$occd80 <- as.numeric(with(df_14_16, (occ3 >= 7700 & occ3 <= 8965)))
df_14_16$occd90 <- as.numeric(with(df_14_16, (occ3 >= 9000 & occ3 <= 9750 & occ3 != 9130)))
df_14_16$occd91 <- as.numeric(with(df_14_16, (occ3 == 9130)))
df_14_16$occun <- as.numeric(with(df_14_16, (occd11 == 1 | occd12 == 1 | occd21 == 1 | occd22 == 1 | occd23 == 1 | occd24 == 1 | occd25 == 1 | occd30 == 1 | occd40 == 1 | occd41 == 1 | occd42 == 1 | occd50 == 1 | occd60 == 1 | occd70 == 1 | occd80 == 1 | occd90 == 1 | occd91 == 1)))


# Sum of occupation categories
df_14_16$occsum <- rowSums(df_14_16[, grep('occd', names(df_14_16))], na.rm = TRUE)

# Flag for any occupation
df_14_16$occun <- as.numeric(rowSums(df_14_16[, grep('occd', names(df_14_16))] == 1, na.rm = TRUE) >= 1)

# Industry
df_14_16$indd1 <- as.numeric(with(df_14_16, (ind3 >= 170 & ind3 <= 490)))
df_14_16$indd2 <- as.numeric(with(df_14_16, (ind3 == 770)))
df_14_16$indd3 <- as.numeric(with(df_14_16, ((ind3 >= 3360 & ind3 <= 3690) | (ind3 >= 2170 & ind3 <= 2390) | ind3 == 3960 | ind3 == 3180)))
df_14_16$indd4 <- as.numeric(with(df_14_16, ((ind3 >= 2470 & ind3 <= 3170) | (ind3 >= 3190 & ind3 <= 3290) | (ind3 >= 3770 & ind3 <= 3990 & ind3 != 3960) | (ind3 >= 1070 & ind3 <= 2090))))
df_14_16$indd5 <- as.numeric(with(df_14_16, (ind3 >= 4070 & ind3 <= 4590)))
df_14_16$indd6 <- as.numeric(with(df_14_16, (ind3 >= 4670 & ind3 <= 5790)))
df_14_16$indd7 <- as.numeric(with(df_14_16, ((ind3 >= 6070 & ind3 <= 6390) | (ind3 >= 570 & ind3 <= 690))))
df_14_16$indd8 <- as.numeric(with(df_14_16, ((ind3 >= 6470 & ind3 <= 6480) | (ind3 >= 6570 & ind3 <= 6670) | (ind3 >= 6770 & ind3 <= 6780))))
df_14_16$indd9 <- as.numeric(with(df_14_16, (ind3 >= 6870 & ind3 <= 7190)))
df_14_16$indd10 <- as.numeric(with(df_14_16, ((ind3 >= 7290 & ind3 <= 7460) | ind3 == 6490 | (ind3 >= 6675 & ind3 <= 6695))))
df_14_16$indd11 <- as.numeric(with(df_14_16, ((ind3 >= 7270 & ind3 <= 7280) | (ind3 >= 7470 & ind3 <= 7790))))
df_14_16$indd12 <- as.numeric(with(df_14_16, (ind3 >= 7860 & ind3 <= 8470)))
df_14_16$indd13 <- as.numeric(with(df_14_16, (ind3 >= 8560 & ind3 <= 9290)))
df_14_16$indd14 <- as.numeric(with(df_14_16, (ind3 >= 9370 & ind3 <= 9590)))

df_14_16$base <- as.numeric(with(df_14_16, (covered == 0 & nonwhite == 0 & marr == 1 & ed2 == 1 & ex5 == 1 & occd70 == 1 & indd2 == 1)))


# Sum of industry categories
df_14_16$indsum <- rowSums(df_14_16[, grep('indd', names(df_14_16))], na.rm = TRUE)

# Flag for any industry
df_14_16$indun <- as.numeric(rowSums(df_14_16[, grep('indd', names(df_14_16))] == 1, na.rm = TRUE) >= 1)

# Base case
df_14_16$base <- as.numeric(with(df_14_16, (covered == 0 & nonwhite == 0 & marr == 1 & ed2 == 1 & ex5 == 1 & occd70 == 1 & indd2 == 1)))


# Year 14-16: Piketty and Saez (top 1%)
if (nrow(df_14_16) > 0) {
  df_14_16 <- df_14_16 %>%
    mutate(
      ranuni = ifelse(topcode == 1, runif(n()), NA),
      pareto = ifelse(!is.na(ranuni) & year >= 114 & year <= 116, 1 / (ranuni^(1 / 1.87)), NA),
      lwage2 = ifelse(topcode == 0, lwage, lwage + log(pareto))
    )
}

# Combine the two data frames back into one
men8816 <- bind_rows(df_88_90, df_14_16)

# Adjusting to monthly deflated data
## Reading CPI data (customize your data file)
cpi_month <- read.dta13("data-raw/cpimonth.dta")

## Merging the CPI data

men8816 <- merge(men8816, cpi_month, by = c("year", "cmonth"), all = FALSE)

## Deflating wages
men8816$lwage1 <- men8816$lwage1 - log(men8816$cpi / men8816$acpi)
men8816$lwage2 <- men8816$lwage2 - log(men8816$cpi / men8816$acpi)

# Rebase to 2010
men8816$lwage1 <- men8816$lwage1 + log(218.1 / 72.6)
men8816$lwage2 <- men8816$lwage2 + log(218.1 / 72.6)



# Select only variables used in FFL examples
sel <- c("year", "lwage", "covered", "nonwhite", "nmarr", "age", "ed",  "ex", "occd", "indd", "pub", "base", "eweight")

# Create the regular expression pattern
pattern <- paste(sel, collapse="|")

# Subset the data
men8816 <- men8816[, grep(pattern, names(men8816))]


# Recode Year
men8816$year <- factor(men8816$year,
                       levels = c(88, 89, 90, 114, 115, 116),
                       labels = c("88-90", "88-90", "88-90", "14-16", "14-16", "14-16"))


# Create variable labels
var_labels <- c(
  lwage = "Log hourly wage (nominal)",
  lwage1 = "Log hourly wage (real)",
  lwage2 = "Log hourly wage (top-coded real)",
  covered = "Union Covered",
  nonwhite = "Non-White",
  nmarr = "Non-Married",
  age = "Age",
  ed0 = "Primary",
  ed1 = "Some HS",
  ed2 = "High School",
  ed3 = "Some College",
  ed4 = "College",
  ed5 = "Post-Grad",
  ex1 = "Exper < 5 yrs",
  ex2 = "Exper 5-10 yrs",
  ex3 = "Exper 10-15 yrs",
  ex4 = "Exper 15-20 yrs",
  ex5 = "Exper 20-25 yrs",
  ex6 = "Exper 25-30 yrs",
  ex7 = "Exper 30-35 yrs",
  ex8 = "Exper 35-40 yrs",
  ex9 = "Exper >= 40 yrs",
  occd11 = "Upper Management",
  occd12 = "Lower Management",
  occd21 = "Engineers & Computer Occ.",
  occd22 = "Other Scientists",
  occd23 = "Social Support Occ.",
  occd24 = "Lawyers & Doctors",
  occd25 = "Health Treatment Occ.",
  occd30 = "Clerical Occ.",
  occd40 = "Sales Occ.",
  occd41 = "Insur. & Real Estate Sales",
  occd42 = "Financial Sales",
  occd50 = "Service Occ.",
  occd60 = "Primary Occ.",
  occd70 = "Construction & Repair Occ.",
  occd80 = "Production Occ.",
  occd90 = "Transportation Occ.",
  occd91 = "Truckers",
  indd1 = "Agriculture & Mining",
  indd2 = "Construction",
  indd3 = "Hi-Tech Manufac",
  indd4 = "Low-Tech Manufac",
  indd5 = "Wholesale Trade",
  indd6 = "Retail Trade",
  indd7 = "Transportation & Utilities",
  indd8 = "Information except Hi-Tech",
  indd9 = "Financial Activities",
  indd10 = "Hi-Tech Services",
  indd11 = "Business Services",
  indd12 = "Education & Health Services",
  indd13 = "Personal Services",
  indd14 = "Public Admin",
  pub = "Public Sector",
  base = "Base group"
)

attr(men8816, "var.labels") <- var_labels

# Saving the data
usethis::use_data(men8816, overwrite = TRUE)





