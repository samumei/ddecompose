
**************************************************
*** Program to Import CPS Data 1988-1990 (T=0) ***
**************************************************

clear
drop _all
set more off

// random number seed requires older version
version 12
*help whatsnew13to14

// customize your file path
cd "C:\Users\yiged\Dropbox\Research\PRFFL\Programs"
global data "`c(pwd)'\Data"

// customize your log file
capture log close
log using "$data\read_cps_8890.log" , replace

use $data\morgm_all8816.dta , clear
keep if year>=88 & year<=90		//limit to year 1988-1990
keep if lwage1!=. & uhrswk!=.  	//limit to non-missing wage & usual hours per week

*** generate covariates ***

* unmarried
gen nmarr=(1-marr)

* public sector
gen pub=0 if class==1 & year==88
replace pub=1 if class==2 & year==88
replace pub=0 if class==1 & (year>=89 & year<=90)
replace pub=1 if (class>=2 & class<=4) & (year>=89 & year<=90)

* education
gen ed0=(educ<9)
gen ed1=(educ<12 & educ>=9)
gen ed2=(educ>=12 & educ<13)
gen ed3=(educ>=13 & educ<=15)
gen ed4=(educ==16)
gen ed5=(educ>16)

* experience: exper=age-educ-6
gen ex1=(exper<5)
gen ex2=(exper>=5 & exper<10)
gen ex3=(exper>=10 & exper<15)
gen ex4=(exper>=15 & exper<20)
gen ex5=(exper>=20 & exper<25)
gen ex6=(exper>=25 & exper<30)
gen ex7=(exper>=30 & exper<35)
gen ex8=(exper>=35 & exper<40)
gen ex9=(exper>=40)

* occupation
gen occd11=(occ3>=1 & occ3<=13) | (occ3==19)  //upper management
gen occd12=(occ3>=14 & occ3<=18) | (occ3>=20 & occ3<=37) | (occ3>=473 & occ3<=476)  //lower management 
gen occd21=(occ3>=43 & occ3<=68) | (occ3>=213 & occ3<=218) | (occ3==229)  //engineers & computer specialists 
gen occd22=(occ3>=69 & occ3<=83) | (occ3>=166 & occ3<=173) ///
	| (occ3>=223 & occ3<=225) | (occ3==235)  //scientists
gen occd23=(occ3>=113 & occ3<=165) | (occ3>=174 & occ3<=177) ///
	| (occ3>=183 & occ3<=199) | (occ3==234) | (occ3==228)  //education, social support
gen occd24=(occ3>=84 & occ3<=85) | (occ3>=178 & occ3<=179)  //lawyers & doctors
gen occd25=(occ3>=86 & occ3<=106) | (occ3>=203 & occ3<=208)  //health treatment 
gen occd30=(occ3>=303 & occ3<=389)  //clerical occupations
gen occd40=(occ3>=243 & occ3<=252) | (occ3>=256 & occ3<=285)  //sales occupations
gen occd41=(occ3==253 | occ3==254)  //insurance & real estate sales
gen occd42=(occ3==255) //finance sales
gen occd50=(occ3>=403 & occ3<=470)  //service occupations
gen occd60=(occ3>=477 & occ3<=499)  //primary occ
gen occd70=(occ3>=503 & occ3<=617) | (occ3>=863 & occ3<=869)  //construction & repair
gen occd80=(occ3>=633 & occ3<=799) | (occ3==873 | occ3==233)  //production
gen occd90=(occ3==803) | (occ3>=808 & occ3<=859) | (occ3>=875 & occ3<=889) ///
	| (occ3>=226 & occ3<=227)  //other transportation
gen occd91=(occ3>=804 & occ3<=806)  //truck drivers
gen occsum=occd11+occd12+occd21+occd22+occd23+occd24+occd25+occd30+occd40 ///
	+occd41+occd42+occd50+occd60+occd70+occd80+occd90+occd91 
gen occun=(occd11==1 | occd12==1 | occd21==1 | occd22==1 | occd23==1 ///
	| occd24==1 | occd25==1 | occd30==1 | occd40==1 | occd41==1 | occd42==1 ///
	| occd50==1 | occd60==1 | occd70==1 | occd80==1 | occd90==1 | occd91==1)

* industry
gen indd1=(ind3>=10 & ind3<=50)  //agriculture and mining
gen indd2=(ind3==60)  //construction
gen indd3=(ind3==310) | (ind3>=321 & ind3<=322) | (ind3>=340 & ind3<=372) ///
	| (ind3>=180 & ind3<=192) | (ind3>=210 & ind3<=212)  //hi-tech manufac
gen indd4=(ind3>=100 & ind3<=162) | (ind3>=200 & ind3<=201) ///
	| (ind3>=220 & ind3<=301) | (ind3>=311 & ind3<=320 ) ///
	| (ind3>=331 & ind3<=332) | (ind3>=380 & ind3<=392)  //low-tech manufac
gen indd5=(ind3>=500 & ind3<=571)  //wholesale trade
gen indd6=(ind3>=580 & ind3<=691 & ind3!=641)  //retail trade
gen indd7=(ind3>=400 & ind3<=432) | (ind3>=460 & ind3<=472)  //transportation & utilities
gen indd8=(ind3>=171 & ind3<=172) | (ind3==852)  //information except hi-tech
gen indd9=(ind3>=700 & ind3<=712)  //financial activities
gen indd10=(ind3>=440 & ind3<=442) | (ind3>=732 & ind3<=740) | (ind3==882)  //hi-tech services
gen indd11=(ind3>=721 & ind3<=731) | (ind3>=741 & ind3<=760) | (ind3==890 | ind3==892)  //business services
gen indd12=(ind3>=812 & ind3<=872 & ind3!=852) | (ind3==891)  // education & health
gen indd13=(ind3>=761 & ind3<=802) | (ind3>=880 & ind3<=881) | (ind3==641)  //personal services
gen indd14=(ind3>=900 & ind3<=932)  //public admin
gen indsum=indd1+indd2+indd3+indd4+indd5+indd6+indd7+indd8+indd9+indd10 ///
	+indd11+indd12+indd13+indd14
gen indun=(indd1==1 | indd2==1 | indd3==1 | indd4==1 | indd5==1 | indd6==1 ///
	| indd7==1 | indd8==1 | indd9==1 | indd10==1 | indd11==1 ///
	| indd12==1| indd13==1 | indd14==1)
gen base=(covered==0 & nonwhite==0 & marr==1 & ed2==1 & ex5==1 & occd70==1 & indd2==1)

* label variables
run $data\varlabel.do

*** stochatically impute from Pareto for top-coded obs ***

set seed 5783495  					// Piketty and Saez (top 1%)
gen ranuni=runiform() if topcode==1
gen pareto=1/(ranuni^(1/2.05))
gen lwage2=lwage if topcode==0  	// imputed wage
replace lwage2=lwage+log(pareto) if topcode==1

*** adjust to monthly deflated data ***

sort year cmonth
merge m:1 year cmonth using $data\cpimonth.dta , keep(3) nogen
replace lwage1=lwage1-log(cpi/acpi)
replace lwage2=lwage2-log(cpi/acpi)

*** rebase to 2010 ***

replace lwage1=lwage1+log(218.1/72.6)
replace lwage2=lwage2+log(218.1/72.6)

label var lwage "Log hourly wage (nominal)"
label var lwage1 "Log hourly wage (real)"
label var lwage2 "Log hourly wage (top-coded real)"

quietly: compress
save $data\usmen8890_nocc.dta , replace

log close
