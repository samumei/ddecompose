
**************************************************
*** Program to Import CPS Data 2014-2016 (T=1) ***
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
log using "$data\read_cps_1416.log" , replace

use $data\morgm_all8816.dta , clear
keep if year>=114 & year<=116	//limit to year 2014-2016
keep if lwage1!=. & uhrswk!=.  	//limit to non-missing wage & usual hours per week

*** generate covariates ***

* unmarried
gen nmarr=(1-marr)

* public sector
gen pub=class>= 1 & class<=3
replace pub=0 if class>= 4 & class<=5

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
gen occd11=(occ3>=10 & occ3<200) | (occ3==430)  //upper management
gen occd12=(occ3>=200 & occ3<1000) & occ3~=430  //lower management
gen occd21=(occ3>=1000 & occ3<=1560)  //engineers and computer specialists
gen occd22=(occ3>=1600 & occ3<2000)  //scientists
gen occd23=(occ3>=2000 & occ3<2100) | (occ3>=2140 & occ3<3000)  //education, social support 
gen occd24=(occ3>=2100 & occ3<=2110) | (occ3==3010 | occ3==3060)  //lawyers & doctors 
gen occd25=(occ3==3000) | (occ3>=3030 & occ3<=3050)| (occ3>=3110 & occ3<=3540) //health treatment 
gen occd30=(occ3>=5000 & occ3<=5930)  //clerical occupations
gen occd40=(occ3>=4700 & occ3<=4960 & occ3!=4810 & occ3!=4820 & occ3!=4920)  //sales occupations
gen occd41=(occ3==4810 | occ3==4920)  //insurance & real estate sales
gen occd42=(occ3==4820)  //finance sales
gen occd50=(occ3>=3600 & occ3<4700)  //service occupations
gen occd60=(occ3>=6000 & occ3<=6130)  //primary occ
gen occd70=(occ3>=6200 & occ3<=7630)  //construction & repair
gen occd80=(occ3>=7700 & occ3<=8965)  //production
gen occd90=(occ3>=9000 & occ3<=9750 & occ3!=9130)  //other transportation
gen occd91=(occ3==9130)  //truck drivers
gen occsum=occd11+occd12+occd21+occd22+occd23+occd24+occd25+occd30+occd40 ///
	+occd41+occd42+occd50+occd60+occd70+occd80+occd90+occd91 
gen occun=(occd11==1 | occd12==1 | occd21==1 | occd22==1 | occd23==1 ///
	| occd24==1 | occd25==1 | occd30==1 | occd40==1 | occd41==1 | occd42==1 ///
	| occd50==1 | occd60==1 | occd70==1 | occd80==1 | occd90==1 | occd91==1) 

* industry
gen indd1=(ind3>=170 & ind3<=490)  //agriculture and mining
gen indd2=(ind3==770)  //construction
gen indd3=((ind3>=3360 & ind3<=3690) | (ind3>=2170 & ind3<=2390) ///
	| ind3==3960 | ind3==3180)  //hi-tech manufac
gen indd4=((ind3>=2470 & ind3<=3170) | (ind3>=3190 & ind3<=3290) ///
	| (ind3>=3770 & ind3<=3990 & ind3!=3960) | (ind3>=1070 & ind3<=2090))
	// low-tech manufac
gen indd5=(ind3>=4070 & ind3<=4590)  //wholesale trade
gen indd6=(ind3>=4670 & ind3<=5790)  //retail trade
gen indd7=((ind3>=6070 & ind3<=6390) | (ind3>=570 & ind3<=690))  //transportation & utilities
gen indd8=((ind3>=6470 & ind3<=6480) |(ind3>=6570 & ind3<=6670) ///
	| (ind3>=6770 & ind3<=6780))  //information except hi-tech
gen indd9=(ind3>=6870 & ind3<=7190)  //financial activities
gen indd10=((ind3>=7290 & ind3<=7460) | ind3==6490 | (ind3>=6675 & ind3<=6695))  //hi-tech services
gen indd11=((ind3>=7270 & ind3<=7280) | (ind3>=7470 & ind3<=7790))  //business services
gen indd12=(ind3>=7860 & ind3<=8470)  //education & health services
gen indd13=(ind3>=8560 & ind3<=9290)  //personal services
gen indd14=(ind3>=9370 & ind3<=9590)  //public admin

gen indsum=indd1+indd2+indd3+indd4+indd5+indd6+indd7+indd8+indd9+indd10 ///
	+indd11+indd12+indd13+indd14
gen indun=(indd1==1 | indd2==1 | indd3==1 | indd4==1 | indd5==1 | indd6==1 ///
	| indd7==1 | indd8==1 | indd9==1 | indd10==1 | indd11==1 ///
	| indd12==1| indd13==1 | indd14==1)
gen base=(covered==0 & nonwhite==0 & marr==1 & ed2==1 & ex5==1 & occd70==1 & indd2==1)

* label variables
run $data\varlabel.do

quietly: compress
save $data\usmen1416tc_occ.dta , replace

*** stochatically impute from Pareto for top-coded obs ***

set seed 5783495  					// Piketty and Saez (top 1%)
gen ranuni=runiform() if topcode==1
gen pareto=1/(ranuni^(1/1.87)) if year>=114 & year<=116
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
save $data\usmen1416_nocc.dta , replace

log close
