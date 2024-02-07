
***********************************
*** Program to Produce Table A1 ***
***********************************

clear
drop _all
set more off

// customize your file path
cd "C:\Users\yiged\Dropbox\Research\PRFFL\Programs"
global data "`c(pwd)'\Data"
global table "`c(pwd)'\Tables"

// customize your log file
capture log close
log using "$table\tabA1.log" , replace

*** import data ***

use $data\usmen8890_nocc.dta , clear 
append using $data\usmen1416_nocc.dta 

*** produce Table A1 ***

// weighted variable means by period
global expvar covered nonwhite nmarr age ed0-ed5 occd11-occd91 indd1-indd14 pub
global weight [aweight=eweight]

collapse (mean) lwage1 $expvar (sd) sdlwage1 = lwage1 $weight, by(time)

// label variables
run $data\varlabel.do
label var sdlwage1 "Std of Log Wages"

// find differences
set obs 3
replace time=2 in 3
foreach var of varlist lwage1 sdlwage1 $expvar {
	quietly: replace `var'=`var'[2]-`var'[1] in 3
}

gen timestr="1980/90" if time==0
replace timestr="2014/16" if time==1
replace timestr="Difference" if time==2

// add blank rows
gen Education = .
gen Occupations = .
gen Industries = .

global expvar2 covered nonwhite nmarr age Education ed0-ed5 ///
	Occupations occd11-occd91 Industries indd1-indd14 pub

// export
estpost tabstat lwage1 sdlwage1 $expvar2 , by(timestr) stat(mean) col(stat) nototal
esttab using $table\tabA1.csv , replace main(mean) unstack label b(%5.3f) ///
	nomtitles nonumbers nostar noobs nogaps nonote onecell ///
	title("Table A1. Sample Means")

log close
