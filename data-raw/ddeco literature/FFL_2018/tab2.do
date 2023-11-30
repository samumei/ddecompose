
**********************************
*** Program to Produce Table 2 ***
**********************************

clear
drop _all
set more off

// customize your file path
cd "C:\Users\yiged\Dropbox\Research\PRFFL\Programs"
global data "`c(pwd)'\Data"
global table "`c(pwd)'\Tables"

// customize your log file
capture log close
log using "$table\tab2.log" , replace

*** run regression ***

global expvar covered nonwhite nmarr ed0-ed1 ed3-ed5 ex1-ex4 ex6-ex9 ///
	occd11-occd60 occd80-occd91 indd1 indd3-indd14 pub

foreach time in 8890 1416 {
	use $data\usmen`time'_nocc.dta , clear
	gen wage_var=lwage2 if lwage2<=7.4
	gen wage_gini=exp(lwage2) if lwage2<=7.4
		// use imputed wages if top-coded; limit to below $1636 (2010 dollar)
	
	foreach stat in var gini {
		rifreg wage_`stat' $expvar [aweight=eweight], `stat'
		estimate store cf`time'`stat'
	}
}

*** export table ***

esttab cf8890var cf1416var cf8890gini cf1416gini ///
	using $table\tab2.csv , replace unstack label b(%5.3f) se(%5.3f) ar2 ///
	nogaps nonote nonumbers mtitles("1988/90" "2014/16" "1988/90" "2014/16") ///
	mgroups("Variance" "Gini", pattern(1 0 1 0)) ///
	title("Table 2. RIF Regression of Inequality Measures")

log close
