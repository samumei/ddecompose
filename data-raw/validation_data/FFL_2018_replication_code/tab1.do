
**********************************
*** Program to Produce Table 1 ***
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
log using "$table\tab1.log" , replace

*** define program to run unconditional quantile regression using RIF ***

capture program drop rifuqr	
program define rifuqr
	local qt=`1'		//quantile of interest
	local qc=`qt'/100
	local width=`2'		//kernel density bandwidth
	local expvar `3'	//RIF expvar: call with double quotes
		
	* calculate density
	pctile valx=lwage1 [aweight=eweight], nq(100)
	kdensity lwage1 [aweight=eweight], at(valx) gen(evalt denst) ///
		width(`width') nograph
	
	* calculate RIF
	gen rif_`qt'=evalt[`qt']+`qc'/denst[`qt'] if lwage1>=evalt[`qt']  
	replace rif_`qt'=evalt[`qt']-(1-`qc')/denst[`qt'] if lwage1<evalt[`qt'] 
	drop valx evalt denst
	
	* RIF regression
	reg rif_`qt' `expvar' [aweight=eweight]   // HOW TO ADD BOOTSTRAP?
end

*** run regression ***

global expvar covered nonwhite nmarr ed0-ed1 ed3-ed5 ex1-ex4 ex6-ex9 ///
	occd11-occd60 occd80-occd91 indd1 indd3-indd14 pub

foreach time in 8890 1416 {
	use $data\usmen`time'_nocc.dta , clear
		// Table 1 uses lwage1 (not smoothed at top) instead of lwage2
		
	forvalues qt=10(40)90 {
		rifuqr `qt' 0.065 "$expvar"
		estimate store cf`time'`qt'
	}
}

*** export table ***

esttab cf889010 cf889050 cf889090 cf141610 cf141650 cf141690 ///
	using $table\tab1.csv , replace unstack label b(%5.3f) se(%5.3f) ar2 ///
	nogaps nonote nonumbers mtitles(Q10 Q50 Q90 Q10 Q50 Q90) ///
	mgroups("1988/90" "2014/16", pattern(1 0 0 1 0 0)) ///
	title("Table 1. Unconditional Quantile Regression Coefficients on Log Wages")

log close
