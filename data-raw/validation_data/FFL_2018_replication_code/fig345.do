
**************************************
*** Program to Produce Figures 3-5 ***
**************************************

clear
drop _all
set more off

// customize your file path
cd "C:\Users\yiged\Dropbox\Research\PRFFL\Programs"
global data "`c(pwd)'\Data"
global figure "`c(pwd)'\Figures"

// customize your log file
capture log close
log using "$figure\fig345.log" , replace

*** define program to run unconditional quantile regression using RIF ***

capture program drop rifuqr	
program define rifuqr
	local qt=`1'		// quantile of interest
	local qc=`qt'/100
	local width=`2'		// kernel density bandwidth
	local expvar `3'	// RIF expvar: call with double quotes
		
	* calculate density
	pctile valx=lwage2 [aweight=eweight], nq(100)	// use imputed wages if top-coded
	kdensity lwage2 [aweight=eweight], at(valx) gen(evalt denst) ///
		width(`width') nograph	// find density at desired grid
	
	* calculate RIF
	gen rif_`qt'=evalt[`qt']+`qc'/denst[`qt'] if lwage2>=evalt[`qt']  
	replace rif_`qt'=evalt[`qt']-(1-`qc')/denst[`qt'] if lwage2<evalt[`qt'] 
	drop valx evalt denst
	
	* RIF regression
	reg rif_`qt' `expvar' [aweight=eweight]
end

*** run regression & save estimates ***

global expvar covered nonwhite nmarr ed0 ed1 ed3 ed4 ed5 ex1 ex2 ex3 ex4 ///
	ex6 ex7 ex8 ex9 occd11 occd12 occd21 occd22 occd23 occd24 occd25 occd30 ///
	occd40 occd41 occd42 occd50 occd60 occd80 occd90 occd91 indd1 indd3 indd4 ///
	indd5 indd6 indd7 indd8 indd9 indd10 indd11 indd12 indd13 indd14 pub
matrix drop _all

foreach time in 8890 1416 88rw14 {
	
	use $data\usmen`time'_nocc.dta , clear 	
	forvalues qt=5(5)95 {
		rifuqr `qt' 0.06 "$expvar"
		matrix T`time'Q`qt'=(`qt',e(b))   
		matrix B`time'=(nullmat(B`time')\T`time'Q`qt')
			// accumulate quantile and coefficients
	}
	
	clear
	svmat double B`time', name(coef)  // export quantile and coefficients
	save $data\temp`time'.dta , replace
}

*** format ***

* combine results
use $data\temp8890.dta , clear
gen time=0
append using $data\temp1416.dta
replace time=1 if time==.
append using $data\temp88rw14.dta
replace time=2 if time==.

rename coef1 qtau
replace qtau=qtau/100

* rename and label variables
local j=2
foreach var of global expvar {
	rename coef`j' `var'
	local j=`j'+1
}
do $data\varlabel.do

rename coef const
label var const "Intercept"

*** plot ***

set scheme s1color
global option xtitle("Quantile") xlabel(0(.2)1) ytitle("") legend(off)

foreach var in $expvar const {
  local title: var label `var'
  graph twoway (connected `var' qtau if time==1, m(t) clw(medthick) lp(solid) mc(red) lc(red)) ///
	(connected `var' qtau if time==2, m(+) clw(medium) lp(dash) mc(dkgreen) lc(dkgreen)) ///
	(connected `var' qtau if time==0, m(o) clw(medthick) lp(longdash) mc(blue) lc(blue)) ///
	, $option subtitle("`title'") saving(`var',replace) nodraw
}

*** combine & export ***

graph combine covered.gph nonwhite.gph nmarr.gph ed0.gph ed1.gph ed3.gph ///
	ed4.gph ed5.gph ex1.gph ex2.gph ex3.gph ex4.gph ex6.gph ex7.gph ex8.gph ///
	, altshrink col(3) xsize(12) ysize(20) saving($figure\fig3.gph, replace)
graph export $figure\fig3.pdf , replace

graph combine occd11.gph occd12.gph occd21.gph occd22.gph occd23.gph ///
	occd24.gph occd25.gph occd30.gph occd40.gph occd41.gph occd42.gph ///
	occd50.gph occd60.gph occd80.gph occd90.gph ///
	, altshrink col(3) xsize(12) ysize(20) saving($figure\fig4.gph, replace)
graph export $figure\fig4.pdf , replace

graph combine indd1.gph indd3.gph indd4.gph indd5.gph indd6.gph indd7.gph indd8.gph ///
	indd9.gph indd10.gph indd11.gph indd12.gph indd13.gph indd14.gph pub.gph const.gph ///
	, altshrink col(3) xsize(12) ysize(20) saving($figure\fig5.gph, replace)
graph export $figure\fig5.pdf , replace

rm $data\temp8890.dta
rm $data\temp1416.dta
rm $data\temp88rw14.dta

log close
