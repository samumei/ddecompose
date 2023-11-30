
***********************************
*** Program to Produce Figure 2 ***
***********************************

clear
drop _all
set more off

// customize your file path
cd "C:\Users\yiged\Dropbox\Research\PRFFL\Programs"
global data "`c(pwd)'\Data"
global figure "`c(pwd)'\Figures"

// customize your log file
capture log close
log using "$figure\fig2.log" , replace

*** import top-coded data ***

use $data\usmen8890_nocc.dta , clear
append using $data\usmen1416_nocc.dta

// use imputed data if top-coded
replace lwage1=lwage2 if lwage2!=.

// minimum wage
gen mw=3.35 if year >=88 & year<=90
replace mw=7.25 if year >=114 & year<=116
gen rlmw=ln(mw)+log(218.1/acpi)

*** generate quantile data ***

global option at(eval) gauss nograph
pctile eval=lwage1 if time==0 | time==1 [aweight=eweight], nq(500)
kdensity lwage1 if time==0 & base==1 [aweight=eweight], gen(eval0 dens0) width(0.08) $option
kdensity lwage1 if time==1 & base==1 [aweight=eweight], gen(eval1 dens1) width(0.08) $option

keep dens* eval*
duplicates drop
sort eval
drop if eval > 5.26

*** plot ***

set scheme s1color
graph twoway (connected dens0 eval0, m(i) clw(medium) lp(dash) mc(blue) lc(blue)) ///
	(connected dens1 eval1, m(i) clw(medium) lp(solid) mc(red) lc(red)) , ///
	yti("Density", marg(0 2 0 0)) ylabel(0(.2)1.2) xti("Log Wages", marg(0 0 0 2)) ///
	xlab(1.5(.5)5) xscale(r(1.26 5.26)) text(1.1 1.75 "Minimum  Wages") ///
	text(0.1 4.5 "Top Coded") xline(1.77 2.02, lw(thin) lc(eltblue)) ///
	xline(1.89 2.22, lw(thin) lc(erose) lp(solid)) legend(lab(1 "1988-90") ///
	lab(2 "2014-16") pos(6) col(2) region(lstyle(none)) ///
	symxsize(8) keygap(1) textwidth(32)) saving($figure\fig2.gph , replace)
graph export $figure\fig2.pdf , replace

log close
