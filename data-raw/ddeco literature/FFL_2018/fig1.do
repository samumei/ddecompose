
***********************************
*** Program to Produce Figure 1 ***
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
log using "$figure\fig1.log" , replace

*** import top-coded data ***

use $data\usmen1416tc_occ.dta , clear
gen byte time=3 if (year >=114 & year<=116)

// adjust to monthly deflated data & rebase to 2010
sort year cmonth
merge m:1 year cmonth using $data\cpimonth.dta , keep(3) nogen
replace lwage1=lwage1-log(cpi/acpi)
replace lwage1=lwage1+log(218.1/72.6)

*** append other data periods ***

append using $data\usmen8890_nocc.dta
replace time=0 if year >=88 & year<=90
append using $data\usmen1416_nocc.dta
replace time=1 if time==. & (year >=114 & year<=116)
append using $data\usmen88rw14_nocc.dta
replace time=2 if time==. & (year >=88 & year<=90)

tab time
replace lwage1=lwage2 if time!=3	// other periods use imputed data
bysort time: sum lwage1 lwage2 [aweight=eweight]

*** find minimum wages ***

// federal min wage
gen fmw=3.35 if year >=88 & year<=90
replace fmw=7.25 if year >=114 & year<=116

// lowest state min wage
gen mimw=3.35 if year >=88 & year<=90
replace mimw=7.25 if year >=114 & year<=116

// highest state min wage
gen mxmw=4.25 if year >=88 & year<=90
replace mxmw=10 if year >=114 & year<=116

foreach mw in fmw mimw mxmw {
	gen rl`mw'=ln(`mw')+ log(218.1/acpi)
}

*** generate quantile data ***

global option at(eval) gauss nograph
pctile eval=lwage1 if time==0 | time==1 [aweight=eweight], nq(1500)
kdensity lwage1 if time==0 [aweight=eweight], gen(eval0 dens0) width(0.06) $option
kdensity lwage1 if time==1 [aweight=eweight], gen(eval1 dens1) width(0.08) $option
kdensity lwage1 if time==2 [aweight=eweight], gen(eval2 dens2) width(0.06) $option
kdensity lwage1 if time==3 [aweight=eweight], gen(eval3 dens3) width(0.08) $option

keep dens* eval*
duplicates drop
sort eval
drop if eval>5.26

*** plot ***

set scheme s1color
graph twoway (connected dens0 eval0, m(i) clw(medium) lp(longdash) mc(blue) lc(blue)) ///
	(connected dens3 eval3, m(i) clw(medium) lp(dash_dot) mc(red) lc(red)) ///
	(connected dens1 eval1, m(i) clw(medium) lp(solid) mc(red) lc(red)) ///
	(connected dens2 eval2, m(i) clw(medium) lp(dash) mc(green) lc(green)) , /// 
	yti("Density", marg(0 2 0 0)) ylabel(0(.1).7) xti("Log Wages", marg(0 0 0 2)) ///
	xlab(1.5(.5)5) xscale(r(1.26 5.26)) text(0.65 1.75 "Minimum Wages") ///
	text(0.1 4.75 "Top Coded") xline(1.77 2.02, lw(thin) lc(eltblue) lp(solid)) ///
	xline(1.89 2.22, lw(thin) lc(erose) lp(solid)) legend(lab(1 "1988-90") ///
	lab(2 "2014-16 top-coded") lab(3 "2014-16") lab(4 "1988-90 rwgt 2014-16") ///
	order(1 4 2 3) pos(6) col(2) region(lstyle(none)) symxsize(8) keygap(1) ///
	textwidth(32)) saving($figure\fig1.gph , replace)
graph export $figure\fig1.pdf , replace

log close
