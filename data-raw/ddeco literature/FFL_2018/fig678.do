
**************************************
*** Program to Produce Figures 6-8 ***
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
log using "$figure\fig678.log" , replace

*** import data and variables ***

global expvar covered nonwhite nmarr ed0-ed1 ed3-ed5 ex1-ex4 ex6-ex9 ///
	occd11-occd60 occd80-occd91 indd1 indd3-indd14 pub

use $data\usmen8890_nocc.dta , clear
append using $data\usmen1416_nocc.dta
append using $data\usmen88rw14_nocc.dta
keep lwage2 time $expvar eweight	// use imputed wages if top-coded
keep if lwage2<=7.4					// limit to below $1636 (2010 dollar)
compress
save $data\temp012.dta , replace

*** calculate RIF ***

quietly {
  forvalues t=0/2 {
	pctile valx=lwage2 [aweight=eweight] if time==`t', nq(100)  // 100 percentile grids
	kdensity lwage2 [aweight=eweight] if time==`t', at(valx) ///
		gen(evalt denst) width(0.065) nograph	// find density at each grid
	
	forvalues qt=5(5)95 {
	  if `t'==0 {
		gen rif_`qt'=.
	  }
	  local qc=`qt'/100
	  replace rif_`qt'=evalt[`qt']+`qc'/denst[`qt'] ///
		if lwage2>=evalt[`qt'] & time==`t'
	  replace rif_`qt'=evalt[`qt']-(1-`qc')/denst[`qt'] ///
		if lwage2<evalt[`qt'] & time==`t'
	}
	drop valx evalt denst
  }
}

forvalues qt=5(5)95 {
	di "doing quantile `qt'"
	local qc=`qt'/100

	*** RIF-decomposition ***

	* decomposition without reweighing [E(X_1|t=1)- E(X_0|t=0)]B_0
	quietly: oaxaca rif_`qt' $expvar [aweight=eweight] if time==0 | time==1, ///
		detail(groupun:covered, groupoth: nonwhite nmarr ex1-ex4 ex6-ex9, /// 
		grouped: ed0-ed1 ed3-ed5, groupocc: occd11-occd60 occd80-occd91, ///
		groupind: pub indd1 indd3-indd14) weight(0) swap by(time) 
	matrix Ra=e(b)
	di "oaxaca 1/3 completed..."

	* composition effects with reweighing [E(X_0|t=1)- E(X_0|t=0)]B_c
	quietly: oaxaca rif_`qt' $expvar [aweight=eweight] if time==0 | time==2, ///
		detail(groupun:covered, groupoth: nonwhite nmarr ex1-ex4 ex6-ex9, /// 
		grouped: ed0-ed1 ed3-ed5, groupocc: occd11-occd60 occd80-occd91, ///
		groupind: pub indd1 indd3-indd14) weight(0) swap by(time) 
	matrix Rc=e(b)
	di "oaxaca 2/3 completed..."

	* get wage structure effects E(X_1|t=1)*[B_1-B_c]
	quietly: oaxaca rif_`qt' $expvar [aweight=eweight] if time==1 | time==2, ///
		detail(groupun:covered, groupoth: nonwhite nmarr ex1-ex4 ex6-ex9, /// 
		grouped: ed0-ed1 ed3-ed5, groupocc: occd11-occd60 occd80-occd91, ///
		groupind: pub indd1 indd3-indd14) weight(0) swap by(time) 
	matrix Rw=-1*e(b)
	di "oaxaca 3/3 completed..."
	
	*** collect results ***

	* w/ reweighting
	matrix Rt4=Ra[1,3],Rc[1,4],Rw[1,5],Rc[1,6..10],Rc[1,5],Rw[1,11..16],Rw[1,4], ///
		Rc[1,6]+Rw[1,11],Rc[1,7]+Rw[1,12],Rc[1,8]+Rw[1,13],Rc[1,9]+Rw[1,14], ///
		Rc[1,10]+Rw[1,15],`qc'
	matrix Rt4qtau=(nullmat(Rt4qtau)\Rt4)
}

*** export results and rename variables ***

clear
svmat Rt4qtau, name(b)
rename b1 Dy
rename b2 Dx
rename b3 Dw
rename b4 EXcov
rename b5 EXoth
rename b6 EXeduc
rename b7 EXocc
rename b8 EXind
rename b9 EXres
rename b10 UNcov
rename b11 UNoth
rename b12 UNeduc
rename b13 UNocc
rename b14 UNind
rename b15 UNcons
rename b16 UNres
rename b17 TOTcov
rename b18 TOToth
rename b19 TOTed
rename b20 TOTocc
rename b21 TOTind
rename b22 qtau
gen TOTDx=Dx+EXres
gen TOTDw=UNcov+UNed+UNocc+UNoth+UNind
sort qtau

*** plot ***

set scheme s1color
global option1 xlab(0(0.2)1.0) ylab(-0.1(0.05)0.15) yscale(range(-0.12 0.2)) ///
	yline(0.0, lw(thin) lc(black) lp(dash)) ///
	xti("Quantile") yti("Log Wage Change", margin(0 2 0 0))
global option2 pos(11) ring(0) col(1) region(lstyle(none)) symxsize(8) ///
	keygap(2) textwidth(25)

twoway (connected Dy qtau, m(d) clw(medium) mc(red) lc(red)) ///
	(connected Dw qtau, m(t) lp(dash) clw(medium) mc(green) lc(green)) ///
	(connected Dx qtau, m(o) lp(shortdash) clw(medium) mc(blue) lc(blue)) ///
	, $option1 subtitle("A. Change in Log Wages 2014/16-1988/90") ///
	legend(lab(1 "Total Change") lab(2 "Wage Structure") lab(3 "Composition") ///
	$option2) saving(fig6a,replace)

graph twoway (connected TOTcov qtau, m(d) clw(medium) mc(blue) lc(blue)) ///
    (connected TOTed qtau, m(o) clw(medium) mc(green) lc(green)) ///
    (connected TOTocc qtau, m(t) clw(medium) mc(red) lc(red)) ///
    (connected TOTind qtau, m(s) clw(medium) mc(magenta) lc(magenta)) ///
    (connected TOToth qtau, m(+) clw(medium) lp(dash) mc(black) lc(black)) ///
    , $option1 subtitle("B. Detailed Total Effects")  ///
    legend(lab(1 "Union") lab(2 "Education") lab(3 "Occupations") ///
	lab(4 "Industries") lab(5 "Other") $option2) saving(fig6b,replace)

graph twoway (connected Dx qtau, m(d) clw(medium) mc(blue) lc(blue)) ///
    (connected TOTDx qtau, m(o) clw(medium) mc(green) lc(green)) ///
    (connected EXres qtau, m(t) clw(medium) mc(red) lc(red)) ///
    , $option1 subtitle("A. Aggregate Composition Effects") ///
    legend(lab(1 "Total") lab(2 "Explained") lab(3 "Spec. Error") $option2) ///
    saving(fig7a,replace)

graph twoway (connected EXcov qtau, m(d) clw(medium) mc(blue) lc(blue)) ///
    (connected EXed qtau, m(o) clw(medium) mc(green) lc(green)) ///
    (connected EXocc qtau, m(t) clw(medium) mc(red) lc(red)) ///
    (connected EXind qtau, m(s) clw(medium) mc(magenta) lc(magenta)) ///
    (connected EXoth qtau, m(+) clw(medium) lp(dash) mc(black) lc(black)) ///
    , $option1 subtitle("B. Detailed Composition Effects") ///
    legend(lab(1 "Union") lab(2 "Education") lab(3 "Occupations") ///
	lab(4 "Industries") lab(5 "Other") $option2) saving(fig7b,replace)

graph twoway (connected Dw qtau, m(d) clw(medium) mc(blue) lc(blue)) ///
    (connected TOTDw qtau, m(o) clw(medium) mc(green) lc(green)) ///
    (connected UNcons qtau, m(t) clw(medium) mc(red) lc(red)) ///
    , $option1 subtitle("A. Aggregate Wage Structure Effects") ///
    legend(lab(1 "Total") lab(2 "Explained") lab(3 "Residual") $option2) ///
    saving(fig8a,replace)

graph twoway (connected UNcov qtau, m(d) clw(medium) mc(blue) lc(blue)) ///
    (connected UNed qtau, m(o) clw(medium) mc(green) lc(green)) ///
    (connected UNocc qtau, m(t) clw(medium) mc(red) lc(red)) ///
    (connected UNind qtau, m(s) clw(medium) mc(magenta) lc(magenta)) ///
    (connected UNoth qtau, m(+) clw(medium) lp(dash) mc(black) lc(black)) ///
    , $option1 subtitle("B. Detailed Wage Structure Effects") ///
    legend(lab(1 "Union") lab(2 "Education") lab(3 "Occupations") ///
	lab(4 "Industries") lab(5 "Other") $option2) saving(fig8b,replace)

*** combine and export ***

graph combine fig6a.gph fig6b.gph, xsize(10) saving($figure\fig6.gph, replace)
graph export $figure\fig6.pdf, replace

graph combine fig7a.gph fig7b.gph, xsize(10) saving($figure\fig7.gph, replace)
graph export $figure\fig7.pdf, replace

graph combine fig8a.gph fig8b.gph, xsize(10) saving($figure\fig8.gph, replace)
graph export $figure\fig8.pdf, replace

rm $data\temp012.dta

log close
