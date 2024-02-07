
**************************************
*** Program to Produce Table 3 & 4 ***
**************************************

clear
drop _all
set more off

// customize your file path
cd "C:\Users\yiged\Dropbox\Research\PRFFL\Programs"
global data "`c(pwd)'\Data"
global table "`c(pwd)'\Tables"

// customize your log file
capture log close
log using "$table\tab34.log" , replace

*** import data ***

global expvar covered nonwhite nmarr ed0-ed1 ed3-ed5 ex1-ex4 ex6-ex9 ///
	occd11-occd60 occd80-occd91 indd1 indd3-indd14 pub

use $data\usmen8890_nocc.dta , clear
append using $data\usmen1416_nocc.dta
append using $data\usmen88rw14_nocc.dta
keep lwage1 lwage2 time $expvar eweight	  // lwage2: imputed wages if top-coded
keep if lwage2<=7.4					// limit to below $1636 (2010 dollar)

*** calculate RIF for inequality measures ***

gen rif_10=.
gen rif_50=.
gen rif_90=.
gen rif_var=.
gen rif_gini=.

// quantile gaps
forvalues t=0/2 {
	forvalues qt=10(40)90 {
		local qc = `qt'/100
		quietly {
		pctile valx=lwage2 [aweight=eweight] if time==`t', nq(100)  // 100 grids
		kdensity lwage2 [aweight=eweight] if time==`t', at(valx) ///
		  gen(evalt denst) width(0.065) nograph  // find density at desired grid
		replace rif_`qt'=evalt[`qt']+`qc'/denst[`qt'] if lwage2>=evalt[`qt'] & time==`t'
		replace rif_`qt'=evalt[`qt']-(1-`qc')/denst[`qt'] if lwage2<evalt[`qt'] & time==`t'
		drop valx evalt denst
		}
	}
}
gen rif_9010=rif_90-rif_10
gen rif_5010=rif_50-rif_10
gen rif_9050=rif_90-rif_50
	
// variance and gini
gen wage_var=lwage2
gen wage_gini=exp(lwage2)
forvalues t=0/2 {
  foreach stat in var gini {
	rifreg wage_`stat' $expvar [aweight=eweight] if time==`t', ///
	  `stat' retain(rif`t'_`stat')
	replace rif_`stat'=rif`t'_`stat' if time==`t'
  }
}
replace rif_var=rif_var*100
replace rif_gini=rif_gini*100

*** OB decomposition ***

matrix drop _all
foreach stat in 9010 5010 9050 var gini {
	display  "doing stat `stat'"

	// decomposition without reweighing [E(X_1|t=1)- E(X_0|t=0)]B_0
	di "oaxaca 1/3 completed..."
	quietly: oaxaca rif_`stat' $expvar [aweight=eweight] if time==0 | time==1, ///
		detail(groupun:covered, groupoth: nonwhite nmarr ex1-ex4 ex6-ex9, /// 
		grouped: ed0-ed1 ed3-ed5, groupocc: occd11-occd60 occd80-occd91, ///
		groupind: pub indd1 indd3-indd14) weight(0) swap by(time) 
	matrix Ra=e(b)

	// composition effects with reweighing [E(X_0|t=1)- E(X_0|t=0)]B_c
	di "oaxaca 2/3 completed..."
	quietly: oaxaca rif_`stat' $expvar [aweight=eweight] if time==0 | time==2, ///
		detail(groupun:covered, groupoth: nonwhite nmarr ex1-ex4 ex6-ex9, /// 
		grouped: ed0-ed1 ed3-ed5, groupocc: occd11-occd60 occd80-occd91, ///
		groupind: pub indd1 indd3-indd14) weight(0) swap by(time) 
	matrix Rc=e(b)

	// wage structure effects E(X_1|t=1)*[B_1-B_c]
	di "oaxaca 3/3 completed..."
	quietly: oaxaca rif_`stat' $expvar [aweight=eweight] if time==1 | time==2, ///
		detail(groupun:covered, groupoth: nonwhite nmarr ex1-ex4 ex6-ex9, /// 
		grouped: ed0-ed1 ed3-ed5, groupocc: occd11-occd60 occd80-occd91, ///
		groupind: pub indd1 indd3-indd14) weight(0) swap by(time) 
	matrix Rw=-1*e(b)

	// collect results: w/o reweighting
	matrix Rt3=Ra[1,3..16],Ra[1,6]+Ra[1,11],Ra[1,7]+Ra[1,12],Ra[1,8]+Ra[1,13], ///
		Ra[1,9]+Ra[1,14],Ra[1,10]+Ra[1,15]
	matrix Rt3`stat'=(nullmat(Rt3`stat')\Rt3)

	// collect results: w/ reweighting
	matrix Rt4=Ra[1,3],Rc[1,4],Rw[1,5],Rc[1,6..10],Rc[1,5],Rw[1,11..16],Rw[1,4], ///
		Rc[1,6]+Rw[1,11],Rc[1,7]+Rw[1,12],Rc[1,8]+Rw[1,13],Rc[1,9]+Rw[1,14], ///
		Rc[1,10]+Rw[1,15]
	matrix Rt4`stat'=(nullmat(Rt4`stat')\Rt4)
}

*** define program to export and format results ***

capture program drop mat2table
program define mat2table
	local matnum `1'  // matrix to export: 3 or 4
	foreach stat in 9010 5010 9050 var gini {
		clear
		svmat double Rt`matnum'`stat', name(b)  //export coefficients
		keep in 1  // use coef from 1st loop (original sample)
		keep b*
		
		xpose, clear varname  //transpose row to column
		rename _varname varname

		// round coef and s.e. by 0.001
		gen y`stat'=string(round(v1,0.001))
		
		// add leading zeros
		replace y`stat'=regexr(y`stat',"^[.]","0.")
		replace y`stat'=regexr(y`stat',"^-[.]","-0.")

		// add parentheses to s.e.
		replace y`stat'="("+y`stat'+")" if substr(varname,1,2)=="se"

		// accumulate results (columns of table)
		keep varname y`stat'
		merge 1:1 _n using $table\temptable.dta , nogen
		save $table\temptable.dta , replace
	}
end

*** format and export table 3 ***

clear  // empty file to store results
save $table\temptable.dta , emptyok replace
mat2table 3
order varname y9010 y5010 y9050 yvar ygini

// format column titles
label var varname "Variables"
label var y9010 "90-10"
label var y5010 "50-10"
label var y9050 "90-50"
label var yvar "Variance"
label var ygini "Gini"

// format variable labels
replace varname="" if substr(varname,1,2)=="se"
replace varname="Total Change" if varname=="b1"
replace varname="Composition" if varname=="b2"
replace varname="Wage Structure" if varname=="b3"
replace varname="X:Union" if varname=="b4"
replace varname="X:Other" if varname=="b5"
replace varname="X:Education" if varname=="b6"
replace varname="X:Occupation" if varname=="b7"
replace varname="X:Industry" if varname=="b8"
replace varname="W:Union" if varname=="b9"
replace varname="W:Other" if varname=="b10"
replace varname="W:Education" if varname=="b11"
replace varname="W:Occupation" if varname=="b12"
replace varname="W:Industry" if varname=="b13"
replace varname="W:Constant" if varname=="b14"
replace varname="T:Union" if varname=="b15"
replace varname="T:Other" if varname=="b16"
replace varname="T:Education" if varname=="b17"
replace varname="T:Occupation" if varname=="b18"
replace varname="T:Industry" if varname=="b19"

// export table
export excel using $table\tab3.xls , replace first(varlabel)
rm $table\temptable.dta

*** format and export table 4 ***

clear  // empty file to store results
save $table\temptable.dta , emptyok replace
mat2table 4
order varname y9010 y5010 y9050 yvar ygini

// format column titles
label var varname "Variables"
label var y9010 "90-10"
label var y5010 "50-10"
label var y9050 "90-50"
label var yvar "Variance"
label var ygini "Gini"

// format variable labels
replace varname="" if substr(varname,1,2)=="se"
replace varname="Total Change" if varname=="b1"
replace varname="Composition" if varname=="b2"
replace varname="Wage Structure" if varname=="b3"
replace varname="X:Union" if varname=="b4"
replace varname="X:Other" if varname=="b5"
replace varname="X:Education" if varname=="b6"
replace varname="X:Occupation" if varname=="b7"
replace varname="X:Industry" if varname=="b8"
replace varname="X:Specification Error" if varname=="b9"
replace varname="W:Union" if varname=="b10"
replace varname="W:Other" if varname=="b11"
replace varname="W:Education" if varname=="b12"
replace varname="W:Occupation" if varname=="b13"
replace varname="W:Industry" if varname=="b14"
replace varname="W:Constant" if varname=="b15"
replace varname="W:Reweighting Error" if varname=="b16"
replace varname="T:Union" if varname=="b17"
replace varname="T:Other" if varname=="b18"
replace varname="T:Education" if varname=="b19"
replace varname="T:Occupation" if varname=="b20"
replace varname="T:Industry" if varname=="b21"

// export table
export excel using $table\tab4.xls , replace first(varlabel)
rm $table\temptable.dta

log close
