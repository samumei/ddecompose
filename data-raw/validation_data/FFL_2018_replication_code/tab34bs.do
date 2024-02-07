
*********************************************************************
*** Program to Produce Table 3 & 4 with Bootstrap Standard Errors ***
*********************************************************************

clear
drop _all
set more off

// customize your file path
cd "C:\Users\yiged\Dropbox\Research\PRFFL\Programs"
global data "`c(pwd)'\Data"
global table "`c(pwd)'\Tables"

// customize your log file
capture log close
log using "$table\tab34bs.log" , replace

*** import data ***

global expvar covered nonwhite nmarr ed0-ed1 ed3-ed5 ex1-ex4 ex6-ex9 ///
	occd11-occd60 occd80-occd91 indd1 indd3-indd14 pub

use $data\usmen8890_nocc.dta , clear
append using $data\usmen1416_nocc.dta
append using $data\usmen88rw14_nocc.dta
keep lwage2 time $expvar eweight		// lwage2: use imputed wages if top-coded
keep if lwage2<=7.4						// limit to below $1636 (2010 dollar)
compress
save $data\bstemp012.dta , replace

use $data\bstemp012.dta , clear
gen rif_10=.
gen rif_50=.
gen rif_90=.
gen rif_9010=.
gen rif_5010=.
gen rif_9050=.
gen rif_var=.
gen rif_gini=.

*** RIF-decomposition with bootstrap ***

local reps=100
local j=1
matrix drop _all

while `j'<=`reps' {
  preserve
	
  if `j'>1 {
	quietly: bsample, strata(time)  // use original sample in 1st loop
  }

  * calculate RIF for inequality measures
	
  quietly {
  
    // quantile gaps
	forvalues t=0/2 {
	  forvalues qt=10(40)90 {
		local qc = `qt'/100
		pctile valx=lwage2 [aweight=eweight] if time==`t', nq(100)  // 100 percentile grids
		kdensity lwage2 [aweight=eweight] if time==`t', at(valx) ///
		  gen(evalt denst) width(0.065) nograph  // find density at desired grid
		replace rif_`qt'=evalt[`qt']+`qc'/denst[`qt'] if lwage2>=evalt[`qt'] & time==`t'
		replace rif_`qt'=evalt[`qt']-(1-`qc')/denst[`qt'] if lwage2<evalt[`qt'] & time==`t'
		drop valx evalt denst
	  }
    }
	replace rif_9010=rif_90-rif_10
	replace rif_5010=rif_50-rif_10
	replace rif_9050=rif_90-rif_50
	
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
	
  }  // end of quietly

  * OB decomposition

  foreach stat in 9010 5010 9050 var gini {
	display  "doing rep `j', stat `stat'"
	
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
	
	* collect results

	// w/o reweighting
	matrix Rt3=Ra[1,3..16],Ra[1,6]+Ra[1,11],Ra[1,7]+Ra[1,12],Ra[1,8]+Ra[1,13], ///
		Ra[1,9]+Ra[1,14],Ra[1,10]+Ra[1,15]
	matrix Rt3`stat'=(nullmat(Rt3`stat')\Rt3)

	// w/ reweighting
	matrix Rt4=Ra[1,3],Rc[1,4],Rw[1,5],Rc[1,6..10],Rc[1,5],Rw[1,11..16],Rw[1,4], ///
		Rc[1,6]+Rw[1,11],Rc[1,7]+Rw[1,12],Rc[1,8]+Rw[1,13],Rc[1,9]+Rw[1,14], ///
		Rc[1,10]+Rw[1,15]
	matrix Rt4`stat'=(nullmat(Rt4`stat')\Rt4)
  }  // end of inequality measures loop
	
  restore
  local j=`j'+1
}  // end of bootstrap loop


*** define program to export and format results ***

capture program drop mat2table_bs
program define mat2table_bs
  local matnum `1'  // matrix to export: 3 or 4

  foreach stat in 9010 5010 9050 var gini {
	clear
	svmat double Rt`matnum'`stat', name(b)  //export coefficients

	* calculate bootstrap standard errors
	local j=1  // loop over rows of output table
	local reps=100
	foreach var of varlist b* {
		sum `var'
		gen se`j'=`r(sd)' in 1
		replace se`j'=sqrt(`reps'/(`reps'-1))*se`j' in 1  // adj. degree of freedom
		order se`j', after(`var')  // put s.e right "below" coef
		local j=`j'+1
		}
	keep in 1  // use coef from 1st loop (original sample)
	keep b* se*
	
	xpose, clear varname  //transpose row to column
	rename _varname varname
	
	* format results

	// round coef and s.e. by 0.001
	gen y`stat'=string(round(v1,0.001))
	
	// add leading zeros
	replace y`stat'=regexr(y`stat',"^[.]","0.")
	replace y`stat'=regexr(y`stat',"^-[.]","-0.")

	// add parentheses to s.e.
	replace y`stat'="("+y`stat'+")" if substr(varname,1,2)=="se"
	
	// calculate p-values and add aesthetics next to coef
	gen pval=2*normal(-abs(v1/v1[_n+1])) if substr(varname,1,1)=="b"
	replace y`stat'=y`stat'+"*" if pval<=0.1
	replace y`stat'=y`stat'+"*" if pval<=0.05
	replace y`stat'=y`stat'+"*" if pval<=0.01
	keep varname y`stat'
		
	// accumulate results (columns of table)
	merge 1:1 _n using $table\temptablebs.dta , nogen
	save $table\temptablebs.dta , replace		
  }
end

*** format and export table 3 ***
  
clear  // empty file to store results
save $table\temptablebs.dta , emptyok replace
mat2table_bs 3
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
export excel using $table\tab3bs.xls , replace first(varlabel)
rm $table\temptablebs.dta

*** format and export table 4 ***

clear  // empty file to store results
save $table\temptablebs.dta , emptyok replace
mat2table_bs 4
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
export excel using $table\tab4bs.xls , replace first(varlabel)
rm $table\temptablebs.dta
rm $data\bstemp012.dta

log close
