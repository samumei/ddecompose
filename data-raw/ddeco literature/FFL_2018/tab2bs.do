
*****************************************************************
*** Program to Produce Table 2 with Bootstrap Standard Errors ***
*****************************************************************

clear
drop _all
set more off

// customize your file path
cd "C:\Users\yiged\Dropbox\Research\PRFFL\Programs"
global data "`c(pwd)'\Data"
global table "`c(pwd)'\Tables"

// customize your log file
capture log close
log using "$table\tab2bs.log" , replace

clear  // empty file to store results
save $table\tab2bs.dta , emptyok replace

*** preliminaries ***

global expvar covered nonwhite nmarr ed0-ed1 ed3-ed5 ex1-ex4 ex6-ex9 ///
	occd11-occd60 occd80-occd91 indd1 indd3-indd14 pub

foreach t in 8890 1416 {
	foreach stat in var gini {
		display  "doing years `t', `stat'"
	
		use $data\usmen`t'_nocc.dta , clear
		gen wage_var=lwage2 if lwage2<=7.4
		gen wage_gini=exp(lwage2) if lwage2<=7.4
				// use imputed wages if top-coded; limit to below $1636 (2010 dollar)
		matrix drop _all

		*** resample and run regressions ***
		
		local reps=100
		local j=1
		while `j'<=`reps' {
			display  "doing rep `j'"
			quietly {
				preserve
				if `j'>1 {
					bsample  // use original sample in 1st loop
				}
				rifreg wage_`stat' $expvar [aweight=eweight], `stat'
				matrix B=(nullmat(B)\e(b))	// accumulate coefficients
				restore
				local j=`j'+1
			}
		}   // end of bootstrap loop
		clear
		svmat double B, name(b)  //export coefficients
		
		quietly {
		
			*** calculate bootstrap standard errors ***
		
			local j=1  // loop over variables in $expvar
			foreach var of varlist b* {
				sum `var'
				gen se`j'=`r(sd)' in 1
				replace se`j'=sqrt(`reps'/(`reps'-1))*se`j' in 1  // adj. degree of freedom
				order se`j', after(`var')	// put s.e right "below" coef
				local j=`j'+1
			}  
			keep in 1  // use coefficients from 1st loop (original sample)
			keep b* se*
			
			xpose, clear varname  //transpose row to column
			rename _varname varname
			
			*** format results ***

			// round coef and s.e. by 0.001
			gen y`t'_`stat'=string(round(v1,0.001))
			
			// add leading zeros
			replace y`t'_`stat'=regexr(y`t'_`stat',"^[.]","0.")
			replace y`t'_`stat'=regexr(y`t'_`stat',"^-[.]","-0.")
		
			// add parentheses to s.e.
			replace y`t'_`stat'="("+y`t'_`stat'+")" if substr(varname,1,2)=="se"
			
			// calculate p-values and add aesthetics next to coef
			gen pval=2*normal(-abs(v1/v1[_n+1])) if substr(varname,1,1)=="b"
			replace y`t'_`stat'=y`t'_`stat'+"*" if pval<=0.1   // * p<=0.1
			replace y`t'_`stat'=y`t'_`stat'+"*" if pval<=0.05  // ** p<=0.05
			replace y`t'_`stat'=y`t'_`stat'+"*" if pval<=0.01  // *** p<=0.01
			keep varname y`t'_`stat'
						
			// accumulate results (columns of table)
			merge 1:1 _n using $table\tab2bs.dta , nogen
			save $table\tab2bs.dta , replace		
		}
	}   // end of stat loop
}	// end of period loop

*** format and export table ***

use $table\tab2bs.dta , clear
order varname y8890_var y1416_var y8890_gini y1416_gini

// format column titles
label var varname "Variables"
label var y8890_var "Variance 1988/90"
label var y1416_var "Variance 2014/16"
label var y8890_gini "Gini 1988/90"
label var y1416_gini "Gini 2014/16"

// format variable labels
replace varname="" if substr(varname,1,2)=="se"
replace varname="Union Covered" if varname=="b1"
replace varname="Non-White" if varname=="b2"
replace varname="Non-Married" if varname=="b3"
replace varname="Primary" if varname=="b4"
replace varname="Some HS" if varname=="b5"
replace varname="Some College" if varname=="b6"
replace varname="College" if varname=="b7"
replace varname="Post-Grad" if varname=="b8"
replace varname="Exper < 5 yrs" if varname=="b9"
replace varname="Exper 5-10 yrs" if varname=="b10"
replace varname="Exper 10-15 yrs" if varname=="b11"
replace varname="Exper 15-20 yrs" if varname=="b12"
replace varname="Exper 25-30 yrs" if varname=="b13"
replace varname="Exper 30-35 yrs" if varname=="b14"
replace varname="Exper 35-40 yrs" if varname=="b15"
replace varname="Exper >= 40 yrs" if varname=="b16"
replace varname="Upper Management" if varname=="b17"
replace varname="Lower Management" if varname=="b18"
replace varname="Engineers & Computer Occ." if varname=="b19"
replace varname="Other Scientists" if varname=="b20"
replace varname="Social Support Occ." if varname=="b21"
replace varname="Lawyers & Doctors" if varname=="b22"
replace varname="Health Treatment Occ." if varname=="b23"
replace varname="Clerical Occ." if varname=="b24"
replace varname="Sales Occ." if varname=="b25"
replace varname="Insur. & Real Estate Sales" if varname=="b26"
replace varname="Financial Sales" if varname=="b27"
replace varname="Service Occ." if varname=="b28"
replace varname="Primary Occ." if varname=="b29"
replace varname="Production Occ." if varname=="b30"
replace varname="Transportation Occ." if varname=="b31"
replace varname="Truckers" if varname=="b32"
replace varname="Agriculture & Mining" if varname=="b33"
replace varname="Hi-Tech Manufac" if varname=="b34"
replace varname="Low-Tech Manufac" if varname=="b35"
replace varname="Wholesale Trade" if varname=="b36"
replace varname="Retail Trade" if varname=="b37"
replace varname="Transportation & Utilities" if varname=="b38"
replace varname="Information except Hi-Tech" if varname=="b39"
replace varname="Financial Activities" if varname=="b40"
replace varname="Hi-Tech Services" if varname=="b41"
replace varname="Business Services" if varname=="b42"
replace varname="Education & Health Services" if varname=="b43"
replace varname="Personal Services" if varname=="b44"
replace varname="Public Admin" if varname=="b45"
replace varname="Public Sector" if varname=="b46"
replace varname="Constant" if varname=="b47"

// export table
export excel using "$table\tab2bs.xls", replace first(varlabel)
rm $table\tab2bs.dta

log close
