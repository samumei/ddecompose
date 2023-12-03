
*****************************************************
*** Program to Reweight CPS Data 1988-1990 (T=C) ***
*****************************************************

clear
drop _all
set more off

// random number seed requires older version
version 12
*help whatsnew13to14

// customize your file path
cd "C:\Users\yiged\Dropbox\Research\PRFFL\Programs"
global data "`c(pwd)'\Data"

// customize your log file
capture log close
log using "$data\read_cps_88rw14.log" , replace

use $data\usmen8890_nocc.dta , clear
append using $data\usmen1416_nocc.dta

*** make weights comparable across two periods ***

forvalues t=0/1 { 
	su eweight if time==`t'
	replace eweight=`r(N)'*eweight/`r(sum)' if time==`t'
	su eweight if time==`t'
}

preserve
	keep if time==0
	save $data\usmen8890_nocc.dta , replace
restore

preserve
	keep if time==1
	save $data\usmen1416_nocc.dta , replace
restore

*** probit for year effects ***

forvalues d=0(1)5 {
	gen uned`d'=ed`d'*covered 
	gen pubed`d'=ed`d'*pub
	gen marred`d'=ed`d'*marr
}

forvalues x=1(1)9 {
	gen unex`x'=ex`x'*covered 
	gen pubex`x'=ex`x'*pub
}

forvalues d=0(1)4 {
	forvalues x=1(1)9 {
		gen ex`x'ed`d'=ed`d'*ex`x'
	}
}

forvalues in=1(1)14 {
	forvalues x=1(1)9 {
		gen indd`in'ex`x'=ex`x'*indd`in'
	}
}

foreach occup in occd11 occd12 occd21 occd22 occd23 occd24 occd25 ///
	occd30 occd41 occd42 occd50 occd60 occd70 occd80 occd90 occd91 {
	forvalues d=1(1)5 {
		gen `occup'ed`d'=ed`d'*`occup'
	}
}
drop occd24ed1 occd42ed1

logistic time covered nonwhite marr* ed0-ed1 ed3-ed5 ex1-ex4 ex6-ex9 uned* ///
   unex* ex1ed* ex2ed* ex3ed* ex4ed* ex5ed* ex6ed* ex7ed* ex8ed* ex9ed* pub* ///
   indd1 indd1e* indd3* indd4* indd5* indd6* indd7* indd8* indd9* indd10* ///
   indd11* indd13* indd14* occd* [iweight=eweight]
predict py1416, p
replace py1416=0.99 if py1416>0.99

summ time [weight=eweight]
replace eweight=eweight*py1416/(1-py1416)*((1-`r(mean)')/`r(mean)') if time==0

drop uned* unex* pubed* pubex* ex1ed* ex2ed* ex3ed* ex4ed* ex6ed* ex7ed* ///
	ex8ed* ex9ed* marred* occd11e* occd12ed* occd21e* occd22e* occd23e* ///
	occd24e* occd25e* occd30e* occd41e* occd42e* occd50e* occd60e* occd70e* ///
	occd80e* occd90e* occd91e* indd1e* indd2e* indd3e* indd4e* indd5e* indd6e* ///
	indd7e* indd8e* indd9e* indd10e* indd11e* indd12e* indd13e* indd14e*

keep if time==0
replace time=2

compress
save $data\usmen88rw14_nocc.dta , replace

log close
