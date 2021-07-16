//Analysis using prostate dataset available at https://hbiostat.org/data/ 
//Total effects are estimated: 1. Cause specific cumulative incidence functions 
//                             2. Expected loss in life due to a cause of death before time t* (as well as total loss due to both causes)

/* 
Flexible parametric models can be fitted within Stata using the user-written command stpm2 that can be installed as follows:
  ssc install stpm2

To obtain marginal (and non-marginal) estimates using standardisation, the \texttt{standsurv} command must also be installed.
This  can be installed by running,
 net from https://www.pclambert.net/downloads/standsurv  
*/

// first set working directory
// then, load data
use "prostate", clear

// restrict to placebo and high dose estrogen
keep if inlist(rx,1,4)

// update coding (0 for placebo, 1 for DES)
replace rx = cond(rx==1,0,1)
label define lblrx 0 "placebo" 1 "DES"
label values rx lblrx 

// get rid of zero timevar
// replace with half day
replace dtime = 0.5 if dtime==0

// all-cause deaths indicator and event indicator
gen allcause = status != 1
gen eventType = cond(status==1,0,cond(status==2,1,2))
label define causelab 0 "Alive" 1 "Prostate" 2 "Other"
label values eventType causelab

//create categorical variables 
gen hgBinary = hg<12
egen ageCat = cut(age), at(0,60,75,100)
gen normalAct = pf == 1
// create dummy variables for age
tab ageCat, gen(ageCat)

// generate variable timevar for follow-up, which has 121 observations between 0 and 60 months, and  will be used for the predictions (every half month)		
range timevar 0 60 121



//****************   Total effects - 1. Cause specific cumulative incidence functions ****************//
// fit cause specific models (exit at 60)
// for prostate 
stset dtime, failure(eventType==1) exit(time 60)
stpm2 rx normalAct ageCat2 ageCat3 hx hgBinary, scale(hazard) df(4) tvc(rx) dftvc(2) eform
//store model estimates
estimates store prostate

/*
// get HRs for prostate (varies with follow-up as we have assumed a time-dependent effect for treatment)
predict hr_prostate, hrnum(rx 1) timevar(timevar)
li timevar hr_prostate if timevar==12
li timevar hr_prostate if timevar==36
li timevar hr_prostate if timevar==60
*/

// for other causes
stset dtime, failure(eventType==2) exit(time 60)
stpm2 rx normalAct ageCat2 ageCat3 hx hgBinary, scale(hazard) df(3) eform
// store model estimates
estimates store other

/*
// get HRs for other causes (remains constant across follow-up as we have assumed no time-dependent effect for treatment)
predict hr_other, hrnum(rx 1) timevar(timevar)
li timevar hr_other if timevar==12
li timevar hr_other if timevar==36
li timevar hr_other if timevar==60
*/


// obtain predictions of standardised CIFs as well as their difference using standsurv with the option cif and crmodels()
standsurv, crmodels(prostate other) cif at1(rx 0) at2(rx 1) timevar(timevar) ///
            contrast(difference) ci  atvars(CIF0 CIF1) contrastvars(CIF_diff)

// list results at 60 months
// for the CIF from death due to prostate cancer under placebo and under DES with confidence intervals (CIs)
list CIF0_prostate* CIF1_prostate* if timevar==60, noobs abbrev(25)
// for the CIF from death due to other causes under placebo and under DES with confidence intervals (CIs)
list CIF0_other* CIF1_other* if timevar==60, noobs abbrev(25)
// for the differences between treatment arms
list CIF_diff_prostate* if timevar==60, noobs abbrev(25)
list CIF_diff_other*    if timevar==60, noobs abbrev(25)




// plotting the standardized cause-specific CIFs
// for prostate cancer
twoway  (rarea CIF0_prostate_lci CIF0_prostate_uci timevar, color("200 82 0 %20")) ///
         (line CIF0_prostate timevar, color("200 82 0")) ///
         (rarea CIF1_prostate_lci CIF1_prostate_uci timevar, color(black%15)) ///
         (line CIF1_prostate timevar, color(black)) ///
                 , legend(order(2 "Placebo" 4 "DES") cols(1) ring(0) pos(10) size(medium) region(lwidth(none))) ///
                 ylabel(0.0(0.2)1, angle(h) format(%3.1f)) ///
                 xtitle("Time since randomisation (months)", size(medlarge)) ytitle("Cumulative incidence", size(medlarge)) ///
                 title("Prostate cancer", color(black)) ///
                 name(prostate, replace) ///
				 graphregion(color(white)) ///
				 plotregion(margin(b=0 r=0 l=0)) 

// for other causes                 
twoway  (rarea CIF0_other_lci CIF0_other_uci timevar, color("200 82 0 %20")) ///
         (line CIF0_other timevar, color("200 82 0")) ///
         (rarea CIF1_other_lci CIF1_other_uci timevar, color(black%15)) ///
         (line CIF1_other timevar, color(black)) ///
                 , legend(off) ///
                 ylabel(0.0(0.2)1, angle(h) format(%3.1f)) ///
                 xtitle("Time since randomisation (months)", size(medlarge)) ytitle("Cumulative incidence", size(medlarge)) ///
                 title("Other", color(black)) ///
                 name(other, replace) ///
				 graphregion(color(white)) ///
				 plotregion(margin(b=0 r=0 l=0)) 
                 
 graph combine prostate other, nocopies ycommon graphregion(color(white)) name(cif, replace)
 
 graph export Total_Crude.pdf, replace
 
// plotting their difference
// for prostate cancer
twoway  (rarea CIF_diff_prostate_lci CIF_diff_prostate_uci timevar, color(black%15) lpattern(dash dash)) ///
         (line CIF_diff_prostate timevar, color(black) lpattern(dash)) ///
                 , legend(order(2 "Difference") cols(1) ring(0) pos(10) size(medium) region(lwidth(none))) ///
                 ylabel(-0.4(0.2)0.4, angle(h) format(%3.1f)) ///
                 xtitle("Time since randomisation (months)", size(medlarge)) ytitle("Cumulative incidence", size(medlarge)) ///
                 title("Prostate cancer", color(black)) ///
                 name(prostate_dif, replace) ///
				 graphregion(color(white)) ///
				 plotregion(margin(b=0 r=0 l=0)) 

// for other causes                 
twoway  (rarea CIF_diff_other_lci CIF_diff_other_uci timevar, color(black%15) lpattern(dash dash)) ///
         (line CIF_diff_other timevar, color(black) lpattern(dash)) ///
                 , legend(off) ///
                 ylabel(-0.4(0.2)0.4, angle(h) format(%3.1f)) ///
                 xtitle("Time since randomisation (months)", size(medlarge)) ytitle("Cumulative incidence", size(medlarge)) ///
                 title("Other", color(black)) ///
                 name(other_dif, replace) ///
				 graphregion(color(white)) ///
				 plotregion(margin(b=0 r=0 l=0)) 
                 
 graph combine prostate_dif other_dif, nocopies ycommon graphregion(color(white)) name(dif_cif, replace)
 graph export Total_Crude_diff.pdf, replace
 
 //Combined (Figure 2 in the paper)
 graph combine prostate other prostate_dif other_dif, nocopies graphregion(color(white)) cols(2)
 graph export Total_Crude_comb.pdf, replace

 
 
 
//****************   Total effect - 2. Expected loss in life due to a cause of death before time t*  ********************//

// generate time variable at 60 months (5 years) as t*
gen t_rmft60 = 60 in 1

// obtain the standardised expected loss in life due to a cause of death before 60 months under placebo and under DES as well as their difference
standsurv, crmodels(prostate other) cif rmft  ///
    at1(rx 0) at2(rx 1) timevar(t_rmft60) contrast(difference) ci ///
	atvars(RMFT0 RMFT1) contrastvars(RMFT_diff)
	
// list the estimates of the life years lost due to prostate cancer before 60 months
// under placebo
list t_rmft60 RMFT0_prostate* in 1, noobs abb(22) 	
//under DES
list t_rmft60 RMFT1_prostate* in 1, noobs abb(22)
// their difference
list t_rmft60 RMFT_diff_prostate* in 1, noobs abb(22)

// To list the estimates of the life years lost due to other causes before 60 months
// under placebo
list t_rmft60 RMFT0_other* in 1, noobs abb(22)
// under DES	
list t_rmft60 RMFT1_other* in 1, noobs abb(22) 
// their difference	
list t_rmft60 RMFT_diff_other* in 1, noobs abb(22)



// we can also obtain the total months lost due both causes using option lincom(#...#) that calculates a linear combination of atn options. 
// here we have two atn() options and two models (prostate and other); the first two # corresponds to at1() i.e. under placebo and in specific the standardised months lost due to prostate cancer and the standardised months lost due to other causes
// similarly the following two # in lincom() correspond to at2() i.e. under DES and in specific the standardised months lost due to prostate cancer and the standardised months lost due to other causes

// total months lost (from both prostate and other causes) under placebo, by setting the first two # in lincom() to 1 and the following two in 0
standsurv, crmodels(prostate other) cif rmft  ///
    at1(rx 0) at2(rx 1) timevar(t_rmft60) lincom(1 1 0 0) ci ///
    atvar(RMLT0b RMLT1b) lincomvar(RMLT_total0)
	
// list	   
li RMLT_total0* in 1, noobs abb(22)


// total months lost (from both prostate and other causes) under DES, by setting the first two # in lincom() to 0 and the following two in 1
standsurv, crmodels(prostate other) cif rmft  ///
    at1(rx 0) at2(rx 1) timevar(t_rmft60) lincom(0 0 1 1) ci ///
	atvar(RMFT0c RMFTc) lincomvar(RMFT_total1)

// list		   
li RMFT_total1* in 1, noobs abb(22)
