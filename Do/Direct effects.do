//Analysis using prostate cancer dataset available at https://hbiostat.org/data/ 
//Direct effects are estimated: 1. Net probability of death
//                             

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


// fit cause specific model only for the event of interest (prostate cancer deaths) (exit at 60 months)
// for prostate
stset dtime, failure(eventType==1) exit(time 60)
stpm2 normalAct ageCat2 ageCat3 hx hgBinary rx, scale(hazard) df(4) tvc(rx) dftvc(2)
estimates store prostate

// obtain the standardised net probability of prostate cancer death under placebo and under DES as well as their difference
// with standsurv use option failure
standsurv, failure at1(rx 0) at2(rx 1) timevar(timevar) contrast(difference) ci ///
	atvars(F_net_prostate0 F_net_prostate1) contrastvars(F_net_prostate_diff)

// list estimates at 60 months	
// standardised net probabilty of prostate cancer death under placebo and under DES with confidence intervals (CIs)
list F_net_prostate0* F_net_prostate1* if timevar==60, noobs abbrev(25)	
// their difference with CIs   
list F_net_prostate_diff* if timevar==60, noobs abbrev(25)

// plotting the standardized net probabilities of death under each treatment arm
twoway	(rarea F_net_prostate0_lci F_net_prostate0_uci timevar, color("200 82 0 %20")) ///
	(rarea F_net_prostate1_lci F_net_prostate1_uci timevar, color(black%15)) ///
	(line F_net_prostate0 F_net_prostate1  timevar,  color("200 82 0" black)) ///
		, xtitle("Time since randomisation (months)", size(medlarge)) ///
		ytitle("Net probability of death", size(medlarge)) ///
		legend(order(3 "Placebo" 4 "DES") pos(11) ring(0) cols(1) size(medium) region(lwidth(none)))   ///
		name(net, replace) ///
		ylabel(0.0(0.2)1, angle(h) format(%3.1f)) ///
		graphregion(color(white)) ///
		plotregion(margin(b=0 r=0 l=0)) 

// plotting their difference
twoway	(rarea F_net_prostate_diff_lci F_net_prostate_diff_uci timevar, lpattern(dash dash) color(black%15)) ///
	(line F_net_prostate_diff timevar , lpattern(dash) color(black)) ///
		, xtitle("Time since randomisation (months)", size(medlarge)) ///
		ytitle("Net probability of death", size(medlarge)) ///
		legend(order(2 "Difference") pos(11) ring(0) cols(1) size(medium) region(lwidth(none)))   ///
		name(net_dif, replace) ///
		ylabel(-0.4(0.2)0.4, angle(h) format(%3.1f))	///
		graphregion(color(white)) ///
		plotregion(margin(b=0 r=0 l=0)) 

// both plots together (Figure 3)
graph combine net net_dif, graphregion(color(white)) plotregion(margin(b=0 r=0 l=0)) cols(2) commonscheme 		
graph export Direct_comb.pdf, replace
