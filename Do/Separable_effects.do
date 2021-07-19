//Analysis using prostate cancer dataset available at https://hbiostat.org/data/ 
//Separable effects are estimated.


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


// make copy of treatment variable, so can manipulate separately in standsurv
gen rx_c = rx
gen rx_o = rx

// fit cause specific models (exit at 60)
// for prostate
stset dtime, failure(eventType==1) exit(time 60)
stpm2 normalAct ageCat2 ageCat3 hx hgBinary rx_c, scale(hazard) df(4) tvc(rx_c) dftvc(2)
estimates store prostate

// for other causes
stset dtime, failure(eventType==2) exit(time 60)
stpm2 normalAct ageCat2 ageCat3 hx hgBinary rx_o, scale(hazard) df(3) 
estimates store other

// obtain the standardised CIFs under DES, under placebo and if the effect of other causes of death is removed as well as the differences (separable effects)
// use option cif and crmodels()
standsurv, crmodels(prostate other) cif timevar(timevar) contrast(difference) ci ///
    at1(rx_c 1 rx_o 1)  ///
    at2(rx_c 1 rx_o 0)  ///
    at3(rx_c 0 rx_o 0)  ///
	atvars(F_rx11 F_rx10 F_rx00) contrastvars(F_diff_indirect F_diff_total)
	
	
// list the standardised CIFs from death due to prostate cancer at 36 months with confidence intervals
list F_rx11_prostate* if timevar==36
list F_rx10_prostate* if timevar==36
list F_rx00_prostate* if timevar==36			
			
// list the differences (direct and indirect effect)
list F_diff_total_prostate* if timevar==36	
list F_diff_indirect_prostate* if timevar==36	

// plotting the total and indirect separable effect (figure 5 in the paper)
twoway	(line F_diff_total_prostate F_diff_indirect_prostate timevar , lpattern(solid solid  ) color(black "200 82 0")) ///
		(rarea F_diff_indirect_prostate_lci F_diff_indirect_prostate_uci timevar, color("200 82 0 %20")) ///
		(rarea F_diff_total_prostate_lci F_diff_total_prostate_uci timevar, color(black%15)) ///
		, xtitle("Time since randomisation (months)", size(medlarge)) ///
		ytitle("Cumulative incidence", size(medlarge)) ///
		legend(order(1 "Total difference" 2 "Indirect seperable effect") pos(11) ring(0) cols(1) size(medium) region(lwidth(none)))   ///
		name(seperableeffects, replace) ///
		ylabel(-0.2(0.2)1, angle(h) format(%3.1f)) ///
		graphregion(color(white)) ///
		plotregion(margin(b=0 r=0 l=0)) 


graph export Seperable_effects.pdf, replace



// cumulative incidence of death from any cause
gen F_total_rx11= F_rx11_prostate + F_rx11_other
gen F_total_rx10= F_rx10_prostate + F_rx10_other
gen F_total_rx00= F_rx00_prostate + F_rx00_other


// plotting the cumulative incidence of death from prostate cancer (solid lines) and cumulative incidence of death from any cause (dash lines), under DES, under placebo as well as under the hypothetical treatment where the effect of other causes of death is removed (in blue)
// figure 4 in the paper
twoway (line F_rx11_prostate F_rx10_prostate F_rx00_prostate timevar , lpattern(solid solid  ) color(black "95 158 209" "200 82 0")) ///
		(line F_total_rx11 F_total_rx10 F_total_rx00 timevar , lpattern(dash dash dash) color(black "95 158 209" "200 82 0")) ///
		, xtitle("Time since randomisation (months)", size(medlarge)) ///
		ytitle("Cumulative incidence", size(medlarge)) ///
		legend(order(1 "Under rx_c=1 and rx_o=1" 2 "Under rx_c=1 and rx_o=0" 3 "Under rx_c=0 and rx_o=0") pos(11) ring(0) cols(1) size(medium) region(lwidth(none)))   ///
		name(cifs, replace) ///
		ylabel(0.0(0.2)1, angle(h) format(%3.1f)) ///
		graphregion(color(white)) ///
		plotregion(margin(b=0 r=0 l=0)) 
		
graph export Seperable_effect_withtotal.pdf, replace


