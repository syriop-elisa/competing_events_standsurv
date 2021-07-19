//Analysis using prostate dataset available at https://hbiostat.org/data/ 
//Advance modelling details: 1. adding interactions
//                           2. non-linear effects
//							 3. non-marginalised estimates
//					         4. other contrasts e.g. ratio


/* 
Flexible parametric models can be fitted within Stata using the user-written command stpm2 that can be installed as follows:
  ssc install stpm2

To obtain marginal (and non-marginal) estimates using standardisation, the standsurv command must also be installed.
This  can be installed by running,
 net from https://www.pclambert.net/downloads/standsurv  
 
*/

//To generate the restricted cubic spline functions in Stata the user-written command rcsgen should be installed from SSC:
//ssc install rcsgen //uncomment to install

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


// fit cause specific model for other causes same as in previous do files
stset dtime, failure(eventType==2) exit(time 60)
stpm2 normalAct ageCat2 ageCat3 hx hgBinary rx, scale(hazard) df(3) 
estimates store other

// Below we provide some examples for obtaining cause-specific cumulative incidence after fitting more complex FPMs but other estimates of interest could also be obtained in a similar way. 
// For the remaining do file, we keep the same model for other causes as the one fitten above but allow more complex models for prostate cancer.


//************** 1. INCLUDING INTERACTIONS IN THE MODEL ************//

// generate interaction terms between treatment and age
forvalues i = 2/3 {		
	gen ageCat`i'rx=ageCat`i'*rx  				
}

// fit model for prostate cancer including interaction between rx and ageCat	 
stset dtime, failure(eventType==1) exit(time 60)
stpm2 normalAct ageCat2 ageCat3 hx hgBinary rx ageCat?rx, scale(hazard) df(4) tvc(rx) dftvc(2)
estimates store prostate

// obtain standardised CIFs by treatment group and the difference between them
// need to specify values for the interaction terms at the atn() options
standsurv, crmodels(prostate other) cif timevar(timevar) contrast(difference) ci ///
    at1(rx 0 ageCat2rx 0 ageCat3rx 0)  ///
    at2(rx 1 ageCat2rx=ageCat2 ageCat3rx=ageCat3)  ///
    atvars(CIF0b CIF1b) contrastvars(CIF_diffb)


// plot
twoway  (rarea CIF0b_prostate_lci CIF0b_prostate_uci timevar, color("200 82 0 %20")) ///
         (line CIF0b_prostate timevar, color("200 82 0")) ///
         (rarea CIF1b_prostate_lci CIF1b_prostate_uci timevar, color(black%15)) ///
         (line CIF1b_prostate timevar, color(black)) ///
                 , ylabel(0.0(0.2)1, angle(h) format(%3.1f)) ///
				 legend(order(2 "Placebo" 4 "DES") pos(10) ring(0) cols(1) size(medium) region(lwidth(none))) ///
                 xtitle("Time since randomisation (months)", size(medlarge)) ytitle("Cumulative incidence", size(medlarge)) ///
                 title("Prostate cancer", color(black)) ///
                 name(prostate_int, replace) ///
				 graphregion(color(white)) ///
				 plotregion(margin(b=0 r=0 l=0)) 
		
		
//************** 2. INCLUDING NON-LINEAR EFFECTS   ******************//

// To generate restricted cubic splines with 4 knots (3 restricted cubic spline terms) for age (relaxing linearity) we can use command rcsgen.
// for 3 degrees of freedom, 3 new age spline variables are created, agercs1-agercs3
rcsgen age, gen(agercs) df(3) orthog
// we can store the knot locations and the R Matrix, so that we can derive post-estimation predictions for specific ages later on.
// store knot positions in global macro
global ageknots `r(knots)'
// save matrix for orthogonalization
matrix Rage =r(R)	

// interactions involving the age splines can also be included in the model
// For instance, to generate interactions between age splines and treatment 

forvalues i = 1/3 {
    gen agercs`i'rx = agercs`i'*rx
	}   

// fit model for prostate cancer including interaction between treatment and age splines
stset dtime, failure(eventType==1) exit(time 60)
stpm2 normalAct agercs1 agercs2 agercs3 hx hgBinary rx agercs1rx agercs2rx agercs3rx, scale(hazard) df(4) tvc(rx) dftvc(2)
estimates store prostate
	
// obtain standardised CIFs by treatment group and the difference between them
// need to specify values for the interaction terms			
standsurv, crmodels(prostate other) cif timevar(timevar) contrast(difference) ci ///
    at1(rx 0 agercs1rx 0 agercs2rx 0 agercs3rx 0)  ///
    at2(rx 1 agercs1rx=agercs1 agercs2rx=agercs2 agercs3rx=agercs3)  ///
    atvars(CIF0c CIF1c) contrastvars(CIF_diffc)
			
// plot
twoway  (rarea CIF0c_prostate_lci CIF0c_prostate_uci timevar, color("200 82 0 %20")) ///
         (line CIF0c_prostate timevar, color("200 82 0")) ///
         (rarea CIF1c_prostate_lci CIF1c_prostate_uci timevar, color(black%15)) ///
         (line CIF1c_prostate timevar, color(black)) ///
                 , ylabel(0.0(0.2)1, angle(h) format(%3.1f)) ///
				 legend(order(2 "Placebo" 4 "DES") pos(10) ring(0) cols(1) size(medium) region(lwidth(none))) ///
                 xtitle("Time since randomisation (months)", size(medlarge)) ytitle("Cumulative incidence", size(medlarge)) ///
                 title("Prostate cancer", color(black)) ///
                 name(prostate_int_splines, replace) ///
				 graphregion(color(white)) ///
				 plotregion(margin(b=0 r=0 l=0)) 			
				 
				 
  
//************** 3. NON-MARGINALISED ESTIMATES   ***************//  
// Non-marginalised estimates can be obtained by specifying the entire covariate pattern so that the predictions are not averaged over any covariate distribution.
// For instance age-specific predictions, can be derived by calculating the spline variables at that particular age with the same knot locations and projection matrix as before
// An example is given below when interest is CIFs and we focus on individuals with normal daily activity (normalAct=1), no history of cardiovascular disease (hx=0) 
// and hemoglobin level lower than 12 (g/100ml) (hgBinary=1) and compare CIFs under DES  and under placebo for ages 60, 75 and 80 years old
// Below, the spline variables for specific ages are stored in the local macros c1, c2 and c3. 
				 
				 
foreach age in 60 75 80 {

    rcsgen, scalar(`age') knots($ageknots) rmatrix(Rage) gen(c)
    
    standsurv if _n==1, crmodels(prostate other) cif timevar(timevar) ///
        contrast(difference) ci ///
        at1(rx 0 normalAct 1 hx 0 hgBinary 1 ///
            agercs1 `=c1' agercs2 `=c2' agercs3 `=c3' ///
            agercs1rx 0 agercs2rx 0 agercs3rx 0) ///
        at2(rx 1 normalAct 1 hx 0 hgBinary 1 ///
            agercs1 `=c1' agercs2 `=c2' agercs3 `=c3' ///
            agercs1rx `=c1' agercs2rx `=c2' agercs3rx `=c3') ///
        contrastvars(CIF_diff`age') atvars(CIF0`age' CIF1`age') ///
            
}   

// As we do not average over each observation, we use if \_n == 1 to tell standsurv to only take the first observation in the stacked data to calculate non-marginalised predictions. 


// plot age-specific CIF under placebo and under DES
foreach age in  75 80  {
twoway	(rarea CIF0`age'_prostate_lci CIF0`age'_prostate_uci timevar, color("200 82 0 %20")) ///
	(rarea CIF1`age'_prostate_lci CIF1`age'_prostate_uci timevar, color(black%15)) ///
	(line CIF0`age'_prostate CIF1`age'_prostate timevar , color("200 82 0" black)) ///
		, xtitle("Months", size(medlarge)) ///
		ytitle("Cumulative incidence", size(medlarge)) ///
		ylabel(0.0(0.2)1.2, angle(h) format(%3.1f)) ///
		legend(off)   ///
		title("`age' years old", color(black)) ///
		name(CIF_`age', replace) ///
		graphregion(color(white)) ///
		plotregion(margin(b=0 r=0 l=0)) 
}
foreach age in 60 {
twoway	(rarea CIF0`age'_prostate_lci CIF0`age'_prostate_uci timevar, color("200 82 0 %20")) ///
	(rarea CIF1`age'_prostate_lci CIF1`age'_prostate_uci timevar, color(black%15)) ///
	(line CIF0`age'_prostate CIF1`age'_prostate timevar , color("200 82 0" black)) ///
		, xtitle("Months", size(medlarge)) ///
		ytitle("Cumulative incidence", size(medlarge)) ///
		legend(order(3 "Placebo" 4 "DES") pos(10) ring(0) cols(1) size(medium) region(lwidth(none)))   ///
		title("`age' years old", color(black)) ///
		name(CIF_`age', replace) ///
		ylabel(0.0(0.2)1.2, angle(h) format(%3.1f)) ///
		graphregion(color(white)) ///
		plotregion(margin(b=0 r=0 l=0))  
}



// plot their difference
foreach age in 60 {
twoway	(rarea CIF_diff`age'_prostate_lci CIF_diff`age'_prostate_uci timevar, color(black%15) lpattern(dash dash)) ///
	(line CIF_diff`age'_prostate timevar , color(black) lpattern(dash)) ///
		, xtitle("Months", size(medlarge)) ///
		ytitle("Cumulative incidence", size(medlarge)) ///
		legend(order(2 "Difference") pos(10) ring(0) cols(1) size(medium) region(lwidth(none))) ///
		title("`age' years old", color(black)) ///
		name(CIF_`age'_dif, replace) ///
		ylabel(-0.6(0.2)0.6, angle(h) format(%3.1f))	///
		graphregion(color(white)) ///
		plotregion(margin(b=0 r=0 l=0))  
}

foreach age in 75  80 {
twoway	(rarea CIF_diff`age'_prostate_lci CIF_diff`age'_prostate_uci timevar, color(black%15) lpattern(dash dash)) ///
	(line CIF_diff`age'_prostate timevar , color(black) lpattern(dash)) ///
		, xtitle("Months", size(medlarge)) ///
		ytitle("Cause-specific CIF", size(medlarge)) ///
		legend(off)   ///
		ytitle("Cumulative incidence", size(medlarge)) ///
		title("`age' years old", color(black)) ///
		name(CIF_`age'_dif, replace) ///
		ylabel(-0.6(0.2)0.6, angle(h) format(%3.1f))	///
		graphregion(color(white)) ///
		plotregion(margin(b=0 r=0 l=0))
		//ysc(off) 
}

//figure 6
graph combine CIF_60 CIF_75 CIF_80 CIF_60_dif CIF_75_dif CIF_80_dif, graphregion(color(white)) plotregion(margin(b=0 r=0 l=0)) cols(3)  commonscheme
graph export CIF_agespec.pdf, replace



//************** 4. FORMING OTHER CONTRASTS  **************//  
// Instead of the difference, the ratio can also be calculated with the option contrast(ratio). 
// For instance, the ratio in standardised CIFs under DES and under placebo can be obtained by specifying contrast(ratio) within standsurv.
standsurv, crmodels(prostate other) cif timevar(timevar) contrast(ratio) ci ///
    at1(rx 0 agercs1rx 0 agercs2rx 0 agercs3rx 0)  ///
    at2(rx 1 agercs1rx=agercs1 agercs2rx=agercs2 agercs3rx=agercs3)  ///
	atvars(CIF0d CIF1d) contrastvars(CIF_ratio)


// plot standardised CIF under placebo and under DES
twoway (rarea CIF0d_prostate_lci CIF0d_prostate_uci timevar, color("200 82 0 %20")) ///
	(rarea CIF1d_prostate_lci CIF1d_prostate_uci timevar, color(black%15)) ///
    (line CIF0d_prostate CIF1d_prostate timevar, color("200 82 0" black)) ///
		 , legend(order(3 "Placebo" 4 "DES") pos(11) ring(0) cols(1) size(medium) region(lwidth(none)))   ///
		 ylabel(0.0(0.2)1, angle(h) format(%3.1f)) ///
		 xtitle("Time since randomisation (months)", size(medlarge)) ytitle("Cumulative incidence", size(medlarge)) ///
		 name(prostate_b, replace) ///
		 graphregion(color(white)) ///
		 plotregion(margin(b=0 r=0 l=0)) 			

// plot the ratio		 
twoway (rarea CIF_ratio_prostate_lci CIF_ratio_prostate_uci timevar if timevar>1 , color(black%15) lpattern(dash dash)) ///
    (line CIF_ratio_prostate timevar if timevar>1, color(black) lpattern(dash)) ///
		 , legend(order(2 "Ratio" ) pos(11) ring(0) cols(1) size(medium) region(lwidth(none)))  ///
		 ylabel(0(2)12, angle(h)) ///
		 xtitle("Time since randomisation (months)", size(medlarge)) ytitle("Cumulative incidence", size(medlarge))  ///
		 name(prostate_ratio, replace) ///
		 graphregion(color(white)) ///
		 plotregion(margin(b=0 r=0 l=0)) 			

// figure 7
graph combine prostate_b prostate_ratio, graphregion(color(white)) plotregion(margin(b=0 r=0 l=0)) cols(3)  commonscheme  
graph export CIF_ratio.pdf, replace

