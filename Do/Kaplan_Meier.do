//Analysis using prostate dataset available at https://hbiostat.org/data/ 
//Kaplan-Meier for all-cause failure by treatment group (Figure 1 in the paper)

// first need to set working directory

// then, load data
use "prostate", clear

// restrict to placebo and high dose estrogen (DES)
keep if inlist(rx,1,4)

// update coding (0 for placebo, 1 for DES)
replace rx = cond(rx==1,0,1)
label define lblrx 0 "placebo" 1 "DES"
label values rx lblrx 

// get rid of zero timevar
// replace with half day
replace dtime = 0.5 if dtime==0

// generate all-cause death indicator 
gen allcause = status != 1
// generate event indicator
gen eventType = cond(status==1,0,cond(status==2,1,2))
label define causelab 0 "Alive" 1 "Prostate" 2 "Other"
label values eventType causelab


// declare survival data for all-cause deaths
stset dtime, failure(allcause==1) exit(time 60)

// Kaplan-Meier failure curve for DES and placebo
sts graph, by(rx) failure  ///
	xtitle("Time since randomisation (months)", size(medlarge)) ///
	ytitle("Probability of death", size(medlarge)) ///
	title("", color(black))  ///
	plot1(lcolor("200 82 0")) plot2(lcolor(black)) ///
	legend(order(1 "Placebo" 2 "DES") pos(11) ring(0) cols(1) size(medium) region(lwidth(none)))   ///
	ylabel(, angle(h) format(%3.1f)) ///
	graphregion(color(white)) ///
	plotregion(margin(b=0 r=0 l=0))
		
graph export KM_failure.pdf, replace
