# Estimating causal effects in the presence of competing events using regression standardisation with the Stata command `standsurv`

This repository contains the code described in the manuscript titled Estimating causal effects in the presence of competing events using regression standardisation with the Stata command `standsurv` by Syriopoulou et al. (2021).

In this paper, we outline causal effects that might be of interest in the presence of competing events and show how to estimate those using regression standardisation with the Stata command `standsurv`.
Command `standsurv` applies regression standardisation and calculate the estimates as the average over all individual-specific predictions. 
Confidence intervals can also be derived using the delta method.

Causal effects can be defined as the total effect of treatment through all causal pathways between treatment and the event of interest (i.e. cumulative incidence and expected loss in life due to a cause of death) as well as the directs effect of treatment on the event of interest that does not capture the effect of treatment on the competing event (i.e. net probability of death). 
For settings where the  treatment effect can be decomposed into distinct components, separable effects have also been defined, with the separable indirect effect of treatment corresponding to the component of treatment that affects the competing event.
We demonstrate how to obtain estimates for all statistical estimands of interest and the causal effects using an example of publicly available prostate cancer data.
Data include 502 individuals that were randomly assigned estrogen therapy and are available at [`http://biostat.mc.vanderbilt.edu/DataSets`](http://biostat.mc.vanderbilt.edu/DataSets).


For the analysis, we use some user-written Stata commands.
These can be installed within Stata from the Boston College Statistical Software Components (SSC) archive as follows:

To fit the flexible parametric survival models
* `ssc install stpm2 ` 

To generate the restricted cubic spline functions
* `ssc install rcsgen`


`standsurv` command can be installed by running

`net from https://www.pclambert.net/downloads/standsurv`

Files of Stata code for obtaining several estimates of interest are available in folder [`Do`]():

* [`Kaplan Meier.do`](): produces the Kaplan-Meier for all-cause failure by treatment group (Figure 1 in the paper)
* [`Total effects.do`](): shows how to obtain total effects (i.e. cause specific cumulative incidence functions 
and expected loss in life due to a cause of death before time t*)
* [`Direct effects.do`](): shows how to obtain direct effects (i.e. net probability of death)
* [`Separable effects.do`](): shows how to obtain separable effects 
* [`Advanced modelling details.do`](): described how to obtain estimates after fitting complex models e.g. with interactions, non-linear effects as well as how to obtain non-marginalised estimates and other contrasts e.g. ratio.
