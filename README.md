
<!-- README.md is generated from README.Rmd. Please edit that file -->

# causalPAF

<!-- badges: start -->
<!-- badges: end -->

The causalPAF package contains a suite of functions for causal analysis
calculations of population attributable fractions (PAF) given a causal
diagram which apply both:

1.  Pathway-specific population attributable fractions (PS-PAFs)
    O’Connell and Ferguson (2020)
    <https://doi.org/10.1101/2020.10.15.20212845> and
2.  Sequential population attributable fractions Ferguson, O’Connell,
    and O’Donnell (2020) <https://doi.org/10.1186/s13690-020-00442-x>.

Results are presentable in both table and plot format.

## Installation

You can install the released version of causalPAF from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("causalPAF")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("MauriceOConnell/causalPAF")
```

## Example

You start with causalPAF by supplying data and a causal directed acyclic
graph (DAG). In the R code below, the DAG is defined in the variable
‘in\_out’. You can then calculate causal population attributable
fractions (PAFs) using functions such as causalPAFplot() for pathway
specific PAFs (PS-PAF) with bootstrapped confidence intervals; or the
pointEstimate() function if only point estimates (with no confidence
intervals) of the PS-PAF are required. Sequential PAFs are calculated
using the sequential\_PAF() function.

### Stroke Data

The ‘strokedata’ data within the causalPAF package is a fictional case
control dataset \#containing key causal and modifiable risk factors for
stroke. The data is restricted to \#ischemic stroke patients and their
matched controls according to age, gender and region.

``` r
library(causalPAF)
head(strokedata)
#>   regionnn7 case esex eage htnadmbp nevfcur global_stress2 whrs2tert phys
#> 1         1    0    1   61        1       1              1         3    1
#> 2         1    0    1   60        1       1              2         2    1
#> 3         1    0    2   81        1       1              1         2    2
#> 4         1    0    2   75        1       1              1         3    2
#> 5         1    0    1   60        1       2              1         2    1
#> 6         1    0    1   73        1       1              1         1    1
#>   alcohfreqwk dmhba1c2 cardiacrfcat ahei3tert apob_apoatert subeduc moteduc
#> 1           2        1            1         2             3       3       2
#> 2           1        2            1         3             2       3       3
#> 3           2        1            2         1             3       5       3
#> 4           1        1            1         3             1       3       3
#> 5           1        1            1         3             2       5       2
#> 6           2        1            1         3             1       5       5
#>   fatduc subhtn       whr apob_apoa weights
#> 1      2      1 0.8806173 0.6603774  0.9965
#> 2      4      2 0.9051095 0.7500000  0.9965
#> 3      5      2 1.0734177 0.8540146  0.9965
#> 4      4      2 1.0635593 1.3809524  0.9965
#> 5      3      1 0.8173077 0.7132353  0.9965
#> 6      5      1 0.8232044 0.8161765  0.9965
```

Example R code, calculating pathway specific PAFs (PS-PAFs) using
causalPAFplot() and pointEstimate() and ‘strokedata’ are shown below.

``` r
library(causalPAF)

stroke_reduced <- strokedata

# The data should contain a column of weights for case control matching.
# strokedata$weights
# Weigths are not needed for cohort/cross sectional designs.


# Next, define the causal structure or directed acyclic graph (DAG) of the causal Bayesian
# network defined by the data. We list the parents of each exposure or risk factor or outcome
# in a vector as follows:

# Note it is important that the order at which the variables are defined is such that all
# parents of that variable are defined before it. 

in_phys <- c("subeduc","moteduc","fatduc")
in_ahei <- c("subeduc","moteduc","fatduc")
in_nevfcur <- c("subeduc","moteduc","fatduc")
in_alcohfreqwk <- c("subeduc","moteduc","fatduc")
in_global_stress2 <- c("subeduc","moteduc","fatduc")
in_subhtn <- c("subeduc","moteduc","fatduc","phys","ahei3tert","nevfcur","alcohfreqwk",
                "global_stress2")
 in_apob_apoa <- c("subeduc","moteduc","fatduc","phys","ahei3tert","nevfcur","alcohfreqwk",
                   "global_stress2")
 in_whr <- c("subeduc","moteduc","fatduc","phys","ahei3tert","nevfcur","alcohfreqwk",
             "global_stress2")

# Note splines can be fitted within the causal structure as shown below especially if splines
# are to be used in the fitted models.
# It is important that splines of parent variables are "typed" or "spelt" consistently
# (including spaces) throughout as causalPAF can fit models automatically provided variables are
# spelt consistently. Also if a parent variable is a spline it should be defined in spline
# format in all occurences of the parent variable.
 in_cardiacrfcat <- c("subeduc","moteduc","fatduc","phys","ahei3tert","nevfcur","alcohfreqwk",
                      "global_stress2",
"ns(apob_apoa,knots=quantile(apob_apoa,c(.25,.5,.75)),Boundary.knots=quantile(apob_apoa,c(.001,.95)))",
 "ns(whr,df=5)","subhtn")
 in_dmhba1c2 <- c("subeduc","moteduc","fatduc","phys","ahei3tert","nevfcur","alcohfreqwk",
                  "global_stress2",
"ns(apob_apoa,knots=quantile(apob_apoa,c(.25,.5,.75)),Boundary.knots=quantile(apob_apoa,c(.001,.95)))",
 "ns(whr,df=5)","subhtn")
 in_case <- c("subeduc","moteduc","fatduc","phys","ahei3tert","nevfcur","alcohfreqwk",
              "global_stress2",
"ns(apob_apoa,knots=quantile(apob_apoa,c(.25,.5,.75)),Boundary.knots=quantile(apob_apoa,c(.001,.95)))",
 "ns(whr,df=5)","subhtn","cardiacrfcat","dmhba1c2")

# Then we define a two dimensional list consisting of
# 1. inlist i.e. a list of the parents of each variable of interest corresponding to its column
# name in the data. Splines should be included here if they are to be modelled as splines.
# 2. outlist i.e. a list of each variable of interest corresponding to its column name in the
# data. Splines should not be input here, only the column names of the variables of interest in
# the data.
# Again the order is such that each variable is defined after all its parents.

in_out <- list(inlist=list(in_phys,in_ahei,in_nevfcur,in_alcohfreqwk,in_global_stress2,
                in_subhtn,in_apob_apoa,in_whr,in_cardiacrfcat,in_dmhba1c2,in_case),
                outlist=c("phys","ahei3tert","nevfcur","alcohfreqwk","global_stress2","subhtn",
                          "apob_apoa","whr","cardiacrfcat","dmhba1c2","case"))

# If splines are to be used for variables listed in in_out$outlist, then the splines should be
# defined in the same order as variables appear in in_out$outlist as follows. It is necessary to
# list variables in in_out$outlist without splines if no spline is to be applied.
# It is important that Splines_outlist is defined in the following format
# list(c("splinename1","splinename2","splinename3")) for the package to be applied correctly.
# And Splines_outlist should not be an empty list(). If there are no splines it should be
# defined the same as in_out[[2]] and in the same order as variables defined in_out[[2]].
 Splines_outlist = list( c("phys","ahei3tert","nevfcur","alcohfreqwk","global_stress2","subhtn",
"ns(apob_apoa,knots=quantile(apob_apoa,c(.25,.5,.75)),Boundary.knots=quantile(apob_apoa,c(.001,.95)))",
 "ns(whr,df=5)","cardiacrfcat","dmhba1c2","case") )

# To fit these models to case control data, one needs to perform weighted maximum-likelihood
# estimation to imitate estimation using a random sample from the population. We chose weights
# of 0.0035 (for each case) and 0.9965 (for each control), reflective of a yearly incidence of
# first ischemic stroke of 0.35%, or 3.5 strokes per 1,000 individuals. These weights were
# chosen according to average incidences across country, age, group and gender within
# INTERSTROKE according to the global burden of disease.
w <- rep(1,nrow(stroke_reduced))
w[stroke_reduced$case==0] <- 0.9965
w[stroke_reduced$case==1] <- 0.0035

# It is important to assign stroke_reduced$weights to the updated weights defined in w.
# Otherwise if stroke_reduced$weights <- w is not set, the alternative weights supplied in the
#  fictional data will be used. In this case, we want to use weigths as defined in w.
stroke_reduced$weights <- w

#The checkMarkovDAG() function in the causalPAF package should be used before running
# causalPAFplot() to ensure:
#1. The causal Markov condition holds for the causal structure defined in the variable in_out.
#2. The variables in in_out are listed in the order so that no variable is defined before a
# parent or direct cause. Note: if this order does not hold, checkMarkovDAG() will automatically
# reorder the variables in, in_out, provided it is a Markov DAG.

#The causal analysis requires that the causal structure is a Markov DAG. The Causal Markov (CM)
# condition states that, conditional on the set of all its direct causes, a node is independent
# of all variables which are not direct causes or direct effects of that node. In the event that
# the structure of a Bayesian network accurately depicts causality, the two conditions are
# equivalent. However, a network may accurately embody the Markov condition without depicting
# causality, in which case it should not be assumed to embody the causal Markov condition.

# in_out is as defined above and input into this code.
 if(checkMarkovDAG(in_out)$IsMarkovDAG & !checkMarkovDAG(in_out)$Reordered){
   print("Your in_out DAG is a Markov DAG.")
   } else if( checkMarkovDAG(in_out)$IsMarkovDAG & checkMarkovDAG(in_out)$Reordered ) {

       in_out <- checkMarkovDAG(in_out)[[2]]

           print("Your in_out DAG is a Markov DAG.The checkMarkovDAG function has reordered your
                in_out list so that all parent variables come before descendants.")
           } else{ print("Your ``in_out'' list is not a Bayesian Markov DAG so the methods in the
                         causalPAF package cannot be applied for non-Markov DAGs.")}
# The pointEstimate() function evaluates Point Estimates for Total PAF, Direct PAF, Indirect PAF
# and Path Specific PAF for a user inputted number of integral simulations. There is no bootstap
# applied in this fucntion.
# Since bootstraps are not applied, the pointEstimate() function will run quicker than the
# alternative causalPAFplot() function which calculates bootstrap estimates which can take
# longer to run.

          pointEstimate(dataframe = stroke_reduced,
                        exposure="phys",
                        mediator=c("subhtn","apob_apoa","whr"),
                        response="case",
                        response_model_mediators = list(),
                        response_model_exposure = list(),
                        in_outArg = in_out,
                        Splines_outlist = Splines_outlist,
                        splinesDefinedIn_in_outDAG = TRUE,
                        model_listArg = list(),
                        weights = w,
                        NumSimulation = 3,
                        addCustom = TRUE,
                        custom = "regionnn7*ns(eage,df=5)+esex*ns(eage,df=5)")


# The causalPAFplot() function will perform Pathway-Specific Population Attributable Fraction
# (PS-PAF) calculations and output results based on an exposure, mediators and response input
# by the user according to the columns names of these variables defined in the dataframe.

# Setting model_listArg, response_model_mediators and response_model_exposure by default to an
# empty list will instruct the causalPAF package to fit these models automatically based on the
# causal DAG supplied in the in _outArg. Alternatively the user can supply their custom fitted,
# model_listpop, response_model_mediators and response_model_exposure which should be consistent
# with the causal structure.

# Note we fit a custom interaction for the outcome (or case or response) regression
# ( custom = "regionnn7*ns(eage,df=5)+esex*ns(eage,df=5)") ). Care should be taken that the
# customised regression should not contain variables that might affect the causal interpretation
# of the regression e.g. in this case we have used baseline confounders (i.e. regionn, eage and
# esex) with interactions and splines. In general, using baseline confounders in custom should
# not affect any causal interpretations whereas using variables far ``downstream'' might block
# causal pathways. The user is required to apply discretion in using ``addCustom'' and
# ``Custom'' in ensuring a causal interpretation remains. If no customisation is required the
# user can input addCustom = FALSE and custom = "" which is the default setting.

# Finally we call the causalPAFplot function for the pathway-specific PAF calculations as
# follows:
          causalPAFplot(dataframe = stroke_reduced,
                        exposure="phys",
                        mediator=c("subhtn","apob_apoa","whr"),
                        response="case",
                        response_model_mediators = list(),
                        response_model_exposure = list(),
                        in_outArg = in_out,
                        Splines_outlist = Splines_outlist,
                        splinesDefinedIn_in_outDAG = TRUE,
                        model_listArg = list(),
                        weights = w,
                        NumBootstrap = 2,
                        NumSimulation = 2,
                        plot = "bar",
                        fill= "skyblue",
                        colour="orange",
                        addCustom = TRUE,
                        custom = "regionnn7*ns(eage,df=5)+esex*ns(eage,df=5)")
```

Sequential PAFs can be calculated using the sequential\_PAF() function.
Example R code, using ‘strokedata’ is shown below.

``` r
stroke_reduced <- strokedata

 in_phys <- c("subeduc","moteduc","fatduc")
 in_ahei <- c("subeduc","moteduc","fatduc")
in_nevfcur <- c("subeduc","moteduc","fatduc")
in_alcohfreqwk <- c("subeduc","moteduc","fatduc")
in_global_stress2 <- c("subeduc","moteduc","fatduc")
in_htnadmbp <- c("subeduc","moteduc","fatduc","phys","ahei3tert","nevfcur","alcohfreqwk",
                  "global_stress2")
in_apob_apoatert <- c("subeduc","moteduc","fatduc","phys","ahei3tert","nevfcur","alcohfreqwk",
                       "global_stress2")
in_whrs2tert <- c("subeduc","moteduc","fatduc","phys","ahei3tert","nevfcur","alcohfreqwk",
                   "global_stress2")
in_cardiacrfcat <- c("subeduc","moteduc","fatduc","phys","ahei3tert","nevfcur","alcohfreqwk",
                      "global_stress2", "apob_apoatert","whrs2tert","htnadmbp")
in_dmhba1c2 <- c("subeduc","moteduc","fatduc","phys","ahei3tert","nevfcur","alcohfreqwk",
                   "global_stress2", "apob_apoatert","whrs2tert","htnadmbp")
in_case <- c("subeduc","moteduc","fatduc","phys","ahei3tert","nevfcur","alcohfreqwk",
"global_stress2", "apob_apoatert","whrs2tert","htnadmbp","cardiacrfcat","dmhba1c2")

in_out <- list(inlist=list(in_phys,in_ahei,in_nevfcur,in_alcohfreqwk,in_global_stress2,
                in_htnadmbp, in_apob_apoatert,in_whrs2tert,in_cardiacrfcat,
                in_dmhba1c2,in_case),
                outlist=c("phys","ahei3tert","nevfcur","alcohfreqwk","global_stress2",
                          "htnadmbp","apob_apoatert", "whrs2tert","cardiacrfcat",
                          "dmhba1c2","case"))


 if(checkMarkovDAG(in_out)$IsMarkovDAG & !checkMarkovDAG(in_out)$Reordered){
   print("Your in_out DAG is a Markov DAG.")
 } else if( checkMarkovDAG(in_out)$IsMarkovDAG & checkMarkovDAG(in_out)$Reordered ) {

   in_out <- checkMarkovDAG(in_out)[[2]]

   print("Your in_out DAG is a Markov DAG.The checkMarkovDAG function has reordered your
           in_out list so that all parent variables come before descendants.")
 } else{ print("Your ``in_out'' list is not a Bayesian Markov DAG so the methods in the
                causalPAF package cannot be applied for non-Markov DAGs.")}

 w <- rep(1,nrow(stroke_reduced))
 w[stroke_reduced$case==0] <- 0.9965
 w[stroke_reduced$case==1] <- 0.0035

 stroke_reduced$weights <- w

 sequentialPAF <- sequential_PAF( dataframe = stroke_reduced,
                                  model_list_var = list(),
                                  weights = w,
                                  in_outDAG = in_out,
                                  response = "case",
                                  NumOrderRiskFactors = 3,
                                  addCustom = TRUE,
                                  custom = "regionnn7*ns(eage,df=5)+esex*ns(eage,df=5)" )

 sequentialPAF$SAF_summary
          
```
