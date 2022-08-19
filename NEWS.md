# causalPAF 1.2.1

* Added a `NEWS.md` file to track changes to the package.

# causalPAF 1.2.2

* Updated the date and patch number in the description file.
* Released causalPAF package to CRAN for CRAN checks and to see if CRAN are happy to release the package.


# causalPAF 1.2.3

* After submitting 1.2.1 to Cran, I received two errors which were rectified on this version as follows. Firstly, 
the R examples were taking more than 5 seconds to run in pointEstimate.R, causalPAFplot.R and sequential_PAF.R. These require this time due to simulations and bootstraps. As a result I entered \dontrun{ } for the R examples here so that they do not produce an error when submitting to Cran.
* Secondly, the R examples in pointEstimate.R and causalPAFplot.R had a spline as follows:
"ns(apob_apoa,knots=quantile(apob_apoa,c(.25,.5,.75)),Boundary.knots=quantile(apob_apoa,c(.001,.95)))",
This was more than 106 lines long, and needed to be written on 1 line in the R example (but did not have to be written in one line in the main R code and functions, since Roxygen examples added in a, new line break, "  \n "
which was a fault with Roxygen rather than the code). As a result of feedback from CRAN, I had to make sure all occurrences of this line were under 101 characters in the R examples. This was solved by changing the variable name of "apob_apoa" to "apb" in the R examples which ensured the line containing the spline was 100 characters or less as required to remove the CRAN submission error.

# causalPAF 1.2.4
* After submitting 1.2.3 to Cran, I received feedback from Cran. I made changes as follows based on that feedback.
* Added references into Description file. 
* Replaced \dontrun{} with \dontest{} for Roxygen examples, which take more than 5 seconds. 
* Created a smaller package dataset called ‘strokedata_smallSample’,  (5,000 rows of fictional patient data) and this smaller dataset was used in all Roxygen examples so that they take less time to run. The larger ‘strokedata’ (containing 16,623 rows of fictional data is still available for more accurate results). 
* Updated the package number to 1.2.4 and the date to 2021-10-22.  
* When the new package data , ‘strokedata_smallSample’, was put into the package data, it was necessary to remove an error as follows by adding in R  “--resave-data” into the  Build \ Configure Build Tools… (which is accessed via R studio and pressing Build and then pressing Configure Build Tools ) in both of these lines: 1. Install and Restart – R CMD INSTALL additional options and 2. Build Source Package – R CMD additional options.
* Updated doi urls in README.Rmd.
* Updated Description file to 1.2.4 and updated data to 1.2.4
* After adding in a new package dataset called 'strokedata_smallSample' it was necessary to add Depends: R (>= 2.10) in the description file which seemed to have to do with compressing the file.

# causalPAF 1.2.4.9000
* Created a branch called 'devel'.
* Updated the description package number for the 'devel' branch by adding .9000 to show it is the development branch version i.e. 1.2.4.9000
* Updated NEWS.md for devel 1.2.4.9000 as shown here.

# causalPAF 1.2.5.9000
* Updated description to 1.2.5.9000
* Made changes necessary for the new development version of CRAN R being released soon as per an email from CRAN. These changes are listed below are relate to changes to the documentation.
* Changed Roxygen code for the following 5 files: causalPAFplot.R (parameter model_listArg definition removed dash dollar sign) , eval_make_formula.R (parameter model_list definition removed dash dollar sign), pointEstimate.R (parameter model_listArg definition removed dash dollar sign), sim_outnode.R (parameter model_list definition removed dash dollar sign) and sim_outnode2.R (parameter model_list definition removed dash dollar sign). In each of these files the symbols dash dollar sign were removed since these were causing errors for linux and some windows documentation compilation. 
* Changed Roxygen code forNumBootstrap parameter definition in the file causalPAFplot.R to read plus or minus rather than the symbols plus or minus as it may solve issue for development version of R of linux and windows that is coming out soon that does not process documentation when there is a unicode minus character as per email received from CRAN.
* Removed minus sign in roxygen documentation in addInSplinesTo_in_out.R  in definitions of parameters 'count' and 'Subset_adjustmentSet' to read `length( in_outDAG_SplinesRemoved[[2]] ) minus count'.
* In causalPAFplot.R roxygen documentation code removed hyphen in 'maximum-likelihood' and 'non-Markov' and 'Pathway-Specific ' and 'PS-PAF' and 'pathway-specific' and 'pre-fit'.
* In checkMarkovDAG.R roxygen documentation code removed hyphen in 'non Markov'.
* In indirect_PAF_Sjolander_onesimulation.R roxygen documentation code removed hyphen in 'waist hip ratio' and 'reweighted'.
* In overall_direct.R roxygen documentation code removed hyphen in  'waist hip ratio'
* In path_specific_onesimulation.R roxygen documentation code removed hyphen in  'waist hip ratio' and 'reweighted'
* In pointestimate.R roxygen documentation code removed hyphen in  'maximum likelihood' and 'non Markov' and 'pre fit'.
* In sequential_PAF.R roxygen documentation code removed hyphen in 'Monte Carlo' and 'k minus 1' and 'non Markov'.
* In strokedata.R roxygen documentation code removed hyphen in 'preexisting' and 'Preclinical' and 'no education' and 'stage wise' and 'vice versa' and 'waist to hip ratio' and 'Leisure Physical' and '6.5 yes or no' and 'disease yes or no' and '2=1 to 8, 3=9 to 12' and '2=1 to 8, 3=9 to 12' and '2=1 to 8, 3=9 to 12' and 'Waist to hip' and 'maximum likelihood'. Changed citation to Revisiting sequential attributable fractions, J Ferguson, M O'Connell, M O'Donnell, Archives of Public Health, 2020.
* In strokedata_smallSample.R roxygen documentation code removed hyphen in 'preexisting' and 'Preclinical' and 'no education' and 'stage wise' and 'vice versa' and 'waist to hip ratio' and 'Leisure Physical' and '6.5 yes or no' and 'disease yes or no' and '2=1 to 8, 3=9 to 12' and '2=1 to 8, 3=9 to 12' and '2=1 to 8, 3=9 to 12' and 'Waist to hip' and 'maximum likelihood'. Changed citation to Revisiting sequential attributable fractions, J Ferguson, M O'Connell, M O'Donnell, Archives of Public Health, 2020.

# causalPAF 1.2.5.9019
* In indirect_PAF_Sjolander_onesimulation.R and path_specific_onesimulation.R, changed Roxygen wording for parameter weighting. In particular, changed the formula for case control weighting in roxygen documentation wording from latex formula format to word format as this formula was causing errors in the new R devel platforms in Windows and Linux. This addressed the errors as per email from Prof Brian D Ripley and Cran on 11th August 2022.
* Undone the change in the third bullet point above in 1.2.5.9000, by adding in the dollar sign as \eqn{\$} in each of the files causalPAFplot.R, eval_make_formula.R, pointEstimate.R , sim_outnode.R  and sim_outnode2.R .
* Added the published paper Pathway-specific population attributable fractions (PS-PAFs) O’Connell and Ferguson (2022) <https://doi.org/10.1093/ije/dyac079> in the description file. And updated the date and number version of the description file.


# causalPAF 1.2.5
* In Description changed version from 1.2.5.9019 to 1.2.5.

