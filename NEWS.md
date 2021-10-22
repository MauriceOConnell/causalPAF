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
