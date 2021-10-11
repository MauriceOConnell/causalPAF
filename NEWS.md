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

