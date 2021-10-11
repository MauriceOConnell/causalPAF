## Test environments
* local R installation, R 4.1.1
* ubuntu 16.04 (on travis-ci), R 4.1.1
* win-builder (devel)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

This is a resubmission. I had two notes on the first submission which were coming from the
package examples in the functions pointEstimate(), causalPAFplot() and sequential_PAF().
Although I was aware of these two notes, I submitted them as I felt there were notes that
needed to be included in the submission to CRAN. 

However, the CRAN submission required I fix these two notes, so I rectified the two notes as follows:

1. The causalPAF package functions require simulations and bootstraps which take time to run. I had included
these R examples, in the first CRAN submission, with a small number of simulations and bootstraps to reduce run time,
since I thought it would be better to have running examples over 5 seconds rather than running examples not run. But this was still taking over 5 seconds. So now the R examples are present but enclosed with a \dontrun{} so it does not
create the error that the CRAN submission informed me to remove.

2. The second note requested, that in the R examples, I ensure the spline defined in the line below was 100
characters or less. Again, I was aware of this note when I submitted the causalPAF package but I had thought
it was better submitted with a comment showing the characters truncated. This line can be split into multipe lines int the R code, but when using Roxygen, in the R examples, Roxygen adds in a new line symbol when passing to the functions, " \n ", which should not be inserted. This is an issue with Roxygen.

"ns(apob_apoa,knots=quantile(apob_apoa,c(.25,.5,.75)),Boundary.knots=quantile(apob_apoa,c(.001,.95)))",

However, the CRAN submission informed me to rectify this. So I rectified it as follows.
I changed the variable name in the R examples, "apob_apoa" to "apb" which reduced each occurence of the line
to 88 characters as follows, which avoids truncation in the documents.
"ns(apb,knots=quantile(apb,c(.25,.5,.75)),Boundary.knots=quantile(apb,c(.001,.95)))",
