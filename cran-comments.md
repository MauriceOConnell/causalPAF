## Test environments
* local R installation, R 4.1.1
* ubuntu 16.04 (on travis-ci), R 4.1.1
* win-builder (devel)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

I performed the checks on two platforms: Mac OS X and Windows.
The R package was developed and checked on Mac OS X. 
Then the causalPAF package was checked on windows using devtools::check_win_release().

The note is because this line of code below, from the package examples is 106 characters long and is truncated in the documentation rather than 100 characters long, but it needs to be written on one line as 106 characters for the package functions to work. This is because writing it in two lines in the examples, within the inverted commas, creates a space within the inverted commas, which is not recognised by the checkMarkovDag() function. I have added a comment in the example, which gives details of the code that is truncated, so a user can see what was truncated. I believe this note is ok.

Rd file 'pointEstimate.Rd':
    \examples lines wider than 100 characters:
       "ns(apob_apoa,knots=quantile(apob_apoa,c(.25,.5,.75)),Boundary.knots=quantile(apob_apoa,c(.001,.95)))",
