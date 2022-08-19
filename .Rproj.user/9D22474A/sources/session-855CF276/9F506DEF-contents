## Test environments
* R version 4.2.1 (2022-06-23) -- "Funny-Looking Kid" Copyright (C) 2022 The R Foundation for Statistical Computing Platform: x86_64-apple-darwin17.0 (64-bit)
* win-builder (devel)

This is a resubmission.

## R CMD check results
R CMD check succeeded

There was one note i.e. Version contains large components (1.2.5.9019). But this note will be removed when the updates are moved to the master file and the version number is shortened to 1.2.5.

1. I addressed all the errors, warnings and notes as per the email from CRAN Prof Brian D Ripley. This was achieved by making the following change. In indirect_PAF_Sjolander_onesimulation.R and path_specific_onesimulation.R, changed Roxygen wording for parameter weighting. In particular, changed the formula for case control weighting in roxygen documentation wording from latex formula format to word format as this formula was causing errors in the new R devel platforms in Windows and Linux. This addressed the errors as per email from Prof Brian D Ripley and Cran on 11th August 2022. I also noted that we have a draft paper which is being prepared for publication, at present, which will be linked to the documentation in the future. This paper will give examples of how to deal with weighting for different study designs when prevalence is known and unknown and other considerations e.g. rare diseases.
2. I checked the updated causalPAF package on the windows development platform and it passed all these tests and checks.


## Test environments
* local R installation, R 4.1.1
* ubuntu 16.04 (on travis-ci), R 4.1.1
* win-builder (devel)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

This is a resubmission. I received feedback from CRAN and I made changes as follows:
1. I added references into the Description file.
2. I created a second dataset for the 'causalPAF' package called 'strokedata_smallSample.rda'. This consists of 5,000 rows of data whereas the other package dataset contains 16,623 rows of data. The Roxygen examples are now run on the smaller dataset, 'strokedata_smallSample', to reduce runtime. Although the examples also refer to the larger dataset should the user want to run this to obtain more accurate results. It was hoped the smaller dataset would reduce the run time to under 5 seconds. Although it has reduced the run time by about 3 minutes for the entire package examples, the run time is still over 5 seconds for most examples. As a result, I have used \dontest{} as advised by CRAN and remove \dontrun{}.
3. Some code lines in examples are commented out because they were already run in creating the package dataset. The purpose is to inform the user of the formatting of the data for the 'causalPAF' package. Including them again will increase the run time of the examples. But I can uncomment those lines if required.
4.   'causalPAF' was written in single quotes as advise by CRAN: 'causalPAF'.
5. After adding in a new package dataset called 'strokedata_smallSample' it was necessary to add Depends: R (>= 2.10) in the description file. 
