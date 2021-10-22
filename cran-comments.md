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
5. The following warning appears when running devtools::check_win_devel().

checking data for ASCII and uncompressed saves ... WARNING
  
  Note: significantly better compression could be obtained
        by using R CMD build --resave-data
                             old_size new_size compress
  strokedata.rda                269Kb    186Kb       xz
  strokedata_smallSample.rda     89Kb     67Kb       xz

But this warning was removed by adding in --resave-data in R Studio Build/Build Configure Tools and running ==> devtools::check(args = c('--as-cran'), build_args = c('--resave-data'))
