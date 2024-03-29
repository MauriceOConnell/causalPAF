% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/addInSplinesTo_in_out.R
\name{addInSplinesTo_in_out}
\alias{addInSplinesTo_in_out}
\title{Creates a DAG for input into package Dagitty to identify adjustmentSets given exposure and outcome}
\usage{
addInSplinesTo_in_out(
  in_outDAG,
  splinesVariables,
  in_outDAG_SplinesRemoved,
  Splines_outlist_Var,
  splinesDefinedIn_in_outDAG,
  count,
  Subset_adjustmentSet
)
}
\arguments{
\item{in_outDAG}{This defines the causal directed acyclic graph (DAG). A list of length 2. It is defined as a two dimensional list consisting of, firstly, the first list, inlist, i.e. a list of the parents of each variable of interest corresponding to its column name in the data. Splines can be included here if they are to be modelled as splines. Secondly, the second list, outlist, contains a list of a single name of exposure or risk factor or outcome in form of characters i.e. a list of each variable of interest (risk factors, exposures and outcome) corresponding to its column name in the data. Splines should not be input here, only the column names of the variables of interest in the data. The order at which variables are defined must satisfy (i) It is important that variables are defined in the same order in both lists e.g. the first risk factor defined in outlist has its parents listed first in inlist, the second risk factor defined in outlist has its parents listed secondly in inlist and so on. The package assumes this ordering and will not work if this order is violated. (ii) Note it is important also that the order at which the variables are defined is such that all parents of that variable are defined before it. See example in tutorial.}

\item{splinesVariables}{A list of the names of the splines listed in in_outDAG. Note these can only be splines of variables listed in in_out[[2]].}

\item{in_outDAG_SplinesRemoved}{in_outDAG as defined above but with the splines removed.}

\item{Splines_outlist_Var}{A list defined of same size and order of variables as defined in in_outArg[[2]]. If splines are to be used for variables listed in in_outArg[[2]], then the splines should be defined in Splines_outlist in the same order as variables appear in in_outArg[[2]]. It is necessary to list variables in Splines_outlist the same as in in_outArg[[2]] without splines if no spline is to be applied. It should not be input as an empty list, list(), if no splines. A warning will show if input as an empty list requiring the user to populate Splines_outlist either the same as in_outArg[[2]] (if no splines) or in the same order as in_outArg[[2]] with splines (if splines).  See example in tutorial.}

\item{splinesDefinedIn_in_outDAG}{Logical TRUE or FALSE indicating whether the user has defined splines in the causal DAG, in_out, if TRUE. If FALSE and splines are defined in Splines_outlist_Var, then it is necessary for the package to populate the in_out DAG with splines listed in Splines_outlist_Var.}

\item{count}{This variable is calculated within the function make_DAG_AdjustmentSets_in_out.R It is a count variable that should lie somewhere between 0 and length( in_outDAG_SplinesRemoved[[2]] ). If the count is calculated from make_DAG_AdjustmentSets_in_out.R to be equal to length( in_outDAG_SplinesRemoved[[2]] ) then it suggests that all in_outDAG_SplinesRemoved[[1]][[1:length( in_outDAG_SplinesRemoved[[2]] )]] are all valid adjustment sets for each of their outcomes in in_outDAG_SplinesRemoved[[2]] respectively. If count is less than length( in_outDAG_SplinesRemoved[[2]] ) then ( length( in_outDAG_SplinesRemoved[[2]] ) minus count) adjustment sets have been updated in in_outDAG_SplinesRemoved[[1]][[1:length( in_outDAG_SplinesRemoved[[2]] )]] such that they are valid adjustment sets for each of their outcomes in in_outDAG_SplinesRemoved[[2]] respectively.}

\item{Subset_adjustmentSet}{is a list of length length( in_outDAG_SplinesRemoved[[2]] ). It is calculated within the function make_DAG_AdjustmentSets_in_out.R. If all indices of Subset_adjustmentSet[[]] are empty this means that there was no updates to the adjustment sets for each of in_outDAG_SplinesRemoved[[1]][[1:length( in_outDAG_SplinesRemoved[[2]] )]] causal parents of each of in_outDAG_SplinesRemoved[[2]]  and this should coincide with count = length( in_outDAG_SplinesRemoved[[2]] ). If count < length( in_outDAG_SplinesRemoved[[2]] ), then ( length( in_outDAG_SplinesRemoved[[2]] ) minus count) adjustment sets have been updated in in_outDAG_SplinesRemoved[[1]][[1:length( in_outDAG_SplinesRemoved[[2]] )]] such that they are valid adjustment sets for each of their outcomes in in_outDAG_SplinesRemoved[[2]] respectively. And the changes to the adjustments sets are stored in each of Subset_adjustmentSet[[1:length( in_outDAG_SplinesRemoved[[2]] )]], where the non-empty index, say i in Subset_adjustmentSet[[i]],  corresponds to the variable in in_outDAG_SplinesRemoved[[2]][[i]] that has had its adjustment set in in_outDAG_SplinesRemoved[[1]][[i]] updated.}
}
\value{
Returns in_outDAG with spline variables included
}
\description{
Creates a DAG for input into package Dagitty to identify adjustmentSets given exposure and outcome
}
\keyword{internal}
