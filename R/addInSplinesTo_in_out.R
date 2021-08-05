#' @title Creates a DAG for input into package Dagitty to identify adjustmentSets given exposure and outcome
#' @description Creates a DAG for input into package Dagitty to identify adjustmentSets given exposure and outcome
#' @param in_outDAG A list of length 2. The first list contains a list of character vectors of the parents of the exposure or risk factor or outcome which are either causes or confounders of the exposure or risk factor or outcome. The second list contains a list of a single name of exposure or risk factor or outcome in form of characters. See tutorial examples for examples.
#' @param splinesVariables A list of the names of the splines listed in in_outDAG. Note these can only be splines of variables listed in in_out[[2]].
#' @param in_outDAG_SplinesRemoved in_outDAG as defined above but with the splines removed.
#' @param Splines_outlist_Var A list defined of same size and order of variables as defined in in_outArg[[2]]. If splines are to be used for variables listed in in_outArg[[2]], then the splines should be defined in the same order as variables appear in in_outArg[[2]]. It is necessary to list variables in in_outArg[[2]] without splines if no spline is to be applied. Should be input as list() if no splines.
#' @param splinesDefinedIn_in_outDAG Logical TRUE or FALSE indicating whether the user has defined splines in the causal DAG, in_out, if TRUE. If FALSE and splines are defined in Splines_outlist_Var, then it is necessary for the package to populate the in_out DAG with splines listed in Splines_outlist_Var.
#' @param count This variable is calculated within the function make_DAG_AdjustmentSets_in_out.R It is a count variable that should lie somewhere between 0 and length( in_outDAG_SplinesRemoved[[2]] ). If the count is calculated from make_DAG_AdjustmentSets_in_out.R to be equal to length( in_outDAG_SplinesRemoved[[2]] ) then it suggests that all in_outDAG_SplinesRemoved[[1]][[1:length( in_outDAG_SplinesRemoved[[2]] )]] are all valid adjustment sets for each of their outcomes in in_outDAG_SplinesRemoved[[2]] respectively. If count is less than length( in_outDAG_SplinesRemoved[[2]] ) then ( length( in_outDAG_SplinesRemoved[[2]] ) - count) adjustment sets have been updated in in_outDAG_SplinesRemoved[[1]][[1:length( in_outDAG_SplinesRemoved[[2]] )]] such that they are valid adjustment sets for each of their outcomes in in_outDAG_SplinesRemoved[[2]] respectively.
#' @param Subset_adjustmentSet is a list of length length( in_outDAG_SplinesRemoved[[2]] ). It is calculated within the function make_DAG_AdjustmentSets_in_out.R. If all indices of Subset_adjustmentSet[[]] are empty this means that there was no updates to the adjustment sets for each of in_outDAG_SplinesRemoved[[1]][[1:length( in_outDAG_SplinesRemoved[[2]] )]] causal parents of each of in_outDAG_SplinesRemoved[[2]]  and this should coincide with count = length( in_outDAG_SplinesRemoved[[2]] ). If count < length( in_outDAG_SplinesRemoved[[2]] ), then ( length( in_outDAG_SplinesRemoved[[2]] ) - count) adjustment sets have been updated in in_outDAG_SplinesRemoved[[1]][[1:length( in_outDAG_SplinesRemoved[[2]] )]] such that they are valid adjustment sets for each of their outcomes in in_outDAG_SplinesRemoved[[2]] respectively. And the changes to the adjustments sets are stored in each of Subset_adjustmentSet[[1:length( in_outDAG_SplinesRemoved[[2]] )]], where the non-empty index, say i in Subset_adjustmentSet[[i]],  corresponds to the variable in in_outDAG_SplinesRemoved[[2]][[i]] that has had its adjustment set in in_outDAG_SplinesRemoved[[1]][[i]] updated.
#' @export
#' @import MASS stats utils
#' @keywords causal DAG with splines added in
#' @return Returns in_outDAG with spline variables included

addInSplinesTo_in_out <- function(in_outDAG , splinesVariables, in_outDAG_SplinesRemoved, Splines_outlist_Var,
                                  splinesDefinedIn_in_outDAG,
                                  count,
                                  Subset_adjustmentSet){



      # splinesDefinedIn_in_outDAG
      # count
      # Subset_adjustmentSet[[i]]

#########################################################################
#########################################################################
      # 4 Cases
      #         isAdjustmentSet()      splinesDefinedIn_in_outDAG
      # Case 1.    All True                 Yes
      # Case 2.    All True                 No
      # Case 3.    Not All True             Yes
      # Case 4.    Not All True             No
#########################################################################
#########################################################################

      # Used at the end of the function, but defined at the start of the function.
      if( splinesDefinedIn_in_outDAG ){

              if( count == length( in_outDAG[[2]] ) ){
                    # Case 1.
                    in_outDAG_SplinesRemovedWithSplinesReturn <- in_outDAG
              } else if(count < length( in_outDAG[[2]] ) ){
                      # Case 3.
                      in_outDAG_SplinesRemovedWithSplinesReturn <- in_outDAG
                      Subset_adjustmentSetWithSplinesReturn <- Subset_adjustmentSet
                      for( i in 1:length( in_outDAG[[2]] ) ){

                              splinesInSubset_adjustmentSetIndices <- vector(mode = "list", length = length( in_outDAG_SplinesRemoved[[2]] ) )
                              Subset_adjustmentSetIndicesToSplines <- vector(mode = "list", length = length( in_outDAG_SplinesRemoved[[2]] ) )

                              splinesInSubset_adjustmentSetIndicesInin_outDAG2 <- vector(mode = "list", length = length( in_outDAG_SplinesRemoved[[2]] ) )

                              if( length(Subset_adjustmentSet[[i]]) > 0 ){ # This i in in_outDAG_SplinesRemoved[[1]][[i]] was updated to be an adjustment set in function make_DAG_AdjustmentSets_in_out.R

                                    # Need to add in splines to Subset_adjustmentSet[[i]]
                                    splinesInSubset_adjustmentSetIndices[[i]] <- which( splinesVariables %in% Subset_adjustmentSet[[i]] )

                                    Subset_adjustmentSetIndicesToSplines[[i]] <- which(  Subset_adjustmentSet[[i]] %in% splinesVariables  )

                                    if( length( splinesInSubset_adjustmentSetIndices[[i]] ) > 0 ){

                                            splinesInSubset_adjustmentSet <- splinesVariables[ splinesInSubset_adjustmentSetIndices[[i]] ]

                                            splinesInSubset_adjustmentSetIndicesInin_outDAG2[[i]] <- which(lapply(in_outDAG[[2]], function(data_input) all(data_input %in% splinesInSubset_adjustmentSet )  ) > 0 )

                                            ### WILL ONLY WORK IF Splines_outlist_Var DEFINED IN FORMAT list(c("splinename1","splinename2","splinename3")) since need Splines_outlist_Var[[1]][splinesInMediatorIndicesInin_outDAG2] in code
                                            Subset_adjustmentSetWithSplinesReturn[[i]][Subset_adjustmentSetIndicesToSplines[[i]] ] <- Splines_outlist_Var[[1]][splinesInSubset_adjustmentSetIndicesInin_outDAG2[[i]] ]

                                            # in_outDAG[[1]][[i]] <- list( c( in_outDAG[[1]][[i]],
                                            #                                 Subset_adjustmentSetWithSplinesReturn[[i]] ) )
                                            in_outDAG_SplinesRemovedWithSplinesReturn[[1]][[i]] <- c( in_outDAG[[1]][[i]],
                                                                                                      Subset_adjustmentSetWithSplinesReturn[[i]] )
                                      }

                              }
                      }
              } else{
                stop("Error possibly with count in functions make_DAG_AdjustmentSets_in_out.R and addinSplinesTo_in_out.R .")
              }


      } else{
              # if( count == length( in_outDAG[[2]] ) ) {

                    # Case 2. and Case 4. (since define in_outDAG_SplinesRemovedWithSplinesReturn <-  in_outDAG_SplinesRemoved wehere in_outDAG_SplinesRemoved  has had additional adjustment sets added in from function make_DAG_AdjustmentSets_in_out.R)
                    ####
                    ## NB ONLY WORKS IF in_outDAG_SplinesRemoved  has had additional adjustment sets added in from function make_DAG_AdjustmentSets_in_out.R
                    ####
                    in_outDAG_SplinesRemovedWithSplinesReturn <-  in_outDAG_SplinesRemoved # in_outDAG

                    splinesInin_outDAG_SplinesRemovedIndices <- vector(mode = "list", length = length( in_outDAG_SplinesRemoved[[2]] ) )
                    in_outDAG_SplinesRemovedIndicesToSplines <- vector(mode = "list", length = length( in_outDAG_SplinesRemoved[[2]] ) )

                    splinesInin_outDAG_SplinesRemovedIndicesInin_outDAG2 <- vector(mode = "list", length = length( in_outDAG_SplinesRemoved[[2]] ) )
                    for( i in 1:length( in_outDAG[[2]] ) ){

                          # splinesInMediatorIndices <- which( splinesVariables %in% mediator )
                          splinesInin_outDAG_SplinesRemovedIndices[[i]] <- which( splinesVariables %in% in_outDAG_SplinesRemoved[[1]][[i]] )

                          # MediatorIndicesToSplines <- which(  mediator %in% splinesVariables  )
                          in_outDAG_SplinesRemovedIndicesToSplines[[i]] <- which(  in_outDAG_SplinesRemoved[[1]][[i]] %in% splinesVariables  )
                          #### NEED TO CHECK IF WORKS FOR ANY ORDER OF mediator=c("subhtn","apob_apoa","whr") e.g. mediator=c("whr","subhtn","apob_apoa")
                          # if( length(splinesInMediatorIndices) > 0 ){
                          if( length( splinesInin_outDAG_SplinesRemovedIndices[[i]] ) > 0 ){


                                # splinesInMediator <- splinesVariables[ splinesInMediatorIndices ]
                                splinesInin_outDAG_SplinesRemoved <- splinesVariables[ splinesInin_outDAG_SplinesRemovedIndices[[i]] ]

                                # splinesInMediatorIndicesInin_outDAG2 <- which(lapply(in_outDAG[[2]], function(data_input) all(data_input %in% splinesInMediator )  ) > 0 )
                                splinesInin_outDAG_SplinesRemovedIndicesInin_outDAG2[[i]] <- which(lapply(in_outDAG[[2]], function(data_input) all(data_input %in% splinesInin_outDAG_SplinesRemoved )  ) > 0 )

                                # mediatorsWithSplinesReturn <- mediator
                                #in_outDAG_SplinesRemovedWithSplinesReturn <- in_outDAG_SplinesRemoved[[i]]

                                ### WILL ONLY WORK IF Splines_outlist_Var DEFINED IN FORMAT list(c("splinename1","splinename2","splinename3")) since need Splines_outlist_Var[[1]][splinesInMediatorIndicesInin_outDAG2] in code
                                # mediatorsWithSplinesReturn[MediatorIndicesToSplines] <- Splines_outlist_Var[[1]][splinesInMediatorIndicesInin_outDAG2]
                                in_outDAG_SplinesRemovedWithSplinesReturn[[1]][[i]][in_outDAG_SplinesRemovedIndicesToSplines[[i]] ] <- Splines_outlist_Var[[1]][splinesInin_outDAG_SplinesRemovedIndicesInin_outDAG2[[i]] ]

                          }

                    }

              # } else if(count < length( in_outDAG[[2]] ) ){
              #
              #              # Case 4.
              #              # Case 4 can be joined above with case 2.
              #
              # } else{
              #   stop("Error possibly with count in functions make_DAG_AdjustmentSets_in_out.R and addinSplinesTo_in_out.R .")
              # }

      }

       list_causal <- list("in_outDAG_WithSplines" = in_outDAG_SplinesRemovedWithSplinesReturn)

    return(list_causal)

}
