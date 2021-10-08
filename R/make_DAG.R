#' @title Creates a DAG for input into package Dagitty to identify adjustmentSets given exposure and outcome
#' @description Creates a DAG for input into package Dagitty to identify adjustmentSets given exposure and outcome
#' @param in_outDAG This defines the causal directed acyclic graph (DAG). A list of length 2. It is defined as a two dimensional list consisting of, firstly, the first list, inlist, i.e. a list of the parents of each variable of interest corresponding to its column name in the data. Splines can be included here if they are to be modelled as splines. Secondly, the second list, outlist, contains a list of a single name of exposure or risk factor or outcome in form of characters i.e. a list of each variable of interest (risk factors, exposures and outcome) corresponding to its column name in the data. Splines should not be input here, only the column names of the variables of interest in the data. The order at which variables are defined must satisfy (i) It is important that variables are defined in the same order in both lists e.g. the first risk factor defined in outlist has its parents listed first in inlist, the second risk factor defined in outlist has its parents listed secondly in inlist and so on. The package assumes this ordering and will not work if this order is violated. (ii) Note it is important also that the order at which the variables are defined is such that all parents of that variable are defined before it. See example in tutorial.
#' @param exposure The name of the exposure column variable within dataframe in text format e.g. "phys".
#' @param response The name of the response column variable within dataframe in text format e.g. "case". The cases should be coded as 1 and the controls as 0.
#' @param mediator The name of the mediator column variables within dataframe in text format. There can be more than one mediator of interest. It can be a vector of mediators names within the dataframe e.g. c("subhtn","apob_apoa","whr").
#' @param Splines_outlist_Var A list defined of same size and order of variables as defined in in_outArg[[2]]. If splines are to be used for variables listed in in_outArg[[2]], then the splines should be defined in Splines_outlist in the same order as variables appear in in_outArg[[2]]. It is necessary to list variables in Splines_outlist the same as in in_outArg[[2]] without splines if no spline is to be applied. It should not be input as an empty list, list(), if no splines. A warning will show if input as an empty list requiring the user to populate Splines_outlist either the same as in_outArg[[2]] (if no splines) or in the same order as in_outArg[[2]] with splines (if splines).  See example in tutorial.
#' @param splinesDefinedIn_in_outDAG Logical TRUE or FALSE indicating whether the user has defined splines in the causal DAG, in_out, if TRUE. If FALSE and splines are defined in Splines_outlist_Var, then it is necessary for the package to populate the in_out DAG with splines listed in Splines_outlist_Var.
#' @param addCustomExposureAdjustmentSet Logical TRUE or FALSE indicating whether a customised interaction term is to be added to the each regression for the Exposure adjustment set. The interaction term can include splines. NOTE variables in addCustom are input by the User and the User should decide whether they are to be included in the adjustmentSet and whether these variables are not already present. The DAG will be fit without nodes in addCustom since causal DAG nodes will not recognise customised nodes e.g. with a customised formula in the node name. The addCustom is added in additively as adjustment set variables after the causal diagram has been defined and adjustment sets have been identified. Care should be taken when using addCustom not to add in variables eg. twice if already included as variables in the causal dag.
#' @param customExposureAdjustmentSet text containing the customised interaction term to be added to each regression for the Exposure adjustment set. The text should be enclosed in inverted commas. Splines can be included within the interaction terms. If there are more than one exposure and the custom is to differ per exposure, custom can be defined in a vector of form c(), where each customised adjustment is in the same order as the variables defined in exposure if more than one exposure.  See tutorial for examples. NOTE variables in addCustom are input by the User and the User should decide whether they are to be included in the adjustmentSet and whether these variables are not already present. The DAG will be fit without nodes in addCustom since causal DAG nodes will not recognise customised nodes e.g. with a customised formula in the node name. The addCustom is added in additively as adjustment set variables after the causal diagram has been defined and adjustment sets have been identified. Care should be taken when using addCustom not to add in variables eg. twice if already included as variables in the causal dag.
#' @param addCustomMediatorAdjustmentSet Logical TRUE or FALSE indicating whether a customised interaction term is to be added to the each regression for the Mediator adjustment set. The interaction term can include splines. NOTE variables in addCustom are input by the User and the User should decide whether they are to be included in the adjustmentSet and whether these variables are not already present. The DAG will be fit without nodes in addCustom since causal DAG nodes will not recognise customised nodes e.g. with a customised formula in the node name. The addCustom is added in additively as adjustment set variables after the causal diagram has been defined and adjustment sets have been identified. Care should be taken when using addCustom not to add in variables eg. twice if already included as variables in the causal dag.
#' @param customMediatorAdjustmentSet text containing the customised interaction term to be added to each regression for the Mediator adjustment set. The text should be enclosed in inverted commas. Splines can be included within the interaction terms.  If there are more than one mediator and the custom is to differ per mediator, custom can be defined in a vector of form c(), where each customised adjustment is in the same order as the variables defined in mediator if more than one mediator. See tutorial for examples. NOTE variables in addCustom are input by the User and the User should decide whether they are to be included in the adjustmentSet and whether these variables are not already present. The DAG will be fit without nodes in addCustom since causal DAG nodes will not recognise customised nodes e.g. with a customised formula in the node name. The addCustom is added in additively as adjustment set variables after the causal diagram has been defined and adjustment sets have been identified. Care should be taken when using addCustom not to add in variables eg. twice if already included as variables in the causal dag.
#' @export
#' @importFrom ggdag dagify tidy_dagitty ggdag theme_dag
#' @import splines MASS stats utils dagitty ggplot2
#' @keywords internal
#' @return Returns DAG and adjustment set given exposure and outcome.

make_DAG <- function(in_outDAG , exposure, response, mediator, Splines_outlist_Var = list(), splinesDefinedIn_in_outDAG = list(), addCustomExposureAdjustmentSet = FALSE, customExposureAdjustmentSet, addCustomMediatorAdjustmentSet = FALSE, customMediatorAdjustmentSet ){


  # ggproto <- ggplot2::ggproto

######################
######################
### NEED TO RUN dagitty on in_out[[2]] and see if all have adjustment sets
######################
######################

# CREATE A FUNCTION AND CALL IT WITHIN THIS FUNCTION USING dagified

# isAdjustmentSet(dagified, in_outDAG_SplinesRemoved[[1]][[1]], exposure = in_outDAG_SplinesRemoved[[1]][[1]], outcome = in_outDAG_SplinesRemoved[[2]][[1]])

# adjustmentSets(dagified, exposure = in_outDAG_SplinesRemoved[[1]][[1]], outcome = in_outDAG_SplinesRemoved[[2]][[1]] , type = "canonical", effect = "total")
#
# adjustmentSets(dagified, exposure = in_outDAG_SplinesRemoved[[1]][[1]][1], outcome = in_outDAG_SplinesRemoved[[2]][[1]] , type = "canonical", effect = "total")
#
# adjustmentSets(dagified, exposure = in_outDAG_SplinesRemoved[[1]][[1]][2], outcome = in_outDAG_SplinesRemoved[[2]][[1]] , type = "canonical", effect = "total")
#
# adjustmentSets(dagified, exposure = in_outDAG_SplinesRemoved[[1]][[1]][3], outcome = in_outDAG_SplinesRemoved[[2]][[1]] , type = "canonical", effect = "total")

#######################
#######################


if( (length(splinesDefinedIn_in_outDAG) == 0) || (length(Splines_outlist_Var) == 0) ){

  stop("Please ensure the logical variable splinesDefinedIn_in_outDAG and list variable Splines_outlist_Var are defined. splinesDefinedIn_in_outDAG is a logical TRUE or FALSE indicating whether the user has defined splines in the causal DAG, in_out, if TRUE. If splinesDefinedIn_in_outDAG is set to FALSE and splines are defined in Splines_outlist_Var, then this informs the package to populate the in_out DAG with splines listed in Splines_outlist_Var.
       Splines_outlist_Var is a list defined of same size and order of variables as defined in in_outArg[[2]]. If splines are to be used for variables listed in in_outArg[[2]], then the splines should be defined in Splines_outlist_Var in the same order as variables appear in in_outArg[[2]]. It is necessary to list variables in Splines_outlist_Var as in in_outArg[[2]] without splines if no spline is to be applied. If Splines_outlist_Var is input as a list(),
       this error message will appear prompting the user to populate Splines_outlist_Var similar to in_outArg[[2]], with or without splines as appropriate.")

}

## NB NEED TO ADD IN IF STATEMENT CODE AS MIGHT NOT BE ANY SPLINES

# Splines_outlist_Var MUST BE DEFINED AS LIST OF SIZE 1 for Splines_outlist_Var[[1]] to work below
# splinesVariables_splines <- <- Splines_outlist_Var[[1]][ splinesVariablesIndices ]

# If models below are listed as empty list(), then expect them to be populated otherwise do not populate. Identify this in causalPAFplot.r
# response_model_mediators = response_vs_mediator,
# response_model_exposure = list(),


  # make_DAG( in_outDAG = in_out,
  #            exposure = "phys" ,
  #            response = "case",
  #            mediator = c("subhtn","apob_apoa","whr") ,
  #            Splines_outlist_Var = Splines_outlist,
  # addCustomExposureAdjustmentSet = TRUE,
  # customExposureAdjustmentSet =  "regionnn7*ns(eage,df=5)+esex*ns(eage,df=5)",
  # addCustomMediatorAdjustmentSet = TRUE,
  # customMediatorAdjustmentSet = c("regionnn7*ns(eage,df=5)+esex*ns(eage,df=5)",
  #                                 "regionnn7*ns(eage,df=5)+esex*ns(eage,df=5)",
  #                                 "regionnn7*ns(eage,df=5)+esex*ns(eage,df=5)"),
  # splinesDefinedIn_in_outDAG = TRUE)



  #            addCustom = TRUE,
  #            custom = "regionnn7*ns(eage,df=5)+esex*ns(eage,df=5) + ")$exposuresAdjustmentSetCanonical

  ######
  ######
  # DAG will not include custom in DAG for adjustmentSet (as dagitty cannot recognise e.g. splines in custom etc..) as assumes user has to input custom themselves. So user needs to be aware whether custom is
  # in adjustment set and whether that variable has been included or excluded already
  #######
  #######

  # Decided to add this in the end after DAG as DAG will not reconize it likley e.g. if slines etc...
  # This can be done in other function to add in custom.Maybe with make_formula function
  # addCustom = FALSE, custom = "~ regionnn7*ns(eage,df=5)+esex*ns(eage,df=5) + "

  ####
  # MIGTH NEED TO CHANGE FORMAT OF Splines_outlist_Var AND in_outDAG[[2]] IN CASE THE USER DEFINES IT SOMEOTHER WAY IE. LIST, DATAFRAME ETC
  ####

  # in_outDAG = in_out
  # exposure="phys"
  # response="case"
  # mediator=c("subhtn","apob_apoa","whr")
  # Splines_outlist_Var = list( c("phys","ahei3tert","nevfcur","alcohfreqwk","global_stress2","subhtn","ns(apob_apoa, knots = quantile(apob_apoa,c(.25,0.5,0.75)), Boundary.knots = quantile(apob_apoa,c(.001,0.95)))","ns(whr,df=5)","cardiacrfcat","dmhba1c2","case") )
  # addCustomExposureAdjustmentSet = TRUE
  # customExposureAdjustmentSet =  "regionnn7*ns(eage,df=5)+esex*ns(eage,df=5)"
  # addCustomMediatorAdjustmentSet = TRUE
  # customMediatorAdjustmentSet = c("regionnn7*ns(eage,df=5)+esex*ns(eage,df=5)",
  #                                 "regionnn7*ns(eage,df=5)+esex*ns(eage,df=5)",
  #                                 "regionnn7*ns(eage,df=5)+esex*ns(eage,df=5)")

# addCustom = FALSE, custom = "~ regionnn7*ns(eage,df=5)+esex*ns(eage,df=5) + "

################################################################################################################################################
################################################################################################################################################
### To create Causal DAG, Need to remove splines from in_outDAG since since dagitty will not recognised spline nodes in order to create adjustment set
################################################################################################################################################
################################################################################################################################################

# Will only work if in same order
# finds exact match
splinesVariablesIndices <- which( as.data.frame( in_outDAG[[2]] ) != as.data.frame( Splines_outlist_Var ) )



######################
######################
######################
# # if there are splines
## NB NEED TO ADD IN IF STATEMENT CODE AS MIGHT NOT BE ANY SPLINES
########################
######################
######################
# if there are splines and these splines are defined in in_outDAG as per the logical splinesDefinedIn_in_outDAG, then we need to remove those splines for daggity to operate first.
if( (length( splinesVariablesIndices) > 0) & (splinesDefinedIn_in_outDAG == TRUE) ){

          #splinesVariables <- in_outDAG[[2]][ which( as.data.frame( in_outDAG[[2]] ) != as.data.frame( Splines_outlist_Var ) ) ]
          splinesVariables <- in_outDAG[[2]][ splinesVariablesIndices ]

          # MUST BE DEFINED AS LIST OF SIZE 1
          splinesVariables_splines <- Splines_outlist_Var[[1]][ splinesVariablesIndices ]

          # #############
          # #############
          # #############
          # # Rows
          # which(lapply(in_outDAG[[1]] , function(data_input) any(data_input %in% splinesVariables_splines[1] )  ) > 0 )
          # which(lapply(in_outDAG[[1]] , function(data_input) any(data_input %in% splinesVariables_splines[2] )  ) > 0 )
          #
          # # Columns
          # which(lapply(in_outDAG[[1]][[11]] , function(data_input) data_input %in% splinesVariables_splines[1] )  > 0 )
          # which(lapply(in_outDAG[[1]][[11]] , function(data_input) data_input %in% splinesVariables_splines[2] )  > 0 )
          #
          # ############
          # ############
          # ############


          # ARE THERE VARIABLES IN as.data.frame( in_outDAG[[2]] ) THAT CONTAIN all of words in splinesVariables or part of words in splinesVariables (which could resutls in mixing up variables)

          # v[order(nchar(v), v, decreasing = TRUE)]
          # Order with longest variables first and then in decreasing alphabetical order, this avoids mixin up eg. apoa and apoah since apoah is taken out first.
          # This approach is not needed now if searching for exact match of spline instead
          # splinesVariables <- splinesVariables[order(nchar(splinesVariables), splinesVariables, decreasing = TRUE)]
          # Next need to remove index if used for a variabe so e.g. apoa is not mixed with with apoah rather than apoa

          #############
          ## Rows
          #############
          splinesIndices <- list()
          # splinesIndicesRows <- list()
          #splinesIndices_splines <- list()
          # Rows
          for( i in 1:length( splinesVariables ) ){

            # THIS WILL NOT WORK IF VARIABLES WITH WORDS SHARED E.G. apoa and apoah and apoa33
            # finds any word containing it
            #splinesIndices[[i]] <- grep(paste( splinesVariables[i] ,sep=''), in_outDAG[[1]] ,perl=TRUE)

            # SPLINES NEED TO BE WRITTEN IN EXACT SAME WAY THROUGHOUT I.E. SAME SPACES ETC OTHERWISE WILL NOT WORK
            #splinesIndicesRows[[i]] <- which(lapply(in_outDAG[[1]] , function(data_input) any(data_input %in% splinesVariables_splines[i] )  ) > 0 )
            splinesIndices[[i]] <- which(lapply(in_outDAG[[1]] , function(data_input) any(data_input %in% splinesVariables_splines[i] )  ) > 0 )

          }


          # grep(paste( splinesVariables[1] ,sep=''), in_outDAG[[1]][splinesIndices[[1]] ][[1]] ,perl=TRUE)
          # grep(paste( splinesVariables[1] ,sep=''), in_outDAG[[1]][splinesIndices[[1]] ][[2]] ,perl=TRUE)
          # grep(paste( splinesVariables[1] ,sep=''), in_outDAG[[1]][splinesIndices[[1]] ][[3]] ,perl=TRUE)
          #
          # grep(paste( splinesVariables[2] ,sep=''), in_outDAG[[1]][splinesIndices[[2]] ][[1]] ,perl=TRUE)
          # grep(paste( splinesVariables[2] ,sep=''), in_outDAG[[1]][splinesIndices[[2]] ][[2]] ,perl=TRUE)
          # grep(paste( splinesVariables[2] ,sep=''), in_outDAG[[1]][splinesIndices[[2]] ][[3]] ,perl=TRUE)
          # Will not work if variable text has multiple variables with same text in splines
          ########
          # COLUMNS
          #########
          storeSplines <- vector(mode = "list", length = length(splinesVariables) )
          for(i in 1:length(splinesVariables ) ){
                storeSplinesIndices <- vector(mode = "list", length = length( splinesIndices[[i]] ) )
                for(j in 1:length( splinesIndices[[i]] ) ){
                      # THIS WILL NOT WORK IF VARIABLES WITH WORDS SHARED E.G. apoa and apoah and apoa33
                      # finds any word containing it
                      #storeSplinesIndices[[i]][[j]] <- grep(paste( splinesVariables[i] ,sep=''), in_outDAG[[1]][splinesIndices[[i]] ][[j]] ,perl=TRUE)

                      # which(lapply(in_outDAG[[1]][[11]] , function(data_input) data_input %in% splinesVariables_splines[1] )  > 0 )
                      #storeSplinesIndicesTest[[i]][[j]] <- which(lapply(in_outDAG[[1]][splinesIndices[[i]] ][[j]] , function(data_input) data_input %in% splinesVariables_splines[i] )  > 0 )
                      storeSplinesIndices[[i]][[j]] <- which(lapply(in_outDAG[[1]][splinesIndices[[i]] ][[j]] , function(data_input) data_input %in% splinesVariables_splines[i] )  > 0 )
                }
            storeSplines[[i]] <- storeSplinesIndices[[i]]
            storeSplinesIndices <- list()

          }

          # in_outDAG[[1]][[ splinesIndices[[1]][1]  ]][[ storeSplines[[1]][[1]] ]]

          in_outDAG_SplinesRemoved <- in_outDAG

          for( i in 1:length(splinesVariables) ){
                for(j in 1:length( splinesIndices[[i]] ) ){
                                                        ## Rows                    ## Columns
                      in_outDAG_SplinesRemoved[[1]][[ splinesIndices[[i]][j]  ]][[ storeSplines[[i]][[j]] ]] <- splinesVariables[[i]]

                }
          }


} else{

          # There are no splines. To reduce need for extra code just set  in_outDAG_SplinesRemoved <- in_outDAG
          # could just use in_outDAG but easier as follows as reduces need for extra code.
          in_outDAG_SplinesRemoved <- in_outDAG
}

################################################################
################################################################



# assumes in_outArg[[1]] and in_outArg[[2]] have all variables and all variables are defined in same order and of same length
for( i in 1:length( in_outDAG_SplinesRemoved[[2]] ) ){
    # in_outDAG[[1]]  abd in_outDAG[[2]] must be of same lenght and order




# in_outArg[[1]]
#
# in_outArg[[2]]
#
# if(addCustom){
#         # result <- paste(outvar,"~ regionnn7*ns(eage,df=5)+esex*ns(eage,df=5) + ",in_vars[1])
#         # result <- paste(outvar, custom, in_vars[1])
#         result <- paste(outvar, custom, in_vars[1])
#   }else{
#         # result <- paste(outvar,"~ regionnn7*ns(eage,df=5)+esex*ns(eage,df=5) + ",in_vars[1])
#         result <- paste(outvar,"~ ",in_vars[1])
#   }
#         if(length(in_vars)>=2){
#
#                 for(i in 2:length(in_vars)){
#
#                         result <- paste(result,"+ ",in_vars[i],sep='')
#
#                 }
#         }
#
#         result

  if( i == 1 ){

                #################
                #################
                ## NB NOTE NB
                ## addCustom is removed here since Dagitty may bot recognise e.g. splines
                ## will add in addCustom at the end, but User will have to know whether the values should be in the adjustmentSet and whether they were already included and hence double counted if added in again in custom?
                #################
                #################
                #   ## ONLY NEED dagify() at start on i=1
                #   #  NB if variables included in custom then they are not to be included in in_out as would be included in custom
                #     if(addCustom){
                #       # result <- paste(outvar,"~ regionnn7*ns(eage,df=5)+esex*ns(eage,df=5) + ",in_vars[1])
                #       # result <- paste(outvar, custom, in_vars[1])
                #       DAG <- paste("dagify(", in_outDAG_SplinesRemoved[[2]][[1]], custom, in_outDAG_SplinesRemoved[[1]][[1]][1])
                # }else{
                # result <- paste(outvar,"~ regionnn7*ns(eage,df=5)+esex*ns(eage,df=5) + ",in_vars[1])

                        DAG <- paste("dagify(", in_outDAG_SplinesRemoved[[2]][1],"~ ",in_outDAG_SplinesRemoved[[1]][[1]][1] )
                # }


                      if(length( in_outDAG_SplinesRemoved[[1]][[1]] )>=2){

                              for(m in 2:length( in_outDAG_SplinesRemoved[[1]][[1]] )){

                                      DAG <- paste(DAG,"+ ",in_outDAG_SplinesRemoved[[1]][[1]][m],sep='')

                              }
                      }

  }else{
          ## for i>1 ONLY NEED dagify() at start on i=1
          #  NB if variables included in custom then they are not to be included in in_out as would be included in custom
            # if(addCustom){
            #   # result <- paste(outvar,"~ regionnn7*ns(eage,df=5)+esex*ns(eage,df=5) + ",in_vars[1])
            #   # result <- paste(outvar, custom, in_vars[1])
            #   DAG <- paste("dagify(", in_outDAG[[2]][i], custom, in_outDAG[[1]][i])
            # }else{
                  # result <- paste(outvar,"~ regionnn7*ns(eage,df=5)+esex*ns(eage,df=5) + ",in_vars[1])

                  DAG <- paste(",", in_outDAG_SplinesRemoved[[2]][i],"~ ",in_outDAG_SplinesRemoved[[1]][[i]][1] )
            # }
                  if(length( in_outDAG_SplinesRemoved[[1]][[i]] )>=2){

                        for(m in 2:length( in_outDAG_SplinesRemoved[[1]][[i]] )){

                                DAG <- paste(DAG,"+ ",in_outDAG_SplinesRemoved[[1]][[i]][m],sep='')

                        }
                  }

   }



        # add at end
        # , exposure = "x", outcome = "y")

        if( i == length( in_outDAG_SplinesRemoved[[2]] ) ){
           # only need this on last i
            DAG <- paste(DAG,", exposure =", deparse(substitute(exposure)) , ", outcome =", deparse(substitute(response)) , ")",sep='')
        }

   if(i == 1){
          to_execute <-  DAG
   }else{
          to_execute <-  paste(to_execute, DAG ,sep='')
   }

}





    # library(dagitty)
    # library(ggdag)
    # library(ggplot2)
    dagified <- eval(parse(text = to_execute ) )

    tidy_ggdag <- tidy_dagitty(dagified)

    # https://ggdag.malco.io
    # quartz()
    #### TRY AND ADD IN CODE TO HAVE PARENTS IN LEFT AND TIDY UP GRAPH
    #CausalDagPlot <- ggdag::ggdag(tidy_ggdag) + theme_dag()
    #CausalDagPlot <- ggdag(tidy_ggdag) + theme_dag()

  #  quartz()
  #   ggdag_adjustment_set(tidy_ggdag,exposure = exposure, outcome = response, node_size = 14) +
  # theme(legend.position = "bottom")

    # adjustmentSets(x, exposure = NULL, outcome = NULL, type = canonical", effect = c("total", "direct"))
    # Might be more than 1 exposure required
    exposureAdjustmentSetCanonical <- vector(mode = "list", length = length(exposure) )
    for( i in 1:length(exposure) ){
          exposureAdjustmentSetCanonical[i] <- adjustmentSets(dagified, exposure = exposure, outcome = response , type = "canonical", effect = c("total", "direct"))

    }



    mediatorAdjustmentSetCanonical <- vector(mode = "list", length = length(mediator) )
    for( i in 1:length(mediator) ){

      mediatorAdjustmentSetCanonical[i] <- adjustmentSets(dagified, exposure = mediator[i], outcome = response , type = "canonical", effect = c("total", "direct"))

    }

#################################################################################################################
##############################################
## ADD IN CHECKS REQUESTED BY John Ferguson
##############################################
#############################################

Full_MediatorAdjustmentSet <- vector(mode = "list", length = length( mediator ) )
Subset_MediatorAdjustmentSet <- vector(mode = "list", length = length( mediator ) )

updateMediatorAdjustSetCanon_noResponse <- vector(mode = "list", length = length( mediator ) )
updateMediatorAdjustSetCanon_noResponse_All <- vector(mode = "list", length = length( mediator ) )

indexMediator  <- which(lapply(in_outDAG_SplinesRemoved[[2]] , function(data_input) data_input %in% mediator )  > 0 )

in_outDAG_SplinesRemoved_IndexUpdated <- vector(mode = "list", length = length( mediator ) )

mediatorAdjustmentSetAll  <- vector(mode = "list", length = length(mediator) )

Full_MediatorAdjustmentSet_All <- vector(mode = "list", length = length( mediator ) )
Subset_MediatorAdjustmentSet_ALL <- vector(mode = "list", length = length( mediator ) )


# # The elements of setdiff(x,y) are those elements in x but not in y.
# # setdiff( c("phys", "ahei3tert" , "nevfcur", "alcohfreqwk","global_stress2"), exposure )
# # setdiff( c( "ahei3tert" , "nevfcur", "alcohfreqwk","global_stress2"), exposure )
# setdiff( in_outDAG_SplinesRemoved[[1]][[indexMediator[i] ]], exposure)

#######
## dagitty
## adjustmentSets()
## For type="canonical", a single adjustment set is returned that consists of all (possible) ancestors of exposures and outcomes, minus (possible) descendants of nodes on proper causal paths. This canonical adjustment set is always valid if any valid set exists at all.
#######

###########
###########
### Testing
# create test cases for each if else loop below
###########
###########
## Scenario 1: when mediator[1] is "subhtn" and the in_outDAG_SplinesRemoved[[1]][[indexMediator[1] ]] is set incorrectly as follows:
### It should be as follows:
### in_outDAG_SplinesRemoved[[1]][[indexMediator[1] ]] <- c("subeduc" , "moteduc" , "fatduc" , "phys" , "ahei3tert", "nevfcur", "alcohfreqwk","global_stress2")
### But is changed to this
# in_outDAG_SplinesRemoved[[1]][[indexMediator[1] ]] <- c("subeduc" , "phys" )

## Scenario 2: Same as scenario 1 but skip onto next part of loop.
## Scenario 3: Same as scenario 1 but skip onto next part of loop.


for( i in 1:length( mediator ) ){

        ## for loop after then
        if ( all(in_outDAG_SplinesRemoved[[1]][[indexMediator[i] ]] %in% mediatorAdjustmentSetCanonical[[i]]) & isAdjustmentSet(dagified, setdiff( in_outDAG_SplinesRemoved[[1]][[indexMediator[i] ]], exposure ), exposure = exposure, outcome = mediator[i]) ){

           ## OK MAYBE DO NOTHING, IT PASSES OUR CHECK
           next

        } else{

                ## need to for loop to see if an adjustment set exists
                ###############################################################################
                # Try canonical set first.Before looping through all possible adjustment sets.
                ###############################################################################
                #updateMediatorAdjustSetCanon_noResponse[[i]] <- adjustmentSets(dagified, exposure = in_outDAG_SplinesRemoved[[1]][[indexMediator[i] ]], outcome = mediator[i] , type = "canonical", effect = c("total", "direct"))
                updateMediatorAdjustSetCanon_noResponse[[i]] <- adjustmentSets(dagified, exposure = exposure, outcome = mediator[i] , type = "canonical", effect = c("total", "direct"))

                ##############
                ##############
                ### 1. >0
                ##############
                ##############
                ###########

                if( length(updateMediatorAdjustSetCanon_noResponse[[i]][[1]]) > 0 ){
                ###########
                ### Start
                ###########
                ###########

                        Subset_MediatorAdjustmentSet[[i]]  <- updateMediatorAdjustSetCanon_noResponse[[i]][[1]]

                        Full_MediatorAdjustmentSet[[i]] <- c( exposure, Subset_MediatorAdjustmentSet[[i]] )

                        if( all( Full_MediatorAdjustmentSet[[i]] %in% mediatorAdjustmentSetCanonical[[i]]) & isAdjustmentSet(dagified,setdiff(Full_MediatorAdjustmentSet[[i]] , exposure ) , exposure = exposure, outcome = mediator[i]) ){
                                in_outDAG_SplinesRemoved[[1]][[ indexMediator[i] ]] <- Full_MediatorAdjustmentSet[[i]]
                                #1. store that this has been updated
                                in_outDAG_SplinesRemoved_IndexUpdated[[i]] <- indexMediator[i]
                                next
                        } else{
##################################################################
##################################################################
                        ## **Loop over all adjustment sets. Keeping canonical mediator outcome adjustment set fixed ## SAME CODE AS BELOW
##################################################################
##################################################################
                               updateMediatorAdjustSetCanon_noResponse_All[[i]] <- adjustmentSets(dagified, exposure = exposure, outcome = mediator[i] , type = "all", effect = c("total", "direct"))

                               stopCounter <- FALSE
                               for( j in 1:length( updateMediatorAdjustSetCanon_noResponse_All[[i]] ) ){

                                        # Then already dealth with this case above, so just want to populate list rather than leaving it empty.
                                        if( length(updateMediatorAdjustSetCanon_noResponse_All[[i]][[j]]) == 0 ){

                                                Subset_MediatorAdjustmentSet[[i]][[j]]  <- updateMediatorAdjustSetCanon_noResponse_All[[i]][[j]]

                                                Full_MediatorAdjustmentSet[[i]][[j]] <- c( exposure )

                                                if( all( Full_MediatorAdjustmentSet[[i]][[j]] %in% mediatorAdjustmentSetCanonical[[i]]) & isAdjustmentSet(dagified, c(), exposure = exposure, outcome = mediator[i]) ){
                                                        in_outDAG_SplinesRemoved[[1]][[ indexMediator[i] ]] <- Full_MediatorAdjustmentSet[[i]][[j]]
                                                        #2. store that this has been updated
                                                        in_outDAG_SplinesRemoved_IndexUpdated[[i]] <- indexMediator[i]
                                                        stopCounter <- TRUE
                                                        break
                                                }

                                        } else if( length(updateMediatorAdjustSetCanon_noResponse_All[[i]][[j]]) > 0 ){
                                              # Check if issue if Subset_MediatorAdjustmentSet[[i]] in one iteration of for loop Subset_MediatorAdjustmentSet[[i]][[j]] in the other iteration of the for loop
                                              # ONLY issue if [[i]] is replace with a call of form [[i]][[j]] and the i is the same.
                                              Subset_MediatorAdjustmentSet[[i]][[j]]  <- updateMediatorAdjustSetCanon_noResponse_All[[i]][[j]]

                                              Full_MediatorAdjustmentSet[[i]][[j]] <- c( exposure,
                                                                                         Subset_MediatorAdjustmentSet[[i]][[j]] )

                                              if( all( Full_MediatorAdjustmentSet[[i]][[j]] %in% mediatorAdjustmentSetCanonical[[i]]) & isAdjustmentSet(dagified, setdiff( Full_MediatorAdjustmentSet[[i]][[j]] , exposure ) , exposure = exposure, outcome = mediator[i]) ){
                                                      in_outDAG_SplinesRemoved[[1]][[ indexMediator[i] ]] <- Full_MediatorAdjustmentSet[[i]][[j]]
                                                      #3. store that this has been updated
                                                      in_outDAG_SplinesRemoved_IndexUpdated[[i]] <- indexMediator[i]
                                                      stopCounter <- TRUE
                                                      break
                                              }
                                        }
                            }
      ##################################################################
      ##################################################################

                                    ###############################################################
                                    ###############################################################
                                    ###############################################################
                                    ##### Loop over both (1) mediator exposure adjustment sets and (2) mediator outcome adjustment sets
                                    ###############################################################
                                    ###############################################################
                                    ###############################################################

                                    stop = FALSE
                                    if( stopCounter == FALSE ){


                                        # for( p in 1:length(mediator) ){

                                          # mediatorAdjustmentSetCanonical[p] <- adjustmentSets(dagified, exposure = mediator[p], outcome = response , type = "canonical", effect = c("total", "direct"))
                                          mediatorAdjustmentSetAll[[i]] <- adjustmentSets(dagified, exposure = mediator[i], outcome = response , type = "all", effect = c("total", "direct"))

                                        # }

                                      for( m in 1:length( mediatorAdjustmentSetAll[[i]] ) ){

                                            for(j in 1:length( updateMediatorAdjustSetCanon_noResponse_All[[i]] ) ){

                                                                      # Then already dealth with this case above, so just want to populate list rather than leaving it empty.
                                                                      if( length(updateMediatorAdjustSetCanon_noResponse_All[[i]][[j]]) == 0 ){

                                                                              Subset_MediatorAdjustmentSet_ALL[[i]][[j]]  <- updateMediatorAdjustSetCanon_noResponse_All[[i]][[j]]

                                                                              Full_MediatorAdjustmentSet_All[[i]][[j]] <- c( exposure )

                                                                              if( all( Full_MediatorAdjustmentSet_All[[i]][[j]] %in% mediatorAdjustmentSetAll[[i]][[m]]) & isAdjustmentSet(dagified, c(), exposure = exposure, outcome = mediator[i]) ){
                                                                                      in_outDAG_SplinesRemoved[[1]][[ indexMediator[i] ]] <- Full_MediatorAdjustmentSet_All[[i]][[j]]
                                                                                      #2. store that this has been updated
                                                                                      in_outDAG_SplinesRemoved_IndexUpdated[[i]] <- indexMediator[i]
                                                                                      # NEED TO CHECK IF THIS CHANGE TO mediatorAdjustmentSetCanonical[i] MAKES A DIFFERENCE
                                                                                       mediatorAdjustmentSetCanonical[[i]] <- mediatorAdjustmentSetAll[[i]][[m]]
                                                                                      stop <- TRUE
                                                                                      break
                                                                              }

                                                                      } else if( length(updateMediatorAdjustSetCanon_noResponse_All[[i]][[j]]) > 0 ){
                                                                            # Check if issue if Subset_MediatorAdjustmentSet[[i]] in one iteration of for loop Subset_MediatorAdjustmentSet[[i]][[j]] in the other iteration of the for loop
                                                                            # ONLY issue if [[i]] is replace with a call of form [[i]][[j]] and the i is the same.
                                                                            Subset_MediatorAdjustmentSet_ALL[[i]][[j]]  <- updateMediatorAdjustSetCanon_noResponse_All[[i]][[j]]

                                                                            Full_MediatorAdjustmentSet_All[[i]][[j]] <- c( exposure,
                                                                                                                       Subset_MediatorAdjustmentSet_ALL[[i]][[j]] )

                                                                            if( all( Full_MediatorAdjustmentSet_All[[i]][[j]] %in% mediatorAdjustmentSetAll[[i]][[m]]) & isAdjustmentSet(dagified, setdiff( Full_MediatorAdjustmentSet_All[[i]][[j]] , exposure ) , exposure = exposure, outcome = mediator[i]) ){
                                                                                    in_outDAG_SplinesRemoved[[1]][[ indexMediator[i] ]] <- Full_MediatorAdjustmentSet_All[[i]][[j]]
                                                                                    #3. store that this has been updated
                                                                                    in_outDAG_SplinesRemoved_IndexUpdated[[i]] <- indexMediator[i]
                                                                                    mediatorAdjustmentSetCanonical[[i]] <- mediatorAdjustmentSetAll[[i]][[m]]
                                                                                    stop <- TRUE
                                                                                    break
                                                                            }
                                                                      }
                                            }
                                            if(stop){ break }
                                      }


                                    }

                                    ###############################################################
                                    ###############################################################
                                    ###############################################################
                                    ## END
                                    ###############################################################
                                    ###############################################################
                                    ###############################################################




                        }
                ###########
                ###########
                ### End
                ###########
                ###########
                }


                ################
                ################
                ## 2. == 0
                ################
                ################
                if( (length(updateMediatorAdjustSetCanon_noResponse[[i]]) == 0) & isAdjustmentSet(dagified, c(), exposure = exposure, outcome = mediator[i]) & all( exposure %in% mediatorAdjustmentSetCanonical[[i]])  ){

                        in_outDAG_SplinesRemoved[[1]][[ indexMediator[i] ]] <- exposure
                        #4. store that this has been updated
                        in_outDAG_SplinesRemoved_IndexUpdated[[i]] <- indexMediator[i]
                        next

                        } else if( length(updateMediatorAdjustSetCanon_noResponse[[i]]) == 0 ){

      ##################################################################
      ##################################################################
                              ## **Loop over all adjustment sets.Keeping canonical mediator outcome adjustment set fixed ## SAME CODE AS ABOVE
      ##################################################################
      ##################################################################
                              # Then need to do a biggier lopp
                               # updateMediatorAdjustSetCanon_noResponse_All[[i]] <- adjustmentSets(dagified, exposure = in_outDAG_SplinesRemoved[[1]][[indexMediator[i] ]], outcome = mediator[i] , type = "all", effect = c("total", "direct"))
                               updateMediatorAdjustSetCanon_noResponse_All[[i]] <- adjustmentSets(dagified, exposure = exposure, outcome = mediator[i] , type = "all", effect = c("total", "direct"))

                               stopCounter <- FALSE
                               for( j in 1:length( updateMediatorAdjustSetCanon_noResponse_All[[i]] ) ){

                                        # Then already dealth with this case above, so just want to populate list rather than leaving it empty.
                                        if( length(updateMediatorAdjustSetCanon_noResponse_All[[i]][[j]]) == 0 ){

                                                Subset_MediatorAdjustmentSet[[i]][[j]]  <- updateMediatorAdjustSetCanon_noResponse_All[[i]][[j]]

                                                # c( in_outDAG_SplinesRemoved[[1]][[i]], Subset_MediatorAdjustmentSet[[i]][[j]] )
                                                Full_MediatorAdjustmentSet[[i]][[j]] <- c( exposure )

                                                if( all( Full_MediatorAdjustmentSet[[i]][[j]] %in% mediatorAdjustmentSetCanonical[[i]]) & isAdjustmentSet(dagified, c(), exposure = exposure, outcome = mediator[i]) ){
                                                        in_outDAG_SplinesRemoved[[1]][[ indexMediator[i] ]] <- Full_MediatorAdjustmentSet[[i]]
                                                        #5. store that this has been updated
                                                        in_outDAG_SplinesRemoved_IndexUpdated[[i]] <- indexMediator[i]
                                                        stopCounter <- TRUE
                                                        break
                                                        ## NOTE MIGHT NEED ANOTHER BREAK ALSO TO EXIT OTHER FOR LOOPS
                                                }

                                        } else if( length(updateMediatorAdjustSetCanon_noResponse_All[[i]][[j]]) > 0 ){

                                              Subset_MediatorAdjustmentSet[[i]][[j]]  <- updateMediatorAdjustSetCanon_noResponse_All[[i]][[j]]

                                              Full_MediatorAdjustmentSet[[i]][[j]] <- c( exposure,
                                                                                         Subset_MediatorAdjustmentSet[[i]][[j]] )

                                              if( all( Full_MediatorAdjustmentSet[[i]][[j]] %in% mediatorAdjustmentSetCanonical[[i]]) & isAdjustmentSet(dagified, setdiff(Full_MediatorAdjustmentSet[[i]][[j]] , exposure ) , exposure = exposure, outcome = mediator[i]) ){
                                                      in_outDAG_SplinesRemoved[[1]][[ indexMediator[i] ]] <- Full_MediatorAdjustmentSet[[i]]
                                                      #6. store that this has been updated
                                                      in_outDAG_SplinesRemoved_IndexUpdated[[i]] <- indexMediator[i]
                                                      stopCounter <- TRUE
                                                      break
                                                        ## NOTE MIGHT NEED ANOTHER BREAK ALSO TO EXIT OTHER FOR LOOPS
                                              }
                                      }
                               }


                                  ###############################################################
                                  ###############################################################
                                  ###############################################################
                                  ##### Loop over both (1) mediator exposure adjustment sets and (2) mediator outcome adjustment sets
                                  ###############################################################
                                  ###############################################################
                                  ###############################################################

                                  stop = FALSE
                                  if( stopCounter == FALSE ){


                                      # for( p in 1:length(mediator) ){

                                        # mediatorAdjustmentSetCanonical[p] <- adjustmentSets(dagified, exposure = mediator[p], outcome = response , type = "canonical", effect = c("total", "direct"))
                                        mediatorAdjustmentSetAll[[i]] <- adjustmentSets(dagified, exposure = mediator[i], outcome = response , type = "all", effect = c("total", "direct"))

                                      # }

                                    for( m in 1:length( mediatorAdjustmentSetAll[[i]] ) ){

                                          for(j in 1:length( updateMediatorAdjustSetCanon_noResponse_All[[i]] ) ){

                                                                    # Then already dealth with this case above, so just want to populate list rather than leaving it empty.
                                                                    if( length(updateMediatorAdjustSetCanon_noResponse_All[[i]][[j]]) == 0 ){

                                                                            Subset_MediatorAdjustmentSet_ALL[[i]][[j]]  <- updateMediatorAdjustSetCanon_noResponse_All[[i]][[j]]

                                                                            Full_MediatorAdjustmentSet_All[[i]][[j]] <- c( exposure )

                                                                            if( all( Full_MediatorAdjustmentSet_All[[i]][[j]] %in% mediatorAdjustmentSetAll[[i]][[m]]) & isAdjustmentSet(dagified, c(), exposure = exposure, outcome = mediator[i]) ){
                                                                                    in_outDAG_SplinesRemoved[[1]][[ indexMediator[i] ]] <- Full_MediatorAdjustmentSet_All[[i]][[j]]
                                                                                    #2. store that this has been updated
                                                                                    in_outDAG_SplinesRemoved_IndexUpdated[[i]] <- indexMediator[i]
                                                                                    # NEED TO CHECK IF THIS CHANGE TO mediatorAdjustmentSetCanonical[i] MAKES A DIFFERENCE
                                                                                     mediatorAdjustmentSetCanonical[[i]] <- mediatorAdjustmentSetAll[[i]][[m]]
                                                                                    stop <- TRUE
                                                                                    break
                                                                            }

                                                                    } else if( length(updateMediatorAdjustSetCanon_noResponse_All[[i]][[j]]) > 0 ){
                                                                          # Check if issue if Subset_MediatorAdjustmentSet[[i]] in one iteration of for loop Subset_MediatorAdjustmentSet[[i]][[j]] in the other iteration of the for loop
                                                                          # ONLY issue if [[i]] is replace with a call of form [[i]][[j]] and the i is the same.
                                                                          Subset_MediatorAdjustmentSet_ALL[[i]][[j]]  <- updateMediatorAdjustSetCanon_noResponse_All[[i]][[j]]

                                                                          Full_MediatorAdjustmentSet_All[[i]][[j]] <- c( exposure,
                                                                                                                     Subset_MediatorAdjustmentSet_ALL[[i]][[j]] )

                                                                          if( all( Full_MediatorAdjustmentSet_All[[i]][[j]] %in% mediatorAdjustmentSetAll[[i]][[m]]) & isAdjustmentSet(dagified, setdiff( Full_MediatorAdjustmentSet_All[[i]][[j]] , exposure ) , exposure = exposure, outcome = mediator[i]) ){
                                                                                  in_outDAG_SplinesRemoved[[1]][[ indexMediator[i] ]] <- Full_MediatorAdjustmentSet_All[[i]][[j]]
                                                                                  #3. store that this has been updated
                                                                                  in_outDAG_SplinesRemoved_IndexUpdated[[i]] <- indexMediator[i]
                                                                                  mediatorAdjustmentSetCanonical[[i]] <- mediatorAdjustmentSetAll[[i]][[m]]
                                                                                  stop <- TRUE
                                                                                  break
                                                                          }
                                                                    }
                                          }
                                          if(stop){ break }
                                    }


                                  }

                                  ###############################################################
                                  ###############################################################
                                  ###############################################################
                                  ## END
                                  ###############################################################
                                  ###############################################################
                                  ###############################################################

                  }
        }

}

#############################################
#############################################
#############################################
#################################################################################################################



######################
######################
######################
# # if there are splines
## NB NEED TO ADD IN IF STATEMENT CODE AS MIGHT NOT BE ANY SPLINES
########################
######################
######################
# if there are splines
if( length( splinesVariablesIndices) > 0 ){

          ##########################################################
          ##########################################################
          ### Need to undo splines then i.e. get splines back again in two places and to return them from function with splines, maybe plot casual DAG also returned
              # 1. Add splines back into exposureAdjustmentSetCanonical
              # 2. Add splines back into mediatorAdjustmentSetCanonical
              # 3. Return plot of causal DAG
              # 4. NO need to add into in_outDAG_SplinesRemoved within function, as outside function in_outDAG should still contain the splines
              # 5. If in_outDAG fails the check above and is updated with a new adjustment set for at least one of the mediators, then we need to update in_outDAG and ensure splines are added in.
          ##########################################################
          ##########################################################

          # # if there are splines
          # if( length( splinesVariablesIndices) > 0 ){
          #
          # }

          # 1. Add splines back into exposureAdjustmentSetCanonical
          for ( i in 1:length(exposureAdjustmentSetCanonical) ) {

                    if( any(splinesVariables %in% exposureAdjustmentSetCanonical[[i]] ) ){

                                   for( j in 1:length( splinesVariables ) ){

                                          indicesAddSplinesInExposureAdjSet <- which( exposureAdjustmentSetCanonical[[i]] == in_outDAG[[2]][ splinesVariablesIndices[1]:splinesVariablesIndices[length(splinesVariablesIndices)] ][j]  )

                                        if( length(indicesAddSplinesInExposureAdjSet) == 0 ){

                                                # Do nothing since no splines in adjustment set.

                                        } else if( length(indicesAddSplinesInExposureAdjSet) > 0 ){

                                                 exposureAdjustmentSetCanonical[[ i]][[indicesAddSplinesInExposureAdjSet ]] <- Splines_outlist_Var[[1]][ splinesVariablesIndices[j]  ]

                                        } else{

                                                stop("Error in processing Splines_outlist_Var. Splines_outlist_Var variable may be defined incorrectly.")

                                        }


                                       }

                    }
          }

          # 2. Add splines back into mediatorAdjustmentSetCanonical
          #if( any(splinesVariables %in% unlist(mediatorAdjustmentSetCanonical) ) ){

                    # indicesAddSplinesInMediatorAdjSet <- vector(mode = "list", length = length(splinesVariables) )

          for ( i in 1:length(mediatorAdjustmentSetCanonical) ) {

                        if( any(splinesVariables %in% mediatorAdjustmentSetCanonical[[i]] ) ){

                                   for( j in 1:length( splinesVariables ) ){

                                        # indicesAddSplinesInMediatorAdjSet <- which( mediatorAdjustmentSetCanonical[[1]] == in_outDAG[[2]][ splinesVariablesIndices[1]:splinesVariablesIndices[length(splinesVariablesIndices)] ][1]  )
                                          indicesAddSplinesInMediatorAdjSet <- which( mediatorAdjustmentSetCanonical[[i]] == in_outDAG[[2]][ splinesVariablesIndices[1]:splinesVariablesIndices[length(splinesVariablesIndices)] ][j]  )

                                        if( length(indicesAddSplinesInMediatorAdjSet) == 0 ){

                                                # Do nothing since no splines in adjustment set.

                                        } else if( length(indicesAddSplinesInMediatorAdjSet) > 0 ){

                                                 mediatorAdjustmentSetCanonical[[ i]][[indicesAddSplinesInMediatorAdjSet ]] <- Splines_outlist_Var[[1]][ splinesVariablesIndices[j]  ]

                                        } else{

                                                stop("Error in processing Splines_outlist_Var. Splines_outlist_Var variable may be defined incorrectly. E.G. Is Splines_outlist_Var defined in form list(c(text1_invertedcommas,text2_invertedcommas,text3_invertedcommas)) with splines names in text frmat with inverted commas.")

                                        }


                                   }
                        }

           }


          # 3. Return plot of causal DAG

             # CausalDagPlot

          # 4.(A) Are any of mediator=c("subhtn","apob_apoa","whr") splines? If so return them as splines

              # mediator=c("subhtn","apob_apoa","whr")

               # splinesInMediator <- splinesVariables[which( splinesVariables %in% mediator )]
               splinesInMediatorIndices <- which( splinesVariables %in% mediator )

               MediatorIndicesToSplines <- which(  mediator %in% splinesVariables  )
               #### NEED TO CHECK IF WORKS FOR ANY ORDER OF mediator=c("subhtn","apob_apoa","whr") e.g. mediator=c("whr","subhtn","apob_apoa")
               if( length(splinesInMediatorIndices) > 0 ){

                    # splinesInMediator <- splinesVariables[which( splinesVariables %in% mediator )]
                      splinesInMediator <- splinesVariables[ splinesInMediatorIndices ]

                      splinesInMediatorIndicesInin_outDAG2 <- which(lapply(in_outDAG[[2]], function(data_input) all(data_input %in% splinesInMediator )  ) > 0 )



                    mediatorsWithSplinesReturn <- mediator

                    ### WILL ONLY WORK IF Splines_outlist_Var DEFINED IN FORMAT list(c("splinename1","splinename2","splinename3")) since need Splines_outlist_Var[[1]][splinesInMediatorIndicesInin_outDAG2] in code
                    mediatorsWithSplinesReturn[MediatorIndicesToSplines] <- Splines_outlist_Var[[1]][splinesInMediatorIndicesInin_outDAG2]


               } else{

                 mediatorsWithSplinesReturn <- mediator
               }
          #### NEED TO CHECK IF WORKS FOR ANY ORDER OF mediator=c("subhtn","apob_apoa","whr") e.g. mediator=c("whr","subhtn","apob_apoa
          #### NEED TO CHECK IF WORKS FOR ANY ORDER OF mediator=c("subhtn","apob_apoa","whr") e.g. mediator=c("whr","subhtn","apob_apoa
               #### NEED TO CHECK IF WORKS FOR ANY ORDER OF mediator=c("subhtn","apob_apoa","whr") e.g. mediator=c("whr","subhtn","apob_apoa
               #### NEED TO CHECK IF WORKS FOR ANY ORDER OF mediator=c("subhtn","apob_apoa","whr") e.g. mediator=c("whr","subhtn","apob_apoa
               #### NEED TO CHECK IF WORKS FOR ANY ORDER OF mediator=c("subhtn","apob_apoa","whr") e.g. mediator=c("whr","subhtn","apob_apoa


               # 4.(B) Are any of exposures=c(, , ) splines? (could be more than one expsoure) If so return them as splines.

              # exposure="phys"   or exposure=c("phys", "...")

                   # splinesInMediator <- splinesVariables[which( splinesVariables %in% mediator )]
               splinesInExposureIndices <- which( splinesVariables %in% exposure )

               ExposureIndicesToSplines <- which(  exposure %in% splinesVariables  )
               #### NEED TO CHECK IF WORKS FOR ANY ORDER OF exposure=c("subhtn","apob_apoa","whr") e.g. exposure=c("whr","subhtn","apob_apoa")
               if( length(splinesInExposureIndices) > 0 ){

                    # splinesInExposure <- splinesVariables[which( splinesVariables %in% exposure )]
                      splinesInExposure <- splinesVariables[ splinesInExposureIndices ]

                      splinesInExposureIndicesInin_outDAG2 <- which(lapply(in_outDAG[[2]], function(data_input) all(data_input %in% splinesInExposure )  ) > 0 )



                    exposuresWithSplinesReturn <- exposure

                    ### WILL ONLY WORK IF Splines_outlist_Var DEFINED IN FORMAT list(c("splinename1","splinename2","splinename3")) since need Splines_outlist_Var[[1]][splinesInExposureIndicesInin_outDAG2] in code
                    exposuresWithSplinesReturn[ExposureIndicesToSplines] <- Splines_outlist_Var[[1]][splinesInExposureIndicesInin_outDAG2]


               } else{

                  exposuresWithSplinesReturn <- exposure

               }
          #### NEED TO CHECK IF WORKS FOR ANY ORDER OF exposure=c("subhtn","apob_apoa","whr") e.g. exposure=c("whr","subhtn","apob_apoa
          #### NEED TO CHECK IF WORKS FOR ANY ORDER OF exposure=c("subhtn","apob_apoa","whr") e.g. exposure=c("whr","subhtn","apob_apoa
               #### NEED TO CHECK IF WORKS FOR ANY ORDER OF exposure=c("subhtn","apob_apoa","whr") e.g. exposure=c("whr","subhtn","apob_apoa
               #### NEED TO CHECK IF WORKS FOR ANY ORDER OF exposure=c("subhtn","apob_apoa","whr") e.g. exposure=c("whr","subhtn","apob_apoa
               #### NEED TO CHECK IF WORKS FOR ANY ORDER OF exposure=c("subhtn","apob_apoa","whr") e.g. exposure=c("whr","subhtn","apob_apoa

########################################################
########################################################
          # 5. If in_outDAG fails the check above and is updated with a new adjustment set for at least one of the mediators, then we need to update in_outDAG and ensure splines are added in.

          indicesUpdated <- which(lapply(in_outDAG_SplinesRemoved_IndexUpdated , function(data_input) !is.null(data_input)  ) > 0 )

          if( length( indicesUpdated ) > 0 ){

                in_outDAG_SplinesAddIn <- in_outDAG

                for( w in 1:length(indicesUpdated) ){

                      # This is the mediator adjustment set that was updated and we need to identify if there are variables which will be analysed as splines
                      in_outDAG_SplinesAddIn[[1]][[  in_outDAG_SplinesRemoved_IndexUpdated[[ indicesUpdated[w] ]]   ]] <- in_outDAG_SplinesRemoved[[1]][[  in_outDAG_SplinesRemoved_IndexUpdated[[ indicesUpdated[w] ]]   ]]



                      # for ( i in 1:length(exposureAdjustmentSetCanonical) ) {
                      # for ( i in 1:length(in_outDAG_SplinesAddIn[[1]][[  in_outDAG_SplinesRemoved_IndexUpdated[[ indicesUpdated[w] ]]   ]]) ) {

                            # if( any(splinesVariables %in% exposureAdjustmentSetCanonical[[i]] ) ){
                              if( any(splinesVariables %in% in_outDAG_SplinesAddIn[[1]][[  in_outDAG_SplinesRemoved_IndexUpdated[[ indicesUpdated[w] ]]   ]] ) ){

                                  for( j in 1:length( splinesVariables ) ){

                                        # indicesAddSplinesInExposureAdjSet <- which( exposureAdjustmentSetCanonical[[i]] == in_outDAG[[2]][ splinesVariablesIndices[1]:splinesVariablesIndices[length(splinesVariablesIndices)] ][j]  )
                                        indicesAddSplinesInMediatorCheckAdjSet <- which( in_outDAG_SplinesAddIn[[1]][[  in_outDAG_SplinesRemoved_IndexUpdated[[ indicesUpdated[w] ]]   ]] == in_outDAG[[2]][ splinesVariablesIndices[1]:splinesVariablesIndices[length(splinesVariablesIndices)] ][j]  )

                                        if( length(indicesAddSplinesInMediatorCheckAdjSet) == 0 ){

                                                # Do nothing since no splines in adjustment set.

                                        } else if( length(indicesAddSplinesInMediatorCheckAdjSet) > 0 ){

                                                 # exposureAdjustmentSetCanonical[[ i]][[indicesAddSplinesInMediatorCheckAdjSet ]] <- Splines_outlist_Var[[1]][ splinesVariablesIndices[j]  ]
                                                 in_outDAG_SplinesAddIn[[1]][[  in_outDAG_SplinesRemoved_IndexUpdated[[ indicesUpdated[w] ]]   ]][[indicesAddSplinesInMediatorCheckAdjSet ]] <- Splines_outlist_Var[[1]][ splinesVariablesIndices[j]  ]

                                        } else{

                                                stop("Error in processing Splines_outlist_Var. Splines_outlist_Var variable may be defined incorrectly.")

                                        }
                                 }
                          }
                     # }
                }
           }

########################################################
########################################################


} else{

          # There are no splines. To avoid longer code just code this as follows if there are no splines.
          mediatorsWithSplinesReturn <- mediator

          exposuresWithSplinesReturn <- exposure

}

# # 4.(B) Are any of exposures=c(, , ) splines? (could be more than one expsoure) If so return them as splines.
#
#               # exposure="phys"   or exposure=c("phys", "...")
#
#                    # splinesInMediator <- splinesVariables[which( splinesVariables %in% mediator )]
#                splinesInExposureIndices <- which( splinesVariables %in% exposure )
#
#                ExposureIndicesToSplines <- which(  exposure %in% splinesVariables  )
#                #### NEED TO CHECK IF WORKS FOR ANY ORDER OF exposure=c("subhtn","apob_apoa","whr") e.g. exposure=c("whr","subhtn","apob_apoa")
#                if( length(splinesInExposureIndices) > 0 ){
#
#                     # splinesInExposure <- splinesVariables[which( splinesVariables %in% exposure )]
#                       splinesInExposure <- splinesVariables[ splinesInExposureIndices ]
#
#                       splinesInExposureIndicesInin_outDAG2 <- which(lapply(in_outDAG[[2]], function(data_input) all(data_input %in% splinesInExposure )  ) > 0 )
#
#
#
#                     exposuresWithSplinesReturn <- exposure
#
#                     ### WILL ONLY WORK IF Splines_outlist_Var DEFINED IN FORMAT list(c("splinename1","splinename2","splinename3")) since need Splines_outlist_Var[[1]][splinesInExposureIndicesInin_outDAG2] in code
#                     exposuresWithSplinesReturn[ExposureIndicesToSplines] <- Splines_outlist_Var[[1]][splinesInExposureIndicesInin_outDAG2]
#
#
#                }
#           #### NEED TO CHECK IF WORKS FOR ANY ORDER OF exposure=c("subhtn","apob_apoa","whr") e.g. exposure=c("whr","subhtn","apob_apoa
#           #### NEED TO CHECK IF WORKS FOR ANY ORDER OF exposure=c("subhtn","apob_apoa","whr") e.g. exposure=c("whr","subhtn","apob_apoa
#                #### NEED TO CHECK IF WORKS FOR ANY ORDER OF exposure=c("subhtn","apob_apoa","whr") e.g. exposure=c("whr","subhtn","apob_apoa
#                #### NEED TO CHECK IF WORKS FOR ANY ORDER OF exposure=c("subhtn","apob_apoa","whr") e.g. exposure=c("whr","subhtn","apob_apoa
#                #### NEED TO CHECK IF WORKS FOR ANY ORDER OF exposure=c("subhtn","apob_apoa","whr") e.g. exposure=c("whr","subhtn","apob_apoa
#
#
#
# } else{
#
#           # There are no splines. To avoid longer code just code this as follows if there are no splines.
#           exposuresWithSplinesReturn <- exposure
#
# }



##########################################
##########################################
##########################################

# This can be done in other function to add in custom.Maybe with make_formula function
# mediatoraddCustom = FALSE, custom = "~ regionnn7*ns(eage,df=5)+esex*ns(eage,df=5) + "

#################################
#################################

     if( addCustomExposureAdjustmentSet ){

       # 1. How to input customExposureAdjustmentSet of same length as length(exposure) but only splines for some exposures and not others
       # Require it input in the form c(FALSE, "custom input here for 2nd exposure only", FALSE)
        if( (length(customExposureAdjustmentSet) == length(exposure)) & ( FALSE %in% customExposureAdjustmentSet )  ){

                     indicesNotFALSE_exposure <- setdiff(1:length(customExposureAdjustmentSet), which(lapply(customExposureAdjustmentSet, function(data_input) FALSE %in% data_input   ) > 0 ) )

                     for( i in indicesNotFALSE_exposure ){

                      # append(x, values, after = length(x))
                      exposureAdjustmentSetCanonical[[i]]  <- append( exposureAdjustmentSetCanonical[[i]], customExposureAdjustmentSet[i], after = length( exposureAdjustmentSetCanonical[[i]] ))

                      }
        }
       # Consider 4 scenarios (other scenarios possible)
       # 2. length(customExposureAdjustmentSet) == length(exposure)
            else if( length(customExposureAdjustmentSet) == length(exposure) ){

                      for( i in 1:length(customExposureAdjustmentSet) ){

                      # append(x, values, after = length(x))
                      exposureAdjustmentSetCanonical[[i]]  <- append( exposureAdjustmentSetCanonical[[i]], customExposureAdjustmentSet[i], after = length( exposureAdjustmentSetCanonical[[i]] ))

                      }
        # 3. length(customExposureAdjustmentSet) == 1 with no FALSE inputs as in 2.   ASSUME IT IS TO BE APPLIED TO ALL exposure regressions.
            } else if( length(customExposureAdjustmentSet) == 1 ){

                    for( i in 1:length(exposureAdjustmentSetCanonical) ){

                      # append(x, values, after = length(x))
                      exposureAdjustmentSetCanonical[[i]]  <- append( exposureAdjustmentSetCanonical[[i]], customExposureAdjustmentSet[1], after = length( exposureAdjustmentSetCanonical[[i]] ))

                      }

            }
       # 4. length(customExposureAdjustmentSet) < length(exposure) and not equal to 1 with no FALSE in it
        else if( (length(customExposureAdjustmentSet) < length(exposure)) & ( length(customExposureAdjustmentSet) != 1) ){

              stop("customExposureAdjustmentSet is less than the length(exposure). Please ensure customExposureAdjustmentSet is either the same length as the number of exposures of interest (with FALSE input if no custom applied to that exposure) or alternatively of length 1 if the same custom customExposureAdjustmentSet is to be applied to all exposures of interst.")

            } else{
              stop("Error in if statment in make_DAG.R.")
          }




     }

#################################
#################################


     if( addCustomMediatorAdjustmentSet ){

       # 1. How to input customMediatorAdjustmentSet of same length as length(mediator) but only splines for some mediators and not others
       # Require it input in the form c(FALSE, "custom input here for 2nd mediator only", FALSE)
        if( (length(customMediatorAdjustmentSet) == length(mediator)) & ( FALSE %in% customMediatorAdjustmentSet )  ){

                     indicesNotFALSE <- setdiff(1:length(customMediatorAdjustmentSet), which(lapply(customMediatorAdjustmentSet, function(data_input) FALSE %in% data_input   ) > 0 ) )

                     for( i in indicesNotFALSE ){

                      # append(x, values, after = length(x))
                      mediatorAdjustmentSetCanonical[[i]]  <- append( mediatorAdjustmentSetCanonical[[i]], customMediatorAdjustmentSet[i], after = length( mediatorAdjustmentSetCanonical[[i]] ))

                      }
        }
       # Consider 4 scenarios (other scenarios possible)
       # 2. length(customMediatorAdjustmentSet) == length(mediator)
            else if( length(customMediatorAdjustmentSet) == length(mediator) ){

                      for( i in 1:length(customMediatorAdjustmentSet) ){

                      # append(x, values, after = length(x))
                      mediatorAdjustmentSetCanonical[[i]]  <- append( mediatorAdjustmentSetCanonical[[i]], customMediatorAdjustmentSet[i], after = length( mediatorAdjustmentSetCanonical[[i]] ))

                      }
        # 3. length(customMediatorAdjustmentSet) == 1 with no FALSE inputs as in 2.   ASSUME IT IS TO BE APPLIED TO ALL mediator regressions.
            } else if( length(customMediatorAdjustmentSet) == 1 ){

                    for( i in 1:length(mediatorAdjustmentSetCanonical) ){

                      # append(x, values, after = length(x))
                      mediatorAdjustmentSetCanonical[[i]]  <- append( mediatorAdjustmentSetCanonical[[i]], customMediatorAdjustmentSet[1], after = length( mediatorAdjustmentSetCanonical[[i]] ))

                      }

            }
       # 4. length(customMediatorAdjustmentSet) < length(mediator) and not equal to 1 with no FALSE in it
        else if( (length(customMediatorAdjustmentSet) < length(mediator)) & ( length(customMediatorAdjustmentSet) != 1) ){

              stop("customMediatorAdjustmentSet is less than the length(mediator). Please ensure customMediatorAdjustmentSet is either the same length as the number of mediators of interest (with FALSE input if no custom applied to that exposure) or alternatively of length 1 if the same custom customMediatorAdjustmentSet is to be applied to all mediators of interst.")

            } else{
              stop("Error in if statment in make_DAG.R.")
          }


     }

####################
####################
# create regressions
####################
####################
     # 1. exposure
    resultExposure <- vector(mode = "list", length = length(exposureAdjustmentSetCanonical) )

    for( i in 1:length(exposureAdjustmentSetCanonical) ){
              # result <- paste(outvar,"~ regionnn7*ns(eage,df=5)+esex*ns(eage,df=5) + ",in_vars[1])
              resultExposure[[i]] <- paste(response,"~ ", exposuresWithSplinesReturn[i]," + ",  exposureAdjustmentSetCanonical[[i]][1] )

              if( length( exposureAdjustmentSetCanonical[[i]] ) >= 2 ){

                    for(j in 2:length(exposureAdjustmentSetCanonical[[i]] )){

                        resultExposure[[i]] <- paste(resultExposure[[i]] ," + ", exposureAdjustmentSetCanonical[[i]][j], sep='')

                      }
              }

    }


     # 2. mediator
      resultMediator <- vector(mode = "list", length = length(mediatorAdjustmentSetCanonical) )

    for( i in 1:length(mediatorAdjustmentSetCanonical) ){
              # result <- paste(outvar,"~ regionnn7*ns(eage,df=5)+esex*ns(eage,df=5) + ",in_vars[1])
              # resultMediator[[i]] <- paste(response,"~ ", mediatorAdjustmentSetCanonical[[i]][1] )
              resultMediator[[i]] <- paste(response,"~ ", mediatorsWithSplinesReturn[i]," + ", mediatorAdjustmentSetCanonical[[i]][1] )

              if( length( mediatorAdjustmentSetCanonical[[i]] ) >= 2 ){

                    for(j in 2:length(mediatorAdjustmentSetCanonical[[i]] )){

                        resultMediator[[i]] <- paste(resultMediator[[i]] ," + ", mediatorAdjustmentSetCanonical[[i]][j], sep='')

                      }
              }

    }

####################
####################

# If models below are listed as empty list(), then expect them to be populated otherwise do not populate. Identify this in causalPAFplot.r
# response_model_mediators = response_vs_mediator,
# response_model_exposure = list(),


      if( length( indicesUpdated ) > 0 ){

            my_list_causal <- list(#"plot" = CausalDagPlot,
                                   "exposuresAdjustmentSetCanonical" = exposureAdjustmentSetCanonical,
                                   "mediatorsAdjustmentSetCanonical" = mediatorAdjustmentSetCanonical,
                                   "exposuresWithSplinesReturn" = exposuresWithSplinesReturn,
                                   "mediatorWithSplines" = mediatorsWithSplinesReturn,
                                   "resultExposure" = resultExposure,
                                   "resultMediator" = resultMediator,
                                   "in_outDAG_updatedAfterCheck" = in_outDAG_SplinesAddIn)  # Only returned if check is applied above.

      }else{
            my_list_causal <- list(#"plot" = CausalDagPlot,
                                   "exposuresAdjustmentSetCanonical" = exposureAdjustmentSetCanonical,
                                   "mediatorsAdjustmentSetCanonical" = mediatorAdjustmentSetCanonical,
                                   "exposuresWithSplinesReturn" = exposuresWithSplinesReturn,
                                   "mediatorWithSplines" = mediatorsWithSplinesReturn,
                                   "resultExposure" = resultExposure,
                                   "resultMediator" = resultMediator ,
                                   "in_outDAG_updatedAfterCheck" = NULL)  # Can check if this is NULL to identify if the check has been applied.
      }




    return(my_list_causal)

      # in_outDAG_SplinesRemoved
      # in_outDAG    needs to be updated based on in_outDAG_SplinesRemoved
      # only mediator indices would have been changed


}


