#' @title Checks if the causal DAG satisfies the Markov condition
#' @description The functions checks if the Markov condition holds for the Directed Acyclic Graph (DAG) defined. Sometimes called the Markov assumption, is an assumption made in Bayesian probability theory, that every node in a Bayesian network is conditionally independent of its nondescendants, given its parents. In other words, it is assumed that a node has no bearing on nodes which do not descend from it. This is equivalent to stating that a node is conditionally independent of the entire network, given its Markov blanket. The related Causal Markov condition states that, conditional on the set of all its direct causes, a node is independent of all variables which are not direct causes or direct effects of that node.
#' @param in_out A list of length 2. The first list contains a list of character vectors of the parents of the exposure or risk factor or outcome which are either causes or confounders of the exposure or risk factor or outcome. The second list contains a list of a single name of exposure or risk factor or outcome in form of characters. See tutorial examples for examples.
#' @export
#' @import splines MASS stats forestplot utils grid magrittr checkmate ggplot2 rlist
#' @return  \item{IsMarkovDAG }{Returns a logical TRUE or FALSE whether it is a Markov DAG provided in_out is input as described in the documentation.}
#' \item{in_out}{The in_out list supplied in the function is returns the same of the input if IsMarkovDAG is returned TRUE. If IsMarkovDAG is returned FALSE the order of the in_out list is updated such that all parent variables come before ancestors in both i_out[[1]] and in_out[[2]]. This corrects any error where variables from a given Markov Bayesian DAG are input to the package in the incorrect order.}
#' \item{Reordered}{Reordered is FALSE if in_out is left in the same order as input. Reordered is FALSE if in_out has been reordered so that parents of variables could before descendants.}
#' @examples
#' # Loads some data (fictional Stroke data from the package 'causalPAF')
#' # In this example, we use a small data set called 'strokedata_smallSample' consisting of 5,000
#' # rows of fictional patient data. For more accurate results, a larger data set is available
#' # called 'strokedata'which contains 16,623 rows of fictional patient data. The methodology
#' # applied in the 'causalPAF' package is more accurate the larger the dataset. To use the larger
#' # 'strokedata' dataset, simply call
#' # stroke_reduced <- strokedata
#' stroke_reduced <- strokedata_smallSample
#'
#' in_phys <- c("subeduc","moteduc","fatduc")
#' in_ahei <- c("subeduc","moteduc","fatduc")
#' in_nevfcur <- c("subeduc","moteduc","fatduc")
#' in_alcohfreqwk <- c("subeduc","moteduc","fatduc")
#' in_global_stress2 <- c("subeduc","moteduc","fatduc")
#' in_htnadmbp <- c("subeduc","moteduc","fatduc","phys","ahei3tert","nevfcur","alcohfreqwk",
#'                  "global_stress2")
#' in_apob_apoatert <- c("subeduc","moteduc","fatduc","phys","ahei3tert","nevfcur","alcohfreqwk",
#'                       "global_stress2")
#' in_whrs2tert <- c("subeduc","moteduc","fatduc","phys","ahei3tert","nevfcur","alcohfreqwk",
#'                   "global_stress2")
#' in_cardiacrfcat <- c("subeduc","moteduc","fatduc","phys","ahei3tert","nevfcur","alcohfreqwk",
#'                      "global_stress2", "apob_apoatert","whrs2tert","htnadmbp")
#' in_dmhba1c2 <- c("subeduc","moteduc","fatduc","phys","ahei3tert","nevfcur","alcohfreqwk",
#'                   "global_stress2", "apob_apoatert","whrs2tert","htnadmbp")
#' in_case <- c("subeduc","moteduc","fatduc","phys","ahei3tert","nevfcur","alcohfreqwk",
#' "global_stress2", "apob_apoatert","whrs2tert","htnadmbp","cardiacrfcat","dmhba1c2")
#'
#' in_out <- list(inlist=list(in_phys,in_ahei,in_nevfcur,in_alcohfreqwk,in_global_stress2,
#'                in_htnadmbp, in_apob_apoatert,in_whrs2tert,in_cardiacrfcat,
#'                in_dmhba1c2,in_case),
#'                outlist=c("phys","ahei3tert","nevfcur","alcohfreqwk","global_stress2",
#'                          "htnadmbp","apob_apoatert", "whrs2tert","cardiacrfcat",
#'                          "dmhba1c2","case"))
#'
#'
#' if(checkMarkovDAG(in_out)$IsMarkovDAG & !checkMarkovDAG(in_out)$Reordered){
#'   print("Your in_out DAG is a Markov DAG.")
#' } else if( checkMarkovDAG(in_out)$IsMarkovDAG & checkMarkovDAG(in_out)$Reordered ) {
#'
#'   in_out <- checkMarkovDAG(in_out)[[2]]
#'
#'   print("Your in_out DAG is a Markov DAG.The checkMarkovDAG function has reordered your
#'           in_out list so that all parent variables come before descendants.")
#' } else{ print("Your ``in_out'' list is not a Bayesian Markov DAG so the methods in the
#'                'causalPAF' package cannot be applied for non-Markov DAGs.")}

checkMarkovDAG <- function(in_out){


Boolean <- vector(mode = "list", length = ( length(in_out[[2]]) - 1) )

  for(i in (1:(length( in_out[[2]] ) - 1)) ){
        m <- 1
        for( j in (i + 1):length( in_out[[2]] )){

              Boolean[[i]][[m]] <- all(in_out[[1]][[i]] %in% in_out[[1]][[ j ]] )

              m = m + 1
        }
  }

CheckAllTrue <- vector(mode = "list", length = length(Boolean) )

  for( w in 1:length(CheckAllTrue)){
        CheckAllTrue[[w]] <- all(Boolean[[w]])
  }


List1 <- in_out[[1]]
List2 <- in_out[[2]]

ListReduce1 <- in_out[[1]]
ListReduce2 <- in_out[[2]]


 if( all(CheckAllTrue) ){
        IsMarkovDAG <- all(CheckAllTrue)

        return(list(IsMarkovDAG=IsMarkovDAG,in_out = in_out, Reordered = FALSE ))

  } else if(!all(CheckAllTrue) ){
      # move in_out[[2]][[i]] and corresponding in_out[[1]][[i]] to index below its first occurrence in in_out[[1]].
      # HERE IT IS REQUIRED AS PER THE PACKAGE THAT THE in_out[[1]] and in_out[[2]] are listed in the same order otherwise this will not work.
      # And unable to check this since user can name variables what they want which is beyond checking.
#################
#################
#################

      while( length(ListReduce2) !=0 ){

        VarInLoop1 <- ListReduce1[[1]]
        VarInLoop2 <- ListReduce2[[1]]

        h = grep(VarInLoop2,List2 )
        # h = grep(paste('^',VarInLoop2,'$',sep=''),List2,perl=TRUE )

         # if(  length(which(lapply(List1, function(data_input) grep(VarInLoop2, data_input)) > 0 )) == 0 ){
        if(  length(which(lapply(List1, function(data_input) grep(VarInLoop2, data_input)) > 0 )) == 0 ){
        # Check if response variable then move to end of list if response variable
        # if(  length(which(lapply(List1, function(data_input) grep(paste('^',VarInLoop2,'$',sep=''),data_input,perl=TRUE)) > 0 )) == 0 ){
            # if length() == 0 implies it is a variable which is a not a parent of any variable e.g. response variable
            # Do nothing as let next if statement move other variables before this.

             # h = grep(VarInLoop2,List2 )

             List1 <- list.remove(List1, h)
             List2 <- list.remove(List2, h)

             List1 <- list.insert(List1, length(in_out[[1]]), VarInLoop1 )
             List2 <- list.insert(List2, length(in_out[[2]]), VarInLoop2 )

        # }else if(h > which(sapply(List1, function(data_input) VarInLoop2 %in% data_input)) ){
        }else if( any(h > which(lapply(List1, function(data_input) grep(VarInLoop2, data_input)) > 0 ) ) ){
        # }else if( any(h > which(lapply(List1, function(data_input) grep(paste('^',VarInLoop2,'$',sep=''),data_input,perl=TRUE )) > 0 ) ) ){

               ToMove1 <- List1[[h]]
               ToMove2 <- List2[[h]]

        IndexToInsertAt <- which(lapply(List1, function(data_input) grep(VarInLoop2, data_input)) > 0 )[1]
        # IndexToInsertAt <- which(lapply(List1, function(data_input) grep(paste('^',VarInLoop2,'$',sep=''),data_input,perl=TRUE )) > 0 )[1]

             List1 <- list.remove(List1, h)
             List2 <- list.remove(List2, h)

             List1 <- list.insert(List1, IndexToInsertAt, ToMove1 )
             List2 <- list.insert(List2, IndexToInsertAt, ToMove2 )

           }


        ListReduce1 <- list.remove(ListReduce1, 1)
        ListReduce2 <- list.remove(ListReduce2, 1)

      }



#################
#################
##################
      # for( h in 1:length(in_out_test[[2]]) ){
      #
      #
      #   # if( length( which(sapply(in_out_test[[1]], function(data_input) in_out_test[[2]][[h]] %in% data_input)) ) == 0 ){
      #   if(  length(which(lapply(in_out_test[[1]], function(data_input) grep(in_out_test[[2]][[h]], data_input)) > 0 )) == 0 ){
      #       # if length() == 0 implies it is a variable which is a not a parent of any variable e.g. response variable
      #       # Do nothing as let next if statement move other variables before this.
      #
      #        List1 <- list.remove(in_out_test[[1]], h)
      #        List2 <- list.remove(in_out_test[[2]], h)
      #
      #        # in_out_test[[1]] <- list.insert(in_out_test[[1]], IndexToInsertAt, ToMove1 )
      #        # in_out_test[[2]] <- list.insert(in_out_test[[2]], IndexToInsertAt, ToMove2 )
      #        List1 <- list.insert(List1, length(in_out_test[[1]]), in_out_test[[1]][[h]] )
      #        List2 <- list.insert(List2, length(in_out_test[[2]]), in_out_test[[2]][[h]] )
      #
      #        # StoreList1 <- List1
      #        # StoreList2 <- List2
      #
      #   # }else if(h > which(sapply(in_out_test[[1]], function(data_input) in_out_test[[2]][[h]] %in% data_input)) ){
      #   }else if( any(h > which(lapply(in_out_test[[1]], function(data_input) grep(in_out_test[[2]][[h]], data_input)) > 0 ) ) ){
      #
      #        # h > which(sapply(in_out[[1]], function(data_input) in_out[[2]][[10]] %in% data_input))
      #       # which(sapply(x, function(y) x %in% y))
      #
      #        # ToMove1 <- in_out_test[[1]][[h]]
      #        # ToMove2 <- in_out_test[[2]][[h]]
      #          ToMove1 <- List1[[h]]
      #          ToMove2 <- List2[[h]]
      #
      # # Insert before this index (but when remove first you need to insert at this index then)
      # # IndexToInsertAt <- which(sapply(in_out_test[[1]], function(data_input) in_out_test[[2]][[h]] %in% data_input))
      #   # IndexToInsertAt <- which(lapply(in_out_test[[1]], function(data_input) grep(in_out_test[[2]][[h]], data_input)) > 0 )[1]
      #   IndexToInsertAt <- which(lapply(List1, function(data_input) grep(List2[[h]], data_input)) > 0 )[1]
      #
      #
      #        # in_out_test[[1]] <- list.remove(in_out_test[[1]], h)
      #        # in_out_test[[2]] <- list.remove(in_out_test[[2]], h)
      #        # List1 <- list.remove(in_out_test[[1]], h)
      #        # List2 <- list.remove(in_out_test[[2]], h)
      #        List1 <- list.remove(List1, h)
      #        List2 <- list.remove(List2, h)
      #
      #
      #        # in_out_test[[1]] <- list.insert(in_out_test[[1]], IndexToInsertAt, ToMove1 )
      #        # in_out_test[[2]] <- list.insert(in_out_test[[2]], IndexToInsertAt, ToMove2 )
      #        List1 <- list.insert(List1, IndexToInsertAt, ToMove1 )
      #        List2 <- list.insert(List2, IndexToInsertAt, ToMove2 )
      #
      #      }
      # }


in_out_updated <- list(List1, List2)


Boolean <- vector(mode = "list", length = ( length(in_out_updated[[2]]) - 1) )

  for(i in (1:(length( in_out_updated[[2]] ) - 1)) ){
        m <- 1
        for( j in (i + 1):length( in_out_updated[[2]] )){

              Boolean[[i]][[m]] <- all(in_out_updated[[1]][[i]] %in% in_out_updated[[1]][[ j ]] )

              m = m + 1
        }
  }

CheckAllTrue <- vector(mode = "list", length = length(Boolean) )

  for( w in 1:length(CheckAllTrue)){
        CheckAllTrue[[w]] <- all(Boolean[[w]])
  }


##############################
##############################
##############################
                # if( all(CheckAllTrue) ){
                #         all(CheckAllTrue)
                # } else if(!all(CheckAllTrue) ){
                #       #  Model is not defined as per the documentation so function will not order the variables correctly so return FALSE
                #         all(CheckAllTrue)
                #       }
               # If statement above not needed since same value returned in both cases but commented out to show logic of code above
               IsMarkovDAG <- all(CheckAllTrue)

              return(list(IsMarkovDAG = IsMarkovDAG,in_out = in_out_updated, Reordered = TRUE))
##############################
##############################
##############################

  }





}
