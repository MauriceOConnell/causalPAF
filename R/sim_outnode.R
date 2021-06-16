#' @title Set a variable defined in col_num to its reference level (or otherwise 0 i.e. symbolising the absence of the risk factor) and then simulate all it's cildren.
#' @description Set variable defined in col_num to its reference level or otherwise to 0. Then simulate each variable if there is a direct arrow from the variable thats set to its reference value (or otherwise 0) to that variable i.e. if the variable set to its reference level (or otherwise 0) is a parent  of a variable simulate the child.
#' @param dataframe A wide format dataframe containing all the risk factors, confounders, exposures and outcomes within the causal DAG Bayesian network.
#' @param col_num Column number of variable in dataframe current_mat that is to be set to its reference level (if a factor) or to 0 otherwise i.e. the absence of the risk factor.
#' @param current_mat Data frame containing the data that is to be manipulated by the fucntion. The data frame has cases in rows and variables in columns.
#' @param in_outArg  A list of length 2. The first list contains a list of character vectors of the parents of the exposure or risk factor or outcome which are either causes or confounders of the exposure or risk factor or outcome. The second list contains a list of a single name of exposure or risk factor or outcome in form of characters. See tutorial examples for examples.
#' @param col_list is a list which gives the column numbers of each column in dataframe and current_mat that are in the same order as in_outArg  and model_list. See tutorial for an example.
#' @param model_list is a list of models fitted for each of the variables in in_out$outlist based on its parents given in in_out$inlist. By default this is set to an empty list. In the default setting, the models are fitted based on the order of the variables input in the parameter in_outArg. See the tutorial for more examples. Alternatively, the user can supply their own fitted models here by populating ``model_listArg'' with their own fitted models for each risk factor, mediator, exposure and response varialble. But the order of these models must be in the same order of the variables in the second list of in_outArg. See tutorial for further examples.
#' @export
#' @import stats
#' @keywords models Regression
#' @return \item{ current_mat }{ Returns the dataframe current_mat updated by setting the exposure or variable in column, col_num, of current_mat to its reference level (if a factor) or to 0 otherwise i.e. the absence of the risk factor. current_mat is then returned with all children of column, col_num, of current_mat simulated given col_num, of current_mat has been set to its reference level (if a factor) or to 0 otherwise i.e. the absence of the risk factor. In summary, current_mat is returned setting the variable defined in col_num to its reference level (or otherwise 0 i.e. symbolising the absence of the risk factor) and then simulates all it's cildren based on this reference value of the variable defined at col_cum of current_mat. This allows for indirect effects and direct effects.   }
#' @examples \dontrun{
#' # I don't want you to run this
#' }
sim_outnode <- function(dataframe, col_num, current_mat, in_outArg, col_list, model_list){

  ## do_sim NOT THIS VERSION
  ## do_sim_sequentialPAF BUT THIS VERSION

        ##  set current_mat[,col_num] to reference, otherwise to 0
        # if(is.factor(current_mat[,col_num])) current_mat[,col_num] <- levels(stroke_reduced[,col_num])[1]
        if(is.factor(current_mat[,col_num])) current_mat[,col_num] <- levels(dataframe[,col_num])[1]
        if(is.numeric(current_mat[,col_num])) current_mat[,col_num] <- 0

        colname <- colnames(current_mat)[col_num]

        #cols_to_simulate <- c()

        # simulate variable if direct arrow from colname to variable i.
        for(i in 1:length(in_outArg[[1]])){
                # MIGHT NOT WORK IF MULTIPLE COLUMNS WITH SAME NAME CONTAINED IN IT
                if(colname %in% in_outArg[[1]][[i]]){
                        if(length(table(current_mat[,col_list[[i]]] ))==1) next ##  don't alter variables that have already been changed

                        if(is.factor(current_mat[,col_list[i]])) current_mat[,col_list[i]] <- factor(do_sim_sequentialPAF(col_list[i],current_mat,model_list[[i]]),levels=levels(current_mat[,col_list[i]]))
                        if(!is.factor(current_mat[,col_list[i]])) current_mat[,col_list[i]] <- do_sim_sequentialPAF(col_list[i],current_mat, model_list[[i]])
                }
        }
        return(current_mat)
}
