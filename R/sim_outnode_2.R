#' @title Set a variable defined in col_num to its reference level (or otherwise 0 i.e. symbolising the absence of the risk factor) and then simulate all it's children.
#' @description Set variable defined in col_num to its reference level or otherwise to 0. Then simulate each variable if there is a direct arrow from the variable thats set to its reference value (or otherwise 0) to that variable i.e. if the variable set to its reference level (or otherwise 0) is a parent  of a variable simulate the child.
#' @param dataframe  A wide format dataframe containing all the risk factors, confounders, exposures and outcomes within the causal DAG Bayesian network.
#' @param col_num Column number of variable in dataframe current_mat that is to be set to its reference level (if a factor) or to 0 otherwise i.e. the absence of the risk factor.
#' @param current_mat Data frame containing the data. The data frame has cases in rows and variables in columns.
#' @param col_list is a list which gives the column numbers of each column in dataframe and current_mat that are in the same order as in_outArg  and model_list. See tutorial for an example.
#' @param model_list is a list of models fitted for each of the variables in in_outArg[[2]] (or in_outArg\eqn{\$}outlist ) based on its parents given in in_outArg[[1]] ( or in_out\eqn{\$}inlist ). By default this is set to an empty list. In the default setting, the models are fitted automatically by the 'causalPAF' package based on the order of the variables input in the parameter in_outArg. See the tutorial for more examples. Alternatively, the user can supply their own fitted models here by populating ``model_listArg'' with their own fitted models for each risk factor, mediator, exposure and response variable. But the order of these models must be in the same order of the variables in the second list of in_outArg ( in_outArg[[2]] ) and these models be defined within a list, list(), of the same length as in_outArg[[2]]. See tutorial for further examples.
#' @param response_col_num The column of current_mat and dataframe that defines the response variable.
#' @param in_outArg This defines the causal directed acyclic graph (DAG). A list of length 2. It is defined as a two dimensional list consisting of, firstly, the first list, inlist, i.e. a list of the parents of each variable of interest corresponding to its column name in the data. Splines can be included here if they are to be modelled as splines. Secondly, the second list, outlist, contains a list of a single name of exposure or risk factor or outcome in form of characters i.e. a list of each variable of interest (risk factors, exposures and outcome) corresponding to its column name in the data. Splines should not be input here, only the column names of the variables of interest in the data. The order at which variables are defined must satisfy (i) It is important that variables are defined in the same order in both lists e.g. the first risk factor defined in outlist has its parents listed first in inlist, the second risk factor defined in outlist has its parents listed secondly in inlist and so on. The package assumes this ordering and will not work if this order is violated. (ii) Note it is important also that the order at which the variables are defined is such that all parents of that variable are defined before it. See example in tutorial.
#' @export
#' @import stats
#' @keywords internal
#' @return \item{ current_mat }{ Returns the dataframe current_mat updated by setting the exposure or variable in column, col_num, of current_mat to its reference level (if a factor) or to 0 otherwise i.e. the absence of the risk factor. current_mat is then returned with only the response model updated given col_num, of current_mat has been set to its reference level (if a factor) or to 0 otherwise i.e. the absence of the risk factor. In summary, current_mat is returned setting the variable defined in col_num to its reference level (or otherwise 0 i.e. symbolising the absence of the risk factor) and then simulates only the response model based on this reference value of the variable defined at col_cum of current_mat. This allows for only direct effects. There are no indirect effects simulated. }

sim_outnode_2 <- function(dataframe, col_num, current_mat, col_list, model_list, response_col_num, in_outArg ){

        ##  set current_mat[,col_num] to reference, otherwise to 0
        if(is.factor(current_mat[,col_num])) current_mat[,col_num] <- levels(dataframe[,col_num])[1]
        if(is.numeric(current_mat[,col_num])) current_mat[,col_num] <- 0

        colname <- colnames(current_mat)[col_num]

        #cols_to_simulate <- c()

        # response_col_num <- (1:length(col_list))[lapply(in_outArg[[2]], function(data_input) data_input == response ) == TRUE ]
        # simulate variable if direct arrow from colname to variable i (just last column) or 11th column
        # i <- 11
        # assumes response is defined in the last column which may not be the case
        # i <- length(col_list)
        i <- response_col_num


        if(is.factor(current_mat[,col_list[i]])) current_mat[,col_list[i]] <- factor(do_sim_sequentialPAF(col_list[i],current_mat,model_list[[i]]),levels=levels(current_mat[,col_list[i]]))
        if(!is.factor(current_mat[,col_list[i]])) current_mat[,col_list[i]] <- do_sim_sequentialPAF(col_list[i],current_mat,model_list[[i]])

        return(current_mat)
}




