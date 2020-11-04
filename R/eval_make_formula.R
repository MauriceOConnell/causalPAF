#' @title Evaluates and Makes Formula for regression of exposure or risk factor or outcome on its parents in a causal Bayesian network directed acyclic graph.
#' @description Evaluates and Makes Formula for regression of exposure or risk factor or outcome on its parents in a causal Bayesian network directed acyclic graph. Given a causal Bayesian network, directed acyclic graph (DAG) where arrows representing
#' causal dependencies between confounders, risk factors/exposure and disease, together with a sensible probability distribution on
#' the graph that respects these causal dependencies. To consistently estimate causal effects that risk factors may have on each
#' other and on disease, we need to make a strong no unmeasured confounding assumption: that is common causes of nodes in the graph,
#' which may be causes of two risk factors or a cause of risk factor and disease, are also included as nodes in the graph.
#' Causal Bayesian networks have a local Markov property that the conditional probability distribution of any node Xj, given values
#' for the other variables in the network, only depends on the values $x_{pa}_{j}$ of the parent nodes.
#' @param dataframe A wide format dataframe containing all the risk factors, confounders, exposures and outcomes within the causal DAG Bayesian network.
#' @param in_out A list of length 2. The first list contains the ``in_vars'' which is a list of character vectors of the parents of the exposure or risk factor or outcome which are either causes or confounders of the exposure or risk factor or outcome. The second list conttains a list of a single name of exposure or risk factor or outcome in form of characters.
#' @param model_list list of a single name of exposure or risk factor or outcome in form of characters.
#' @param w Column of weights for case control matching listing in same order as patients in dataframe.
#' @export
#' @import splines MASS stats dplyr
#' @keywords models Regression
#' @return \item{model_list }{model list}
#' @examples \dontrun{
#' # I don't want you to run this
#' }
#' in_vars = c("subeduc","moteduc","fatduc")
#' outvar = c("phys")
#' make_formula(in_vars,outvar)
eval_make_formula <- function(dataframe,in_out, model_list, w){

  if(length(model_list) == 0){

    dataframe_text <- deparse(substitute(dataframe))

    model_list_text <- deparse(substitute(model_list))

            for(i in 1:length(in_out[[2]])){

                          column <- (1:length(colnames(dataframe)))[colnames(dataframe) %in% in_out[[2]][i]]
                          formula_text <- make_formula(in_out[[1]][[i]],in_out[[2]][i])
                          y <- dataframe[,column]
                          if(length(table(y))==2){
                                  theform <- paste("glm(",formula_text,",data=",dataframe_text,",family='binomial',w=w)",sep='')
                          }
                          if(length(table(y))>2 & is.factor(y)){
                                  theform <- paste("polr(",formula_text,",data=",dataframe_text,",w=w)",sep='')
                          }
                          if(length(table(y))>2 & is.numeric(y)){
                                  theform <- paste("lm(",formula_text,",data=",dataframe_text,",w=w)",sep='')
                          }
                          to_execute <- paste(model_list_text,"[[i]] <-", theform,sep='')
                          eval(parse(text=to_execute))
            }
            model_list

  } else{
            model_list
        }

}
