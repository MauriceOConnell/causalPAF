#' @title Evaluates and Makes Formula for regression of exposure or risk factor or outcome on its parents in a causal Bayesian network directed acyclic graph.
#' @description Evaluates and Makes Formula for regression of exposure or risk factor or outcome on its parents in a causal Bayesian network directed acyclic graph. Given a causal Bayesian network, directed acyclic graph (DAG) where arrows representing
#' causal dependencies between confounders, risk factors, exposure and disease, together with a sensible probability distribution on
#' the graph that respects these causal dependencies. To consistently estimate causal effects that risk factors may have on each
#' other and on disease, we need to make a strong no unmeasured confounding assumption: that is common causes of nodes in the graph,
#' which may be causes of two risk factors or a cause of risk factor and disease, are also included as nodes in the graph.
#' Causal Bayesian networks have a local Markov property that the conditional probability distribution of any node \eqn{X_j}, given values
#' for the other variables in the network, only depends on the values \eqn{x_{pa_{j}} } of the parent nodes.
#' @param data A wide format data containing all the risk factors, confounders, exposures and outcomes within the causal DAG Bayesian network.
#' @param in_out This defines the causal directed acyclic graph (DAG). A list of length \eqn{2}. It is defined as a two dimensional list consisting of, firstly, the first list, inlist, i.e. a list of the parents of each variable of interest corresponding to its column name in the data. Splines can be included here if they are to be modelled as splines. Secondly, the second list, outlist, contains a list of a single name of exposure or risk factor or outcome in form of characters i.e. a list of each variable of interest (risk factors, exposures and outcome) corresponding to its column name in the data. Splines should not be input here, only the column names of the variables of interest in the data. The order at which variables are defined must satisfy (i) It is important that variables are defined in the same order in both lists e.g. the first risk factor defined in outlist has its parents listed first in inlist, the second risk factor defined in outlist has its parents listed secondly in inlist and so on. The package assumes this ordering and will not work if this order is violated. (ii) Note it is important also that the order at which the variables are defined is such that all parents of that variable are defined before it. See example in tutorial.
#' @param model_list is a list of models fitted for each of the variables in in_outArg[[2]] (or in_outArg\eqn{\$}outlist ) based on its parents given in in_outArg[[1]] ( or in_out\eqn{\$}inlist ). By default this is set to an empty list. In the default setting, the models are fitted automatically by the causalPAF package based on the order of the variables input in the parameter in_outArg. See the tutorial for more examples. Alternatively, the user can supply their own fitted models here by populating ``model_listArg'' with their own fitted models for each risk factor, mediator, exposure and response variable. But the order of these models must be in the same order of the variables in the second list of in_outArg ( in_outArg[[2]] ) and these models be defined within a list, list(), of the same length as in_outArg[[2]]. See tutorial for further examples.
#' @param w Column of weights for case control matching listing in same order as patients in data.
#' @param addCustom Logical TRUE or FALSE indicating whether a customised interaction term is to be added to the each regression. The interaction term can include splines.
#' @param custom text containing the customised interaction term to be added to each regression. The text should be enclosed in inverted commas. Splines can be included within the interaction terms. See tutorial for examples.
#' @export
#' @import splines MASS stats forestplot utils grid magrittr checkmate
#' @keywords internal
#' @return \itemize{
#' \item{model_listReturn[[1]] }{model list A}
#' \item{model_listReturn[[2]] }{model list B}
#' \item{model_listReturn[[3]] }{model list C}
#' \item{model_listReturn[[4]] }{model list D}
#' \item{model_listReturn[[5]] }{model list E}
#' \item{model_listReturn[[6]] }{model list F}
#' \item{model_listReturn[[7]] }{model list G}
#' \item{model_listReturn[[8]] }{model list H}
#' \item{model_listReturn[[9]] }{model list I}
#' \item{model_listReturn[[10]] }{model list J}
#' \item{model_listReturn[[11]] }{model list K}
#' }

eval_make_formula <- function(data,in_out, model_list, w, addCustom = FALSE, custom = ""){

  count <- 0
  if( ( length(model_list) == 0 ) & (count == 0) ){

    data_text <- deparse(substitute(data))

     model_list_text <- deparse(substitute(model_list))
    #model_list_text <- deparse(substitute(test))

    w_text <- deparse(substitute(w))


    # model_listReturn <- list()
    model_listReturn <- vector(mode = "list", length = length(in_out[[2]]) )

            for(i in 1:length(in_out[[2]]) ){

                          column <- (1:length(colnames(data)))[colnames(data) %in% in_out[[2]][i]]
                          # formula_text <- make_formula(in_out[[1]][[i]],in_out[[2]][i], addCustom = TRUE, custom = "~ regionnn7*ns(eage,df=5)+esex*ns(eage,df=5) + ")
                          formula_text <- make_formula(in_out[[1]][[i]],in_out[[2]][i], addCustom , custom )
                          y <- data[,column]
                          if(length(table(y))==2){
                                  theform <- paste("glm(",formula_text,",data=",data_text,",family='binomial',w=",w_text,")",sep='')
                          }
                          if(length(table(y))>2 & is.factor(y)){
                                  theform <- paste("polr(",formula_text,",data=",data_text,",w=",w_text,")",sep='')
                          }
                          if(length(table(y))>2 & is.numeric(y)){
                                  theform <- paste("lm(",formula_text,",data=",data_text,",w=",w_text,")",sep='')
                          }
                           to_execute <- paste(model_list_text,"[[i]] <-", theform,sep='')
                          model_listReturn[[i]] <- to_execute
                          # eval(parse(text=to_execute))
                          ## model_listReturn[[i]] <- eval(parse(text=to_execute))
                          #model_listReturn[[i]] <- eval( parse(text= paste(model_list_text,"[[i]]",sep='') ) )

                          ## E.G. TRY model$terms <- eval(model$call$formula)
                          ##model_listReturn[[i]]$terms <- eval(model_listReturn[[i]]$call$formula)
            }
             count <- 1
             # model_listReturn <- model_list
             return(model_listReturn)
  } else if( ( length(model_list) != 0 ) & (count == 0)  ) {
             model_listReturn <- model_list
             return(model_listReturn)
  } else{
        return("error in if statement in function eval_make_formula()")
    }

}
