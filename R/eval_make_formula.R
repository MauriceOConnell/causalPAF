#' @title Evaluates and Makes Formula for regression of exposure or risk factor or outcome on its parents in a causal Bayesian network directed acyclic graph.
#' @description Evaluates and Makes Formula for regression of exposure or risk factor or outcome on its parents in a causal Bayesian network directed acyclic graph. Given a causal Bayesian network, directed acyclic graph (DAG) where arrows representing
#' causal dependencies between confounders, risk factors/exposure and disease, together with a sensible probability distribution on
#' the graph that respects these causal dependencies. To consistently estimate causal effects that risk factors may have on each
#' other and on disease, we need to make a strong no unmeasured confounding assumption: that is common causes of nodes in the graph,
#' which may be causes of two risk factors or a cause of risk factor and disease, are also included as nodes in the graph.
#' Causal Bayesian networks have a local Markov property that the conditional probability distribution of any node Xj, given values
#' for the other variables in the network, only depends on the values $x_{pa}_{j}$ of the parent nodes.
#' @param data A wide format data containing all the risk factors, confounders, exposures and outcomes within the causal DAG Bayesian network.
#' @param in_out A list of length 2. The first list contains the ``in_vars'' which is a list of character vectors of the parents of the exposure or risk factor or outcome which are either causes or confounders of the exposure or risk factor or outcome. The second list conttains a list of a single name of exposure or risk factor or outcome in form of characters.
#' @param model_list list of a single name of exposure or risk factor or outcome in form of characters.
#' @param w Column of weights for case control matching listing in same order as patients in data.
#' @export
#' @import splines MASS stats dplyr forestplot utils grid magrittr checkmate
#' @keywords models Regression
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
#' @examples \dontrun{
#' # I don't want you to run this
#' }
#' in_vars = c("subeduc","moteduc","fatduc")
#' outvar = c("phys")
#' make_formula(in_vars,outvar)
eval_make_formula <- function(data,in_out, model_list, w){

  count <- 0
  if( ( length(model_list) == 0 ) & (count == 0) ){

    #test <- list()

    data_text <- deparse(substitute(data))

     model_list_text <- deparse(substitute(model_list))
    #model_list_text <- deparse(substitute(test))

    w_text <- deparse(substitute(w))


    # model_listReturn <- list()
    model_listReturn <- vector(mode = "list", length = length(in_out[[2]]) )

            for(i in 1:length(in_out[[2]]) ){

                          column <- (1:length(colnames(data)))[colnames(data) %in% in_out[[2]][i]]
                          formula_text <- make_formula(in_out[[1]][[i]],in_out[[2]][i])
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
