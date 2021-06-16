#' @title Makes Formula for regression of exposure or risk factor or outcome on its parents in a causal Bayesian network directed acyclic graph.
#' @description Given a causal Bayesian network, directed acyclic graph (DAG) where arrows representing
#' causal dependencies between confounders, risk factors/exposure and disease, together with a sensible probability distribution on
#' the graph that respects these causal dependencies. To consistently estimate causal effects that risk factors may have on each
#' other and on disease, we need to make a strong no unmeasured confounding assumption: that is common causes of nodes in the graph,
#' which may be causes of two risk factors or a cause of risk factor and disease, are also included as nodes in the graph.
#' Causal Bayesian networks have a local Markov property that the conditional probability distribution of any node Xj, given values
#' for the other variables in the network, only depends on the values $x_{pa}_{j}$ of the parent nodes.
#' @param in_vars a list of character vectors of the parents of the exposure or risk factor or outcome which are either causes or confounders of the exposure or risk factor or outcome
#' @param outvar list of a single name of exposure or risk factor or outcome in form of characters
#' @param addCustom is a logical which is TRUE if a custom regression is to be added in or FALSE otherwise. The customised part is defined in the parameter "custom".
#' @param custom custom must start and end as follows: "~  + " with custom values between ~ and +. Custom allows customised input into regresssion not stated already e.g. allow for interaction of variables not listed yet e.g. "~ regionnn7 x ns(eage,df=5)+esex x ns(eage,df=5) + "
#' @export
#' @keywords models Regression
#' @return \item{result }{result}
#' @examples \dontrun{
#' # I don't want you to run this
#' }
#' in_vars = c("subeduc","moteduc","fatduc")
#' outvar = c("phys")
#' make_formula(in_vars,outvar)
make_formula <- function(in_vars,outvar, addCustom = FALSE, custom = ""){

  # changed from this
  # custom = "~ regionnn7*ns(eage,df=5)+esex*ns(eage,df=5) + "

  # # Before update (which was working)
  # result <- paste(outvar,"~ regionnn7*ns(eage,df=5)+esex*ns(eage,df=5) + ",in_vars[1])
  #       if(length(in_vars)>=2){
  #
  #               for(i in 2:length(in_vars)){
  #
  #                       result <- paste(result,"+ ",in_vars[i],sep='')
  #
  #               }
  #       }
  #
  #       result

  # Need custom variable since e.g. "~ regionnn7*ns(eage,df=5)+esex*ns(eage,df=5) + " is specific to our example and allows interaction
  if(addCustom){
        # result <- paste(outvar,"~ regionnn7*ns(eage,df=5)+esex*ns(eage,df=5) + ",in_vars[1])
        # changed from this
        # result <- paste(outvar, custom, in_vars[1])
        result <- paste(outvar,"~ ", in_vars[1],"+ ", custom, sep='')
  }else{
        # result <- paste(outvar,"~ regionnn7*ns(eage,df=5)+esex*ns(eage,df=5) + ",in_vars[1])
        result <- paste(outvar,"~ ",in_vars[1])
  }
        if(length(in_vars)>=2){

                for(i in 2:length(in_vars)){

                        result <- paste(result,"+ ",in_vars[i],sep='')

                }
        }

        result

}


