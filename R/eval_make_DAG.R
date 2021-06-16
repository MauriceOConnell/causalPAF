#' @title Evaluates and Makes Formula for regression of exposure or risk factor or outcome on its parents in a causal Bayesian network directed acyclic graph.
#' @description Evaluates and Makes Formula for regression of exposure or risk factor or outcome on its parents in a causal Bayesian network directed acyclic graph. Given a causal Bayesian network, directed acyclic graph (DAG) where arrows representing
#' causal dependencies between confounders, risk factors/exposure and disease, together with a sensible probability distribution on
#' the graph that respects these causal dependencies. To consistently estimate causal effects that risk factors may have on each
#' other and on disease, we need to make a strong no unmeasured confounding assumption: that is common causes of nodes in the graph,
#' which may be causes of two risk factors or a cause of risk factor and disease, are also included as nodes in the graph.
#' Causal Bayesian networks have a local Markov property that the conditional probability distribution of any node Xj, given values
#' for the other variables in the network, only depends on the values $x_{pa}_{j}$ of the parent nodes.
#' @param data A wide format data containing all the risk factors, confounders, exposures and outcomes within the causal DAG Bayesian network.
#' @param regressionExposure Regression of respsonse given exposure based on canonical adjustment set output from function make_DAG.R.
#' @param regressionMediator Regression of respsonse given exposure (mediator as exposoure) based on canonical adjustment set output from function make_DAG.R.
#' @param response The name of the response column variable within dataframe in text format e.g. "case". The cases should be coded as 1 and the controls as 0.
#' @param response_model_mediators A model fitted for the response in a causal Bayesian network excluding ``children'' of the mediators in the causal Bayesian network. See example in tutorial.
#' @param response_model_exposure A model fitted for the response in a causal Bayesian network excluding ``children'' of the exposure and risk factors in the causal Bayesian network. See example in tutorial.
#' @param w Column of weights for case control matching listing in same order as patients in data.
#' @export
#' @import splines MASS stats forestplot utils grid magrittr checkmate
#' @keywords models Regression
#' @return \itemize{
#' \item{regressionExposure_listReturn }{model list regressionExposure_listReturn}
#' \item{regressionMediator_listReturn }{model list regressionMediator_listReturn}
#' }
#' @examples \dontrun{
#' # I don't want you to run this
#' }
#' in_vars = c("subeduc","moteduc","fatduc")
#' outvar = c("phys")
#' make_formula(in_vars,outvar,addCustom = FALSE, custom = "~ regionnn7 + ")
eval_make_DAG <- function(data, regressionExposure , regressionMediator, response, response_model_mediators, response_model_exposure, w){


#############
# FUNCTIONS USING OUTPUT RETURNED NEED TO CHECK IF DIMENSION IS DIFFERENT FROM EXPECTED I.E. E.G. MIGHT BE 1 DMENSIONAL AS OUTLINED BELOW.
# If using the models response_model_mediators, response_model_exposure populated by the user instead. Output could be 1 dimensional when more than one dimensional expected so look out for this in functions that are using output.
#############

#####################
## HAD ONLY STARTD THIS FUNCTION ON 30TH APRIL 2021 SO NEED TO START AGAIN FROM HERE
## NEED TO USE make_DAG.R output and evaluaate it here
## also nice in make_DAG.R to have plot look nice and have paretns on left of graph
#####################

# resultExposure, resultMediator

      countVAR <- 0
        # if( ( length(model_list) == 0 ) & (countVAR == 0) ){
        if( ( length(response_model_exposure) == 0 ) & (countVAR == 0) ){

          #test <- list()

          data_text <- deparse(substitute(data))

           # model_list_text <- deparse(substitute(model_list))

           regressionExposure_list_text <- deparse(substitute(regressionExposure))
           regressionMediator_list_text <- deparse(substitute(regressionMediator))


          w_text <- deparse(substitute(w))


          # model_listReturn <- vector(mode = "list", length = length(in_out[[2]]) )

          regressionExposure_listReturn <- vector(mode = "list", length = length(regressionExposure) )
          regressionMediator_listReturn <- vector(mode = "list", length = length(regressionMediator) )

          column <- (1:length(colnames(data)))[colnames(data) %in% response ]
          y <- data[,column]

          for(i in 1:length(regressionExposure) ){

                                # column <- (1:length(colnames(data)))[colnames(data) %in% in_out[[2]][i]]
                                # # formula_text <- make_formula(in_out[[1]][[i]],in_out[[2]][i], addCustom = TRUE, custom = "regionnn7*ns(eage,df=5)+esex*ns(eage,df=5) + ")
                                # formula_text <- make_formula(in_out[[1]][[i]],in_out[[2]][i], addCustom , custom )
                                formula_text <- regressionExposure[[i]]
                                # y <- data[,column]
                                if(length(table(y))==2){
                                        theform <- paste("glm(",formula_text,",data=",data_text,",family='binomial',w=",w_text,")",sep='')
                                }
                                if(length(table(y))>2 & is.factor(y)){
                                        theform <- paste("polr(",formula_text,",data=",data_text,",w=",w_text,")",sep='')
                                }
                                if(length(table(y))>2 & is.numeric(y)){
                                        theform <- paste("lm(",formula_text,",data=",data_text,",w=",w_text,")",sep='')
                                }
                                # to_execute <- paste(model_list_text,"[[i]] <-", theform,sep='')
                                # model_listReturn[[i]] <- to_execute
                                to_execute <- paste(regressionExposure_list_text,"[[i]] <-", theform,sep='')
                                regressionExposure_listReturn[[i]] <- to_execute

                  }
                   countVAR <- 1
                   # model_listReturn <- model_list
                   #return(regressionExposure_listReturn)
        } else if( ( length(response_model_exposure) != 0 ) & (countVAR == 0)  ) {

                   # exposure
                   #
                   # response
                   #
                   # mediator

                   # But could be only 1 dimensional now rather than a list as above
                   regressionExposure_listReturn <- response_model_exposure
                   #return(model_listReturn)
        } else{
              return("error in if statement in function eval_make_DAG()")
        }

##########################################################
##########################################################
      countVAR <- 0
      if( ( length(response_model_mediators) == 0 ) & (countVAR == 0) ){

                                for(i in 1:length(regressionMediator) ){

                                # column <- (1:length(colnames(data)))[colnames(data) %in% in_out[[2]][i]]
                                # # formula_text <- make_formula(in_out[[1]][[i]],in_out[[2]][i], addCustom = TRUE, custom = "regionnn7*ns(eage,df=5)+esex*ns(eage,df=5) + ")
                                # formula_text <- make_formula(in_out[[1]][[i]],in_out[[2]][i], addCustom , custom )
                                formula_text <- regressionMediator[[i]]
                                # y <- data[,column]
                                if(length(table(y))==2){
                                        theform <- paste("glm(",formula_text,",data=",data_text,",family='binomial',w=",w_text,")",sep='')
                                }
                                if(length(table(y))>2 & is.factor(y)){
                                        theform <- paste("polr(",formula_text,",data=",data_text,",w=",w_text,")",sep='')
                                }
                                if(length(table(y))>2 & is.numeric(y)){
                                        theform <- paste("lm(",formula_text,",data=",data_text,",w=",w_text,")",sep='')
                                }
                                # to_execute <- paste(model_list_text,"[[i]] <-", theform,sep='')
                                # model_listReturn[[i]] <- to_execute
                                to_execute <- paste(regressionMediator_list_text,"[[i]] <-", theform,sep='')
                                regressionMediator_listReturn[[i]] <- to_execute

                  }
                   countVAR <- 1
                   # model_listReturn <- model_list
                   # return(regressionMediator_listReturn)

                    # my_list_RegressionCausal <- list("regressionExposure_listReturn" = regressionExposure_listReturn,
                    #                                  "regressionMediator_listReturn" = regressionMediator_listReturn )
                    #
                    # return(my_list_RegressionCausal)

##########################################################
##########################################################

##########################################################
##########################################################
                  # for(i in 1:length(in_out[[2]]) ){
                  #
                  #               column <- (1:length(colnames(data)))[colnames(data) %in% in_out[[2]][i]]
                  #               # formula_text <- make_formula(in_out[[1]][[i]],in_out[[2]][i], addCustom = TRUE, custom = "~ regionnn7*ns(eage,df=5)+esex*ns(eage,df=5) + ")
                  #               formula_text <- make_formula(in_out[[1]][[i]],in_out[[2]][i], addCustom , custom )
                  #               y <- data[,column]
                  #               if(length(table(y))==2){
                  #                       theform <- paste("glm(",formula_text,",data=",data_text,",family='binomial',w=",w_text,")",sep='')
                  #               }
                  #               if(length(table(y))>2 & is.factor(y)){
                  #                       theform <- paste("polr(",formula_text,",data=",data_text,",w=",w_text,")",sep='')
                  #               }
                  #               if(length(table(y))>2 & is.numeric(y)){
                  #                       theform <- paste("lm(",formula_text,",data=",data_text,",w=",w_text,")",sep='')
                  #               }
                  #                to_execute <- paste(model_list_text,"[[i]] <-", theform,sep='')
                  #               model_listReturn[[i]] <- to_execute
                  #               # eval(parse(text=to_execute))
                  #               ## model_listReturn[[i]] <- eval(parse(text=to_execute))
                  #               #model_listReturn[[i]] <- eval( parse(text= paste(model_list_text,"[[i]]",sep='') ) )
                  #
                  #               ## E.G. TRY model$terms <- eval(model$call$formula)
                  #               ##model_listReturn[[i]]$terms <- eval(model_listReturn[[i]]$call$formula)
                  # }
                  #  countVAR <- 1
                  #  # model_listReturn <- model_list
                  #  return(model_listReturn)
        } else if( ( length(response_model_mediators) != 0 ) & (countVAR == 0)  ) {

                   # exposure
                   #
                   # response
                   #
                   # mediator

                   regressionMediator_listReturn <- response_model_mediators
                   # return(model_listReturn)
        } else{
              return("error in if statement in function eval_make_DAG()")
        }


      my_list_RegressionCausal <- list("regressionExposure_listReturn" = regressionExposure_listReturn,
                                       "regressionMediator_listReturn" = regressionMediator_listReturn )

      return(my_list_RegressionCausal)



}
