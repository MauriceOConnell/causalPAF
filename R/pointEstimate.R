#' @title Evaluates Point Estimates for Total PAF, Direct PAF, Indirect PAF and Path Specific PAF for a user inputted number of integral simulations. There is no bootstap applied in this fucntion.
#' @description Evaluates Total PAF, Direct PAF, Indirect PAF and Path Specific PAF for a user inputted number of bootstraps and integral simulations
#' @param dataframe A wide format dataframe containing all the risk factors, confounders, exposures and outcomes within the causal DAG Bayesian network.
#' @param exposure The name of the exposure column variable within dataframe in text format e.g. "phys".
#' @param mediator The name of the mediator column variables within dataframe in text format. There can be more than one mediator of interest. It can be a vector of mediators names within the dataframe e.g. c("subhtn","apob_apoa","whr").
#' @param response The name of the response column variable within dataframe in text format e.g. "case". The cases should be coded as 1 and the controls as 0.
#' @param response_model_mediators A model fitted for the response in a causal Bayesian network excluding ``children'' of the mediators in the causal Bayesian network. See example in tutorial.
#' @param response_model_exposure A model fitted for the response in a causal Bayesian network excluding ``children'' of the exposure and risk factors in the causal Bayesian network. See example in tutorial.
#' @param in_outArg A list of length 2. The first list contains a list of character vectors of the parents of the exposure or risk factor or outcome which are either causes or confounders of the exposure or risk factor or outcome. The second list contains a list of a single name of exposure or risk factor or outcome in form of characters. See tutorial examples for examples.
#' @param Splines_outlist A list defined of same size and order of variables as defined in in_outArg[[2]]. If splines are to be used for variables listed in in_outArg[[2]], then the splines should be defined in the same order as variables appear in in_outArg[[2]]. It is necessary to list variables in in_outArg[[2]] without splines if no spline is to be applied.
#' @param splinesDefinedIn_in_outDAG Logical TRUE or FALSE indicating whether the user has defined splines in the causal DAG, in_out, if TRUE. If FALSE and splines are defined in Splines_outlist_Var, then it is necessary for the package to populate the in_out DAG with splines listed in Splines_outlist_Var.
#' @param model_listArg is a list of models fitted for each of the variables in in_out$outlist based on its parents given in in_out$inlist. By default this is set to an empty list. In the default setting, the models are fitted based on the order of the variables input in the parameter in_outArg. See the tutorial for more examples. Alternatively, the user can supply their own fitted models here by populating ``model_listArg'' with their own fitted models for each risk factor, mediator, exposure and response varialble. But the order of these models must be in the same order of the variables in the second list of in_outArg. See tutorial for further examples.
#' @param weights Column of weights for case control matching listed in the same order as the patients in the data e.g. weights = strokedata$weights.
#' @param NumSimulation This is the number of simulatons requested by the user to estimate integrals. The larger the number of simulations the more accurate the results but the longer the code takes to run. Therefore the user may wish to balance speed with accuracy depending on which is of more value in the specific context of interest. The integrals for continuous variables are estimated using simulation methods.
#' @param addCustom Logical TRUE or FALSE indicating whether a customised interaction term is to be added to the each regression. The interaction term can include splines.
#' @param custom text containing the customised interaction term to be added to each regression. The text should be enclosed in inverted commas. Splines can be included within the interactin terms. See tutorial for examples.
#' @export
#' @import splines MASS stats utils
#' @keywords models Regression Population Attributable Fraction
#' @return Estimates point estimates for 5 results that are:(1)Total Population Attributable Fraction (PAF),(2)Direct Effect Population Attributable Fraction (PAF) using  Sjolanders definition, (3)Indirect Effect Population Attributable Fraction (PAF) using  Sjolanders definition, (4)Path Specific Population Attributable Fraction (PAF), (5)Overall Direct Population Attributable Fraction (PAF)
#' @examples \dontrun{
#' # I don't want you to run this
#' }
#' in_vars = c("subeduc","moteduc","fatduc")
#' outvar = c("phys")
#' make_formula(in_vars,outvar)
pointEstimate <- function(dataframe,
                          exposure="phys",
                          mediator=c("subhtn","apob_apoa","whr"),
                          response="case",
                          # response_model_mediators=response_vs_mediator,
                          # response_model_exposure=response_vs_phys,
                          response_model_mediators = list(),
                          response_model_exposure = list(),
                          in_outArg,
                          Splines_outlist,
                          splinesDefinedIn_in_outDAG,
                          model_listArg,
                          # weights = stroke_reduced$weights,
                          # weights = strokedata$weights,
                          weights = 1,  ## At moment, weigths are populated in the dataframe before calling function causalPAFplot.r. The default is weights = 1 which is not a case-control data format but rather a cohort weighting by default.
                          # prevalence ADD IN
                          # CHECK IS WEIGHT HARD CODED
                          NumSimulation,
                          addCustom = FALSE,
                          custom = ""
                          ){

  # ggproto <- ggplot2::ggproto

  ##################################
  ## Test run entries
  ##################################
                          # dataframe = stroke_reduced
                          # exposure="phys"
                          # mediator=c("subhtn","apob_apoa","whr")
                          # response="case"
                          # response_model_mediators = list()
                          # response_model_exposure = list()
                          # in_outArg = in_out
                          # Splines_outlist = Splines_outlist  # needs to be input as list() if no splines. Assumes if appears as spline once needs to appear as a spline in all occurences
                          # splinesDefinedIn_in_outDAG = TRUE
                          # model_listArg = list()
                          # weights = stroke_reduced$weights
                          # NumSimulation = 2
                          # addCustom = TRUE
                          # custom = "regionnn7*ns(eage,df=5)+esex*ns(eage,df=5)" # need to update make_formula() so that it runs without needing ~ first
                          # # custom = "+ regionnn7*ns(eage,df=5)+esex*ns(eage,df=5) "


# TestpointEstimate  <- pointEstimate(dataframe = stroke_reduced,
#                                     exposure="phys",
#                                     mediator=c("subhtn","apob_apoa","whr"),
#                                     response="case",
#                                     response_model_mediators = list(),
#                                     response_model_exposure = list(),
#                                     in_outArg = in_out,
#                                     Splines_outlist = Splines_outlist,  # needs to be input as list() if no splines. Assumes if appears as spline once needs to appear as a spline in all occurences
#                                     splinesDefinedIn_in_outDAG = TRUE,
#                                     model_listArg = list(),
#                                     weights = stroke_reduced$weights,
#                                     NumSimulation = 1000,
#                                     addCustom = TRUE,
#                                     custom = "regionnn7*ns(eage,df=5)+esex*ns(eage,df=5)" # need to update make_formula() so that it runs without needing ~ first
#                                     # custom = "+ regionnn7*ns(eage,df=5)+esex*ns(eage,df=5) "
#                                     )


  # TestpointEstimate$results_mediatorPointEstimate
  # TestpointEstimate$mediators
  # TestpointEstimate$response_vs_physPointEstimate
  # TestpointEstimate$response_vs_mediatorPointEstimate

# #####################
# #####################
#  NumSimulation = 1000
# #####################
# #####################

####################
####################
## ADD IN POINT ESTIMATE WITHOUT BOOTSTRAP
##
####################
####################
# library(splines)


results_subhtn <- matrix(nrow= 1, ncol=5)
colnames(results_subhtn) <- c("overall","direct Sjolander","indirect Sjolander","path specific","overall Direct")


results_whr <- matrix(nrow= 1, ncol=5)
colnames(results_whr) <- c("overall","direct Sjolander","indirect Sjolander","path specific","overall Direct")


results_apob_apoa <- matrix(nrow= 1, ncol=5 )
colnames(results_apob_apoa) <- c("overall","direct Sjolander","indirect Sjolander","path specific","overall Direct")


#####################
#######################
#######################
# ADDED IN FROM causalPAFplot
results_mediatorPointEstimate <- list()
for(i in 1:length(mediator)){
  # results_mediatorPointEstimate[[i]] <- matrix(nrow= dataframe,ncol=5)
  results_mediatorPointEstimate[[i]] <- matrix(nrow = 1,ncol=5)
  colnames(results_mediatorPointEstimate[[i]]) <- c("overall","direct Sjolander","indirect Sjolander","path specific","overall Direct")
}

model_list_usePointEstimate <- vector(mode = "list", length = length(in_outArg[[2]]) )
model_list_evalPointEstimate <- vector(mode = "list", length = length(in_outArg[[2]]) )

##########################################################################
##########################################################################
# Perform check and dagitty check
##########################################################################
##########################################################################
## ADD IN make_DAG() here to perform check and dagitty on in_outArg
# MOVE make_DAG() below up here
# This works since in_outArg is used to define the DAG based on the original input of the user, but is then updated based on the check in the following if statement.
make_DAG_output <- make_DAG(in_outDAG = in_outArg ,
                            exposure = exposure,
                            response = response,
                            mediator = mediator,
                            Splines_outlist_Var = Splines_outlist,
                            # splinesDefinedIn_in_outDAG = TRUE,
                            splinesDefinedIn_in_outDAG = splinesDefinedIn_in_outDAG,
                            addCustomExposureAdjustmentSet = addCustom,
                            customExposureAdjustmentSet = custom,
                            addCustomMediatorAdjustmentSet = addCustom,
                            customMediatorAdjustmentSet = custom )   # custom = "regionnn7*ns(eage,df=5)+esex*ns(eage,df=5)"

## If is NULL, it means there was no update applied and no change is required.
## So if not Null, then in_outArg needs to be updated based on (1) dagitty check and (2) check required based on derivation.
if( !is.null(make_DAG_output$my_list_causal$in_outDAG_updatedAfterCheck) ){

        ##in_outArg is a parameter of the function but what to update it based on make_DAG (1) check requested by John Ferguson and (2) Dagitty check
        # MAYBE WE DO NOT WANT TO CHANGE THE INPUT PARAMETER in_outArg and rather we want to change
        # eval_make_formula(...., in_out = make_DAG_output$my_list_causal$in_outDAG_updatedAfterCheck, ....)
        ## check whether changing the parameter input in_outArg upsets anything else within the function?
         in_outArg <- make_DAG_output$my_list_causal$in_outDAG_updatedAfterCheck
}

###########################################################################
###########################################################################
###########################################################################

# If user has not fitted their own model_listArg then fit as follows
if( length(model_listArg) == 0 ){
      model_list_inputPointEstimate <- list()

      # model_list_usePointEstimate <- eval_make_formula(data = dataframe, in_out=in_outArg,model_list=model_list_inputPointEstimate, w=dataframe$weights, addCustom, custom)
      model_list_usePointEstimate <- eval_make_formula(data = dataframe, in_out = in_outArg, model_list=model_list_inputPointEstimate, w = weights, addCustom, custom)
      #   # TRY TO GET IN FUNCTION ENVIRONMENT E.G. model$terms <- eval(model$call$formula)
    #   # Return these from function and evaluate in this evironment!! Hopefully solve it.

      ## N.B. FOR LOOP HAS TO BE i AS MODELS DEFINED IN [[I]]
      for( i in 1:length(model_list_usePointEstimate)){
        eval(parse(text=model_list_usePointEstimate[[i]]))
         model_list_evalPointEstimate[[i]] <- model_list_inputPointEstimate[[i]]
      }
      # model_list_usePointEstimate <- eval(parse(text=model_list_usePointEstimate))
      # to_execute <- paste(model_list_text,"[[i]] <-", theform,sep='')
      # eval(parse(text=to_execute))
} else{
       model_list_evalPointEstimate <- model_listArg
}

# Select models for mediators input by user from model_listArg
  mediator_orderPointEstimate <- c()
  for(i in 1:length(model_list_evalPointEstimate)){
    mediator_orderPointEstimate[i] <-  as.character(formula(model_list_evalPointEstimate[[i]]))[2]
  }
  model_list_mediator_orderPointEstimate <- (1:length( mediator_orderPointEstimate ))[ mediator_orderPointEstimate %in% mediator]

  mediator_modelPointEstimate <- list()
  for(i in 1:length(model_list_mediator_orderPointEstimate)){
    mediator_modelPointEstimate[[i]] = model_list_evalPointEstimate[[ model_list_mediator_orderPointEstimate[i] ]]
  }

### NO NEED TO UPDATE AS SHOULD BE RUN ON ORIGINAL DATA AS IN causalPAFplot
# response_vs_mediator <- update(response_model_mediators, data = Bootstrap, weights = Bootstrap$weights )
#
# response_vs_phys <- update(response_model_exposure, data = Bootstrap, weights = Bootstrap$weights )


#######################
#######################
#######################

# model_list <- list()
# w <- rep(1,nrow(stroke_reduced))
# w[stroke_reduced$case==0] <- 0.9965
# w[stroke_reduced$case==1] <- 0.0035

#######################
#######################
# weights is input into function instead as some users may prefer a more detailed weighting calculation. So instead this is prepopulated as input into the functions.
#######################
# weightControl = 0.9965 # POPULATE THESE FROM FUNCTION
# weightCase = 0.0035  # POPULATE THESE FROM FUNCTION
# w <- rep(1,nrow(dataframe))
#
# column_response <- (1:length(colnames(dataframe)))[colnames(dataframe) == response ]
# w[ dataframe[,column_response] == 0 ] = weightControl
# w[ dataframe[,column_response] == 1 ] = weightCase
# dataframe <- mutate(dataframe, weigths_use = w)
#######################
#######################
#######################

# # response
# column <- (1:length(colnames(data)))[colnames(data) %in% response ]
#           y <- data[,column]
# deparse(substitute(response))
# w <- rep(1,nrow(dataframe))
# w[dataframe$case==0] <- 0.9965
# w[dataframe$case==1] <- 0.0035


# ##################################
# ##################################
# # If user has not fitted their own model_listArg then fit as follows
# if( length(model_listArg) == 0 ){
#       model_list_input <- list()
#       # model_list_use <- eval_make_formula(data = Bootstrap, in_out=in_outArg,model_list=model_listArg, w=Bootstrap$weights)
#       model_list_use <- eval_make_formula(data = dataframe, in_out=in_outArg,model_list=model_list_input, w=dataframe$weights, addCustom, custom)
#       #   # TRY TO GET IN FUNCTION ENVIRONMENT E.G. model$terms <- eval(model$call$formula)
#     #   # model$terms <- eval(model$call$formula)
#     #   model_list_use[[1]]$terms <- eval(model_list_use[[1]]$call$formula)
#     #   eval(model_list_use[[1]]$call$data )
#     #   eval(model_list_use[[1]]$call$family )
#     #   eval(model_list_use[[1]]$call$weights )
#     # #   glm(formula = phys ~ regionnn7 * ns(eage, df = 5) + esex * ns(eage,
#     # # df = 5) + subeduc + moteduc + fatduc, family = "binomial",
#     # # data = Bootstrap, weights = Bootstrap$weights)
#     #   glm(formula = model_list_use[[1]]$call$formula, family = model_list_use[[1]]$call$family,
#     # data = model_list_use[[1]]$call$data, weights = model_list_use[[1]]$call$weights )
#     #   # Return these from function and evaluate in this evironment!! Hopefully solve it.
#       # eval(parse(text=to_execute))
#                           ## model_listReturn[[i]] <- eval(parse(text=to_execute))
#                           #model_listReturn[[i]] <- eval( parse(text= paste(model_list_text,"[[i]]",sep='') ) )
#
#       ## N.B. FOR LOOP HAS TO BE i AS MODELS DEFINED IN [[I]]
#       for( i in 1:length(model_list_use)){
#         eval(parse(text=model_list_use[[i]]))
#          model_list_eval[[i]] <- model_list_input[[i]]
#       }
#       # model_list_use <- eval(parse(text=model_list_use))
#       # to_execute <- paste(model_list_text,"[[i]] <-", theform,sep='')
#       # eval(parse(text=to_execute))
# } else{
#        model_list_eval <- model_listArg
# }
#
# # Select models for mediators input by user from model_listArg
#   mediator_order <- c()
#   for(i in 1:length(model_list_eval)){
#     mediator_order[i] <-  as.character(formula(model_list_eval[[i]]))[2]
#   }
#   model_list_mediator_order <- (1:length( mediator_order ))[ mediator_order %in% mediator]
#
#   mediator_model <- list()
#   for(i in 1:length(model_list_mediator_order)){
#     mediator_model[[i]] = model_list_eval[[ model_list_mediator_order[i] ]]
#   }
#
# ###################################
# ###################################


################
## Maybe fit from make_DAG.R both below based on similar code to for loop below
##  resultExposure
##  resultMediator
################
################
### NEED TO RUN HERE WITH custom = "regionnn7*ns(eage,df=5)+esex*ns(eage,df=5)" RATHER THAN custom = "~ regionnn7*ns(eage,df=5)+esex*ns(eage,df=5) + "
# NEED TO UPDATE CUSTOM AS CURRENTLY make_formula() requires ~ but should change so runs without ~ needed
# custom = "~ regionnn7*ns(eage,df=5)+esex*ns(eage,df=5) + "
# custom = "regionnn7*ns(eage,df=5)+esex*ns(eage,df=5)"
# ############
# ## Moved make_DAG_output <- make_DAG(...) call up as needs to be used earlier
# ############
# make_DAG_output <- make_DAG(in_outDAG = in_outArg ,
#                             exposure = exposure,
#                             response = response,
#                             mediator = mediator,
#                             Splines_outlist_Var = Splines_outlist,
#                             splinesDefinedIn_in_outDAG = TRUE,
#                             addCustomExposureAdjustmentSet = addCustom,
#                             customExposureAdjustmentSet = custom,
#                             addCustomMediatorAdjustmentSet = addCustom,
#                             customMediatorAdjustmentSet = custom )   # custom = "regionnn7*ns(eage,df=5)+esex*ns(eage,df=5)"


# make_DAG_output_resultExposure <- make_DAG(in_outDAG = in_outArg ,
#                                            exposure = exposure,
#                                            response = response,
#                                            mediator = mediator,
#                                            Splines_outlist_Var = Splines_outlist,
#                                            splinesDefinedIn_in_outDAG = TRUE,
#                                            addCustomExposureAdjustmentSet = addCustom,
#                                            customExposureAdjustmentSet = custom,
#                                            addCustomMediatorAdjustmentSet = addCustom,
#                                            customMediatorAdjustmentSet = custom )$resultExposure   # custom = "regionnn7*ns(eage,df=5)+esex*ns(eage,df=5)"
#
# make_DAG_output_resultMediator <- make_DAG(in_outDAG = in_outArg ,
#                                            exposure = exposure,
#                                            response = response,
#                                            mediator = mediator,
#                                            Splines_outlist_Var = Splines_outlist,
#                                            splinesDefinedIn_in_outDAG = TRUE,
#                                            addCustomExposureAdjustmentSet = addCustom,
#                                            customExposureAdjustmentSet = custom,
#                                            addCustomMediatorAdjustmentSet = addCustom,
#                                            customMediatorAdjustmentSet = custom )$resultMediator   # custom = "regionnn7*ns(eage,df=5)+esex*ns(eage,df=5)"


# test$resultExposure
#
# test$resultMediator

eval_make_DAG_output <- eval_make_DAG(data = dataframe,
                                      regressionExposure = make_DAG_output$resultExposure,
                                      regressionMediator = make_DAG_output$resultMediator,
                                      response = response,
                                      response_model_mediators = response_model_mediators,
                                      response_model_exposure = response_model_exposure,
                                      w = weights )

# eval_make_DAG_output_regressionExposure_listReturn <- eval_make_DAG(data = dataframe,
#                                                                     regressionExposure = make_DAG_output_resultExposure,
#                                                                     regressionMediator = make_DAG_output_resultMediator,
#                                                                     response = response,
#                                                                     response_model_mediators = response_model_mediators,
#                                                                     response_model_exposure = response_model_exposure,
#                                                                     w = weights )$regressionExposure_listReturn
#
# eval_make_DAG_output_regressionMediator_listReturn <- eval_make_DAG(data = dataframe,
#                                                                     regressionExposure = make_DAG_output_resultExposure,
#                                                                     regressionMediator = make_DAG_output_resultMediator,
#                                                                     response = response,
#                                                                     response_model_mediators = response_model_mediators,
#                                                                     response_model_exposure = response_model_exposure,
#                                                                     w = weights )$regressionMediator_listReturn


#######
### NEED TO EVALUATE THEM NEXT
#######

#################################
#################################
#################################
# eval_make_DAG_output
#
# $regressionExposure_listReturn
# $regressionExposure_listReturn[[1]]
# [1] "make_DAG_output$resultExposure[[i]] <-glm(case ~  phys  +  ahei3tert + alcohfreqwk + fatduc + global_stress2 + moteduc + nevfcur + subeduc + regionnn7*ns(eage,df=5)+esex*ns(eage,df=5),data=dataframe,family='binomial',w=weights)"
#
#
# $regressionMediator_listReturn
# $regressionMediator_listReturn[[1]]
# [1] "make_DAG_output$resultMediator[[i]] <-glm(case ~  subhtn  +  ahei3tert + alcohfreqwk + ns(apob_apoa, knots = quantile(apob_apoa,c(.25,0.5,0.75)), Boundary.knots = quantile(apob_apoa,c(.001,0.95))) + fatduc + global_stress2 + moteduc + nevfcur + phys + subeduc + ns(whr,df=5) + regionnn7*ns(eage,df=5)+esex*ns(eage,df=5),data=dataframe,family='binomial',w=weights)"
#
# $regressionMediator_listReturn[[2]]
# [1] "make_DAG_output$resultMediator[[i]] <-glm(case ~  ns(apob_apoa, knots = quantile(apob_apoa,c(.25,0.5,0.75)), Boundary.knots = quantile(apob_apoa,c(.001,0.95)))  +  ahei3tert + alcohfreqwk + fatduc + global_stress2 + moteduc + nevfcur + phys + subeduc + subhtn + ns(whr,df=5) + regionnn7*ns(eage,df=5)+esex*ns(eage,df=5),data=dataframe,family='binomial',w=weights)"
#
# $regressionMediator_listReturn[[3]]
# [1] "make_DAG_output$resultMediator[[i]] <-glm(case ~  ns(whr,df=5)  +  ahei3tert + alcohfreqwk + ns(apob_apoa, knots = quantile(apob_apoa,c(.25,0.5,0.75)), Boundary.knots = quantile(apob_apoa,c(.001,0.95))) + fatduc + global_stress2 + moteduc + nevfcur + phys + subeduc + subhtn + regionnn7*ns(eage,df=5)+esex*ns(eage,df=5),data=dataframe,family='binomial',w=weights)"
#################################
#################################
#################################

## N.B. FOR LOOP HAS TO BE i AS MODELS DEFINED IN [[I]]
      # for( i in 1:length(model_list_usePointEstimate)){
      #   eval(parse(text=model_list_usePointEstimate[[i]]))
      #    model_list_evalPointEstimate[[i]] <- model_list_inputPointEstimate[[i]]
      # }

# Eval_regressionExposure_listReturn <- eval_make_DAG_output_regressionExposure_listReturn

for(i in 1:length(eval_make_DAG_output$regressionExposure_listReturn) ){

                  # eval(parse(text=model_list_usePointEstimate[[i]]))
                  eval(parse(text = eval_make_DAG_output$regressionExposure_listReturn[[i]]))
                  # eval(parse(text = eval_make_DAG_output_regressionExposure_listReturn[[i]]))
                  # model_list_evalPointEstimate[[i]] <- model_list_inputPointEstimate[[i]]
                  # THIS IS HOW THE REGRESSION IS CALLED WITHIN eval_make_DAG.R
                  make_DAG_output$resultExposure[[i]] <- make_DAG_output$resultExposure[[i]]

            }


# Eval_regressionMediator_listReturn <- eval_make_DAG_output_regressionMediator_listReturn

for(i in 1:length(eval_make_DAG_output$regressionMediator_listReturn) ){

                  # eval(parse(text=model_list_usePointEstimate[[i]]))
                  eval(parse(text = eval_make_DAG_output$regressionMediator_listReturn[[i]]))
                  # eval(parse(text = eval_make_DAG_output_regressionMediator_listReturn[[i]]))
                  # model_list_evalPointEstimate[[i]] <- model_list_inputPointEstimate[[i]]
                  make_DAG_output$resultMediator[[i]] <- make_DAG_output$resultMediator[[i]]

            }




###############
################
################

# #####
# ## MOC NOTE: index in for loop here may need to be "i" sunce make-formula function uses "i" in it and passes it into this below as "i".
# #####
# for(i in 1:length(in_out[[2]])){
#
#
#         column <- (1:length(colnames(stroke_reduced)))[colnames(stroke_reduced) %in% in_out[[2]][i]]
#         formula_text <- make_formula(in_out[[1]][[i]],in_out[[2]][i])
#         y <- stroke_reduced[,column]
#         if(length(table(y))==2){
#                 theform <- paste("glm(",formula_text,",data=stroke_reduced,family='binomial',w=w)",sep='')
#         }
#         if(length(table(y))>2 & is.factor(y)){
#                 theform <- paste("polr(",formula_text,",data=stroke_reduced,w=w)",sep='')
#         }
#         if(length(table(y))>2 & is.numeric(y)){
#                 theform <- paste("lm(",formula_text,",data=stroke_reduced,w=w)",sep='')
#         }
#         to_execute <- paste("model_list[[i]] <-", theform,sep='')
#         eval(parse(text=to_execute))
# }



#### MOC: Added in all 3 mediators in one model. Separate models are listed underneath as used previously in earlier versions of model
# response_vs_mediator <-  glm("case ~ regionnn7*ns(eage,df=5)+esex*ns(eage,df=5) +  subeduc+ moteduc+ fatduc+ phys+ ahei3tert+ nevfcur+ alcohfreqwk+ global_stress2+ subhtn + ns(apob_apoa, knots = quantile(apob_apoa,c(.25,0.5,0.75)), Boundary.knots = quantile(apob_apoa,c(.001,0.95)))+ns(whr,df=5)",data = stroke_reduced,family='binomial',w=w)
# could be more than 1 dimensional
# response_vs_mediator <- eval_make_DAG_output$regressionMediator_listReturn
##### for loop above evaluates make_DAG_output$resultMediator
response_vs_mediatorPointEstimate <- make_DAG_output$resultMediator


# ###  Model that estimates causal effect of high blood pressure on stroke
#
# # MOC: USED IN PREVIOUS VERSIONS OF MODEL, BUT UPDATED TO INCLUDE ALL 3 MEDIATORS IN "Response_vs_mediator"
#  response_vs_HBP <- glm("case ~ regionnn7*ns(eage,df=5)+esex*ns(eage,df=5) +  subeduc+ moteduc+ fatduc+ phys+ ahei3tert+ nevfcur+ alcohfreqwk+ global_stress2+ subhtn",data = stroke_reduced,family='binomial',w=w)
#
# # MOC: USED IN PREVIOUS VERSIONS OF MODEL, BUT UPDATED TO INCLUDE ALL 3 MEDIATORS IN "Response_vs_mediator"
#   response_vs_apob_apoa <- glm("case ~ regionnn7*ns(eage,df=5)+esex*ns(eage,df=5) +  subeduc+ moteduc+ fatduc+ phys+ ahei3tert+ nevfcur+ alcohfreqwk+ global_stress2+ ns(apob_apoa, knots = quantile(stroke_reduced$apob_apoa,c(.25,0.5,0.75)), Boundary.knots = quantile(stroke_reduced$apob_apoa,c(.001,0.95)))", data = stroke_reduced,family = 'binomial', w = w )
#
#
# # MOC: USED IN PREVIOUS VERSIONS OF MODEL, BUT UPDATED TO INCLUDE ALL 3 MEDIATORS IN "Response_vs_mediator"
#   response_vs_whr <- glm("case ~ regionnn7*ns(eage,df=5)+esex*ns(eage,df=5) +  subeduc+ moteduc+ fatduc+ phys+ ahei3tert+ nevfcur+ alcohfreqwk+ global_stress2+ ns(whr , df = 5)",data=stroke_reduced,family='binomial',w=w)

#######
#######

# MOC: leaving this without mediators
# response_vs_phys <- glm("case ~ regionnn7*ns(eage,df=5)+esex*ns(eage,df=5) +  subeduc+ moteduc+ fatduc+ phys+ ahei3tert+ nevfcur+ alcohfreqwk+ global_stress2",data=stroke_reduced,family='binomial',w=w)
# could be more than 1 dimensional
# response_vs_phys <- eval_make_DAG_output$regressionExposure_listReturn
# for loop above evaluates make_DAG_output$resultExposure
response_vs_physPointEstimate <- make_DAG_output$resultExposure

#####################################
#####################################
# add in text "PointEstimate"
#####################################

  results_mediator_simulationStorePointEstimate <- list()
  for(i in 1:length(mediator)){
    results_mediator_simulationStorePointEstimate[[i]] <- matrix(nrow = NumSimulation,ncol=5)
    colnames(results_mediator_simulationStorePointEstimate[[i]]) <- c("overall","direct Sjolander","indirect Sjolander","path specific","overall Direct")
  }

  for(med in 1:length(results_mediatorPointEstimate) ){

          for(i in 1:NumSimulation){

# EXPOSURE COULD BE MORE THAN 1 DIMENSIONAL BUT ONLY CODED HERE FRO 1 DEMSIONAL RESULT I.E. MIGHT NEED FOR LOOP FOR EACH EXPOSURE

              results_mediator_simulationStorePointEstimate[[med]][i,1:3] <- indirect_PAF_Sjolander_onesimulation(
                data_frame = dataframe,
                exposure=exposure, # is it an issue that exposure = exposure?
                mediator=mediator[med],
                response = response, # is it an issue that response = response?
                # mediator_model = mediator_modelPointEstimate[[med]], # CHECK IF THIS CHANGE WORKS i.e. ADDING ON [[med]]. # Is it an issue that mediator_model = mediator_model? AND THAT THIS NOW CAN BE MORE THAN 1 DIMENSIONAL
                mediator_model = mediator_modelPointEstimate,
                # response_model=response_vs_mediator,
                response_model = response_vs_mediatorPointEstimate[[med]],  # CHECK IF THIS CHANGE WORKS i.e. ADDING ON [[med]]
                # response_model_2=response_vs_phys,
                response_model_2=response_vs_physPointEstimate[[1]], # EXPOSURE COULD BE MORE THAN 1 DIMENSIONAL BUT ONLY CODED HERE FRO 1 DEMSIONAL RESULT I.E. MIGHT NEED FOR LOOP FOR EACH EXPOSURE
                weights= weights)  # SHOULD BE = weights

              results_mediator_simulationStorePointEstimate[[med]][i,4] <- path_specific_onesimulation(
                data_frame = dataframe,
                exposure = exposure, # is it an issue that exposure = exposure?
                mediator = mediator[med],
                response = response, # is it an issue that response = response?
                # mediator_model = mediator_modelPointEstimate[[med]], # CHECK IF THIS CHANGE WORKS i.e. ADDING ON [[med]]. # is it an issue that mediator_model = mediator_model?
                mediator_model = mediator_modelPointEstimate,
                response_model=response_vs_mediatorPointEstimate[[med]],  # CHECK IF THIS CHANGE WORKS i.e. ADDING ON [[med]]
                response_model_2=response_vs_physPointEstimate[[1]], # EXPOSURE COULD BE MORE THAN 1 DIMENSIONAL BUT ONLY CODED HERE FRO 1 DEMSIONAL RESULT I.E. MIGHT NEED FOR LOOP FOR EACH EXPOSURE
                weights= weights)  # SHOULD BE = weights

              results_mediator_simulationStorePointEstimate[[med]][i,5] <- overall_direct(
                data_frame = dataframe,
                exposure=exposure, # is it an issue that exposure = exposure?
                mediator = mediator[med],
                response = response, # is it an issue that response = response?
                # mediator_model = mediator_modelPointEstimate[[med]], # CHECK IF THIS CHANGE WORKS i.e. ADDING ON [[med]]. # is it an issue that mediator_model = mediator_model?
                mediator_model = mediator_modelPointEstimate,
                response_model=response_vs_mediatorPointEstimate[[med]],  # CHECK IF THIS CHANGE WORKS i.e. ADDING ON [[med]]
                response_model_2=response_vs_physPointEstimate[[1]], # EXPOSURE COULD BE MORE THAN 1 DIMENSIONAL BUT ONLY CODED HERE FRO 1 DEMSIONAL RESULT I.E. MIGHT NEED FOR LOOP FOR EACH EXPOSURE
                weights= weights)  # SHOULD BE = weights

              flush.console()
              print(i)
          }

      # v IS ONLY NEEDED FOR BOOTSTRAP
      # results_mediatorPointEstimate[[med]][v,] = apply(results_mediator_simulationStorePointEstimate[[med]],2,mean)
        results_mediatorPointEstimate[[med]] = apply(results_mediator_simulationStorePointEstimate[[med]],2,mean)

  }

###########
###########
##### ONLY FOR BOOTSTRAP
###########
###########
# results_mediator_tablePointEstimate <- list()
# for(i in 1:length(mediator)){
#     results_mediator_tablePointEstimate[[i]] <- matrix(nrow = 3,ncol=5)
#     results_mediator_tablePointEstimate[[i]][1,] <- apply(results_mediatorPointEstimate[[i]],2,mean)
#     results_mediator_tablePointEstimate[[i]][2,] <- apply(results_mediatorPointEstimate[[i]],2,mean) - 1.96*apply(results_mediatorPointEstimate[[i]],2,sd)
#     results_mediator_tablePointEstimate[[i]][3,] <- apply(results_mediatorPointEstimate[[i]],2,mean) + 1.96*apply(results_mediatorPointEstimate[[i]],2,sd)
#     colnames(results_mediator_tablePointEstimate[[i]]) <- c("overall","direct Sjolander","indirect Sjolander","path specific","overall Direct")
#     rownames(results_mediator_tablePointEstimate[[i]]) <- c("Mean", "Lower 95% C.I.","Upper 95% C.I.")
#     results_mediator_tablePointEstimate[[i]]
# }



# my_listPointEstimate <- list("plot" = plotList, "mediators" = mediator, "table" = results_mediator_table )
# return("results_mediator_tablePointEstimate" = results_mediator_tablePointEstimate)
my_listPointEstimate <- list("results_mediatorPointEstimate" = results_mediatorPointEstimate, "mediators" = mediator, "response_vs_physPointEstimate" = response_vs_physPointEstimate , "response_vs_mediatorPointEstimate" = response_vs_mediatorPointEstimate )


return( my_listPointEstimate )

# make_DAG_output$resultExposure
# make_DAG_output$resultMediator

#####################################
#####################################
#####################################

#   results_subhtn_simulationStore = matrix(nrow = NumSimulation,ncol=5)
#   colnames(results_subhtn_simulationStore) <- c("overall","direct Sjolander","indirect Sjolander","path specific","overall Direct")
#
# for(i in 1:NumSimulation){
#   results_subhtn_simulationStore[i,1:3] <- indirect_PAF_Sjolander_onesimulation(data_frame = stroke_reduced, mediator="subhtn")
#   # data_frame, exposure, mediator, response, mediator_model, response_model, response_model_2, weights
#   results_subhtn_simulationStore[i,4] <- path_specific_onesimulation(data_frame = stroke_reduced, mediator="subhtn")
#   results_subhtn_simulationStore[i,5] <- overall_direct(data_frame = stroke_reduced, mediator="subhtn")
#   flush.console()
#   print(i)
# }
# results_subhtn[1,] = apply(results_subhtn_simulationStore,2,mean)
#
#
#   results_whr_simulationStore = matrix(nrow = NumSimulation,ncol=5)
#   colnames(results_whr_simulationStore) <- c("overall","direct Sjolander","indirect Sjolander","path specific","overall Direct")
#
# for(i in 1:NumSimulation){
#   results_whr_simulationStore[i,1:3] <- indirect_PAF_Sjolander_onesimulation(data_frame = stroke_reduced, mediator="whr")
#   results_whr_simulationStore[i,4] <- path_specific_onesimulation(data_frame = stroke_reduced, mediator="whr")
#   results_whr_simulationStore[i,5] <- overall_direct(data_frame = stroke_reduced, mediator="whr")
#   flush.console()
#   print(i)
#  }
# results_whr[1,] = apply(results_whr_simulationStore,2,mean)
#
#   results_apob_apoa_simulationStore = matrix(nrow = NumSimulation,ncol=5)
#   colnames(results_apob_apoa_simulationStore) <- c("overall","direct Sjolander","indirect Sjolander","path specific","overall Direct")
#
#   for(i in 1:NumSimulation){
#   results_apob_apoa_simulationStore[i,1:3] <- indirect_PAF_Sjolander_onesimulation(data_frame = stroke_reduced, mediator="apob_apoa")
#   results_apob_apoa_simulationStore[i,4] <- path_specific_onesimulation(data_frame = stroke_reduced, mediator="apob_apoa")
#   results_apob_apoa_simulationStore[i,5] <- overall_direct(data_frame = stroke_reduced, mediator="apob_apoa")
#   flush.console()
#   print(i)
#   }
#
#   results_apob_apoa[1,] = apply(results_apob_apoa_simulationStore,2,mean)
#
#   results_subhtn
#   results_whr
#   results_apob_apoa
#
#  #############################
#   #############################
#   #############################
#
# ####################
# ####################
# ####################
#
#
#
# # ##########################################
# # ##########################################
# # ##########################################
# # ## MOC CHECK
# # ##########################################
# # ##########################################
# # ##########################################
# # indirect_PAF_Sjolander_onesimulation(data_frame = stroke_reduced, mediator="subhtn")
# # path_specific_onesimulation(data_frame = stroke_reduced, mediator="subhtn")
# # overall_direct(data_frame = stroke_reduced, mediator="subhtn")
# #
# # indirect_PAF_Sjolander_onesimulation(data_frame = stroke_reduced, mediator="whr")
# # path_specific_onesimulation(data_frame = stroke_reduced, mediator="whr")
# # overall_direct(data_frame = stroke_reduced,mediator="whr")
# #
# # indirect_PAF_Sjolander_onesimulation(data_frame = stroke_reduced, mediator="apob_apoa")
# # path_specific_onesimulation(data_frame = stroke_reduced, mediator="apob_apoa")
# # overall_direct(data_frame = stroke_reduced,mediator="apob_apoa")
#


}
