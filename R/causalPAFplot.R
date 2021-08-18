#' @title Evaluates Total PAF, Direct PAF, Indirect PAF and Path Specific PAF for a user inputted number of bootstraps and integral simulations
#' @description Evaluates Total PAF, Direct PAF, Indirect PAF and Path Specific PAF for a user inputted number of bootstraps and integral simulations
#' @param dataframe A wide format dataframe containing all the risk factors, confounders, exposures and outcomes within the causal DAG Bayesian network.
#' @param exposure The name of the exposure column variable within dataframe in text format e.g. "phys".
#' @param mediator The name of the mediator column variables within dataframe in text format. There can be more than one mediator of interest. It can be a vector of mediators names within the dataframe e.g. c("subhtn","apob_apoa","whr").
#' @param response The name of the response column variable within dataframe in text format e.g. "case". The cases should be coded as 1 and the controls as 0.
#' @param response_model_mediators A regression model fitted for the response in a causal Bayesian network excluding ``children'' of the mediators in the causal Bayesian network. See example in tutorial.This model can be listed either as (1) an empty list ( response_model_mediators = list() ) or (2) the user can specify their own customised causal regression model(s) to use. When it is listed as an empty list the causalPAF package will fit the response_model_mediators regression model automatically based on the causal DAG supplied by the user in in_outArg. Alternatively, the user can specify the exact model(s) that the user wishes to use, these model(s) must be in list format (list() where length(response_model_mediators) == length(mediator) ), the same length as the parameter, mediator, with the user customised model for each mediator listed in the same order as in the parmeter, mediator, and if there is only one model, it must be listed each time within the list() so that length(response_model_mediators) == length(mediator).
#' @param response_model_exposure A regression model fitted for the response in a causal Bayesian network excluding ``children'' of the exposure in the causal Bayesian network. This regression model will not adjust for mediators (exclude mediators) of the exposure in the regression model so that the total effect of the exposure on the response can be modelled. This model can be listed either as (1) an empty list ( response_model_exposure = list() ) or (2) the user can specify their own customised causal regression model to use. If specified as an empty list, list(), then the causalPAF function will define and fit the model automatically based on the causal DAG defined by the in_outArg parameter. Alternatively, the user can specify the exact model that the user wishes to use, this model must be in list format (list() where length(response_model_exposure) == 1 ), of length 1, assuming only one exposure of interest (other exposures can be risk factors) and the model must be defined within a list() since the package assumes a list() format is supplied. See example in tutorial. E.G. If physical exercise ("exer") in the example given in the diagram is the exposure. Then the regression would include all parents of "exer" (i.e. sex, region, educ, age) as well as risk factors at the same level of the causal Bayesian network (i.e. stress, smoke, diet, alcoh).
#' @param in_outArg This defines the causal directed acyclic graph (DAG). A list of length 2. It is defined as a two dimensional list consisting of, firstly, the first list, inlist, i.e. a list of the parents of each variable of interest corresponding to its column name in the data. Splines can be included here if they are to be modelled as splines. Secondly, the second list, outlist, contains a list of a single name of exposure or risk factor or outcome in form of characters i.e. a list of each variable of interest (risk factors, exposures and outcome) corresponding to its column name in the data. Splines should not be input here, only the column names of the variables of interest in the data. The order at which variables are defined must satisfy (i) It is important that variables are defined in the same order in both lists e.g. the first risk factor defined in outlist has its parents listed first in inlist, the second risk factor defined in outlist has its parents listed secondly in inlist and so on. The package assumes this ordering and will not work if this order is violated. (ii) Note it is important also that the order at which the variables are defined is such that all parents of that variable are defined before it. See example in tutorial.
#' @param Splines_outlist A list defined of same size and order of variables as defined in in_outArg[[2]]. If splines are to be used for variables listed in in_outArg[[2]], then the splines should be defined in Splines_outlist in the same order as variables appear in in_outArg[[2]]. It is necessary to list variables in Splines_outlist the same as in in_outArg[[2]] without splines if no spline is to be applied. It should not be input as an empty list, list(), if no splines. A warning will show if input as an empty list requiring the user to populate Splines_outlist either the same as in_outArg[[2]] (if no splines) or in the same order as in_outArg[[2]] with splines (if splines).  See example in tutorial.
#' @param splinesDefinedIn_in_outDAG Logical TRUE or FALSE indicating whether the user has defined splines in the causal DAG, in_out, if TRUE. If FALSE and splines are defined in Splines_outlist_Var, then it is necessary for the package to populate the in_out DAG with splines listed in Splines_outlist_Var.
#' @param model_listArg is a list of models fitted for each of the variables in in_outArg[[2]] (or in_outArg\$outlist ) based on its parents given in in_outArg[[1]] ( or in_out\$inlist ). By default this is set to an empty list. In the default setting, the models are fitted automatically by the causalPAF package based on the order of the variables input in the parameter in_outArg. See the tutorial for more examples. Alternatively, the user can supply their own fitted models here by populating ``model_listArg'' with their own fitted models for each risk factor, mediator, exposure and response varialble. But the order of these models must be in the same order of the variables in the second list of in_outArg ( in_outArg[[2]] ) and these models be defined within a list, list(), of the same length as in_outArg[[2]]. See tutorial for further examples.
#' @param weights Column of weights for case control matching listed in the same order as the patients in the data e.g. weights = strokedata$weights.
#' @param NumBootstrap The number of bootstraps the user wants to use to calculate confidence intervals for the effect. A minimum of 200 bootstrap repilcations (Efron (2016), Computer Age Statistical Inference, page 162) are recommended to calculate standard errors (for intervals of the form: estimate +/-1.96*(standard error of boostrap estimate. However increasing the number of bootstraps can result in the package taking a long time to run. So the user may decide to balance speed with accuracy depedning on which is of more value in the specific context.
#' @param NumSimulation This is the number of simulatons requested by the user to estimate integrals. The larger the number of simulations the more accurate the results but the longer the code takes to run. Therefore the user may wish to balance speed with accuracy depending on which is of more value in the specific context of interest. The integrals for continuous variables are estimated using simulation methods.
#' @param plot plot can be text inputs "forestplot" or "bar" where:"forestplot" plots a forest plot."bar" plots a bar chart with error bars.
#' @param fill The colour for the fill in the bar chart is set here in text format. The default is fill= "skyblue".
#' @param colour The colour for the error bar in the bar chart is set here in text format. The default is colour = "orange".
#' @param addCustom Logical TRUE or FALSE indicating whether a customised interaction term is to be added to the each regression. The interaction term can include splines.
#' @param custom text containing the customised interaction term to be added to each regression. The text should be enclosed in inverted commas. Splines can be included within the interactin terms. See tutorial for examples.
#' @export
#' @importFrom dplyr bind_rows mutate
#' @importFrom ggdag dagify tidy_dagitty ggdag theme_dag
#' @import splines MASS stats forestplot utils grid magrittr checkmate ggplot2
#' @keywords models Regression Population Attributable Fraction
#' @return Prints a forest plot or a bar chart with error bars of the 5 results for each mediator. The 5 results are:(1)Total Population Attributable Fraction (PAF),(2)Direct Effect Population Attributable Fraction (PAF) using alternative definition, (3)Indirect Effect Population Attributable Fraction (PAF) using alternative definition, (4)Path Specific Population Attributable Fraction (PAF), (5)Overall Direct Population Attributable Fraction (PAF)
#' @examples \dontrun{
#' # I don't want you to run this
#' }
#' in_vars = c("subeduc","moteduc","fatduc")
#' outvar = c("phys")
#' make_formula(in_vars,outvar)
causalPAFplot <- function(dataframe,
                          exposure="phys",
                          mediator=c("subhtn","apob_apoa","whr"),
                          response="case",
                          # response_model_mediators=response_vs_mediator, # WILL CAUSE ERROR IF CALL = response_vs_mediator SINCE response_vs_mediator IS A FUNCTION AND ASSIGNS IT TO THE FUNCTION
                          # response_model_exposure=response_vs_phys,
                          response_model_mediators = list(),
                          response_model_exposure = list(),
                          in_outArg,
                          Splines_outlist = list(),   # needs to be input as in_outArg[[2]] if no splines. Assumes exposure is not written in spline format
                          splinesDefinedIn_in_outDAG = list(),
                          model_listArg = list(),
                          # weights = strokedata$weights,
                          weights = 1,  # weights = stroke_reduced$weights
                          # prevalence ADD IN
                          # CHECK IS WEIGHT HARD CODED
                          NumBootstrap,
                          NumSimulation,
                          plot = "bar",
                          # errorbar = "errorbar",
                          fill = "skyblue",
                          colour="orange",
                          # boxCol="royalblue", # @param boxCol The colour for the box in the forest plot is set here in text format. The default is box = "royalblue".
                          # lineCol="darkblue", # @param lineCol The colour for the lines in the forest plot is set here in text format. The default is line = "darkblue".
                          # summaryCol="royalblue", # summaryCol The colour for a summary in the forest plot is set here in text format. The default is summary = "royalblue"
                          addCustom = FALSE,
                          custom = ""
                          ){


  # @param boxCol The colour for the box in the forest plot is set here in text format. The default is box = "royalblue".
  # @param lineCol The colour for the lines in the forest plot is set here in text format. The default is line = "darkblue".
  # @param summaryCol The colour for a summary in the forest plot is set here in text format. The default is summary = "royalblue"

  # ggproto <- ggplot2::ggproto


  # library(dplyr)
  # library(splines)
  # library(MASS)
  # library(stats)
  # library(forestplot)
  # library(utils)
  # library(grid)
  # library(magrittr)
  # library(checkmate)
  # library(ggplot2)
  # ####
  # library(ggdag)
  # library(dagitty)

                          # dataframe = stroke_reduced
                          # exposure="phys"
                          # mediator=c("subhtn","apob_apoa","whr")
                          # response="case"
                          # # response_model_mediators = response_vs_mediator # THIS CAUSES ERROR SINCE A FUNCTION IS called response_vs_mediator AND IT ASSIGNS THE FUNCTION RATHER THAN THE MEDIATOR
                          # response_model_mediators = list()
                          # response_model_exposure = list()
                          # in_outArg = in_out
                          # Splines_outlist = Splines_outlist  # needs to be input as list() if no splines. Assumes if appears as spline once needs to appear as a spline in all occurences
                          # splinesDefinedIn_in_outDAG = TRUE
                          # model_listArg = list()
                          # weights = stroke_reduced$weights
                          # NumBootstrap = 2
                          # NumSimulation = 2
                          # plot = "bar"
                          # # errorbar = "NA"
                          # fill= "skyblue"
                          # colour="orange"
                          # boxCol="royalblue"
                          # lineCol="darkblue"
                          # summaryCol="royalblue"
                          # addCustom = TRUE
                          # custom = "regionnn7*ns(eage,df=5)+esex*ns(eage,df=5)"


                          # # If USER WANTS TO FIT OWN MODELS
                          # response_model_mediators = response_vs_mediator
                          # response_model_exposure = response_vs_phys
                          # model_listArg = model_listArgFit



## At moment, weigths are populated in the dataframe before calling function causalPAFplot.r. The default is weights = 1 which is not a case-control data format but rather a cohort weighting by default.
## CHANGE TO MAKE: NEED TO UPDATE causalPAFplot to link weights to bootstrap$weights below as not defined properly
## CHANGE TO MAKE: Build in so response_model_mediators can be built into causalPAFplot similar to how response_model_exposure has been built in


dataframe <- dplyr::mutate(dataframe, weightsAppendColumnToData = weights)

# addCustom = TRUE, custom = "~ regionnn7*ns(eage,df=5)+esex*ns(eage,df=5) + "

# https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
# MOC Note: Try this to get rid of error: "causalPAFplot: no visible binding for global variable..." and "Undefined global functions or variables:..."
# Bootstrap <- NULL
# globalVariables(c("Bootstrap"))


# checkMarkovDAG(in_outArg)[1]
# checkMarkovDAG(in_outArg)[2]

##############################################
##############################################

# pointEstimate$results_mediatorPointEstimate
# pointEstimate$regressionExposure_listReturn
# pointEstimate$regressionMediator_listReturn

results_mediatorPointEstimate <- pointEstimate(dataframe = dataframe,
                                               exposure = exposure,
                                               mediator = mediator,
                                               response = response,
                                               response_model_mediators = response_model_mediators,
                                               response_model_exposure = response_model_exposure,
                                               in_outArg = in_outArg,
                                               Splines_outlist = Splines_outlist,  # needs to be input as list() if no splines. Assumes if appears as spline once needs to appear as a spline in all occurences
                                               # splinesDefinedIn_in_outDAG = TRUE,
                                               splinesDefinedIn_in_outDAG = splinesDefinedIn_in_outDAG,
                                               model_listArg = model_listArg,
                                               weights = weights,
                                               NumSimulation = NumSimulation,
                                               addCustom = addCustom,
                                               custom = custom # need to update make_formula() so that it runs without needing ~ first
                                               )
#### OUTPUT FROM pointEstimate()
# my_listPointEstimate <- list("results_mediatorPointEstimate" = results_mediatorPointEstimate, "mediators" = mediator, "response_vs_physPointEstimate" = response_vs_physPointEstimate , "response_vs_mediatorPointEstimate" = response_vs_mediatorPointEstimate )

# results_mediatorPointEstimate$results_mediatorPointEstimate
# results_mediatorPointEstimate$mediators
# results_mediatorPointEstimate$response_vs_physPointEstimate
# results_mediatorPointEstimate$response_vs_mediatorPointEstimate

###############################################
###############################################

# ########################
# #########################
# # to be called in causalPAFplot.R
# #########################
#   if(length(response_model_exposure) == 0 ){
#
#         response_model_exposure_text <- response_vs_exposure(data = dataframe ,
#                                                           exposure=exposure,
#                                                           response=response,
#                                                           # in_out = in_out, # needs to be in_out = in_outArg,
#                                                           in_out = in_outArg,
#                                                           w = weights,
#                                                           Splines_outlist = Splines_outlist) ## CHANGE TO MAKE: NEED TO UPDATE causalPAFplot to link weights to bootstrap$weights below as not defined properly
# # In meantime not running function use this...then delete
# # response_model_exposure_text <- to_execute
# # NOTE this eval() will create the model variable named response_model_exposure_text
# eval(parse(text = response_model_exposure_text ) )
#
# response_model_exposure <-response_model_exposure_text
# ########################
# ########################
# ########################
#
#   }
# # else{ response_model_exposure IS DEFINED BY THE USER IN THE FUNCTION}

response_model_exposure <- results_mediatorPointEstimate$response_vs_physPointEstimate

response_model_mediators <- results_mediatorPointEstimate$response_vs_mediatorPointEstimate

# Controls = dataframe[dataframe$response_name == 0, ]
Controls = dataframe[ dataframe[ ,grep(paste('^',response,'$',sep=''),colnames(dataframe),perl=TRUE)] == 0, ]


# Cases = dataframe[dataframe$response_name == 1, ]
Cases = dataframe[ dataframe[ ,grep(paste('^',response,'$',sep=''),colnames(dataframe),perl=TRUE)] == 1, ]


results_mediator <- list()
for(i in 1:length(mediator)){
  results_mediator[[i]] <- matrix(nrow= NumBootstrap,ncol=5)
  # colnames(results_mediator[[i]]) <- c("overall","direct Sjolander","indirect Sjolander","path specific","overall Direct")
  colnames(results_mediator[[i]]) <- c("Total PAF","PAF_{Direct,M^j}","PAF_{Indirect,M^j}","PS-PAF_{A->M^j=>Y}","Direct PAF_{A->Y}")
}

model_list_use <- vector(mode = "list", length = length(in_outArg[[2]]) )
model_list_eval <- vector(mode = "list", length = length(in_outArg[[2]]) )

response_vs_mediatorBootstrap <- list()

## This works since in_outArg is used to define the DAG based on the original input of the user, but is then updated based on the check in the following if statement.
make_DAG_output_check <- make_DAG(in_outDAG = in_outArg ,
                                  exposure = exposure,
                                  response = response,
                                  mediator = mediator,
                                  Splines_outlist_Var = Splines_outlist,
                                  splinesDefinedIn_in_outDAG = TRUE,
                                  addCustomExposureAdjustmentSet = addCustom,
                                  customExposureAdjustmentSet = custom,
                                  addCustomMediatorAdjustmentSet = addCustom,
                                  customMediatorAdjustmentSet = custom )   # custom = "regionnn7*ns(eage,df=5)+esex*ns(eage,df=5)"


for ( v in 1:NumBootstrap ){


  BootstrapControls = as.data.frame( Controls[sample(nrow(Controls), replace = T), ] )

  BootstrapCases = as.data.frame( Cases[sample(nrow(Cases), replace = T), ] )

  # MOC NOTE: need to add in this line "Bootstrap <- as.data.frame(dplyr::bind_rows(BootstrapCases, BootstrapControls) )"
  # with "<-" just before it is defined globally with "<<-" in order to remove error with checking CausalPAF package
  # i.e. remove error "no visible binding for '<<-' assignment to "
  Bootstrap <- as.data.frame(dplyr::bind_rows(BootstrapCases, BootstrapControls) )
  # Bootstrap = as.data.frame(bind_rows(BootstrapCases, BootstrapControls) )
  # TRY GLOBAL NAMING OF BOOTSTRAP
  Bootstrap <<- as.data.frame(dplyr::bind_rows(BootstrapCases, BootstrapControls) )


# If user has not fitted their own model_listArg then fit as follows
if( length(model_listArg) == 0 ){
      model_list_input <- list()
      # model_list_use <- eval_make_formula(data = Bootstrap, in_out=in_outArg,model_list=model_listArg, w=Bootstrap$weights)
      # model_list_use <- eval_make_formula(data = Bootstrap, in_out=in_outArg,model_list=model_list_input, w=Bootstrap$weights, addCustom, custom)

      ## If is NULL, it means there was no update applied and no change is required.
      ## So if not Null, then in_outArg needs to be updated based on (1) dagitty check and (2) check required based on derivation.
      if( !is.null(make_DAG_output_check$my_list_causal$in_outDAG_updatedAfterCheck) ){
          ##in_outArg is a parameter of the function but what to update it based on make_DAG (1) check requested by John Ferguson and (2) Dagitty check
          # in_outArg <- make_DAG_output_check$my_list_causal$in_outDAG_updatedAfterCheck
          model_list_use <- eval_make_formula(data = Bootstrap,
                                              # in_out = in_outArg,
                                              in_out = make_DAG_output_check$my_list_causal$in_outDAG_updatedAfterCheck,
                                              model_list=model_list_input, w=Bootstrap$weightsAppendColumnToData, addCustom, custom)
      }else{
            model_list_use <- eval_make_formula(data = Bootstrap, in_out = in_outArg,model_list=model_list_input, w=Bootstrap$weightsAppendColumnToData, addCustom, custom)
       }
      #   # TRY TO GET IN FUNCTION ENVIRONMENT E.G. model$terms <- eval(model$call$formula)
    #   # model$terms <- eval(model$call$formula)
    #   model_list_use[[1]]$terms <- eval(model_list_use[[1]]$call$formula)
    #   eval(model_list_use[[1]]$call$data )
    #   eval(model_list_use[[1]]$call$family )
    #   eval(model_list_use[[1]]$call$weights )
    # #   glm(formula = phys ~ regionnn7 * ns(eage, df = 5) + esex * ns(eage,
    # # df = 5) + subeduc + moteduc + fatduc, family = "binomial",
    # # data = Bootstrap, weights = Bootstrap$weights)
    #   glm(formula = model_list_use[[1]]$call$formula, family = model_list_use[[1]]$call$family,
    # data = model_list_use[[1]]$call$data, weights = model_list_use[[1]]$call$weights )
    #   # Return these from function and evaluate in this evironment!! Hopefully solve it.
      # eval(parse(text=to_execute))
                          ## model_listReturn[[i]] <- eval(parse(text=to_execute))
                          #model_listReturn[[i]] <- eval( parse(text= paste(model_list_text,"[[i]]",sep='') ) )

      ## N.B. FOR LOOP HAS TO BE i AS MODELS DEFINED IN [[I]]
      for( i in 1:length(model_list_use)){
        eval(parse(text=model_list_use[[i]]))
         model_list_eval[[i]] <- model_list_input[[i]]
      }
      # model_list_use <- eval(parse(text=model_list_use))
      # to_execute <- paste(model_list_text,"[[i]] <-", theform,sep='')
      # eval(parse(text=to_execute))
} else{
       model_list_eval <- model_listArg
}

# Select models for mediators input by user from model_listArg
  mediator_order <- c()
  for(i in 1:length(model_list_eval)){
    mediator_order[i] <-  as.character(formula(model_list_eval[[i]]))[2]
  }
  model_list_mediator_order <- (1:length( mediator_order ))[ mediator_order %in% mediator]

  mediator_model <- list()
  for(i in 1:length(model_list_mediator_order)){
    mediator_model[[i]] = model_list_eval[[ model_list_mediator_order[i] ]]
  }




# # response_vs_mediator <- update(response_model_mediators, data = Bootstrap, weights = Bootstrap$weights )
#########
## response_vs_mediator CAN BE MORE THAN 1 DIMENSIONAL IF MORE THAN 1 MEDIATOR OF INTEREST. SO THIS NEEDS TO BE MOVED WITHIN FOR LOOP FOR med i.e. move it below
#########
# response_vs_mediator <- update(response_model_mediators, data = Bootstrap, weights = Bootstrap$weightsAppendColumnToData )

# response_vs_phys <- update(response_model_exposure, data = Bootstrap, weights = Bootstrap$weights )
##############
## Keeping as one dimensional so leaving here. MIGHT NEED TO BE UPDATED IF WE WANT MORE THAN 1 EXPOSURE RUN AT A TIME. IN THIS CASE IT WOULD NEED TO BE MOVED BELOW INTO THE FOR LOOP FOR med
##############
## NB NB NB NBNB NB NB NB  NB NB NB NB  NB NB NB NB
###### NOTE USED [[1]] [[1]] [[1]] [[1]] [[1]] AS ASSUMING ONLY 1 EXPOSURE
response_vs_phys <- update(response_model_exposure[[1]], data = Bootstrap, weights = Bootstrap$weightsAppendColumnToData )
## NB NB NB NBNB NB NB NB  NB NB NB NB  NB NB NB NB


  results_mediator_simulationStore <- list()
  for(i in 1:length(mediator)){
    results_mediator_simulationStore[[i]] <- matrix(nrow = NumSimulation,ncol=5)
    # colnames(results_mediator_simulationStore[[i]]) <- c("overall","direct Sjolander","indirect Sjolander","path specific","overall Direct")
    colnames(results_mediator_simulationStore[[i]]) <- c("Total PAF","PAF_{Direct,M^j}","PAF_{Indirect,M^j}","PS-PAF_{A->M^j=>Y}","Direct PAF_{A->Y}")
  }

  for(med in 1:length(results_mediator) ){

    # response_vs_mediator <- update(response_model_mediators, data = Bootstrap, weights = Bootstrap$weights )
    response_vs_mediatorBootstrap[[med]] <- update(response_model_mediators[[med]], data = Bootstrap, weights = Bootstrap$weightsAppendColumnToData )

          for(i in 1:NumSimulation){

              results_mediator_simulationStore[[med]][i,1:3] <- indirect_PAF_Sjolander_onesimulation(
                data_frame = Bootstrap,
                exposure=exposure, # is it an issue that exposure = exposure?
                mediator=mediator[med],
                response = response, # is it an issue that response = response?
                mediator_model = mediator_model, # is it an issue that mediator_model = mediator_model?
                response_model=response_vs_mediatorBootstrap[[med]],
                response_model_2=response_vs_phys,
                #weights= Bootstrap$weights
                weights= Bootstrap$weightsAppendColumnToData)

              results_mediator_simulationStore[[med]][i,4] <- path_specific_onesimulation(
                data_frame = Bootstrap,
                exposure=exposure, # is it an issue that exposure = exposure?
                mediator = mediator[med],
                response = response, # is it an issue that response = response?
                mediator_model = mediator_model, # is it an issue that mediator_model = mediator_model?
                response_model=response_vs_mediatorBootstrap[[med]],
                response_model_2=response_vs_phys,
                # weights= Bootstrap$weights
                weights= Bootstrap$weightsAppendColumnToData)

              results_mediator_simulationStore[[med]][i,5] <- overall_direct(
                data_frame = Bootstrap,
                exposure=exposure, # is it an issue that exposure = exposure?
                mediator = mediator[med],
                response = response, # is it an issue that response = response?
                mediator_model = mediator_model, # is it an issue that mediator_model = mediator_model?
                response_model=response_vs_mediatorBootstrap[[med]],
                response_model_2=response_vs_phys,
                # weights= Bootstrap$weights
                weights= Bootstrap$weightsAppendColumnToData)

              flush.console()
              print(i)
          }

      results_mediator[[med]][v,] = apply(results_mediator_simulationStore[[med]],2,mean)

  }


}


#######################################
#######################################
### call pointEstimate() function
#######################################
#######################################

            # pointEstimate(dataframe = stroke_reduced,
            #               exposure="phys",
            #               mediator=c("subhtn","apob_apoa","whr"),
            #               response="case",
            #               response_model_mediators = list(),
            #               response_model_exposure = list(),
            #               in_outArg = in_out,
            #               Splines_outlist = Splines_outlist,  # needs to be input as list() if no splines. Assumes if appears as spline once needs to appear as a spline in all occurences
            #               splinesDefinedIn_in_outDAG = TRUE,
            #               model_listArg = list(),
            #               weights = stroke_reduced$weights,
            #               NumSimulation = 2,
            #               plot = "bar",
            #               fill= "skyblue",
            #               colour="orange",
            #               boxCol="royalblue",
            #               lineCol="darkblue",
            #               summaryCol="royalblue",
            #               addCustom = TRUE,
            #               custom = "regionnn7*ns(eage,df=5)+esex*ns(eage,df=5)" # need to update make_formula() so that it runs without needing ~ first
            #               )

############################
############################

############################################################
############################################################
############################################################

results_mediator_table <- list()
for(i in 1:length(mediator)){
    results_mediator_table[[i]] <- matrix(nrow = 3,ncol=5)
    #########
    ## apply(results_mediator[[i]],2,mean) THIS IS BAGGING AND WE DON NOT WANT THE BOOTSTRAPPED MEAN BUT WE WANT THE POINT ESTIMATE WITHOUT BOOTSTRAPING FROM THE ORGINAL DATA
    #########
    # results_mediator_table[[i]][1,] <- apply(results_mediator[[i]],2,mean)
    ## results_mediatorPointEstimate IS THE POINT ESTIMATE FROM THE ORINGAL UNBOOTSTRAPPED DATA GOT FROM pointEstimate() FUNCTION
    results_mediator_table[[i]][1,] <- results_mediatorPointEstimate$results_mediatorPointEstimate[[i]]
    # THESE ARE THE BOOTSTRAPPED CONFIDENCE INTERVALS AND WE WANT TO BOOTSTRAP FOR THE CONFIDENE INTERVALS.
    results_mediator_table[[i]][2,] <- apply(results_mediator[[i]],2,mean) - 1.96*apply(results_mediator[[i]],2,sd)
    results_mediator_table[[i]][3,] <- apply(results_mediator[[i]],2,mean) + 1.96*apply(results_mediator[[i]],2,sd)
    # colnames(results_mediator_table[[i]]) <- c("overall","direct Sjolander","indirect Sjolander","path specific","overall Direct")
    colnames(results_mediator_table[[i]]) <- c("Total PAF","PAF_{Direct,M^j}","PAF_{Indirect,M^j}","PS-PAF_{A->M^j=>Y}","Direct PAF_{A->Y}")
    rownames(results_mediator_table[[i]]) <- c("Mean", "Lower 95% C.I.","Upper 95% C.I.")
    results_mediator_table[[i]]
}

############################################################
############################################################
############################################################




plotList <- list()
cochrane_from_rmeta <- list()
tabletext <- list()
data <- list()
mean <- list()
lower95CI <- list()
upper95CI <- list()

     if( plot == "forestplot" ){
         # for(Nmed in 1:length(mediator)){
         #       cochrane_from_rmeta[[Nmed]] <- structure(list(
         #              mean  = c(NA, results_mediator_table[[Nmed]][1,1], results_mediator_table[[Nmed]][1,2], results_mediator_table[[Nmed]][1,3], results_mediator_table[[Nmed]][1,4], results_mediator_table[[Nmed]][1,5] ),
         #              lower = c(NA, results_mediator_table[[Nmed]][2,1], results_mediator_table[[Nmed]][2,2], results_mediator_table[[Nmed]][2,3], results_mediator_table[[Nmed]][2,4], results_mediator_table[[Nmed]][2,5]),
         #              upper = c(NA, results_mediator_table[[Nmed]][3,1], results_mediator_table[[Nmed]][3,2], results_mediator_table[[Nmed]][3,3], results_mediator_table[[Nmed]][3,4], results_mediator_table[[Nmed]][3,5])),
         #              .Names = c("mean", "lower 95% C.I.", "upper 95% C.I."),
         #              row.names = c(NA, -6L),
         #              class = "data.frame")
         #
         #
         #            tabletext[[Nmed]]<-cbind(
         #            # c("", "overall", "direct Sjolander", "indirect Sjolander","path specific", "overall Direct"),
         #            c("", "Total PAF", "PAF_{Direct,M^j}", "PAF_{Indirect,M^j}","PS-PAF_{A->M^j=>Y}", "Direct PAF_{A->Y}"),
         #            c("Mean", paste(round(results_mediator_table[[Nmed]][1,1],4)) , paste(round(results_mediator_table[[Nmed]][1,2],4)), paste(round(results_mediator_table[[Nmed]][1,3],4)), paste(round(results_mediator_table[[Nmed]][1,4],4)) ,paste(round(results_mediator_table[[Nmed]][1,5],4)) ),
         #            c("lower 95% C.I.", paste(round(results_mediator_table[[Nmed]][2,1],4)) , paste(round(results_mediator_table[[Nmed]][2,2],4)), paste(round(results_mediator_table[[Nmed]][2,3],4)), paste(round(results_mediator_table[[Nmed]][2,4],4)), paste(round(results_mediator_table[[Nmed]][2,5],4)) ),
         #            c("Upper 95% C.I.", paste(round(results_mediator_table[[Nmed]][3,1],4)), paste(round(results_mediator_table[[Nmed]][3,2],4)), paste(round(results_mediator_table[[Nmed]][3,3],4)), paste(round(results_mediator_table[[Nmed]][3,4],4)), paste(round(results_mediator_table[[Nmed]][3,5],4))))
         #
         #
         #             plotList[[Nmed]] <- forestplot(tabletext[[Nmed]],
         #                        cochrane_from_rmeta[[Nmed]],new_page = TRUE,
         #                        is.summary=c(TRUE,rep(FALSE,5)),
         #                        clip=c(min(results_mediator_table[[Nmed]]) -0.1,max(results_mediator_table[[Nmed]]) +0.1),
         #                        xlog=FALSE,
         #                        col=fpColors(box= boxCol,lines = lineCol, summary = summaryCol))
         # }
     } else if(plot == "bar"){
                # if(errorbar == "errorbar" ){
                #    for(Nmed in 1:length(mediator)){
                #         data[[Nmed]] <- data.frame(
                #         name=c("overall", "direct Sjolander", "indirect Sjolander","path specific", "overall Direct"),
                #         mean = results_mediator_table[[Nmed]][1,],
                #         lower95CI = results_mediator_table[[Nmed]][2,],
                #         upper95CI = results_mediator_table[[Nmed]][3,] )
                #                 # Most basic error bar
                #                 # barchart <-
                #                  plotList[[Nmed]] <-ggplot(data[[med]]) +
                #                     geom_bar( aes(x = name, y = mean), stat="identity", fill = fill, alpha=0.7) +
                #                     geom_errorbar( aes(x = name, ymin = lower95CI, ymax = upper95CI), width=0.4, colour = colour, alpha=0.9, size=1.3)
                #    }
                # } else if(errorbar == "rectangle") {
                #    for(Nmed in 1:length(mediator)){
                #         data[[Nmed]] <- data.frame(
                #         name=c("overall", "direct Sjolander", "indirect Sjolander","path specific", "overall Direct"),
                #         mean = results_mediator_table[[Nmed]][1,],
                #         lower95CI = results_mediator_table[[Nmed]][2,],
                #         upper95CI = results_mediator_table[[Nmed]][3,] )
                #            # rectangle
                #            # barchart <-
                #            plotList[[Nmed]] <- ggplot(data[[med]]) +
                #             geom_bar( aes(x = name, y = mean), stat="identity", fill = fill, alpha=0.7) +
                #             geom_crossbar( aes(x=name, y = mean, ymin = lower95CI, ymax = upper95CI), width=0.4, colour = colour, alpha=0.9, size=1.3)
                #    }
                # }else if( errorbar == "line"){
                #    for(Nmed in 1:length(mediator)){
                #         data[[Nmed]] <- data.frame(
                #         name=c("overall", "direct Sjolander", "indirect Sjolander","path specific", "overall Direct"),
                #         mean = results_mediator_table[[Nmed]][1,],
                #         lower95CI = results_mediator_table[[Nmed]][2,],
                #         upper95CI = results_mediator_table[[Nmed]][3,] )
                #          # line
                #            # barchart <-
                #             plotList[[Nmed]] <- ggplot(data[[med]]) +
                #             geom_bar( aes(x = name, y = mean), stat="identity", fill = fill, alpha=0.7) +
                #             geom_linerange( aes(x = name, ymin = lower95CI, ymax = upper95CI), colour = colour, alpha=0.9, size=1.3)
                #    }
                # } else if ( errorbar == "linedot"){
                #    for(Nmed in 1:length(mediator)){
                #         data[[Nmed]] <- data.frame(
                #         name=c("overall", "direct Sjolander", "indirect Sjolander","path specific", "overall Direct"),
                #         mean = results_mediator_table[[Nmed]][1,],
                #         lower95CI = results_mediator_table[[Nmed]][2,],
                #         upper95CI = results_mediator_table[[Nmed]][3,] )
                #           # line + dot
                #            #barchart <-
                #              plotList[[Nmed]] <- ggplot(data[[med]]) +
                #             geom_bar( aes(x = name, y = mean), stat="identity", fill = fill, alpha=0.7) +
                #             geom_pointrange( aes(x=name, y = mean, ymin=lower95CI, ymax = upper95CI), colour = colour, alpha=0.9, size=1.3)
                #    }
                # } else if(errorbar == "horizontal"){
                #    for(Nmed in 1:length(mediator)){
                #         data[[Nmed]] <- data.frame(
                #         name=c("overall", "direct Sjolander", "indirect Sjolander","path specific", "overall Direct"),
                #         mean = results_mediator_table[[Nmed]][1,],
                #         lower95CI = results_mediator_table[[Nmed]][2,],
                #         upper95CI = results_mediator_table[[Nmed]][3,] )
                #           # horizontal
                #            #barchart <-
                #             plotList[[Nmed]] <- ggplot(data[[med]]) +
                #             geom_bar( aes(x = name, y = mean), stat="identity", fill = fill, alpha=0.7) +
                #             geom_errorbar( aes(x=name, ymin=lower95CI, ymax = upper95CI), width=0.4, colour = colour, alpha=0.9, size=1.3) +
                #             coord_flip()
                #    }
                # } else{
                #   print("Input argument for errorbar.")
                # }
               for( Nmed in 1:length(mediator)){

                      data[[Nmed]] <- data.frame(
                      # name=c("overall", "direct Sjolander", "indirect Sjolander","path specific", "overall Direct"),
                      name=c("Total PAF", "PAF_{Direct,M^j}", "PAF_{Indirect,M^j}","PS-PAF_{A->M^j=>Y}", "Direct PAF_{A->Y}"),
                      mean = results_mediator_table[[Nmed]][1,],
                      lower95CI = results_mediator_table[[Nmed]][2,],
                      upper95CI = results_mediator_table[[Nmed]][3,] )

                      plotList[[Nmed]] <-  ggplot(data[[Nmed]]) +
                      geom_bar( aes(x = name, y = mean), stat="identity", fill = fill, alpha=0.7) +
                      geom_errorbar( aes(x=name, ymin=lower95CI, ymax = upper95CI), width=0.4, colour = colour, alpha=0.9, size=1.3) +
                      coord_flip()
                      }
     } else{ print("Enter value for plot argument.") }



     # plotList
     #
     # results_mediator_table

     my_list <- list("plot" = plotList, "mediators" = mediator, "table" = results_mediator_table )
     return(my_list)

}
