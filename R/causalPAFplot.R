#' @title Evaluates Total PAF, Direct PAF, Indirect PAF and Path Specific PAF for a user inputted number of bootstraps and integral simulations
#' @description Evaluates Total PAF, Direct PAF, Indirect PAF and Path Specific PAF for a user inputted number of bootstraps and integral simulations
#' @param dataframe A wide format dataframe containing all the risk factors, confounders, exposures and outcomes within the causal DAG Bayesian network.
#' @param exposure The name of the exposure column variable within dataframe in text format e.g. "phys".
#' @param mediator The name of the mediator column variables within dataframe in text format. There can be more than one mediator of interest. It can be a vector of mediators names within the dataframe e.g. c("subhtn","apob_apoa","whr").
#' @param response The name of the response column variable within dataframe in text format e.g. "case". The cases should be coded as 1 and the controls as 0.
#' @param mediator_model A vector of mediator models fitted which must correspond to the same order of the mediator names in the variable ``mediator''.
#' @param response_model_mediators A model fitted for the response in a causal Bayesian network excluding ``children'' of the mediators in the causal Bayesian network. See example in tutorial.
#' @param response_model_exposure A model fitted for the response in a causal Bayesian network excluding ``children'' of the exposure and risk factors in the causal Bayesian network. See example in tutorial.
#' @param in_out A list of length 2. The first list contains a list of character vectors of the parents of the exposure or risk factor or outcome which are either causes or confounders of the exposure or risk factor or outcome. The second list conttains a list of a single name of exposure or risk factor or outcome in form of characters. See tutorial examples for examples.
#' @param model_list By default this is set to an empty list. In the default setting, the models are fitted based on the order of the variables input in the parameter in_out. See the tutorial for more examples. Alternatively, the user can supply their own fitted models here by populating ``model_list'' with their own fitted models for each risk factor, mediator, exposure and response varialble. But the order of these models must be in the same order of the variables in the second list of in_out. See tutorial for further examples.
#' @param w Column of weights for case control matching listing in same order as patients in dataframe.
#' @param NumBootstrap The number of bootstraps the user wants to use to calculate confidence intervals for the effect. A minimum of 200 bootstrap repilcations (Efron (2016), Computer Age Statistical Inference, page 162) are recommended to calculate standard errors (for intervals of the form: estimate +/-1.96*(standard error of boostrap estimate. However increasing the number of bootstraps can result in the package taking a long time to run. So the user may make to balance speed with accuracy depedning on which is of more value in context.
#' @param NumSimulation This is the number of simulatons requested by the user to estimate integrals. The larger the number of simulations the more accurate the results but the longer the code takes to run. Therefore the user may wish to balance speed with accuracy depedning on which is of more value in the specific context of interest. The integrals for continuous variables are estimated using simulation methods.
#' @export
#' @import splines MASS stats dplyr forestplot utils
#' @keywords models Regression Population Attributable Fraction
#' @return Prints a forest plot of the 5 results for each mediator. The 5 results are:(1)Total Population Attributable Fraction (PAF),(2)Direct Effect Population Attributable Fraction (PAF) using  Sjolanders definition, (3)Indirect Effect Population Attributable Fraction (PAF) using  Sjolanders definition}, (4)Path Specific Population Attributable Fraction (PAF)}, (5)Overall Direct Population Attributable Fraction (PAF)
#' @examples \dontrun{
#' # I don't want you to run this
#' }
#' in_vars = c("subeduc","moteduc","fatduc")
#' outvar = c("phys")
#' make_formula(in_vars,outvar)
causalPAFplot <- function(dataframe,exposure="phys",mediator=c("subhtn","apob_apoa","whr"),response="case", mediator_model=list(model_list[[6]],model_list[[7]],model_list[[8]]),  response_model_mediators=response_vs_mediator, response_model_exposure=response_vs_phys,in_out, model_list = list(), w, NumBootstrap, NumSimulation ){

Controls = dataframe[dataframe$response == 0, ]

Cases = dataframe[dataframe$response == 1, ]


results_mediator <- list()
for(i in 1:length(mediator)){
  results_mediator[[i]] <- matrix(nrow= NumBootstrap,ncol=5)
  colnames(results_mediator[[i]]) <- c("overall","direct Sjolander","indirect Sjolander","path specific","overall Direct")
}


for ( v in 1:NumBootstrap ){


  BootstrapControls = as.data.frame( Controls[sample(nrow(Controls), replace = T), ] )

  BootstrapCases = as.data.frame( Cases[sample(nrow(Cases), replace = T), ] )

  Bootstrap = as.data.frame(bind_rows(BootstrapCases, BootstrapControls) )


model_list <- eval_make_formula(Bootstrap,in_out,model_list, Bootstrap$weights)

response_vs_mediator <- update(response_model_mediators, data = Bootstrap)

response_vs_phys <- update(response_model_exposure, data = Bootstrap)



  results_mediator_simulationStore <- list()
  for(i in 1:length(mediator)){
    results_mediator_simulationStore[[i]] <- matrix(nrow = NumSimulation,ncol=5)
    results_mediator_simulationStore[[i]] <- c("overall","direct Sjolander","indirect Sjolander","path specific","overall Direct")
  }

  for(med in 1:length(results_mediator) ){

          for(i in 1:NumSimulation){

              results_mediator_simulationStore[[med]][i,1:3] <- indirect_PAF_Sjolander_onesimulation(
                data_frame = Bootstrap,
                exposure=exposure, # is it an issue that exposure = exposure?
                mediator=mediator[med],
                response = response, # is it an issue that response = response?
                mediator_model = mediator_model, # is it an issue that mediator_model = mediator_model?
                response_model=response_vs_mediator,
                response_model_2=response_vs_phys,
                weights=w)

              results_mediator_simulationStore[[med]][i,4] <- path_specific_onesimulation(
                data_frame = Bootstrap,
                exposure=exposure, # is it an issue that exposure = exposure?
                mediator = mediator[med],
                response = response, # is it an issue that response = response?
                mediator_model = mediator_model, # is it an issue that mediator_model = mediator_model?
                response_model=response_vs_mediator,
                response_model_2=response_vs_phys,
                weights=w)

              results_mediator_simulationStore[[med]][i,5] <- overall_direct(
                data_frame = Bootstrap,
                exposure=exposure, # is it an issue that exposure = exposure?
                mediator = mediator[med],
                response = response, # is it an issue that response = response?
                mediator_model = mediator_model, # is it an issue that mediator_model = mediator_model?
                response_model=response_vs_mediator,
                response_model_2=response_vs_phys,
                weights=w)

              flush.console()
              print(i)
          }

      results_mediator[[med]][v,] = apply(results_mediator_simulationStore[[med]],2,mean)

  }


}


results_mediator_table <- list()
for(i in 1:length(mediator)){
    results_mediator_table[[i]] <- matrix(nrow = 3,ncol=5)
    results_mediator_table[[i]][1,] <- apply(results_mediator[[i]],2,mean)
    results_mediator_table[[i]][2,] <- apply(results_mediator[[i]],2,mean) - 1.96*apply(results_mediator[[i]],2,sd)
    results_mediator_table[[i]][3,] <- apply(results_mediator[[i]],2,mean) + 1.96*apply(results_mediator[[i]],2,sd)
    colnames(results_mediator_table[[i]]) <- c("overall","direct Sjolander","indirect Sjolander","path specific","overall Direct")
    rownames(results_mediator_table[[i]]) <- c("Mean", "Lower 95% C.I.","Upper 95% C.I.")
    results_mediator_table[[i]]
  }


cochrane_from_rmeta <- list()
tabletext <- list()
for(med in 1:length(mediator)){
            cochrane_from_rmeta[[med]] <-
            structure(list(
              mean  = c(NA, results_mediator_table[[med]][1,1], results_mediator_table[[med]][1,2], results_mediator_table[[med]][1,3], results_mediator_table[[med]][1,4], results_mediator_table[[med]][1,5] ),
              lower = c(NA, results_mediator_table[[med]][2,1], results_mediator_table[[med]][2,2], results_mediator_table[[med]][2,3], results_mediator_table[[med]][2,4], results_mediator_table[[med]][2,5]),
              upper = c(NA, results_mediator_table[[med]][3,1], results_mediator_table[[med]][3,2], results_mediator_table[[med]][3,3], results_mediator_table[[med]][3,4], results_mediator_table[[med]][3,5])),
              .Names = c("mean", "lower 95% C.I.", "upper 95% C.I."),
              row.names = c(NA, -4L),
              class = "data.frame")



          tabletext[[med]]<-cbind(
            c("", "overall", "direct Sjolander", "indirect Sjolander",
              "path specific", "overall Direct"),
            c("Mean", paste(round(results_mediator_table[[med]][1,1],4)) , paste(round(results_mediator_table[[med]][1,2],4)), paste(round(results_mediator_table[[med]][1,3],4)), paste(round(results_mediator_table[[med]][1,4],4)) ,paste(round(results_mediator_table[[med]][1,5],4)) ),
            c("lower 95% C.I.", paste(round(results_mediator_table[[med]][2,1],4)) , paste(round(results_mediator_table[[med]][2,2],4)), paste(round(results_mediator_table[[med]][2,3],4)), paste(round(results_mediator_table[[med]][2,4],4)), paste(round(results_mediator_table[[med]][2,5],4)) ),
            c("Upper 95% C.I.", paste(round(results_mediator_table[[med]][3,1],4)), paste(round(results_mediator_table[[med]][3,2],4)), paste(round(results_mediator_table[[med]][3,3],4)), paste(round(results_mediator_table[[med]][3,4],4)), paste(round(results_mediator_table[[med]][3,5],4))))

          forestplot(tabletext[[med]],
                     cochrane_from_rmeta[[med]],new_page = TRUE,
                     is.summary=c(TRUE,rep(FALSE,3)),
                     clip=c(min(results_mediator_table[[med]]) -0.1,max(results_mediator_table[[med]]) +0.1),
                     xlog=FALSE,
                     col=fpColors(box="royalblue",line="darkblue", summary="royalblue"))

  }



}
