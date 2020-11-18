#' @title Evaluates Total PAF, Direct PAF, Indirect PAF and Path Specific PAF for a user inputted number of bootstraps and integral simulations
#' @description Evaluates Total PAF, Direct PAF, Indirect PAF and Path Specific PAF for a user inputted number of bootstraps and integral simulations
#' @param dataframe A wide format dataframe containing all the risk factors, confounders, exposures and outcomes within the causal DAG Bayesian network.
#' @param exposure The name of the exposure column variable within dataframe in text format e.g. "phys".
#' @param mediator The name of the mediator column variables within dataframe in text format. There can be more than one mediator of interest. It can be a vector of mediators names within the dataframe e.g. c("subhtn","apob_apoa","whr").
#' @param response The name of the response column variable within dataframe in text format e.g. "case". The cases should be coded as 1 and the controls as 0.
#' @param response_model_mediators A model fitted for the response in a causal Bayesian network excluding ``children'' of the mediators in the causal Bayesian network. See example in tutorial.
#' @param response_model_exposure A model fitted for the response in a causal Bayesian network excluding ``children'' of the exposure and risk factors in the causal Bayesian network. See example in tutorial.
#' @param in_outArg A list of length 2. The first list contains a list of character vectors of the parents of the exposure or risk factor or outcome which are either causes or confounders of the exposure or risk factor or outcome. The second list conttains a list of a single name of exposure or risk factor or outcome in form of characters. See tutorial examples for examples.
#' @param model_listArg By default this is set to an empty list. In the default setting, the models are fitted based on the order of the variables input in the parameter in_outArg. See the tutorial for more examples. Alternatively, the user can supply their own fitted models here by populating ``model_listArg'' with their own fitted models for each risk factor, mediator, exposure and response varialble. But the order of these models must be in the same order of the variables in the second list of in_outArg. See tutorial for further examples.
#' @param NumBootstrap The number of bootstraps the user wants to use to calculate confidence intervals for the effect. A minimum of 200 bootstrap repilcations (Efron (2016), Computer Age Statistical Inference, page 162) are recommended to calculate standard errors (for intervals of the form: estimate +/-1.96*(standard error of boostrap estimate. However increasing the number of bootstraps can result in the package taking a long time to run. So the user may make to balance speed with accuracy depedning on which is of more value in context.
#' @param NumSimulation This is the number of simulatons requested by the user to estimate integrals. The larger the number of simulations the more accurate the results but the longer the code takes to run. Therefore the user may wish to balance speed with accuracy depedning on which is of more value in the specific context of interest. The integrals for continuous variables are estimated using simulation methods.
#' @param plot plot can be "forestplot" or "bar" are text inputs where:"forestplot" plots a forest plot."bar" plots a bar chart with error bars.
#' @param fill The colour for the fill in the bar chart is set here in text format. The default is fill= "skyblue".
#' @param colour The colour for the error bar in teh bar chart is set here in text format. The default is colour = "orange".
#' @param boxCol The colour for the box in the forest plot is set here in text format. The default is box = "royalblue".
#' @param lineCol The colour for the lines in the forest plot is set here in text format. The default is line="darkblue".
#' @param summaryCol The colour for a summary in the forest plot is set here in text format. The default is summary="royalblue"
#' @export
#' @import splines MASS stats dplyr forestplot utils grid magrittr checkmate ggplot2
#' @keywords models Regression Population Attributable Fraction
#' @return Prints a forest plot or a bar chart with error bars of the 5 results for each mediator. The 5 results are:(1)Total Population Attributable Fraction (PAF),(2)Direct Effect Population Attributable Fraction (PAF) using  Sjolanders definition, (3)Indirect Effect Population Attributable Fraction (PAF) using  Sjolanders definition, (4)Path Specific Population Attributable Fraction (PAF), (5)Overall Direct Population Attributable Fraction (PAF)
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
                          response_model_mediators=response_vs_mediator,
                          response_model_exposure=response_vs_phys,
                          in_outArg,
                          model_listArg,
                          NumBootstrap,
                          NumSimulation,
                          plot = "bar",
                          # errorbar = "errorbar",
                          fill= "skyblue",
                          colour="orange",
                          boxCol="royalblue",
                          lineCol="darkblue",
                          summaryCol="royalblue"
                          ){


# Controls = dataframe[dataframe$response_name == 0, ]
Controls = dataframe[ dataframe[ ,grep(paste('^',response,'$',sep=''),colnames(dataframe),perl=TRUE)] == 0, ]


# Cases = dataframe[dataframe$response_name == 1, ]
Cases = dataframe[ dataframe[ ,grep(paste('^',response,'$',sep=''),colnames(dataframe),perl=TRUE)] == 1, ]


results_mediator <- list()
for(i in 1:length(mediator)){
  results_mediator[[i]] <- matrix(nrow= NumBootstrap,ncol=5)
  colnames(results_mediator[[i]]) <- c("overall","direct Sjolander","indirect Sjolander","path specific","overall Direct")
}

model_list_use <- vector(mode = "list", length = length(in_outArg[[2]]) )
model_list_eval <- vector(mode = "list", length = length(in_outArg[[2]]) )

for ( v in 1:NumBootstrap ){


  BootstrapControls = as.data.frame( Controls[sample(nrow(Controls), replace = T), ] )

  BootstrapCases = as.data.frame( Cases[sample(nrow(Cases), replace = T), ] )

  # Bootstrap = as.data.frame(bind_rows(BootstrapCases, BootstrapControls) )
  # TRY GLOBAL NAMING OF BOOTSTRAP
  Bootstrap <<- as.data.frame(bind_rows(BootstrapCases, BootstrapControls) )


# If user has not fitted their own model_listArg then fit as follows
if( length(model_listArg) == 0 ){
      model_list_input <- list()
      # model_list_use <- eval_make_formula(data = Bootstrap, in_out=in_outArg,model_list=model_listArg, w=Bootstrap$weights)
      model_list_use <- eval_make_formula(data = Bootstrap, in_out=in_outArg,model_list=model_list_input, w=Bootstrap$weights)
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




response_vs_mediator <- update(response_model_mediators, data = Bootstrap, weights = Bootstrap$weights )

response_vs_phys <- update(response_model_exposure, data = Bootstrap, weights = Bootstrap$weights )



  results_mediator_simulationStore <- list()
  for(i in 1:length(mediator)){
    results_mediator_simulationStore[[i]] <- matrix(nrow = NumSimulation,ncol=5)
    colnames(results_mediator_simulationStore[[i]]) <- c("overall","direct Sjolander","indirect Sjolander","path specific","overall Direct")
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
                weights= Bootstrap$weights)

              results_mediator_simulationStore[[med]][i,4] <- path_specific_onesimulation(
                data_frame = Bootstrap,
                exposure=exposure, # is it an issue that exposure = exposure?
                mediator = mediator[med],
                response = response, # is it an issue that response = response?
                mediator_model = mediator_model, # is it an issue that mediator_model = mediator_model?
                response_model=response_vs_mediator,
                response_model_2=response_vs_phys,
                weights= Bootstrap$weights)

              results_mediator_simulationStore[[med]][i,5] <- overall_direct(
                data_frame = Bootstrap,
                exposure=exposure, # is it an issue that exposure = exposure?
                mediator = mediator[med],
                response = response, # is it an issue that response = response?
                mediator_model = mediator_model, # is it an issue that mediator_model = mediator_model?
                response_model=response_vs_mediator,
                response_model_2=response_vs_phys,
                weights= Bootstrap$weights)

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

plotList <- list()
cochrane_from_rmeta <- list()
tabletext <- list()
data <- list()
mean <- list()
lower95CI <- list()
upper95CI <- list()

     if( plot == "forestplot" ){
         for(Nmed in 1:length(mediator)){
               cochrane_from_rmeta[[Nmed]] <- structure(list(
                      mean  = c(NA, results_mediator_table[[Nmed]][1,1], results_mediator_table[[Nmed]][1,2], results_mediator_table[[Nmed]][1,3], results_mediator_table[[Nmed]][1,4], results_mediator_table[[Nmed]][1,5] ),
                      lower = c(NA, results_mediator_table[[Nmed]][2,1], results_mediator_table[[Nmed]][2,2], results_mediator_table[[Nmed]][2,3], results_mediator_table[[Nmed]][2,4], results_mediator_table[[Nmed]][2,5]),
                      upper = c(NA, results_mediator_table[[Nmed]][3,1], results_mediator_table[[Nmed]][3,2], results_mediator_table[[Nmed]][3,3], results_mediator_table[[Nmed]][3,4], results_mediator_table[[Nmed]][3,5])),
                      .Names = c("mean", "lower 95% C.I.", "upper 95% C.I."),
                      row.names = c(NA, -6L),
                      class = "data.frame")



                    tabletext[[Nmed]]<-cbind(
                    c("", "overall", "direct Sjolander", "indirect Sjolander","path specific", "overall Direct"),
                    c("Mean", paste(round(results_mediator_table[[Nmed]][1,1],4)) , paste(round(results_mediator_table[[Nmed]][1,2],4)), paste(round(results_mediator_table[[Nmed]][1,3],4)), paste(round(results_mediator_table[[Nmed]][1,4],4)) ,paste(round(results_mediator_table[[Nmed]][1,5],4)) ),
                    c("lower 95% C.I.", paste(round(results_mediator_table[[Nmed]][2,1],4)) , paste(round(results_mediator_table[[Nmed]][2,2],4)), paste(round(results_mediator_table[[Nmed]][2,3],4)), paste(round(results_mediator_table[[Nmed]][2,4],4)), paste(round(results_mediator_table[[Nmed]][2,5],4)) ),
                    c("Upper 95% C.I.", paste(round(results_mediator_table[[Nmed]][3,1],4)), paste(round(results_mediator_table[[Nmed]][3,2],4)), paste(round(results_mediator_table[[Nmed]][3,3],4)), paste(round(results_mediator_table[[Nmed]][3,4],4)), paste(round(results_mediator_table[[Nmed]][3,5],4))))


                     plotList[[Nmed]] <- forestplot(tabletext[[Nmed]],
                                cochrane_from_rmeta[[Nmed]],new_page = TRUE,
                                is.summary=c(TRUE,rep(FALSE,5)),
                                clip=c(min(results_mediator_table[[Nmed]]) -0.1,max(results_mediator_table[[Nmed]]) +0.1),
                                xlog=FALSE,
                                col=fpColors(box= boxCol,line = lineCol, summary = summaryCol))
         }
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
                      name=c("overall", "direct Sjolander", "indirect Sjolander","path specific", "overall Direct"),
                      mean = results_mediator_table[[Nmed]][1,],
                      lower95CI = results_mediator_table[[Nmed]][2,],
                      upper95CI = results_mediator_table[[Nmed]][3,] )

                      plotList[[Nmed]] <-  ggplot(data[[Nmed]]) +
                      geom_bar( aes(x = name, y = mean), stat="identity", fill = fill, alpha=0.7) +
                      geom_errorbar( aes(x=name, ymin=lower95CI, ymax = upper95CI), width=0.4, colour = colour, alpha=0.9, size=1.3) +
                      coord_flip()
                      }
     } else{ print("Enter value for plot argument.") }


     plotList

}
