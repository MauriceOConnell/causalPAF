#' @title Path specific population attributable fraction for a mediating pathway.
#' @description Path specific population attributable fraction for a mediating pathway.
#' This metric has several desirable properties. This is a kind of sequential PAF on pathways,
#' but now corresponding to eliminating the mediating pathway first.
#' \deqn{PAF_{A->M->Y} = (P(Y=1)-P(Y_{A,M_{0}}=1)/P(Y=1) }
#' \if{html}{\figure{CausalDAG.jpg} {options: width=100 alt="Causal Bayesian network DAG"} }
#' \if{latex}{\figure{CausalDAG.jpg}{options: width=1in}}
#' @param data_frame Data frame containing the data. The data frame has cases in rows and variables in columns.
#' @param exposure The exposure name in the form of character string e.g. "phys" for physical exercise.
#' @param mediator The mediator name in the form of character string e.g. "whr" for waist-hip ratio.
#' @param response The outcome name in the form of character string e.g. "case" for a stroke case.
#' @param mediator_model A list containing each of the fitted mediator regression models e.g.
#' \code{mediator_model=list(model_list[[6]],model_list[[7]],model_list[[8]])}.
#' @param response_model is a regression for the outcome on all mediators together with all parents and confounders of the mediators in a Markov causal Bayesian network DAG. A regression model fitted for the response in a causal Bayesian network excluding ``children'' of the mediators in the causal Bayesian network. See example in tutorial.This model can be listed either as (1) an empty list ( response_model_mediators = list() ) or (2) the user can specify their own customised causal regression model(s) to use. When it is listed as an empty list the causalPAF package will fit the response_model_mediators regression model automatically based on the causal DAG supplied by the user in in_outArg. Alternatively, the user can specify the exact model(s) that the user wishes to use, these model(s) must be in list format (list() where length(response_model_mediators) == length(mediator) ), the same length as the parameter, mediator, with the user customised model for each mediator listed in the same order as in the parameter, mediator, and if there is only one model, it must be listed each time within the list() so that length(response_model_mediators) == length(mediator).
#' @param response_model_2 A regression model fitted for the response in a causal Bayesian network excluding ``children'' of the exposure in the causal Bayesian network. This regression model will not adjust for mediators (exclude mediators) of the exposure in the regression model so that the total effect of the exposure on the response can be modelled. This model can be listed either as (1) an empty list ( response_model_exposure = list() ) or (2) the user can specify their own customised causal regression model to use. If specified as an empty list, list(), then the causalPAF function will define and fit the model automatically based on the causal DAG defined by the in_outArg parameter. Alternatively, the user can specify the exact model that the user wishes to use, this model must be in list format (list() where length(response_model_exposure) == 1 ), of length 1, assuming only one exposure of interest (other exposures can be risk factors) and the model must be defined within a list() since the package assumes a list() format is supplied. See example in tutorial. E.G. If physical exercise ("exer") in the example given in the diagram is the exposure. Then the regression would include all parents of "exer" (i.e. sex, region, educ, age) as well as risk factors at the same level of the causal Bayesian network (i.e. stress, smoke, diet, alcoh).
#' @param weights A numeric \eqn{n x 1} vector where n is the number of patients in the case control data frame.
#' For case control studies, a reweighting approach is used which assumes the prevalence of disease,
#' \eqn{\pi} is known, and the sampled disease cases and controls are randomly selected from their respective populations.
#' We assume for simplicity that the case:control matching ratio is \eqn{1:r}, for some \eqn{r \ge 1}.
#' Under these assumptions, the components of the PAF can be found as the corresponding empirical expectations
#' and distributions in the re-weighted dataset where cases are assigned weights 1, and controls are assigned weights
#'  \eqn{(1\pi−1)/r}. Effectively then we can think of the reweighted population as a random sample.
#' @export
#' @import stats
# #' @importFrom stats formula predict filter lag
#' @keywords internal
#' @return \item{path_specific_PAF }{path specific PAF}

path_specific_onesimulation <- function(data_frame, exposure, mediator, response, mediator_model, response_model, response_model_2, weights ){

  mediator_outcomes <- c()
  for(i in 1:length(mediator_model)) mediator_outcomes[i] <-  as.character(formula(mediator_model[[i]]))[2]
  which.model <- grep(paste('^',mediator,'$',sep=''),mediator_outcomes,perl=TRUE)
  data_frame_pathspecific=data_frame
  data_frame_pathspecific[,grep(paste('^',exposure,'$',sep=''),colnames(data_frame_pathspecific),perl=TRUE)] <- levels(data_frame_pathspecific[,grep(paste('^',exposure,'$',sep=''),colnames(data_frame),perl=TRUE)])[1]
  for(i in which.model){

    ### simulate mediators given exposure at reference (except for mediator in question)
    thecol <- grep(paste('^',mediator_outcomes[i],'$',sep=''),colnames(data_frame_pathspecific),perl=TRUE)
    data_frame_pathspecific <- do_sim(mediator_model[[i]],data_frame_pathspecific)

  }
  #  make sure exposure is set at natural value
  data_frame_pathspecific[,grep(paste('^',exposure,'$',sep=''),colnames(data_frame_pathspecific),perl=TRUE)] <- data_frame[,grep(paste('^',exposure,'$',sep=''),colnames(data_frame_pathspecific),perl=TRUE)]
  predicted_response <- predict(response_model,newdata=data_frame_pathspecific,type="response")

  path_specific_PAF <- sum(weights*(predict(response_model,type="response")-predicted_response))/sum(weights*predict(response_model,type="response"))

  return(path_specific_PAF)
}
