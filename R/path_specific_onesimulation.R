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
#' @param response_model is a regression for the outcome on all mediators together with all parents and confounders of the mediators in
#' a Markov causal Bayesian network DAG e.g.
#' @param response_model_2 is a regression for the outcome on the exposure together with  all parents and confounders of the exposure in
#' a Markov causal Bayesian network DAG along with other risk factors at the same level of the causal Bayesian network DAG. E.G. If
#' physical exercise ("exer") in the example given in the diagram is the exposure. Then the regression would include all
#' parents of "exer" (i.e. sex, region, educ, age) as well as risk factors at the same level of the causal Bayesian network
#' (i.e. stress, smoke, diet, alcoh).
#' @param weights A numeric \eqn{n x 1} vector where n is the number of patients in the case control data frame.
#' For case control studies, a reweighting approach is used which assumes the prevalence of disease,
#' \eqn{\pi} is known, and the sampled disease cases and controls are randomly selected from their respective populations.
#' We assume for simplicity that the case:control matching ratio is \eqn{1:r}, for some \eqn{r \ge 1}.
#' Under these assumptions, the components of the PAF can be found as the corresponding empirical expectations
#' and distributions in the re-weighted dataset where cases are assigned weights 1, and controls are assigned weights
#'  \eqn{(1\pi−1)/r}. Effectively then we can think of the reweighted population as a random sample.
#' @export
#' @import stats dplyr
# #' @importFrom stats formula predict filter lag
#' @keywords models Regression
#' @return \item{path_specific_PAF }{path specific PAF}
#' @examples \dontrun{
#' # I don't want you to run this
#' }
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

  # browser()
  path_specific_PAF <- sum(weights*(predict(response_model,type="response")-predicted_response))/sum(weights*predict(response_model,type="response"))

  return(path_specific_PAF)
}
