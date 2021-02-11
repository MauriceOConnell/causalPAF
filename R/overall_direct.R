#' @title Calculation of Population Attributable Fraction (PAF), with a decomposition of the total PAF into direct
#' and indirect components.
#' @description Total PAF
#' @param data_frame Data frame containing the data. The data frame has cases in rows and variables in columns.
#' @param exposure The exposure name in the form of character string e.g. "phys" for physical exercise.
#' @param mediator The mediator name in the form of character string e.g. "whr" for waist-hip ratio.
#' @param response The outcome name in the form of character string e.g. "case" for a stroke case.
#' @param mediator_model A list containing each of the mediator regression models e.g. \code{mediator_model=list(model_list[[6]],model_list[[7]],model_list[[8]])}.
#' @param response_model is a regression model fitted for the outcome on all mediators together with all parents and confounders of the mediators in
#' a Markov causal Bayesian network DAG.
#' @param response_model_2 is a regression model fitted for the outcome
#' @param weights A numeric
#' @export
#' @import stats
# #' @importsFrom stats formula predict filter lag
#' @keywords models Regression
#' @return \item{directPAF}{direct PAF}
#' @examples \dontrun{
#' # I don't want you to run this
#' }
overall_direct <- function(data_frame,exposure,mediator,response,mediator_model,response_model,response_model_2,weights){
  mediator_outcomes <- c()
  for(i in 1:length(mediator_model)) mediator_outcomes[i] <-  as.character(formula(mediator_model[[i]]))[2]
  which.model <- grep(paste('^',mediator,'$',sep=''),mediator_outcomes,perl=TRUE)
  data_frame_direct=data_frame
  data_frame_direct[,grep(paste('^',exposure,'$',sep=''),colnames(data_frame_direct),perl=TRUE)] <- levels(data_frame_direct[,grep(paste('^',exposure,'$',sep=''),colnames(data_frame),perl=TRUE)])[1]
  predicted_response <- predict(response_model,newdata=data_frame_direct,type="response")

  # browser()
  directPAF <- sum(weights*(predict(response_model,type="response")-predicted_response))/sum(weights*predict(response_model,type="response"))

  return(directPAF)

}
