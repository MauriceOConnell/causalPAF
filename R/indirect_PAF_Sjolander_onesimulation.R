#' @title Calculation of Population Attributable Fraction (PAF), with a decomposition of the total PAF into direct
#' and indirect components.
#' @description Calculation of Population Attributable Fraction (PAF), with a decomposition of the total PAF into direct
#' and indirect components. It performs one simulation which can be combined with a bootstrap approach to
#' perform multiple simulations. If we think of \eqn{Y_0} as the potential outcome for
#' an individual if they were never exposed to the risk factor, can be directly interpreted as the relative change in
#' disease prevalence if an exposure was absent from the population.
#' Sjolander introduced the ideas of mediation into the literature for PAF, defining a decomposition of the total PAF
#' into direct and indirect components:  and
#' \if{html}{\figure{CausalDAG.jpg} {options: width=100 alt="Causal Bayesian network DAG"} }
#' \if{latex}{\figure{CausalDAG.jpg}{options: width=1in}}
#' @param data_frame Data frame containing the data. The data frame has cases in rows and variables in columns.
#' @param exposure The exposure name in the form of character string e.g. ``phys'' for physical exercise.
#' @param mediator The mediator name in the form of character string e.g. ``whr'' for waist hip ratio.
#' @param response The outcome name in the form of character string e.g. ``case'' for a stroke case.
#' @param mediator_model A list containing each of the fitted mediator regression models e.g.
#' \code{mediator_model=list(model_list[[6]],model_list[[7]],model_list[[8]])}.
#' @param response_model is a regression for the outcome on all mediators together with all parents and confounders of the mediators in
#' a Markov causal Bayesian network DAG e.g.
#' @param response_model_2 is a regression for the outcome on the exposure together with  all parents and confounders of the exposure in
#' a Markov causal Bayesian network DAG along with other risk factors at the same level of the causal Bayesian network DAG. E.G. If
#' physical exercise (``exer'') in the example given in the diagram is the exposure. Then the regression would include all
#' parents of ``exer'' (i.e. sex, region, educ, age) as well as risk factors at the same level of the causal Bayesian network
#' (i.e. stress, smoke, diet, alcoh).
#' @param weights A numeric n x \eqn{1} vector where n is the number of patients in the case control data frame.
#' Different weighting approaches can be applied as per the literature, Pathway-specific population attributable
#' fractions (PS-PAF) O'Connell and Ferguson (\eqn{2022}), Revisiting sequential population attributable fractions Ferguson, O'Connell and O'Donnell (\eqn{2020}).
#' For more information on weighting, a tutorial paper will be published and linked here when it is published.
#' For example, O'Connell and Ferguson (\eqn{2022}), for a case-control study design, when prevalence pi is known, and
#' the sampled disease cases and controls are randomly selected from their respective populations. We assume for
#' simplicity that the case to control matching ration is \eqn{1} to r, for some r greater than or equal to \eqn{1}.
#' Under assumptions outlined in O'Connell and Ferguson (\eqn{2022}), the components of the PAF can be found as the
#' corresponding empirical expectations and distributions in the re-weighted dataset where cases are assigned the weights,
#' \eqn{1}, and controls are assigned the weights of (((\eqn{1} divided by pi) minus \eqn{1}) all divided by r). Effectively, then we can think of the re-weighted population as a random
#' sample. A tutorial paper will be linked here when published for more information which will show how to apply the weighting for different study designs.
#' @export
#' @import stats
#' @keywords internal
#' @return \item{totalPAF }{total PAF}
#' \item{directPAF}{direct PAF}
#' \item{indirectPAF}{indirect PAF}

indirect_PAF_Sjolander_onesimulation <- function(data_frame, exposure, mediator,response,mediator_model,response_model,response_model_2,weights){
        mediator_outcomes <- c()
        for(i in 1:length(mediator_model)) mediator_outcomes[i] <-  as.character(formula(mediator_model[[i]]))[2]
        which.model <- grep(paste('^',mediator,'$',sep=''),mediator_outcomes,perl=TRUE)
        data_frame_direct=data_frame
        data_frame_direct[,grep(paste('^',exposure,'$',sep=''),colnames(data_frame_direct),perl=TRUE)] <- levels(data_frame_direct[,grep(paste('^',exposure,'$',sep=''),colnames(data_frame),perl=TRUE)])[1]
        for(i in setdiff(1:length(mediator_model),which.model)){
          #browser()
          ### simulate mediators given exposure at reference (except for mediator in question)
          thecol <- grep(paste('^',mediator_outcomes[i],'$',sep=''),colnames(data_frame_direct),perl=TRUE)
          data_frame_direct <- do_sim(mediator_model[[i]],data_frame_direct)

        }
        predicted_response <- predict(response_model,newdata=data_frame_direct,type="response")

       # browser()
        directPAF <- sum(weights*(predict(response_model,type="response")-predicted_response))/sum(weights*predict(response_model,type="response"))


        ## calculate totalPAF
        data_frame_total <- data_frame
        data_frame_total[,grep(paste('^',exposure,'$',sep=''),colnames(data_frame),perl=TRUE)] <- levels(data_frame[,grep(paste('^',exposure,'$',sep=''),colnames(data_frame),perl=TRUE)])[1]

        totalPAF <- sum(weights*(predict(response_model_2,type="response")-predict(response_model_2,newdata=data_frame_total,type="response")))/sum(weights*predict(response_model_2,type="response"))
        return(c(totalPAF=totalPAF,directPAF=directPAF,indirectPAF=totalPAF-directPAF))

}


