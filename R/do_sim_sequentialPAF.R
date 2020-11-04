#' @title Simulates a Fitted Model for a Mediator or Exposure or Risk Factor Allowing for Potential Outcomes in Causal Analysis
#' @description A fitted model for a mediator or exposure or risk factor can be simulated given values of the other risk
#' factors or exposure saved in the data frame \code{current_mat}. This allows for potential outcomes to be measured for
#' causal analysis. For example, for an outcome \eqn{Y_{A,M}} with exposure A and mediators \eqn{M_{1}, M_{3}, \dots M_{K}}
#' the function can measure potential outcomes such as \eqn{Y_{A=0,M_{1},M_{2},M_{3}}} or \eqn{Y_{A=0,M_{1},M_{2}=0,M_{3}=0}} when there are three mediators.
#' The model can be either a binary, continuous or an ordered factor response model.
#' @param colnum Column number of expsoure or risk factor of interest within the data frame. The data frame has cases in rows and variables in columns.
#' @param current_mat The data frame containing the data for which the model can be simulated with. For
#' potential outcomes for example such as \eqn{Y_{A=0,M_{1},M_{2},M_{3}}} requires the exposure in this case
#' to be pre set to zero i.e. \code{current_mat} should have the exposure \eqn{Y_{A=0}} set to zero if simulating
#' e.g. \eqn{M_{1}}.
#' @param model A fitted causal regression model for either a binary, continuous or an ordered factor response.
#' @export
#' @import stats
# #' @importFrom stats formula predict filter lag
#' @keywords models Regression
#' @return \item{simulation }{simulation}
#' @examples \dontrun{
#' # I don't want you to run this
#' }
do_sim_sequentialPAF <- function(colnum,current_mat, model){
        ## polr
        if(names(model)[2]=='zeta'){

                probs <- predict(model,newdata=current_mat,type="probs")
                mynames <- colnames(probs)
                simulation <- apply(probs,1,function(x){base::sample(mynames,size=1,prob=x)})
                return(simulation)
        }
        # glm
        if(grep("glm",model$call)){

                probs <- predict(model,newdata=current_mat,type="response")
                if(is.null(levels(current_mat[,colnum]))) return(apply(cbind(1-probs,probs),1,function(x){base::sample(c(0,1),size=1,prob=x)}))
                simulation <-  apply(cbind(1-probs,probs),1,function(x){base::sample(levels(current_mat[,colnum]),size=1,prob=x)})
                return(simulation)
        }
        # regression
        if(grep("lm",model$call)){

                pred <- predict(model,newdata=current_mat,type="response")
                s_d <- sd(model$residuals)
                simulation <-  pred + rnorm(length(pred),mean=0,sd=s_d)
                return(simulation)
        }
}
