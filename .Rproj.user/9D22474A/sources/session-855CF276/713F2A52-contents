#' @title Makes Formula
#' @description Given a causal
#' @param model A fitted mediator regression model
#' @param dataframe Data frame containing data to be analysed. The data frame has cases in rows and variables in columns.
#' @export
#' @import stats
# #' @importsFrom stats formula predict filter lag
#' @keywords internal
#' @return \item{dataframe }{dataframe}
#' @examples \dontrun{
#' # I don't want you to run this
#' }

do_sim <- function(model,dataframe){

  y_name <-  as.character(formula(model)[2])
  y <- dataframe[,colnames(dataframe)==y_name]
  if(length(table(y))==2){
    probs <- predict(model,newdata=dataframe,type='response')
    predictions <- levels(y)[1+rbinom(n=length(probs),size=1,prob=probs)]
    #browser()
    dataframe[,colnames(dataframe)==y_name] <- predictions
    return(dataframe)
  }
  if(length(table(y))>2 & is.factor(y)){
    mediator_probs <- predict(model,newdata=dataframe,type='probs')
    a <- apply(mediator_probs,1,function(x){sample(1:ncol(mediator_probs), size=1, prob=x)})
    predictions <-  levels(y)[a]
    dataframe[,colnames(dataframe)==y_name] <- predictions
    return(dataframe)
  }
  if(length(table(y))>2 & is.numeric(y)){
    mean <- predict(model,newdata=dataframe)
    dataframe[,colnames(dataframe)==y_name] <- mean + sample(summary(model)$resid,size=length(mean),replace=TRUE)
    return(dataframe)
  }
}



