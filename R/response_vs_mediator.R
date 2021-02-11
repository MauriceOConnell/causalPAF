#' @title Checks if the causal DAG satisfies the Markov condition
#' @description The functions checks if the Markov condition holds for the Directed Acyclic Graph (DAG) defined. Sometimes called the Markov assumption, is an assumption made in Bayesian probability theory, that every node in a Bayesian network is conditionally independent of its nondescendents, given its parents. In other words, it is assumed that a node has no bearing on nodes which do not descend from it. This is equivalent to stating that a node is conditionally independent of the entire network, given its Markov blanket. The related Causal Markov condition states that, conditional on the set of all its direct causes, a node is independent of all variables which are not direct causes or direct effects of that node.
#' @param data data stored as a dataframe with columns as variables and rows as cases and controls in case control study e.g. data = strokedata.
#' @param exposure exposure name in text format e.g. exposure="phys"
#' @param mediator vector containing names of mediators in text format e.g. mediator=c("subhtn","apob_apoa","whr")
#' @param response response name in text format e.g. response="case"
#' @param in_out A list of length 2. The first list contains a list of character vectors of the parents of the exposure or risk factor or outcome which are either causes or confounders of the exposure or risk factor or outcome. The second list conttains a list of a single name of exposure or risk factor or outcome in form of characters. See tutorial examples for examples.
#' @param w Column of weights for case control matching listing in same order as patients in data e.g. w = strokedata$weights
#' @export
#' @import splines MASS stats utils grid magrittr checkmate rlist
#' @keywords models Markov Bayesian Directed Acyclic Graph Population Attributable Fraction
#' @return  \item{to_execute }{to_execute is in text format and represents the response_vs_exposure model to be evaluated including all covariates below and at the same level of the exposure in the causal Bayesian DAG.}
#' @examples \dontrun{
#' # I don't want you to run this
#' }
#' in_vars = c("subeduc","moteduc","fatduc")
#' outvar = c("phys")
#' make_formula(in_vars,outvar)
response_vs_mediator <- function(data, exposure="phys",mediator=c("subhtn","apob_apoa","whr"), response="case",in_out, w ){
# data = strokedata
# w = strokedata$weights

# ONLY MAKES SENSE IF MEDIATORS ARE AT THE SAME LEVEL? CONSIDER IF TRUE

  # what to create this from within causalPAFplot
  # response_vs_phys <- glm("case ~ regionnn7*ns(eage,df=5)+esex*ns(eage,df=5) +  subeduc+ moteduc+ fatduc+ phys+ ahei3tert+ nevfcur+ alcohfreqwk+ global_stress2",data = stroke_reduced,family='binomial',w= stroke_reduced$weights)

# Create function for this also
# response_vs_mediator <-  glm("case ~ regionnn7*ns(eage,df=5)+esex*ns(eage,df=5) +  subeduc+ moteduc+ fatduc+ phys+ ahei3tert+ nevfcur+ alcohfreqwk+ global_stress2+ subhtn + ns(apob_apoa, knots = quantile(apob_apoa,c(.25,0.5,0.75)), Boundary.knots = quantile(apob_apoa,c(.001,0.95)))+ns(whr,df=5)",data = stroke_reduced,family='binomial',w = stroke_reduced$weights )


  # mediator=c("subhtn","apob_apoa","whr")

  findMax <- vector()
  for(i in 1:length(mediator)){
      findMax[i] <- grep(mediator[i], in_out[[2]])
  }
  MaxMediatorIndex <- max(findMax)

  # 1. Include all variables below exposureIndex
  # exposureIndex <- grep(exposure, in_out[[2]])
  # Try this
  exposureIndex <-  MaxMediatorIndex

  Subset1 <- in_out[[1]][1:exposureIndex]
  Subset2 <- in_out[[2]][1:exposureIndex]

  # 2. Include all variable above exposureIndex that are at the same level i.e. combine all variables below exposureIndex
  # and then see which ones above are a subset of it them which includes the ``full'' set.

  CombineSubset1 <- unique(unlist(Subset1, use.names=FALSE))

  if(exposureIndex < length(in_out[[2]])){

        indices <- (exposureIndex + 1):(length( in_out[[2]]))

        # This assumes that e.g. if splines of variables that they are defined in same way before and after exposureIndex, if different it will not work
        indicesAddIn <- which(lapply(in_out[[1]][indices], function(data_input) all(data_input %in% CombineSubset1 )  ) > 0 )

        if( length(indicesAddIn) != 0 ){
           indicesAddIn <- indicesAddIn  + exposureIndex

           Subset1 <- in_out[[1]][c(1:exposureIndex, indicesAddIn) ]
           Subset2 <- in_out[[2]][c(1:exposureIndex, indicesAddIn) ]
        }

  }

  # ISSUE HERE IF variable above to be included but need to be in spline format
  # c(unique(unlist( Subset1, use.names=FALSE)), unique(unlist( Subset2, use.names=FALSE)) )

  #3.Find index of response
  # responseIndex <- grep(response, in_out[[2]])
  # Assumes respones appears as variable only with no e.g. spline etc which is expect in in_out[[2]]
  responseIndex <- grep(paste('^',response,'$',sep=''),in_out[[2]],perl=TRUE)


  #4.Find in_out[[1]][responseIndex] and remove indices not in indicesAddIn

  # ########################
  # #########################
  # # to be called in causalPAFplot.R
  # #########################
  # response_model_exposure <- response_vs_exposure(data = stroke_reduced ,
  #                                                 exposure="phys",
  #                                                 response="case",
  #                                                 in_out = in_out,
  #                                                 w = stroke_reduced$weights )
  # #########################
  # #########################

  column <- responseIndex

  data_text <- deparse(substitute(data))

  w_text <- deparse(substitute(w))

  # model_listReturn <- vector(mode = "list", length = length(in_out[[2]]) )

            # for(i in 1:length(in_out[[2]]) ){

                          # column <- (1:length(colnames(data)))[colnames(data) %in% in_out[[2]][i]]
                            column <- (1:length(colnames(data)))[colnames(data) %in% in_out[[2]][responseIndex]]
                          # formula_text <- make_formula(in_out[[1]][[i]],in_out[[2]][i], addCustom = TRUE, custom = "~ regionnn7*ns(eage,df=5)+esex*ns(eage,df=5) + ")
                          formula_text <- make_formula( c(unique(unlist( Subset1, use.names=FALSE)), unique(unlist( Subset2, use.names=FALSE)) ),
                                    in_out[[2]][responseIndex],
                                    addCustom = TRUE,
                                    custom = "~ regionnn7*ns(eage,df=5)+esex*ns(eage,df=5) + ")
                          y <- data[,column]
                          if(length(table(y))==2){
                                  theform <- paste("glm(",formula_text,",data=",data_text,",family='binomial',w=",w_text,")",sep='')
                          }
                          if(length(table(y))>2 & is.factor(y)){
                                  theform <- paste("polr(",formula_text,",data=",data_text,",w=",w_text,")",sep='')
                          }
                          if(length(table(y))>2 & is.numeric(y)){
                                  theform <- paste("lm(",formula_text,",data=",data_text,",w=",w_text,")",sep='')
                          }
                           to_execute <- paste("response_model_exposure_text <-", theform,sep='')
                          # model_listReturn[[i]] <- to_execute
                          # # eval(parse(text=to_execute))
                          # ## model_listReturn[[i]] <- eval(parse(text=to_execute))
                          # #model_listReturn[[i]] <- eval( parse(text= paste(model_list_text,"[[i]]",sep='') ) )
                          #
                          # ## E.G. TRY model$terms <- eval(model$call$formula)
                          # ##model_listReturn[[i]]$terms <- eval(model_listReturn[[i]]$call$formula)
            # }



  # ########################
  # #########################
  # # to be called in causalPAFplot.R
  # #########################
  # response_model_exposure_text <- response_vs_exposure(data = stroke_reduced ,
  #                                                      exposure="phys",
  #                                                      response="case",
  #                                                      in_out = in_out,
  #                                                      w = stroke_reduced$weights )
  # # In meantime not running function use this...then delete
  # response_model_exposure_text <- to_execute
  # eval(parse(text = response_model_exposure_text ) )
  # return(response_model_exposure_text)
  #########################
  #########################
  #########################

  return(to_execute)


}





