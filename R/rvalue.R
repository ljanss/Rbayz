#' Compute r-values (p-values obtained from random effects) from a bayz model output
#'
#' Rbayz::rvalue() computes significances from random effects based on procedures from Gualdron-Duarte et al. (2013)
#' and Bernal-Rubio et al. (2016) in the Rbayz Bayesian framework. The frequentist G-D & B-R approach is extended to
#' also obtain these r-values from a full Bayesian model fit that has also learned the variance hyper-parameters 
#' (see our publications for full background). The r-value is a backtransformation based on posterior meand and variance to
#' a p-value as if every effect in turn was fitted as fixed effect in a inear mixed model
#' with all other effects remaining random. The Bayesian model can extract these r-values from a single model fit, while
#' the mixed model approach would require refitting the model for every effect. 
#' The reported significance is based on a z-test.
#' Note that the model-term for which r-values are to be computed should have set a 'save' or 'rvalue' option.
#'
#' @param object        A bayz model output
#' @param param         Parameter name for which to extract r-values (if in doubt about naming, check $Parameters from the output)
#' @param splitLabels   Whether labels that contain % (in interactions) such as a%b should be split
#'                      in multiple columns (default TRUE).
#' @param ...           Additional parameters.
#'
#' @return a data frame with parameter estimates (posterior mean and SD), the backtransfored z-statistic and its r (or p)-value.
#' @export
rvalue <- function(object, param=NULL, splitLabels=TRUE, ...){
    par = object$Parameters
    est = object$Estimates
    output_object = data.frame()

    for(i in 1:nrow(par)) {
        estim_select = est[par$EstStart[i]:par$EstEnd[i],]
        if (!grepl("%",rownames(estim_select)[1],fixed=TRUE)) splitLabels = FALSE
        if(splitLabels) {
            splitlabels = t(as.data.frame(strsplit(rownames(estim_select),"%")))
            colnm = splitlabels[1,1]
            splitlabels = splitlabels[,-1]
            estim_table = data.frame(splitlabels, estim_select)
            colnames(estim_table) = c(unlist(strsplit(colnm,":")),colnames(estim_select))
            rownames(estim_table) = rownames(estim_select)
        }
        else {
            estim_table = estim_select
        }
        output_object[[rownames(par)[i]]] = estim_table
    }
    return(output_object)
}
