#' Compute pr-values (p-values from random effects) from a bayz model output
#'
#' Rbayz::prval() evaluates significances from random effects based on
#' procedures from Gualdron-Duarte et al. (2013)
#' and Bernal-Rubio et al. (2016), extended to allow use in a full Bayesian
#' framework where the variance hyper-parameters are also learned from the data.
#' The G-D & B-R procedure applies to fitting high-dimensional predictor sets
#' as random effects, as done in bayz with rr(),
#' and evaluates the significance of individual predictors
#' as if they were fitted one by one in turn as fixed in a linear mixed model
#' with all other effects remaining random.
#' The reported pr-value is based on a z-test.
#' Note that the model-term for which pr-values are to be computed should have
#' set a 'save' or 'prval' option.
#'
#' @param object        A bayz model output
#' @param param         Parameter name for which to extract r-values (if in
#'                      doubt about naming, check $Parameters from the output)
#' @param splitlabels   If the estimates come from an interaction, whether
#'                      labels should be split to allow easy matching to levels
#'                      in the single variables (default FALSE).
#' @param ...           Additional parameters.
#' @return a data frame with parameter estimates (posterior mean and SD), the
#'                      backtransformed z-statistic and its pr (p) -value.
#' @export
#'
prval <- function(object, param = NULL, splitlabels = FALSE, ...) {
  par = object$Parameters
    est = object$Estimates
    output_object = data.frame()

    for(i in 1:nrow(par)) {
        estim_select = est[par$EstStart[i]:par$EstEnd[i],]
        if (!grepl("%",rownames(estim_select)[1],fixed=TRUE)) splitlabels = FALSE
        if(splitlabels) {
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
