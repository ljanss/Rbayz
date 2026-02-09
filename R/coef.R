#' Extract posterior means and SDs for a selection of model parameters from a bayz model output
#'
#' Rbayz::coef() extracts coefficients (posterior mean and SD of fixed or random effects) from the Bayz output for one or more model-terms.
#' Common use is to extract estimates from a single model-term, for example using param="year:location" to retrieve the estimates for a
#' modelled year:location interaction, but param can also take a vector of strings to combine estimates from several model-terms in one output.
#' To extract all coefficients, use fixef() and ranef() functions to retrieve all fixed and random effects, respectively.
#' For interactions, and when retrieving estimates from a single model-term, 
#' the labels of the coefficients are by default split in multiple columns; this can be toggled off setting splitLabels=FALSE.
#' Note that coef() with a single model-term is nearly equivalent to directly accessing $Estimates$`model-term` from the output, but with the default feature
#' to splits labels from interaction terms. For multiple model-terms, coef() formats all estimates in a single data.table, while it is stored in 
#' separate data.tables in the output $Estimates. 
#' 
#' @param object        A bayz model output
#' @param param         A string, or vector of strings, to indicate for which parameter(s) estimates should be extracted.
#' @param splitLabels   Whether labels in interaction terms that are pasted together should be split
#'                      in multiple columns. This is default TRUE when retrieving estimates for a single parameter, but cannot be done
#'                      when retrieving estimates from multiple sets of parameters.
#' @param ...           Additional parameters.
#'
#' @return A data.table object with labels or split labels, posterior mean and SD, and for multiple model-terms a column ModelTerm.
#' @import stats
#' @export
#'
coef.bayz <- function(bayz_output, param=NULL, splitLabels=TRUE, ...){
    par = bayz_output$Parameters
    est = bayz_output$Estimates
    if(is.null(param)) {
        stop("param= must be specified. To retrieve all fixed or random effects use fixef() and ranef().")
    }
    if(length(param)==1) {
        search_term = param
    }
    else {
        search_term = param[1]
    }
    if( !(search_term %in% names(est))) {
        stop(paste0("Parameter '",search_term,"' not found in bayz output. Check the Param column in the output$Parameters for available parameter names."))
    }
    return_object = est[[search_term]]
    if(length(param)>1) {
        for(i in 2:length(param)) {
            search_term = param[i]
            if( !(search_term %in% names(est))) { 
                stop(paste0("Parameter '",search_term,"' not found in bayz output. Check the Param column in the output$Parameters for available parameter names."))
            }
            return_object = rbind(return_object, est[[search_term]])
    par = par[which,]
    new_object = list()
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
        new_object[[rownames(par)[i]]] = estim_table
    }
    return(new_object)
}
