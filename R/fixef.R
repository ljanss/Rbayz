#' Extract posterior means and SDs for "fixed" effects from a bayz model output
#'
#' The fixef() (fully fixef.bayz()) function extracts estimates (posterior mean
#' and SD) for all
#' fixed effects from the bayz output. Fixed effects are here defined as the
#' coefficients for the intercept and fx() and rg() model terms.
#' The fixef() function is a wrapper around estim() and 'unlists' the output to
#' return a single data frame.
#' To extract estimates for one, or a subset, of the fixed effects, use the
#' more generic estim(). Use ranef() to extract estimates for all random
#' effects. fixef.bayz() is an implementation of the generic R fixef()
#' function.
#'
#' @param object        A bayz model output
#' @param ...           Additional parameters passed onto the Model function.
#' @return a data frame with estimates for all fixed effects
#' @importFrom nlme fixef
#' @export
#' 
fixef.bayz <- function(object, ...){
  fixed_model_terms <-
    grepl("mn|fx|rg", object$Parameters$ModelTerm)
  fixed_parameters <- object$Parameters$Param[fixed_model_terms]
  estim(object, param = fixed_parameters, splitlabels = FALSE,
        unlist = TRUE)
}
