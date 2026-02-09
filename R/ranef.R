#' Extract posterior means and SDs for "random" effects from a bayz model output
#'
#' The ranef() (fully ranef.bayz()) function extracts estimates (posterior mean
#' and SD) for all
#' random effects from the bayz output. Random effects are here defined as the
#' coefficients for the rn() and rr() model terms.
#' The ranef() function is a wrapper around estim() and 'unlists' the output to
#' return a single data frame.
#' To extract estimates for one, or a subset, of the random effects, use the
#' more generic estim(). Use fixef() to extract estimates for all fixed
#' effects. ranef.bayz() is an implementation of the generic R ranef()
#' function.
#'
#' @param object        A bayz model output
#' @param ...           Additional parameters passed onto the Model function.
#' @return a data frame with estimates for all random effects
#' @importFrom nlme ranef
#' @export
#'
ranef.bayz <- function(object, ...){
  random_model_terms <-
    grepl("rr|rn", object$Parameters$ModelTerm)
  random_parameters <- object$Parameters$Param[random_model_terms]
  estim(object, param = random_parameters, splitlabels = FALSE,
        unlist = TRUE)
}
