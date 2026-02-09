#' Extract posterior means and SDs for all coefficients from a bayz model
#' output
#'
#' Rbayz::coef() extracts estimates (posterior mean and SD) for all
#' coefficients from the bayz output. Coefficients is here defined as all fixed
#' and random effects, i.e. "regression coefficients" in a broad meaning.
#' The coef() function is a wrapper around estim() and 'unlists' the output to
#' return a single data frame.
#' Other functions to extract estimates from the bayz output are fixef(),
#' ranef(), var(), and the generic estim().
#'
#' @param object        A bayz model output
#' @param ...           Additional parameters.
#' @return A data.frame with estimates for all model coefficients
#' @export
#'
coef.bayz <- function(bayz_output, ...) {
  coef_model_terms <-
    grepl("mn|fx|rg|rr|rn", bayz_output$Parameters$ModelTerm) &
    (bayz_output$Parameters$Variance == "-")
  coef_parameters <- bayz_output$Parameters$Param[coef_model_terms]
  estim(bayz_output, param = coef_parameters, splitlabels = FALSE,
        unlist = TRUE)
}
