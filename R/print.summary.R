#' Print function for 'summarybayz' class
#'
#' Prints the bayz class
#'
#' @param x    bayz object
#' @param ...  additional parameters
#'
#' @return Nothing, only prints
#' @import stats
#' @export
print.summarybayz <- function(object, ...) {
  if (object$Runinfo["Nerror"] > 0) {
    cat("Bayz encountered errors while running:\n")
    for (errormsg in object$Errors){
      cat("  ", errormsg, "\n")
    }
  } else {
    # No errors
    cat("Bayz model run completed successfully.\n")
    cat("\n")
    cat("Model run information:\n")
    print(object$Runinfo)
    cat("\n")
    cat("Model parameters and sizes:\n")
    print(noquote(object$Parameters))
    cat("  *Traced and summarized below.\n")
    cat("   Tracing can be toggled - see HelpIndex#tracing-parameters\n")
    cat("\n")
    cat("Estimates, HPD ", object$HPDprob * 100,
        "% intervals and convergence diagnostics for traced parameters:\n")
    print(object$summarystats)
    if (object$WarnFewSamples) {
      cat("  Warning: Output has few samples",
          ", convergence diagnostics may be unreliable.\n")
    }
    cat("\n")
  }
}
