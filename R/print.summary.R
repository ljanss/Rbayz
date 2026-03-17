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
print.summarybayz <- function(object, ...){
  cat("Summary of bayz model fit\n\n")
  cat("model formula:", deparse(object$terms), "\n\n")
  if (object$nError > 0) {
    cat("Bayz encountered errors while running:\n")
    for (errormsg in object$Errors){
      cat("  ", errormsg, "\n")
    }
  } else {
    # No errors
    cat("Parameters and sizes in the output; * is traced and summarized below:\n")
    print(object$Pameters)
    cat("\n")
    cat("Estimates and HPD intervals for 'traced' parameters:\n")
    print(object$summarystats)
    cat("\n")
  }
}
