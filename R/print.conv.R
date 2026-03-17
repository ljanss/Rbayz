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
  if (object$nError > 0){
    cat("Bayz encountered errors while running:\n")
    for (errormsg in object$Errors){
      cat("  ",errormsg,"\n")
    }
  } else {
    # No errors
    cat("Estimates and HPD intervals for 'traced' parameters:\n")
        if(x$ConvergenceStatus == 1) {
            cat("*** This table is not printed because there are fewer than 10 output samples ***\n")
        }
        else if (x$ConvergenceStatus == 0) {
            print(object$Convergence[,c("postMean","postSD","HPDleft","HPDright")])
        }
        cat("\n")

        cat("Convergence diagnostics on 'traced' parameters:\n")
        if (x$ConvergenceStatus == 2) {
            cat("*** This table is not printed because the coda package is not installed ***\n")
        }
        else if (x$ConvergenceStatus == 0) {
            print(object$Convergence[,c("effSize","GewekeZ","MCSE","MCCV%")])
        }
        cat("\n")

    }
}
