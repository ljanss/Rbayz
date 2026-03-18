#' Summary of model run and main model parameters from bayz model fit.
#'
#' Summary() presents an overview of the model fitted, any errors and warnings,
#' and posterior summary statistics, including HPD intervals, for the main model
#' parameters in a bayz output. Summary() collects the 'traced' parameters,
#' for which samples are saved in the output, allowing to compute HPD intervals
#' (as well as trace plots, density plots and convergence diagnostics as done by
#' plot() and conv()). The 'traced' parameters are by default intercept,
#' fixed effects and regressions with up to 4 levels, and variances
#' and hyper-paramters. Tracing can be user-toggled
#' with a 'trace' option on model terms. Computation of HPD intervals accounts
#' for parameter boundaries, allowing variance intervals starting at zero.
#' The conv() function computes convergence diagnostics on the same set of
#' traced parameters.
#'
#' @param object    bayz output object
#' @param HPDprob   probability for the Highest Posterior Density intervals
#'                  (default 0.95)
#' @param ...       additional parameters
#'
#' @return summarybayz object
#' @import stats coda
#' @export
summary.bayz <- function(object, HPDprob=0.95, ...) {

  output <- list()
  class(output) <- "summarybayz"

  if (object$Runinfo["Nerror"] > 0) {
    output[["Errors"]] <- object[["Errors"]]
    output[["Runinfo"]] <- object[["Runinfo"]]
    return(output)
  }

  output[["Errors"]] <- NULL
  output[["Runinfo"]] <- object[["Runinfo"]]

  # It could be nice with a very compact overview of parameter names and sizes,
  # with a * if traced. Something like:
  par_sizes <- as.character(object$Parameters$Size)
  par_sizes <- paste0(par_sizes, ifelse(object$Parameters$Traced == 1, "*", ""))
  names(par_sizes) <- object$Parameters$Param
  par_sizes <- par_sizes[-1] # remove the 'fitval'
  output[["Parameters"]] <- par_sizes

  # Another useful addition could be various run-statistics:
  # Ndata, Nmissing, Nparameters, Chain-length, Burn-in, Chain-skip, Errors,
  # Notes, Warnings, time taken, ...

  # This summary now only lists the "traced" parameters that are in the Samples
  # table.
  output_cycles <- as.numeric(rownames(object$Samples))
  mcmc_samples <- coda::mcmc(object$Samples,
                             start = output_cycles[1],
                             end = output_cycles[length(output_cycles)],
                             thin = output_cycles[2] - output_cycles[1])
  postMeans <- apply(object$Samples, 2, mean)
  postSDs <- apply(object$Samples, 2, sd)
  HPDbounds <- rep("none", ncol(object$Samples))
  HPDbounds[substr(colnames(object$Samples), 0, 3) == "var"] <- "var"
  HPDs <- HPDbayz(object$Samples, prob = HPDprob, bound = HPDbounds)
  effSizes <- coda::effectiveSize(mcmc_samples)
  MCSEs <- postSDs / sqrt(effSizes)
  MCCVpct <- 100 * MCSEs / abs(postMeans)
  GewekeZ <- abs(coda::geweke.diag(mcmc_samples)$z)
  summary_table <- data.frame(postMeans, postSDs, HPDs,
                              effSizes, GewekeZ, MCSEs, MCCVpct)
  colnames(summary_table) <- c("postMean", "postSD", "HPDleft", "HPDright",
                               "effSize", "GewekeZ", "MCSE", "MCCV%")
  rownames(summary_table) <- colnames(object$Samples)
  output[["summarystats"]] <- summary_table
  output[["HPDprob"]] <- HPDprob
  output[["WarnFewSamples"]] <- ifelse(length(output_cycles) < 10, TRUE, FALSE)
  return(output)

}
