#' Bayesian mixed models, shrinkage and interaction kernel regression
#'
#' The bayz() function fits various mixed-linear, Bayesian shrinkage, (sparse)
#' kernel-regression, (kernel) interaction and multi-trait models with complex
#' covariance structures using an extended R-formula syntax.
#' The basic R formula response ~ terms is extended by wrapping
#' every 'term' by a function that controls the type of variable(s) fitted, how
#' they are fitted, to add variance structures or shrinkage distribtions,
#' and options. As an example:
#' Yield ~ fx(Year:Site) + rn(Variety:Year:Site, V=Kv*Ky*Ks)
#' specifies a fixed Year:site interaction, and a random Variety:Year:Site
#' interaction with variance structure as a 3-D kronecker kernel interaction.
#' Kernels and kernel-interactions are fitted using reduced rank ('factor
#' analytic') approaches allowing to efficiently handle multi-dimensional
#' tensors, with options to control sparsity. Kernels can be replaced by
#' estimated (reduced rank) covariance structures,
#' allowing combinations of multi-trait, multi-environment, spatial,
#' temporal and other models. An rr() model term allows to fit random/ridge
#' regression on large sets of features, with options for homogeneous or
#' heterogeneous shrinkage and significances of individual features as
#' $p_r$ values (p-values from random effects). The Rbayz package supports many
#' common R model summary and extraction methods, such as summary(), coef(),
#' fixef(), ranef(), estim() (retrieval of any set(s) of model paramteres),
#' vcomp() (variance and hyper paramters and functions thereof), predict(),
#' conv() and plot(). The latter two use the R Coda package to diagnose MCMC
#' convergence and plot posterior distributions.
#' Full details are available in the github site ljanss.github.io/Rbayz.
#'
#' @param model   An R formula describing the model to be fitted using an
#'                extended R formula syntax (see details).
#' @param data    Data frame with variables (responses, explanatory variables)
#'                to build the model. Kernel matrices and matrices with large
#'                sets of features are typically not in this data frame.
#' @param Ve      Model for the residual variance as formula or string (see
#'                details).
#' @param chain   Vector c(length, burn-in, skip) with total chain length to
#'                run, burn-in, and skip-interval for saving samples and
#'                collecting posterior statistics.
#' @param method  String to indicate analysis method: "Bayes" (full Bayesian,
#'                default) or "BLUPMC" (BLUE/BLUP solutions with Monte Carlo
#'                to get SD/SE).
#' @param verbose Integer to regulate printing to R console: 0 (quiet), 1
#'                (some), >=2 (more). Default verbose=1.
#' @param workdir  Optional string with path to a directory where bayz
#'                can write output files (currently only used when the 'save'
#'                flag is added on a model term to save samples). When omitted
#'                the R working directory as obtained with getwd() will be
#'                used.
#' @param init    An object of class "bayz", which is output from a previous
#'                bayz run, to supply initialisation values to start a new
#'                chain.
#'
#' @return A list of class "bayz" containing results from the fitted model.
#'         Several methods are available to summarize, extract, plot or compute
#'         contrasts, predictions, etc. on the output.
#'
#' @import stats
#' @export
#'
#' @useDynLib Rbayz, .registration = TRUE
#' @importFrom Rcpp sourceCpp
bayz <- function(model, Ve = "", data = NULL, chain = c(0, 0, 0), method = "",
                 verbose = 1, workdir = NULL, init = NULL) {
  if (!inherits(model, "formula")) {
    stop("The first argument is not a valid formula")
  }
  if (!(class(Ve) == "character" || class(Ve) == "formula")) {
    stop("Ve must be given as a string or formula")
  }
  if (class(method) != "character") {
    stop("method must be given as a string")
  }
  if (is.null(data)) {
    stop("The data= argument is missing")
  }
  if (class(Ve) == "formula") {
    Ve <- deparse(Ve)
  }
  if (method == "") {
    method <- "Bayes"
  }
  if (!is.null(workdir)) {
    currdir <- getwd()
    on.exit(setwd(currdir))
    tryCatch(
      setwd(workdir),
      error = function(err) {
        stop("Failed to set working directory to '", workdir, "': ",
             err$message)
      }
    )
  }
  chain <- as.integer(chain)
  result <- rbayz_cpp(model, Ve, data, chain, method, verbose, init)
  result[["workdir"]] <- getwd()
  class(result) <- "bayz"
  return(result)
}
