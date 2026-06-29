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
#' set a 'save' option.
#' If there is only one rr() term in the model, prval() will automatically use
#' that one to compute pr-values, otherwise the user will need to specify which
#' parameter to use with the 'param' (or second) argument.
#'
#' @param object        A bayz model output
#' @param param         Parameter name for which to extract r-values. Can be
#'                      omitted if there is only one rr() term in the model.
#'                      (for parameter naming, check $Parameters from the
#'                      model output).
#' @param splitlabels   If the estimates come from an interaction, whether
#'                      labels should be split to allow easy matching to levels
#'                      in the single variables (default FALSE).
#' @param ...           Additional parameters.
#' @return a data frame with parameter estimates (posterior mean and SD), the
#'                      backtransformed z-statistic and its pr (p) -value.
#' @export
prval <- function(object, param = NULL, splitlabels = FALSE, ...) {
  if (is.null(param)) {
    # try to find a parameter to work on
    candidate_params <- (object$Parameters$ModelTerm == "rr" &
                           object$Parameters$Variance == "-" &
                           object$Parameters$Saved == "1")
    if (sum(candidate_params) == 0) {
      rr_term_params <- (object$Parameters$ModelTerm == "rr" &
                           object$Parameters$Variance == "-")
      if (sum(rr_term_params) > 0) {
        cat("There is/are candidate rr() term(s) (",
            paste(object$Parameters$Param[rr_term_params],
                  collapse = ", "),
            ")\n, but none with the 'save' option to allow prval() to work.\n",
            "You will need to re-fit the model with the 'save' option ",
            "on the rr() term(s) of interest.\n")
        return(invisible(NULL))
      } else {
        cat("Cannot find any rr() terms in the model output to ",
            "compute p_r values from.\n")
        return(invisible(NULL))
      }
    }
    if (sum(candidate_params) > 1) {
      cat("There are multiple candidate parameters in the output (",
          paste(object$Parameters$Param[candidate_params],
                collapse = ", "),
          ").\n",
          "Please specify which one to use with the 'param' argument.\n")
      return(invisible(NULL))
    }
    # If we get here, there is exactly one candidate parameter
    param <- object$Parameters$Param[candidate_params]
    param_row <- which(object$Parameters$Param == param)
    cat("Using rr() term with '", param, "' to compute p_r-values.\n", sep = "")
  } else {
    # User has given a 'param' to use
    if (!param %in% object$Parameters$Param) {
      cat("Specified 'param' not found in model output parameters.\n",
          "Check $Parameters$Param for valid names.\n")
      return(invisible(NULL))
    }
    # check that the param is from a rr() term and rightly indicated the random
    # effect, not its variance parameter (variance should be indicated as '-').
    param_row <- which(object$Parameters$Param == param)
    if (!(grepl("rr|rn", object$Parameters$ModelTerm[param_row]) &&
            object$Parameters$Variance[param_row] == "-")) {
      cat("Specified 'param' is not the random effects from an rr() term.\n",
          "Check $Parameters for valid parameters to use.\n")
      return(invisible(NULL))
    }
    if (object$Parameters$Saved[param_row] != "1") {
      cat("Specified 'param' did not have the 'save' option set, ",
          "which is required for prval() to work.\n",
          "You will need to re-fit the model with the 'save' option ",
          "on this rr() term.\n")
      return(invisible(NULL))
    }
  }
  param_var <- paste("var.", param, sep = "")
  if (!param_var %in% object$Parameters$Param) {
    cat("Needed variance parameter ", param_var,
        " is not in the model output parameters.\n",
        "(Maybe rr() was run with a non-Gaussian distribution?)\n",
        "See the requirements for computing p_r-values.\n")
    return(invisible(NULL))
  }
  param_var_row <- which(object$Parameters$Param == param_var)
  if (object$Parameters$Variance[param_var_row] != "IDEN" ||
        object$Parameters$Size[param_var_row] != 1) {
    cat("The ", param_var, " structure (",
         object$Parameters$Variance[param_var_row], ", size=",
         object$Parameters$Size[param_var_row], ") is too complex.\n",
         "See the requirements for computing p_r-values.\n")
    return(invisible(NULL))
  }
  if (object$Parameters$Trace[param_var_row] != "1") {
    cat("The variance ", param_var, " is not 'traced' ",
        "but this should be default.\n",
        "You may need to consult developers how this is possible.\n")
    return(invisible(NULL))
  }
  # Get the variance samples from the Samples table
  var_samples <- object$Samples[, param_var]
  if (is.null(var_samples)) {
    cat("Variance samples for ", param_var, " not found in model output.\n")
    return(invisible(NULL))
  }
  # Get the samples for the random effects from the saved samples file
  wd <- object$workdir
  samples_file <- file.path(wd, paste0("samples.", param, ".txt"))
  if (!file.exists(samples_file)) {
    cat("Saved samples file ", samples_file, " not found (moved?).\n")
    return(invisible(NULL))
  }
  n_samples <- nrow(object$Samples) # or in RunInfo["SamplesSaved"]
  n_effects <- object$Parameters$Size[param_row]
  effects_samples <- scan(samples_file, what = double(), quiet = TRUE)
  effects_labels <- object$Estimates[[param]]$Label
  if (length(effects_samples) != n_samples * (n_effects + 1)) {
    cat("Number of samples in ", samples_file, " does not match expected ",
        "number based on model output parameters\n.",
        "Files from different model fits may be mixed up?\n")
    return(invisible(NULL))
  }
  effects_samples <- matrix(effects_samples, nrow = n_samples,
                            ncol = (n_effects + 1), byrow = TRUE)
  effects_samples <- effects_samples[, -1, drop = FALSE] # drop cycle numbers
  # Now for every 'beta' regress out variance-sampling effects (2 times) to
  # obtain conditional posterior variance
  lhs <- matrix(c(n_samples, sum(var_samples), sum(var_samples),
                  sum(var_samples^2)), 2, 2)
  lhs_inv <- solve(lhs)
  rhs <- matrix(0, 2, 1)
  var_pm <- mean(var_samples)
  post_mean <- colMeans(effects_samples)
  post_sd <- apply(effects_samples, 2, sd)
  post_var_cond <- numeric(n_effects)
  for (i in 1:n_effects) {
    rhs[1, 1] <- sum(effects_samples[, i])
    rhs[2, 1] <- sum(effects_samples[, i] * var_samples)
    beta <- lhs_inv %*% rhs
    residuals <- effects_samples[, i] - beta[1] - beta[2] * var_samples
    rhs[1, 1] <- sum(residuals^2)
    rhs[2, 1] <- sum(residuals^2 * var_samples)
    beta <- lhs_inv %*% rhs
    post_var_cond[i] <- beta[1] + beta[2] * var_pm
  }
  z <- post_mean / sqrt(var_pm - post_var_cond)
  pval <- dnorm(abs(z))
  return_data <- data.frame(Label = effects_labels,
    PostMean = post_mean,
    PostSD = post_sd,
    PostSDCond = sqrt(post_var_cond),
    Zstat = z,
    Pr.value = pval
  )

#  if (splitlabels) {
#    # Try to split labels into variable and level, if possible
#    split_labels <- strsplit(return_data$Effect, split = ":")
#    if (all(sapply(split_labels, length) == 2)) {
#      return_data$Variable <- sapply(split_labels, `[`, 1)
#      return_data$Level <- sapply(split_labels, `[`, 2)
#    }
#  }
  return(return_data)
}
