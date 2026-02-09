#' Extract parameter estimates from a bayz model output
#'
#' Rbayz::estim() extracts one, several, or all parameter estimates from a
#' bayz output, reporting posterior means and standard deviations. Parameter
#' estimates also can be accessed directly in the $Estimates list in the
#' output, but the estim() function can split labels of interaction terms,
#' combine several parameters in one output, and can unlist the output, merging
#' estimates from several parameters in a single data frame.
#' Main uses are:
#' - estim(output) returns all parameter estimates as a list of data
#'   frames. The estimates for one set of parameters could then be accessed as
#'   estim(output)$`parameter-name`, for example estim(output)$"year.location".
#' - estim(output, param=name) with param set to a string returns a single data
#'   frame with
#'   estimates for the requested parameter. Example: estim(output, 
#'   param="year.location"), or shortly estim(output, "year.location").
#' - estim(output, param=c(name1, name2, ...)) with param a vector of strings
#'   returns a list of data frames with estimates for the requested parameters.
#' For the first and third usage, adding unlist=TRUE will reformat the default
#' list output into a single data frame.
#' For the second usage, adding splitlabels=TRUE will split labels of
#' interaction terms in multiple columns, which can be easier
#' to read and match with the original data. The second and third usage has
#' better error handling because estim() checks exsistence of the requested
#' parameter(s), and returns NULL if not available; with the first usage, and
#' when directly accessing $Estimates in the output, one risks attempting to
#' access non-existing parameters.
#' Parameter names are based on the variable names in the model-term, but can
#' have small modifications, notably colons in interactions
#' are replaced by dots. Random effect terms will create two slots in the
#' parameter list, one for the random effects themselves
#' and one for the variance (or other) hyper-parameters. If in doubt about
#' parameter names, check the $Parameters data frame in the output. 
#' Apart from estim(), Rbayz also provides fixef() and ranef() to retrieve
#' fixed and random effects, respectively, coef() to retrieve the combined
#' fixed and random effects, and
#' var() to retrieve all variance (and other) hyper-parameters.
#'
#' @param object        A bayz model output
#' @param param         A string with the name of the parameter for which to
#'                      retrieve estimates, or a vector of strings with the
#'                      names of multiple parameters. This argument can be
#'                      omitted to retrieve all parameter estimates.
#' @param splitlabels   Whether labels in interactions should be split in
#'                      multiple columns (default FALSE). Splitting labels is
#'                      only possible when retrieving estimates for a single
#'                      parameter.
#' @param unlist        Whether to return multiple sets of parameter estimates
#'                      in a single data frame
#'                      with an additional column Parameter to indicate from
#'                      which set of parameters each estimate comes. If FALSE
#'                      (default) multiple parameters are returned in a list of
#'                      data frames.
#' @param ...           Additional parameters.
#' @return Dependent on usage, a single data frame with parameter estimates,
#'         or a list of data frames with estimates
#'         for multiple parameters.
#' @export
#' 
estim <- function(object, param = NULL, splitlabels = FALSE,
                  unlist = FALSE, ...) {
  # First prepare estim_return as a list with the parameter(s) to be returned.
  # If no param given, this is simply the whole $Estimates list.
  if (is.null(param)) {
    estim_return <- object$Estimates
  } else {
    param <- as.vector(as.character(param))
    estim_return <- list()
    for (i in 1:seq_along(param)) {
      search_term <- param[i]
      # Repair common mistake: requesting A:B, which exists in the
      # $Parameters$Variables, but should be A.B as parameter-name.
      if (grepl(":", search_term, fixed = TRUE) &&
            search_term %in% object$Parameters$Variables) {
        search_term <- gsub(":", ".", search_term, fixed = TRUE)
      }
      if (!(search_term %in% names(object$Estimates))) {
        cat("Parameter '", search_term, "' not found in bayz output.",
            " Check output$Parameters for available parameter names.\n")
        return(NULL)
      }
      estim_return[[search_term]] <- object$Estimates[[search_term]]
    }
  }
  # Now estim_return is a list with one or more sets of parameters; if only one,
  # it is reformatted as a data frame, and it can optionally split labels.
  # If multiple, it can optionally 'unlist' to a single data frame.
  if (length(estim_return) == 1) {
    estim_return <- estim_return[[1]]
    if (splitlabels) {
      clmnames <- unlist(strsplit(search_term, ".", fixed = TRUE))
      split_labels_list <- strsplit(estim_return$Level, ".", fixed = TRUE)
      split_labels_lengths <- sapply(split_labels_list, length)
      if (any(split_labels_lengths != length(clmnames))) {
        cat("Error splitting labels: maybe there are extra dots",
            " in the original levels? The estimates will be",
            " returned without splitting labels.\n")
        return(estim_return)
      } else {
        split_labels_dframe <- t(as.data.frame(split_labels_list))
        colnames(split_labels_dframe) <- clmnames
        rownames(split_labels_dframe) <- NULL
        estim_return <- data.frame(split_labels_dframe,
                                   PostMean = estim_return$PostMean,
                                   PostSD = estim_return$PostSD)
        return(estim_return)
      }
  } else { # length(estim_return) > 1
    if (unlist) {
      estim_return_dframe = data.frame()
      for (i in 1:length(estim_return)) {
        estim_return_dframe <- rbind(estim_return_dframe, 
          cbind(Parameter=names(estim_return)[i], estim_return[[i]]))
      }
      return(estim_return_dframe)
    } else {
      return(estim_return)
    }
  }
}
