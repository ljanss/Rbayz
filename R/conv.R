#' Convergence diagnostics on traced parameters from bayz model fit
#'
#' Rbayz::conv() computes convergence diagnostics on the 'traced' (samples-saved) parameters from a bayz model fit:
#' effective sample size, Geweke Z-value, Monte Carlo Standard Error (MCSE) and Monte Carlo Coefficient of
#' Variation (MCCV%) using functions from the coda package. Note that the MCCV% can remain high when the parameter's
#' posterior mean is close to zero and the posterior mean is therefore also reported to help interpret the MCCV%.
#' Parameter 'burnin' can be used to increase the burnin setting from the orginal chain and discard additional samples
#' from the start of the chain.  
#'
#' @param object    bayz output object
#' @param burnin    Discard additional samples with cycle number less or equal to this (updated) burnin setting.
#'                  Setting burnin here lower than the original chain burnin has no effect.
#' @param ...       additional parameters
#'
#' @return bayzconv object; if bayz finished without errors, a list with a convergence table and convergence status,
#'         otherwise a list with error messages.
#' @import stats coda
#' @export
summary.bayz <- function(object, burnin=0, ...){

    return_object <- list()
    class(return_object) <- "bayzconv"
    return_object[['Errors']] <- object[['Errors']]
    return_object[['Chain']] <- object[['Chain']]
    if(return_object$nError>0){
        return(return_object)
    }


    # This summary now only lists the "traced" parameters with HPD intervals and convergence diagnostics, but the
    # tracing has been somewhat expanded to also include low-level (fixed) factors, and can be user-toggled.
    output_cycles = as.numeric(rownames(object$Samples))
    if(length(output_cycles)<10) convergence_status = 1       # fail because too few output
    convergence_table = data.frame()
    samp = coda::mcmc(object$Samples, start=output_cycles[1], end=output_cycles[length(output_cycles)],
                 thin=output_cycles[2]-output_cycles[1])
    effSizes = coda::effectiveSize(samp)
    postMeans =  apply(object$Samples,2,mean)
    postSDs  = apply(object$Samples,2,sd)
    MCSEs = postSDs / sqrt(effSizes)
    MCCVpct = 100 * MCSEs / abs(postMeans)
    GewekeZ = abs(coda::geweke.diag(object$Samples)$z)
    convergence_table = data.frame(postMeans,postSDs,effSizes,GewekeZ,MCSEs,MCCVpct,HPDs)
    colnames(convergence_table) = c("postMean","postSD","effSize","GewekeZ","MCSE","MCCV%","HPDleft","HPDright")
    rownames(convergence_table) = colnames(object$Samples)
    convergence_status = 0                                # success
    return_object[['Convergence']] = convergence_table
    return_object[['ConvergenceStatus']] = convergence_status

    return_object[['HPDprob']] = HPDprob

    return(return_object)

}
