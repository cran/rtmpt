
#' Set constants for probability parameters and suppress process times in a \code{rtmpt_model} list
#' 
#' By using \code{parameter = "probs"} you can specify which of the probability parameters should be set to a constant
#'   by using values between zero and one. If you use \code{NA} the probability will be estimated. By using 
#'   \code{parameter = "tau_minus"} or \code{parameter = "tau_plus"} you can suppress process times/rates. Here \code{0} will
#'   suppress the named process and \code{NA} allows the process time/rate to be estimated.
#'
#' @param model A list of the class \code{rtmpt_model}.
#' @param parameter Character of length one indicating the parameter to change. Allowed characters:
#'   \itemize{
#'    \item "probs": probability parameters
#'    \item "tau_minus": rate parameters of the exponential distribution of the process times that lead to a negative outcome
#'    \item "tau_plus": rate parameters of the exponential distribution of the process times that lead to a positive outcome
#'   }
#' @param names Character vector with process names.
#' @param values Numerical vector of length \code{length(names)}. By using \code{parameter = "probs"} you have the following options
#'   \itemize{
#'    \item \code{NA}: estimate the named probability
#'    \item \code{0 < values < 1}: set the named probability to a constant value between zero and one
#'   }
#'   Example: \code{set_params(model = model, parameter = "probs", names = c("do", "dn", "g"), values = c(NA, NA, .5))} will set
#'   the guessing "old" (g) to the constant \code{0.5} in the 2HT model.
#'   By using \code{parameter = "tau_minus"} or \code{parameter = "tau_plus"} you have two options:
#'   \itemize{
#'    \item \code{NA}: estimate the process time/rate
#'    \item \code{0}: suppress the process time/rate
#'   }
#'   Example: \code{set_params(model = model, parameter = "tau_minus", names = c("do", "dn", "g"), values = c(NA, NA, 0))} will suppress
#'   the process-completion time for guessing "new" in the 2HT model. This of course does not make sense here, but for some models it might be useful if you assume 
#'   that a time-consuming process is not associated with certain process-outcome pairs (e.g., for technical parameters not corresponding to a psychological process).
#' @return A list of the class \code{rtmpt_model}.
#' @examples
#' ####################################################################################
#' # Detect-Guess variant of the Two-High Threshold model.
#' # The encoding and motor execution times are assumed to be equal for each category.
#' # The process completion times for both failed detections will be suppressed.
#' ####################################################################################
#' 
#' mdl_2HTM <- "
#' # targets
#' do+(1-do)*g
#' (1-do)*(1-g)
#'
#' # lures
#' (1-dn)*g
#' dn+(1-dn)*(1-g)
#' 
#' # do: detect old; dn: detect new; g: guess
#' "
#' 
#' model <- to_rtmpt_model(mdl_file = mdl_2HTM)
#' 
#' ## removing the process times for the failed detection ("tau_minus") 
#' ## of the detection parameters ("dn", "do")
#' model <- set_params(model = model, parameter = "tau_minus", 
#'                     names = c("dn", "do"), values = c(0,0))
#'                     
#' @seealso \code{\link{set_resps}}
#' @author Raphael Hartmann
#' @export
set_params <- function(model, parameter, names, values = NA) {
  
  if (!is.list(model)) stop("model must be a list.")
  if (!("lines" %in% names(model)) || !("params" %in% names(model)) || !("responses" %in% names(model))) stop("No valid model file.")
  params_list <- model$params
  if (!is.list(params_list)) stop("params must be a list.")
  if ( !("probs" %in% names(params_list)) || 
       !("taus" %in% names(params_list)) ) stop("\"params\" must contain \"probs\" and \"taus\".")
  if (!is.character(parameter)) stop("\"parameters\" must be a character.")
  if (length(parameter) != 1) stop("\"parameters\" must be of length one.")
  if (!is.vector(names)) stop("\"names\" must be a vector")
  if (!is.vector(values)) stop("\"values\" must be a vector")
  if ( !(parameter %in% c("probs", "tau_minus", "tau_plus")) ) stop("Allowed \"parameters\" are \"probs\", \"tau_minus\" or \"tau_plus\".")
  if (!is.character(names)) stop("\"names\" must be characters.")
  if (!is.numeric(values) && !is.logical(values)) stop("\"values\" must be numerical.")
  if (length(values)!=1 & length(names) != length(values)) stop("Length of \"names\" and \"values\" must match.")
  if (length(values)==1 & length(names) != length(values)) values <- rep(values, length(names))
  if (any(!(names %in% names(params_list$probs)))) stop("\"names\" do not match with the names of the processes in \"params\".")
  
  if (parameter == "probs") {
    for (i in 1:length(values)) {
      if (!is.na(values[i])) if ( (values[i] <= 0) || (values[i] >= 1) ) stop("\"values\" for parameter \"probs\" must be larger than zero and lower than one.")
    }
    params_list[[parameter]][names] <- values
  } else {
    uniq_values <- unique(values)
    if (any(!(uniq_values %in% c(0,NA)))) stop("\"values\" for parameter \"tau_plus\" or \"tau_minus\" must either be zero or NA.")
    if (parameter == "tau_minus") {
      params_list[["taus"]]["minus", names] <- values
    }
    if (parameter == "tau_plus") {
      params_list[["taus"]]["plus", names] <- values
    }
  }
  
  model$params <- params_list
  return(model)
  
}
