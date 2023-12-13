

############ THETA TO CONSTANT ############

#' Set process probabilities to constants
#' 
#' Setting process probabilities (thetas) to constants or change it back to be estimated. 
#'
#' @param model An object of the class \code{ertmpt_model}.
#' @param names Character vector with process names.
#' @param constants Numerical vector of length one or \code{length(names)}. You have the following options for the elements of the numeric vector:
#'   \itemize{
#'    \item \code{0 < constants < 1}: set the named probability to a constant value between zero and one
#'    \item \code{NA}: estimate the named probability
#'   }
#' @return An object of the class \code{ertmpt_model}.
#' @examples
#' ####################################################################################
#' # Detect-Guess variant of the Two-High Threshold model.
#' # The encoding and motor execution times are assumed to be equal for each category.
#' # The process probability for guessing (g) will be set to 0.5.
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
#' model <- to_ertmpt_model(mdl_file = mdl_2HTM)
#' 
#' ## setting g to a constant (0.5):
#' new_model <- theta2const(model = model, names = c("g"), constants = c(0.5))
#' new_model
#'                     
#' @seealso \code{\link{delta2delta}}, \code{\link{tau2zero}}, \code{\link{theta2theta}} and \code{\link{tau2tau}}
#' @author Raphael Hartmann
#' @export
theta2const <- function(model, names, constants = NA) {
  
  if (!inherits(model, c("ertmpt_model", "rtmpt_model"))) stop("model must be of class \"ertmpt_model\".")
  if (!("lines" %in% names(model)) || !("params" %in% names(model)) || !("responses" %in% names(model))) stop("No valid model file.")
  
  params_list <- model$params
  if (!is.list(params_list)) stop("params must be a list.")
  if ( !("probs" %in% names(params_list)) || 
       !("taus" %in% names(params_list)) ) stop("\"params\" must contain \"probs\" and \"taus\".")
  if (!is.vector(names)) stop("\"names\" must be a vector")
  if (!is.vector(constants)) stop("\"constants\" must be a vector")
  if (!is.character(names)) stop("\"names\" must be characters.")
  if (!is.numeric(constants) && !is.logical(constants)) stop("\"constants\" must be numerical.")
  if (length(constants)!=1 & length(names) != length(constants)) stop("Length of \"names\" and \"constants\" must match.")
  if (length(constants)==1 & length(names) != 1) constants <- rep(constants, length(names))
  if (any(!(names %in% names(params_list$probs)))) stop("\"names\" do not match with the names of the processes in the \"model\".")
  
  for (i in 1:length(constants)) {
    if (!is.na(constants[i])) if ( (constants[i] <= 0) || (constants[i] >= 1) ) stop("\"constants\" must be larger than zero and lower than one.")
  }
  # for (i in 1:length(names)) {
  #   if (params_list[["probs"]][names[i]] %in% names(params_list[["probs"]])) {
  #     tmp_name <- as.character(params_list[["probs"]][names[i]])
  #     if(!tmp_name %in% names) {
  #       names <- c(tmp_name, names)
  #       constants <- c(constants[i], constants)
  #     }
  #   }
  # }
  for (i in 1:length(names)) {
    if (names[i] %in% params_list[["probs"]] & !is.na(constants[i])) {
      ind <- which(params_list[["probs"]] == names[i])
      tmp_names <- names(params_list[["probs"]])[ind]
      if(any(!tmp_names %in% names)) {
        ind2 <- which(!tmp_names %in% names)
        names <- c(names, tmp_names[ind2])
        constants <- c(constants, rep(constants[i], length(ind2)))
      }
    }
  }
  
  params_list[["probs"]][names] <- constants
  model$params <- params_list
  return(model)
}

#' @rdname theta2const
#' @examples
#' 
#' ## setting g to a constant (0.5):
#' new_model <- set_theta_const(model = model, names = c("g"), constants = c(0.5))
#' new_model
#' @export
set_theta_const <- theta2const



############ TAU TO ZERO ############

#' Set process completion times to zero
#' 
#' Setting process completion times (taus) to zero or change it back to be estimated.
#'
#' @param model A list of the class \code{ertmpt_model}.
#' @param names Character vector with process names.
#' @param outcomes Character vector of length \code{length(names)} indicating for which process outcome the process completion time should 
#'   be zero or changed back to be estimated. Allowed characters are:
#'   \itemize{
#'    \item \code{"minus"}: the negative outcome of the process.
#'    \item \code{"plus"}: the positive outcome of the process.
#'   }
#' @param values Numerical vector of length one or \code{length(names)}. You have the following options for the elements of the numeric vector:
#'   \itemize{
#'    \item \code{0}: suppress the process time/rate, i.e., set the process completion time (tau) with the specified output to zero.
#'    \item \code{NA}: estimate the process time (tau)
#'   }
#' @return A list of the class \code{ertmpt_model}.
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
#' model <- to_ertmpt_model(mdl_file = mdl_2HTM)
#' 
#' ## removing the process times (tau) for the failed (="minus") detection ("do" and "dn") 
#' new_model <- tau2zero(model = model, names = c("dn", "do"),
#'                       outcomes = c("minus", "minus"), values = 0)
#' new_model
#' 
#' @seealso \code{\link{delta2delta}}, \code{\link{theta2const}}, \code{\link{theta2theta}} and \code{\link{tau2tau}}
#' @author Raphael Hartmann
#' @export
tau2zero <- function(model, names, outcomes, values = 0) {
  
  if (!inherits(model, c("ertmpt_model", "rtmpt_model"))) stop("model must be of class \"ertmpt_model\".")
  if (!("lines" %in% names(model)) || !("params" %in% names(model)) || !("responses" %in% names(model))) stop("No valid model file.")

  params_list <- model$params
  if (!is.list(params_list)) stop("params must be a list.")
  if ( !("probs" %in% names(params_list)) || 
       !("taus" %in% names(params_list)) ) stop("\"params\" must contain \"probs\" and \"taus\".")
  if (!is.character(outcomes)) stop("\"outcomes\" must be a character.")
  if (!is.vector(names)) stop("\"names\" must be a vector")
  if (!is.vector(values)) stop("\"values\" must be a vector")
  if ( !(all(outcomes %in% c("minus", "plus"))) ) stop("Allowed \"outcomes\" are \"minus\" or \"plus\".")
  if (!is.character(names)) stop("\"names\" must be characters.")
  if (!is.numeric(values) && !is.logical(values)) stop("\"values\" must be numerical.")
  if (length(outcomes)!=1 & length(names) != length(outcomes)) stop("Length of \"names\" and \"outcomes\" must match.")
  if (length(values)!=1 & length(names) != length(values)) stop("Length of \"names\" and \"values\" must match.")
  if (length(values)==1 & length(names) != length(values)) values <- rep(values, length(names))
  if (any(!(names %in% names(params_list$probs)))) stop("\"names\" do not match with the names of the processes in the \"model\".")
  
  uniq_values <- unique(values)
  if (any(!(uniq_values %in% c(0,NA)))) stop("\"values\" must either be zero or NA.")
  
  ord <- order(names)
  names <- names[ord]
  outcomes <- outcomes[ord]
  values <- values[ord]
  
  combo <- data.frame(names = names, outcomes = outcomes, values = values)
  
  # is there a reference process in names, then make all processes that have a reference to that process equal
  for (i in 1:length(names)) {
    if (names[i] %in% params_list[["taus"]][outcomes[i],] & !is.na(values[i])) {
      ind <- which(params_list[["taus"]][outcomes[i],] == names[i])
      tmp_names <- names(params_list[["taus"]])[ind]
      for(k in 1:length(tmp_names)) {
        if(all(!sapply(1:length(names), function(j) all(c(tmp_names[k], outcomes[i])==combo[j, 1:2])))) {
          names <- c(names, tmp_names[k])
          outcomes <- c(outcomes, outcomes[i])
          values <- c(values, values[i])
          combo <- rbind(combo, c(tmp_names[k], outcomes[i], values[]))
        }
      }
    }
  }
  
  # If one of the "minus" or "plus" outcome values is NA then the other must be NA or 0
  if (anyNA(values)) {
    ind <- which(is.na(values))
    tmp_names <- names[ind]
    for(k in 1:length(tmp_names)) {
      tmp_outc <- outcomes[ind[k]]
      opp_outc <- ifelse(tmp_outc == "minus", "plus", "minus")
      if(params_list[["taus"]][opp_outc, tmp_names[k]] %in% names(params_list[["taus"]])) {
        if(all(!sapply(1:length(names), function(j) all(c(tmp_names[k], opp_outc)==combo[j, 1:2])))) {
          names <- c(names, tmp_names[k])
          outcomes <- c(outcomes, opp_outc)
          values <- c(values, NA)
          combo <- rbind(combo, c(tmp_names[k], opp_outc, NA))
        }
      }
      
    }
  }
  
  # if a process with reference should be set to zero, but the reference process not, then either both outcome values must be 0 or one must be NA
  if (!anyNA(values)) {
    ind <- which(!is.na(values))
    tmp_names <- names[ind]
    for (k in 1:length(tmp_names)) {
      tmp_outc <- outcomes[ind[k]]
      opp_outc <- ifelse(tmp_outc == "minus", "plus", "minus")
      if(params_list[["taus"]][opp_outc, tmp_names[k]] %in% names(params_list[["taus"]])) {
        ref_name <- params_list[["taus"]][opp_outc, tmp_names[k]]
        if(all(!sapply(1:length(names), function(j) all(c(ref_name, outcomes[ind[k]])==combo[j, 1:2])))) {
          names <- c(names, tmp_names[k])
          outcomes <- c(outcomes, opp_outc)
          values <- c(values, NA)
          combo <- rbind(combo, c(tmp_names[k], opp_outc, NA))
        }
      }
    }
  }
  
  for(i in 1:length(outcomes)) {
    params_list[["taus"]][outcomes[i], names[i]] <- values[i]
  }
  for(i in 1:length(params_list[["taus"]][1,])) {
    if(all(!params_list[["taus"]][, i] %in% names(params_list[["taus"]]))) {
      params_list[["taus"]][, i] <- as.numeric(params_list[["taus"]][, i] )
    }
  }
  
  model$params <- params_list
  return(model)
}


#' @rdname tau2zero
#' @examples
#' 
#' ## removing the process times (tau) for the failed (="minus") detection ("do" and "dn") 
#' new_model <- set_tau_zero(model = model, names = c("dn", "do"),
#'                           outcomes = c("minus", "minus"), values = 0)
#' new_model
#' @export
set_tau_zero <- tau2zero



############ MAKE THETAS EQUAL ############

#' Set process probabilities equal
#' 
#' Setting multiple process probabilities (thetas) equal. One of the process probabilities will be estimated and
#'   the other named process(es) will be set to equal the former. The equality can be removed by only using one name of a process. 
#'
#' @param model A list of the class \code{ertmpt_model}.
#' @param names Character vector giving the names of the processes for which the process probabilities should be equal. If 
#'   \code{length(names) = 1} then the corresponding process probability will be estimates (i.e., it will be set to NA)
#' @param keep_consts Can be one of the following
#'   \itemize{
#'    \item logical value: \code{FALSE} (default) means none of the constants for \code{names} in the \code{model} will be kept; The probability of
#'      the reference process (i.e., first of \code{names} in alphabetical order) will be set to \code{NA} (i.e., will be estimated) and the others
#'      will be set to the name of the reference process (i.e., will be set to equal the reference process probability). \code{TRUE} means
#'      the constant of the reference process probability (if specified) is used for all other processes. 
#'    \item numeric value: index for \code{names}. If 1, the constant of the first process in \code{names} (in original order defined by the user) is 
#'      used for all other probabilities of the processes in \code{names}. If 2, the constant of the second process is used. And so on.
#'   }
#' @return A list of the class \code{ertmpt_model}.
#' @note If you use \code{theta2theta()} and \code{tau2tau()} with the same process names you might just change the EQN or MDL file accordingly
#'   by using the same process name for all processes which should have equal process times and probabilities.
#' @examples
#' ####################################################################################
#' # Detect-Guess variant of the Two-High Threshold model.
#' # The encoding and motor execution times are assumed to be equal for each category.
#' # The process probabilities for both detection processes ("do" and "dn") will be
#' # set equal.
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
#' model <- to_ertmpt_model(mdl_file = mdl_2HTM)
#' 
#' ## make do = dn
#' new_model <- theta2theta(model = model, names = c("do", "dn"))
#' new_model
#' 
#' @seealso \code{\link{delta2delta}}, \code{\link{theta2const}}, \code{\link{tau2zero}} and \code{\link{tau2tau}}
#' @author Raphael Hartmann
#' @export
theta2theta <- function(model, names, keep_consts = FALSE) {
  
  if (!inherits(model, c("ertmpt_model", "rtmpt_model"))) stop("model must be of class \"ertmpt_model\".")
  if (!("lines" %in% names(model)) || !("params" %in% names(model)) || !("responses" %in% names(model))) stop("No valid model file.")
  
  params_list <- model$params
  if (!is.list(params_list)) stop("params must be a list.")
  if ( !("probs" %in% names(params_list)) || 
       !("taus" %in% names(params_list)) ) stop("\"params\" must contain \"probs\" and \"taus\".")
  if (!is.vector(names)) stop("\"names\" must be a vector")
  if (!is.character(names)) stop("\"names\" must be characters.")
  if (any(!(names %in% names(params_list$probs)))) stop("\"names\" do not match with the names of the processes in the \"model\".")
  if (!is.logical(keep_consts) & !is.numeric(keep_consts)) stop("\"keep_consts\" must either be a logical value or a numeric value smaller or equal the length of \"names\".")
  if (length(keep_consts)!=1) stop("\"keep_consts\" must be of length one")
  if (length(unique(names)) == 1) stop("\"names\" must be of length two or larger")
  names <- unique(names)
  
  sorted_names <- sort(names)
  if (keep_consts == FALSE) {
    params_list[["probs"]][names] <- sorted_names[1]
    params_list[["probs"]][sorted_names[1]] <- NA
  } else if (keep_consts == TRUE) {
    val <- params_list[["probs"]][sorted_names[1]]
    if (is.na(val)) {
      params_list[["probs"]][names] <- sorted_names[1]
      params_list[["probs"]][sorted_names[1]] <- NA
    } else if (val %in% names(params_list[["probs"]])) {
      params_list[["probs"]][names] <- as.character(val)
      params_list[["probs"]][as.character(val)] <- NA
    } else if (is.numeric(val)) {
      params_list[["probs"]][names] <- as.numeric(val)
    }
  } else {
    val <- as.numeric(params_list[["probs"]][names[keep_consts]])
    if (is.na(val)) {
      params_list[["probs"]][names] <- sorted_names[1]
      params_list[["probs"]][sorted_names[1]] <- NA
    } else if (val %in% names(params_list[["probs"]])) {
      params_list[["probs"]][names] <- sorted_names[1]
      params_list[["probs"]][sorted_names[1]] <- NA
    } else if (is.numeric(val)) {
      params_list[["probs"]][names] <- as.numeric(val)
    }
  }
  
  model$params <- params_list
  return(model)
}


#' @rdname theta2theta
#' @examples
#' 
#' ## make do = dn
#' new_model <- set_thetas_equal(model = model, names = c("do", "dn"))
#' new_model
#' @export
set_thetas_equal <- theta2theta



############ MAKE TAUS EQUAL ############

#' Set process completion times equal
#' 
#' Setting multiple process completion times (taus) equal. This means all process times of negative outcomes will be
#'   set equal and all process times of positive outcomes will be set equal. Only two process times (one for the negative
#'   and one for the positive outcome) of the named processes will be estimated. The equality can be removed by just 
#'   naming only one process name.
#'
#' @param model A list of the class \code{ertmpt_model}.
#' @param names Character vector giving the names of the processes for which the process completion times should be equal. If 
#'   \code{length(names) = 1} then the corresponding process completion times (for negative and positive outcomes) will be 
#'   estimates (i.e., they will be set to NA)
#' @param outcome Character (no vector) indicating for which process outcome the process completion times should 
#'   be set equal. Allowed characters are:
#'   \itemize{
#'    \item \code{"minus"}: the negative outcome of the processes.
#'    \item \code{"plus"}: the positive outcome of the processes.
#'    \item \code{"both"}: the negative and positive outcome of the processes. This will set all process completion times for
#'      the "minus" outcome equal and all process completion times for the "plus" outcome equal.
#'   }
#' @param keep_zeros Can be one of the following
#'   \itemize{
#'    \item logical value: \code{FALSE} (default) means none of the zeros for \code{names} in the \code{model} will be kept; The times of
#'      the reference process (i.e., first of \code{names} in alphabetical order) will be set to \code{NA} (i.e., will be estimated) and the others
#'      will be set to the name of the reference process (i.e., will be set to equal the reference process times). \code{TRUE} means
#'      the zero(s) of the reference process times (if specified) is used for the same outcome of all other processes. 
#'    \item numeric value: index for \code{names}. If 1, the zero(s) of the first process in \code{names} (in original order defined by the user) is 
#'      used for the same outcome of all other processes in \code{names}. If 2, the zero(s) of the second process is used. And so on.
#'   }
#' @return A list of the class \code{ertmpt_model}.
#' @note If you use \code{theta2theta()} and \code{tau2tau()} with the same process names you might just change the EQN or MDL file accordingly
#'   by using the same process name for all processes which should have equal process times and probabilities.
#' @examples
#' ####################################################################################
#' # Detect-Guess variant of the Two-High Threshold model.
#' # The encoding and motor execution times are assumed to be equal for each category.
#' # The process completion times for both detection processes ("do" and "dn") will be
#' # set equal.
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
#' model <- to_ertmpt_model(mdl_file = mdl_2HTM)
#' 
#' ## make do = dn
#' new_model <- tau2tau(model = model, names = c("do", "dn"), outcome = "both")
#' new_model
#' 
#' @seealso \code{\link{delta2delta}}, \code{\link{theta2const}}, \code{\link{tau2zero}} and \code{\link{theta2theta}}
#' @author Raphael Hartmann
#' @export
tau2tau <- function(model, names, outcome, keep_zeros = FALSE) {
  
  if (!inherits(model, c("ertmpt_model", "rtmpt_model"))) stop("model must be of class \"ertmpt_model\".")
  if (!("lines" %in% names(model)) || !("params" %in% names(model)) || !("responses" %in% names(model))) stop("No valid model file.")
  
  params_list <- model$params
  if (!is.list(params_list)) stop("params must be a list.")
  if ( !("probs" %in% names(params_list)) || 
       !("taus" %in% names(params_list)) ) stop("\"params\" must contain \"probs\" and \"taus\".")
  if (!is.vector(names)) stop("\"names\" must be a vector")
  if (!is.character(names)) stop("\"names\" must be characters.")
  if (!is.character(outcome)) stop("\"outcome\" must be a character.")
  if (any(!(names %in% names(params_list$probs)))) stop("\"names\" do not match with the names of the processes in the \"model\".")
  if (!is.logical(keep_zeros) & !is.numeric(keep_zeros)) stop("\"keep_zeros\" must either be a logical value or a numeric value smaller or equal the length of \"names\".")
  if (length(outcome)!=1) stop("Length of \"outcome\" must be 1")
  if (!(outcome %in% c("minus", "plus", "both"))) stop("Allowed \"outcome\" is \"minus\", \"plus\", or \"both\".")
  if (length(keep_zeros)!=1) stop("\"keep_zeros\" must be of length one")
  if (length(unique(names)) == 1) stop("\"names\" must be of length two or larger")
  if (keep_zeros > length(names)) stop("\"keep_zeros\" invalid")
  names <- unique(names)
  
  if (outcome == "both") outcome <- c("minus", "plus")
  
  sorted_names <- sort(names)
  if (keep_zeros == FALSE) {
    
    # default way... ignoring suppressed process times
    params_list[["taus"]][outcome, names] <- sorted_names[1]
    params_list[["taus"]][outcome, sorted_names[1]] <- NA
    
  } else if (keep_zeros == TRUE & is.logical(keep_zeros)) {
    
    val <- params_list[["taus"]][outcome, sorted_names[1]]

    if (any(is.na(val))) {
      if (all(is.na(val))) {
        params_list[["taus"]][outcome, names] <- sorted_names[1]
        params_list[["taus"]][outcome, sorted_names[1]] <- NA
      } else if (any(val==0)) {
        params_list[["taus"]][outcome, names] <- sorted_names[1]
        params_list[["taus"]][outcome, sorted_names[1]] <- NA
        ind_zero <- which(val == 0)
        params_list[["taus"]][ind_zero, names] <- val[ind_zero]
      }
    } else if (any(val %in% names(params_list[["taus"]]))) {
      if (val[1]==val[2]) {
        params_list[["taus"]][outcome, names] <- sorted_names[1]
        params_list[["taus"]][outcome, sorted_names[1]] <- NA
      } else if (any(val == 0)) {
        ind_zero <- which(val == 0)
        ind_proc <- which(as.character(val) %in% names(params_list[["taus"]]))
        params_list[["taus"]][outcome, names] <- sorted_names[1]
        params_list[["taus"]][ind_proc, sorted_names[1]] <- NA
        params_list[["taus"]][ind_zero, sorted_names[1]] <- 0
        params_list[["taus"]][ind_zero, names] <- 0
      } else if (all(val %in% names(params_list[["taus"]])) & val[1]!=val[2]) {
        stop("Something is not right with the model. There cannot be two different process names in one column.")
      }
    } else if (all(val == 0)) {
      params_list[["taus"]][outcome, names] <- 0
    }
    
  } else {
    
    val <- params_list[["taus"]][outcome, names[keep_zeros]]
    
    if (any(is.na(val))) {
      if (all(is.na(val))) {
        params_list[["taus"]][outcome, names] <- sorted_names[1]
        params_list[["taus"]][outcome, sorted_names[1]] <- NA
      } else if (any(val==0)) {
        params_list[["taus"]][outcome, names] <- sorted_names[1]
        params_list[["taus"]][outcome, sorted_names[1]] <- NA
        ind_zero <- which(val == 0)
        params_list[["taus"]][ind_zero, names] <- val[ind_zero]
      }
    } else if (any(val %in% names(params_list[["taus"]]))) {
      if (val[1]==val[2]) {
        params_list[["taus"]][outcome, names] <- sorted_names[1]
        params_list[["taus"]][outcome, sorted_names[1]] <- NA
      } else if (any(val == 0)) {
        ind_zero <- which(val == 0)
        ind_proc <- which(as.character(val) %in% names(params_list[["taus"]]))
        params_list[["taus"]][outcome, names] <- sorted_names[1]
        params_list[["taus"]][ind_proc, sorted_names[1]] <- NA
        params_list[["taus"]][ind_zero, sorted_names[1]] <- 0
        params_list[["taus"]][ind_zero, names] <- 0
      } else if (all(val %in% names(params_list[["taus"]])) & val[1]!=val[2]) {
        stop("Something is not right with the model. There cannot be two different process names in one column.")
      }
    } else if (all(val == 0)) {
      params_list[["taus"]][outcome, names] <- 0
    }
    
  }
  
  for(i in 1:length(params_list[["taus"]][outcome[1],])) {
    if(all(!params_list[["taus"]][outcome, i] %in% names(params_list[["taus"]]))) {
      params_list[["taus"]][outcome, i] <- as.numeric(params_list[["taus"]][outcome, i] )
    }
  }
  
  model$params <- params_list
  return(model)
}


#' @rdname tau2tau
#' @examples
#' 
#' ## make do = dn
#' new_model <- set_taus_equal(model = model, names = c("do", "dn"), outcome = "both")
#' new_model
#' @export
set_taus_equal <- tau2tau














############ A TO CONSTANT ############

#' Set process threshold to constants
#'
#' Setting process thresholds (parameter a) to constants or change it back to be estimated.
#'
#' @param model An object of the class \code{rtmpt_model}.
#' @param names Character vector with process names.
#' @param constants Numerical vector of length one or \code{length(names)}. You have the following options for the elements of the numeric vector:
#'   \itemize{
#'    \item \code{0 < constants}: set the named threshold parameter(s) to constant value(s) larger than zero
#'    \item \code{NA}: estimate the named threshold parameter(s)
#'   }
#' @return An object of the class \code{drtmpt_model}.
#' @examples
#' ####################################################################################
#' # Detect-Guess variant of the Two-High Threshold model.
#' # The encoding and motor execution times are assumed to be equal for each category.
#' # The process threshold for guessing (g) will be set to 1.0.
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
#' model <- to_drtmpt_model(mdl_file = mdl_2HTM)
#'
#' ## setting threshold for g to a constant (1.0):
#' new_model <- a2const(model = model, names = c("g"), constants = c(1.0))
#' new_model
#'
#' @seealso \code{\link{delta2delta}}, \code{\link{a2a}}, \code{\link{nu2const}}, \code{\link{nu2nu}}, \code{\link{omega2const}} and \code{\link{omega2omega}}
#' @author Raphael Hartmann
#' @export
a2const <- function(model, names, constants = NA) {
  
  if (!inherits(model, "drtmpt_model")) stop("model must be of class \"drtmpt_model\".")
  if (!("lines" %in% names(model)) || !("params" %in% names(model)) || !("responses" %in% names(model))) stop("No valid model file.")
  
  params_list <- model$params
  if (!is.list(params_list)) stop("params must be a list.")
  if ( !("threshold" %in% names(params_list)) ||
       !("driftrate" %in% names(params_list)) ||
       !("startpoint" %in% names(params_list)) ) stop("\"params\" must contain \"threshold\", \"driftrate\" and \"startpoint\".")
  if (!is.vector(names)) stop("\"names\" must be a vector")
  if (!is.vector(constants)) stop("\"constants\" must be a vector")
  if (!is.character(names)) stop("\"names\" must be characters.")
  if (!is.numeric(constants) && !is.logical(constants)) stop("\"constants\" must be numerical.")
  if (length(constants)!=1 & length(names) != length(constants)) stop("Length of \"names\" and \"constants\" must match.")
  if (length(constants)==1 & length(names) != 1) constants <- rep(constants, length(names))
  if (any(!(names %in% names(params_list$threshold)))) stop("\"names\" do not match with the names of the processes in the \"model\".")
  
  for (i in 1:length(constants)) {
    if (!is.na(constants[i])) if ( (constants[i] <= 0) ) stop("\"constants\" must be larger than zero.")
  }
  for (i in 1:length(names)) {
    if (names[i] %in% params_list[["threshold"]] & !is.na(constants[i])) {
      ind <- which(params_list[["threshold"]] == names[i])
      tmp_names <- names(params_list[["threshold"]])[ind]
      if(any(!tmp_names %in% names)) {
        ind2 <- which(!tmp_names %in% names)
        names <- c(names, tmp_names[ind2])
        constants <- c(constants, rep(constants[i], length(ind2)))
      }
    }
  }
  
  params_list[["threshold"]][names] <- constants
  model$params <- params_list
  return(model)
}

#' @rdname a2const
#' @examples
#'
#' ## setting threshold of g to a constant (1.0):
#' new_model <- set_a_const(model = model, names = c("g"), constants = c(1.0))
#' new_model
#' @export
set_a_const <- a2const











############ NU TO CONSTANT ############

#' Set process drift rate to constants
#'
#' Setting process drif rate (parameter nu) to constants or change it back to be estimated.
#'
#' @param model An object of the class \code{rtmpt_model}.
#' @param names Character vector with process names.
#' @param constants Numerical vector of length one or \code{length(names)}. You have the following options for the elements of the numeric vector:
#'   \itemize{
#'    \item \code{-Inf < constants < Inf}: set the named drift rate parameter(s) to constant value(s)
#'    \item \code{NA}: estimate the named drift rate parameter(s)
#'   }
#' @return An object of the class \code{drtmpt_model}.
#' @examples
#' ####################################################################################
#' # Detect-Guess variant of the Two-High Threshold model.
#' # The encoding and motor execution times are assumed to be equal for each category.
#' # The process drift rate for guessing (g) will be set to 1.0.
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
#' model <- to_drtmpt_model(mdl_file = mdl_2HTM)
#'
#' ## setting drift rate for g to a constant (1.0):
#' new_model <- nu2const(model = model, names = c("g"), constants = c(1.0))
#' new_model
#'
#' @seealso \code{\link{delta2delta}}, \code{\link{a2const}}, \code{\link{a2a}}, \code{\link{nu2nu}}, \code{\link{omega2const}} and \code{\link{omega2omega}}
#' @author Raphael Hartmann
#' @export
nu2const <- function(model, names, constants = NA) {
  
  if (!inherits(model, "drtmpt_model")) stop("model must be of class \"drtmpt_model\".")
  if (!("lines" %in% names(model)) || !("params" %in% names(model)) || !("responses" %in% names(model))) stop("No valid model file.")
  
  params_list <- model$params
  if (!is.list(params_list)) stop("params must be a list.")
  if ( !("threshold" %in% names(params_list)) ||
       !("driftrate" %in% names(params_list)) ||
       !("startpoint" %in% names(params_list)) ) stop("\"params\" must contain \"threshold\", \"driftrate\" and \"startpoint\".")
  if (!is.vector(names)) stop("\"names\" must be a vector")
  if (!is.vector(constants)) stop("\"constants\" must be a vector")
  if (!is.character(names)) stop("\"names\" must be characters.")
  if (!is.numeric(constants) && !is.logical(constants)) stop("\"constants\" must be numerical.")
  if (length(constants)!=1 & length(names) != length(constants)) stop("Length of \"names\" and \"constants\" must match.")
  if (length(constants)==1 & length(names) != 1) constants <- rep(constants, length(names))
  if (any(!(names %in% names(params_list$driftrate)))) stop("\"names\" do not match with the names of the processes in the \"model\".")
  
  # for (i in 1:length(constants)) {
  #   if (!is.na(constants[i])) if ( (constants[i] <= 0) ) stop("\"constants\" must be larger than zero.")
  # }
  for (i in 1:length(names)) {
    if (names[i] %in% params_list[["driftrate"]] & !is.na(constants[i])) {
      ind <- which(params_list[["driftrate"]] == names[i])
      tmp_names <- names(params_list[["driftrate"]])[ind]
      if(any(!tmp_names %in% names)) {
        ind2 <- which(!tmp_names %in% names)
        names <- c(names, tmp_names[ind2])
        constants <- c(constants, rep(constants[i], length(ind2)))
      }
    }
  }
  
  params_list[["driftrate"]][names] <- constants
  model$params <- params_list
  return(model)
}

#' @rdname nu2const
#' @examples
#'
#' ## setting drift rate of g to a constant (1.0):
#' new_model <- set_nu_const(model = model, names = c("g"), constants = c(1.0))
#' new_model
#' @export
set_nu_const <- nu2const






############ OMEGA TO CONSTANT ############

#' Set process relative starting-point to constants
#'
#' Setting process relative starting-point (parameter omega) to constants or change it back to be estimated.
#'
#' @param model An object of the class \code{rtmpt_model}.
#' @param names Character vector with process names.
#' @param constants Numerical vector of length one or \code{length(names)}. You have the following options for the elements of the numeric vector:
#'   \itemize{
#'    \item \code{0 < constants < 0}: set the named relaitve starting-point parameter(s) to constant value(s) larger than zero and smaller than 1
#'    \item \code{NA}: estimate the named relaitve starting-point parameter(s)
#'   }
#' @return An object of the class \code{drtmpt_model}.
#' @examples
#' ####################################################################################
#' # Detect-Guess variant of the Two-High Threshold model.
#' # The encoding and motor execution times are assumed to be equal for each category.
#' # The process relative starting-point for guessing (g) will be set to 0.5.
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
#' model <- to_drtmpt_model(mdl_file = mdl_2HTM)
#'
#' ## setting relative starting-point for g to a constant (1.0):
#' new_model <- omega2const(model = model, names = c("g"), constants = c(0.5))
#' new_model
#'
#' @seealso \code{\link{delta2delta}}, \code{\link{a2const}}, \code{\link{a2a}}, \code{\link{nu2const}}, \code{\link{nu2nu}}, and \code{\link{omega2omega}}
#' @author Raphael Hartmann
#' @export
omega2const <- function(model, names, constants = NA) {
  
  if (!inherits(model, "drtmpt_model")) stop("model must be of class \"drtmpt_model\".")
  if (!("lines" %in% names(model)) || !("params" %in% names(model)) || !("responses" %in% names(model))) stop("No valid model file.")
  
  params_list <- model$params
  if (!is.list(params_list)) stop("params must be a list.")
  if ( !("threshold" %in% names(params_list)) ||
       !("driftrate" %in% names(params_list)) ||
       !("startpoint" %in% names(params_list)) ) stop("\"params\" must contain \"threshold\", \"driftrate\" and \"startpoint\".")
  if (!is.vector(names)) stop("\"names\" must be a vector")
  if (!is.vector(constants)) stop("\"constants\" must be a vector")
  if (!is.character(names)) stop("\"names\" must be characters.")
  if (!is.numeric(constants) && !is.logical(constants)) stop("\"constants\" must be numerical.")
  if (length(constants)!=1 & length(names) != length(constants)) stop("Length of \"names\" and \"constants\" must match.")
  if (length(constants)==1 & length(names) != 1) constants <- rep(constants, length(names))
  if (any(!(names %in% names(params_list$startpoint)))) stop("\"names\" do not match with the names of the processes in the \"model\".")
  
  for (i in 1:length(constants)) {
    if (!is.na(constants[i])) if ( (constants[i] <= 0) || (constants[i] >= 1) ) stop("\"constants\" must be larger than zero and smaller than 1.")
  }
  for (i in 1:length(names)) {
    if (names[i] %in% params_list[["startpoint"]] & !is.na(constants[i])) {
      ind <- which(params_list[["startpoint"]] == names[i])
      tmp_names <- names(params_list[["startpoint"]])[ind]
      if(any(!tmp_names %in% names)) {
        ind2 <- which(!tmp_names %in% names)
        names <- c(names, tmp_names[ind2])
        constants <- c(constants, rep(constants[i], length(ind2)))
      }
    }
  }
  
  params_list[["startpoint"]][names] <- constants
  model$params <- params_list
  return(model)
}

#' @rdname omega2const
#' @examples
#'
#' ## setting relative starting-point of g to a constant (0.5):
#' new_model <- set_omega_const(model = model, names = c("g"), constants = c(0.5))
#' new_model
#' @export
set_omega_const <- omega2const






############ MAKE A EQUAL ############

#' Set process thresholds equal
#'
#' Setting multiple process thresholds (parameter a) equal. One of the process thresholds will be estimated and the other
#'   named thresholds will be set to equal the former. The equality can be removed by only using one name of a process.
#'
#' @param model A list of the class \code{drtmpt_model}.
#' @param names Character vector giving the names of the processes for which the process thresholds should be equal. If
#'   \code{length(names) = 1} then the corresponding process threshold will be estimated (i.e., it will be set to NA)
#' @param keep_consts Can be one of the following
#'   \itemize{
#'    \item logical value: \code{FALSE} (default) means none of the constants for \code{names} in the \code{model} will be kept; The thresholds of
#'      the reference process (i.e., first of \code{names} in alphabetical order) will be set to \code{NA} (i.e., will be estimated) and the others
#'      will be set to the name of the reference process (i.e., will be set to equal the reference process thresholds). \code{TRUE} means
#'      the constant of the reference process threshold (if specified) is used for all other processes.
#'    \item numeric value: index for \code{names}. If 1, the constant of the first process in \code{names} (in original order defined by the user) is
#'      used for all other thresholds of the processes in \code{names}. If 2, the constant of the second process is used. And so on.
#'   }
#' @return A list of the class \code{drtmpt_model}.
#' @examples
#' ####################################################################################
#' # Detect-Guess variant of the Two-High Threshold model.
#' # The encoding and motor execution times are assumed to be equal for each category.
#' # The process thresholds for both detection processes ("do" and "dn")
#' # will be set equal.
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
#' model <- to_drtmpt_model(mdl_file = mdl_2HTM)
#'
#' ## make do = dn
#' new_model <- a2a(model = model, names = c("do", "dn"))
#' new_model
#'
#' @seealso \code{\link{delta2delta}}, \code{\link{a2const}}, \code{\link{nu2const}}, \code{\link{nu2nu}}, \code{\link{omega2const}}, and \code{\link{omega2omega}}
#' @author Raphael Hartmann
#' @export
a2a <- function(model, names, keep_consts = FALSE) {
  
  if (!inherits(model, "drtmpt_model")) stop("model must be of class \"drtmpt_model\".")
  if (!("lines" %in% names(model)) || !("params" %in% names(model)) || !("responses" %in% names(model))) stop("No valid model file.")
  
  params_list <- model$params
  if (!is.list(params_list)) stop("params must be a list.")
  if ( !("threshold" %in% names(params_list)) ||
       !("driftrate" %in% names(params_list)) ||
       !("startpoint" %in% names(params_list)) ) stop("\"params\" must contain \"threshold\", \"driftrate\" and \"startpoint\".")
  if (!is.vector(names)) stop("\"names\" must be a vector")
  if (!is.character(names)) stop("\"names\" must be characters.")
  if (any(!(names %in% names(params_list$threshold)))) stop("\"names\" do not match with the names of the processes in the \"model\".")
  if (!is.logical(keep_consts) & !is.numeric(keep_consts)) stop("\"keep_consts\" must either be a logical value or a numeric value smaller or equal the length of \"names\".")
  if (length(keep_consts)!=1) stop("\"keep_consts\" must be of length one")
  if (length(unique(names)) == 1) stop("\"names\" must be of length two or larger")
  if (keep_consts > length(names)) stop("\"keep_consts\" invalid")
  names <- unique(names)
  
  sorted_names <- sort(names)
  if (keep_consts == FALSE) {
    params_list[["threshold"]][names] <- sorted_names[1]
    params_list[["threshold"]][sorted_names[1]] <- NA
  } else if (is.logical(keep_consts) && keep_consts == TRUE) {
    val <- unlist(params_list[["threshold"]][sorted_names[1]])
    if (is.na(val)) {
      params_list[["threshold"]][names] <- sorted_names[1]
      params_list[["threshold"]][sorted_names[1]] <- NA
    } else if (val %in% names(params_list[["threshold"]])) {
      params_list[["threshold"]][names] <- as.character(val)
      params_list[["threshold"]][as.character(val)] <- NA
    } else if (is.numeric(val)) {
      params_list[["threshold"]][names] <- as.numeric(val)
    }
  } else if (is.numeric(keep_consts)) {
    val <- as.numeric(params_list[["threshold"]][names[keep_consts]])
    if (is.na(val)) {
      params_list[["threshold"]][names] <- sorted_names[1]
      params_list[["threshold"]][sorted_names[1]] <- NA
    } else if (val %in% names(params_list[["threshold"]])) {
      params_list[["threshold"]][names] <- sorted_names[1]
      params_list[["threshold"]][sorted_names[1]] <- NA
    } else if (is.numeric(val)) {
      params_list[["threshold"]][names] <- as.numeric(val)
    }
  }
  
  model$params <- params_list
  return(model)
}


#' @rdname a2a
#' @examples
#'
#' ## make do = dn
#' new_model <- set_a_equal(model = model, names = c("do", "dn"))
#' new_model
#' @export
set_a_equal <- a2a






############ MAKE NU EQUAL ############

#' Set process drift rates equal
#'
#' Setting multiple process drift rates (nu) equal. One of the process drift rates will be estimated and the other
#'   named drift rates will be set to equal the former. The equality can be removed by only using one name of a process.
#'
#' @param model A list of the class \code{drtmpt_model}.
#' @param names Character vector giving the names of the processes for which the process drift rates should be equal. If
#'   \code{length(names) = 1} then the corresponding process drift rates will be estimated (i.e., it will be set to NA)
#' @param keep_consts Can be one of the following
#'   \itemize{
#'    \item logical value: \code{FALSE} (default) means none of the constants for \code{names} in the \code{model} will be kept; The drift rates of
#'      the reference process (i.e., first of \code{names} in alphabetical order) will be set to \code{NA} (i.e., will be estimated) and the others
#'      will be set to the name of the reference process (i.e., will be set to equal the reference process drift rate). \code{TRUE} means
#'      the constant of the reference process drift rate (if specified) is used for all other processes.
#'    \item numeric value: index for \code{names}. If 1, the constant of the first process in \code{names} (in original order defined by the user) is
#'      used for all other drift rates of the processes in \code{names}. If 2, the constant of the second process is used. And so on.
#'   }
#' @return A list of the class \code{drtmpt_model}.
#' @examples
#' ####################################################################################
#' # Detect-Guess variant of the Two-High Threshold model.
#' # The encoding and motor execution times are assumed to be equal for each category.
#' # The process drift rates for both detection processes ("do" and "dn")
#' # will be set equal.
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
#' model <- to_drtmpt_model(mdl_file = mdl_2HTM)
#'
#' ## make do = dn
#' new_model <- nu2nu(model = model, names = c("do", "dn"))
#' new_model
#'
#' @seealso \code{\link{delta2delta}}, \code{\link{a2const}}, \code{\link{a2a}}, \code{\link{nu2const}}, \code{\link{omega2const}}, and \code{\link{omega2omega}}
#' @author Raphael Hartmann
#' @export
nu2nu <- function(model, names, keep_consts = FALSE) {
  
  if (!inherits(model, "drtmpt_model")) stop("model must be of class \"drtmpt_model\".")
  if (!("lines" %in% names(model)) || !("params" %in% names(model)) || !("responses" %in% names(model))) stop("No valid model file.")
  
  params_list <- model$params
  if (!is.list(params_list)) stop("params must be a list.")
  if ( !("threshold" %in% names(params_list)) ||
       !("driftrate" %in% names(params_list)) ||
       !("startpoint" %in% names(params_list)) ) stop("\"params\" must contain \"threshold\", \"driftrate\" and \"startpoint\".")
  if (!is.vector(names)) stop("\"names\" must be a vector")
  if (!is.character(names)) stop("\"names\" must be characters.")
  if (any(!(names %in% names(params_list$driftrate)))) stop("\"names\" do not match with the names of the processes in the \"model\".")
  if (!is.logical(keep_consts) & !is.numeric(keep_consts)) stop("\"keep_consts\" must either be a logical value or a numeric value smaller or equal the length of \"names\".")
  if (length(keep_consts)!=1) stop("\"keep_consts\" must be of length one")
  if (length(unique(names)) == 1) stop("\"names\" must be of length two or larger")
  if (keep_consts > length(names)) stop("\"keep_consts\" invalid")
  names <- unique(names)
  
  sorted_names <- sort(names)
  if (keep_consts == FALSE) {
    params_list[["driftrate"]][names] <- sorted_names[1]
    params_list[["driftrate"]][sorted_names[1]] <- NA
  } else if (is.logical(keep_consts) && keep_consts == TRUE) {
    val <- unlist(params_list[["driftrate"]][sorted_names[1]])
    if (is.na(val)) {
      params_list[["driftrate"]][names] <- sorted_names[1]
      params_list[["driftrate"]][sorted_names[1]] <- NA
    } else if (val %in% names(params_list[["driftrate"]])) {
      params_list[["driftrate"]][names] <- as.character(val)
      params_list[["driftrate"]][as.character(val)] <- NA
    } else if (is.numeric(val)) {
      params_list[["driftrate"]][names] <- as.numeric(val)
    }
  } else if (is.numeric(keep_consts)) {
    val <- as.numeric(params_list[["driftrate"]][names[keep_consts]])
    if (is.na(val)) {
      params_list[["driftrate"]][names] <- sorted_names[1]
      params_list[["driftrate"]][sorted_names[1]] <- NA
    } else if (val %in% names(params_list[["driftrate"]])) {
      params_list[["driftrate"]][names] <- sorted_names[1]
      params_list[["driftrate"]][sorted_names[1]] <- NA
    } else if (is.numeric(val)) {
      params_list[["driftrate"]][names] <- as.numeric(val)
    }
  }
  
  model$params <- params_list
  return(model)
}


#' @rdname nu2nu
#' @examples
#'
#' ## make do = dn
#' new_model <- set_nu_equal(model = model, names = c("do", "dn"))
#' new_model
#' @export
set_nu_equal <- nu2nu






############ MAKE OMEGA EQUAL ############

#' Set process relaitve starting-point equal
#'
#' Setting multiple process relaitve starting-points (omegas) equal. One of the process relaitve starting-points will be estimated and the
#'   other named relaitve starting-points will be set to equal the former. The equality can be removed by only using one name of a process.
#'
#' @param model A list of the class \code{drtmpt_model}.
#' @param names Character vector giving the names of the processes for which the process relaitve starting-points should be equal. If
#'   \code{length(names) = 1} then the corresponding process relaitve starting-point will be estimated (i.e., it will be set to NA)
#' @param keep_consts Can be one of the following
#'   \itemize{
#'    \item logical value: \code{FALSE} (default) means none of the constants for \code{names} in the \code{model} will be kept; The relaitve starting-points of
#'      the reference process (i.e., first of \code{names} in alphabetical order) will be set to \code{NA} (i.e., will be estimated) and the others
#'      will be set to the name of the reference process (i.e., will be set to equal the reference process relative starting-point). \code{TRUE} means
#'      the constant of the reference process relaitve starting-point (if specified) is used for all other processes.
#'    \item numeric value: index for \code{names}. If 1, the constant of the first process in \code{names} (in original order defined by the user) is
#'      used for all other relaitve starting-points of the processes in \code{names}. If 2, the constant of the second process is used. And so on.
#'   }
#' @return A list of the class \code{drtmpt_model}.
#' @examples
#' ####################################################################################
#' # Detect-Guess variant of the Two-High Threshold model.
#' # The encoding and motor execution times are assumed to be equal for each category.
#' # The process relaitve starting-points for both detection processes ("do" and "dn")
#' # will be set equal.
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
#' model <- to_drtmpt_model(mdl_file = mdl_2HTM)
#'
#' ## make do = dn
#' new_model <- omega2omega(model = model, names = c("do", "dn"))
#' new_model
#'
#' @seealso \code{\link{delta2delta}}, \code{\link{a2const}}, \code{\link{a2a}}, \code{\link{nu2const}}, \code{\link{nu2nu}}, and \code{\link{omega2const}}
#' @author Raphael Hartmann
#' @export
omega2omega <- function(model, names, keep_consts = FALSE) {
  
  if (!inherits(model, "drtmpt_model")) stop("model must be of class \"drtmpt_model\".")
  if (!("lines" %in% names(model)) || !("params" %in% names(model)) || !("responses" %in% names(model))) stop("No valid model file.")
  
  params_list <- model$params
  if (!is.list(params_list)) stop("params must be a list.")
  if ( !("threshold" %in% names(params_list)) ||
       !("driftrate" %in% names(params_list)) ||
       !("startpoint" %in% names(params_list)) ) stop("\"params\" must contain \"threshold\", \"driftrate\" and \"startpoint\".")
  if (!is.vector(names)) stop("\"names\" must be a vector")
  if (!is.character(names)) stop("\"names\" must be characters.")
  if (any(!(names %in% names(params_list$startpoint)))) stop("\"names\" do not match with the names of the processes in the \"model\".")
  if (!is.logical(keep_consts) & !is.numeric(keep_consts)) stop("\"keep_consts\" must either be a logical value or a numeric value smaller or equal the length of \"names\".")
  if (length(keep_consts)!=1) stop("\"keep_consts\" must be of length one")
  if (length(unique(names)) == 1) stop("\"names\" must be of length two or larger")
  if (keep_consts > length(names)) stop("\"keep_consts\" invalid")
  names <- unique(names)
  
  sorted_names <- sort(names)
  if (keep_consts == FALSE) {
    params_list[["startpoint"]][names] <- sorted_names[1]
    params_list[["startpoint"]][sorted_names[1]] <- NA
  } else if (is.logical(keep_consts) && keep_consts == TRUE) {
    val <- unlist(params_list[["startpoint"]][sorted_names[1]])
    if (is.na(val)) {
      params_list[["startpoint"]][names] <- sorted_names[1]
      params_list[["startpoint"]][sorted_names[1]] <- NA
    } else if (val %in% names(params_list[["startpoint"]])) {
      params_list[["startpoint"]][names] <- as.character(val)
      params_list[["startpoint"]][as.character(val)] <- NA
    } else if (is.numeric(val)) {
      params_list[["startpoint"]][names] <- as.numeric(val)
    }
  } else if (is.numeric(keep_consts)) {
    val <- as.numeric(params_list[["startpoint"]][names[keep_consts]])
    if (is.na(val)) {
      params_list[["startpoint"]][names] <- sorted_names[1]
      params_list[["startpoint"]][sorted_names[1]] <- NA
    } else if (val %in% names(params_list[["startpoint"]])) {
      params_list[["startpoint"]][names] <- sorted_names[1]
      params_list[["startpoint"]][sorted_names[1]] <- NA
    } else if (is.numeric(val)) {
      params_list[["startpoint"]][names] <- as.numeric(val)
    }
  }
  
  model$params <- params_list
  return(model)
}


#' @rdname omega2omega
#' @examples
#'
#' ## make do = dn
#' new_model <- set_omegas_equal(model = model, names = c("do", "dn"))
#' new_model
#' @export
set_omegas_equal <- omega2omega

