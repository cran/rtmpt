

############ THETA TO CONSTANT ############

#' Set process probabilities to constants
#' 
#' Setting process probabilities (thetas) to constants or change it back to be estimated. 
#'
#' @param model An object of the class \code{rtmpt_model}.
#' @param names Character vector with process names.
#' @param constants Numerical vector of length one or \code{length(names)}. You have the following options for the elements of the numeric vector:
#'   \itemize{
#'    \item \code{0 < constants < 1}: set the named probability to a constant value between zero and one
#'    \item \code{NA}: estimate the named probability
#'   }
#' @return An object of the class \code{rtmpt_model}.
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
#' model <- to_rtmpt_model(mdl_file = mdl_2HTM)
#' 
#' ## setting g to a constant (0.5):
#' new_model <- theta2const(model = model, names = c("g"), constants = c(0.5))
#' new_model
#'                     
#' @seealso \code{\link{delta2delta}}, \code{\link{tau2zero}}, \code{\link{theta2theta}} and \code{\link{tau2tau}}
#' @author Raphael Hartmann
#' @export
theta2const <- function(model, names, constants = NA) {
  
  if (!inherits(model, "rtmpt_model")) stop("model must be of class \"rtmpt_model\".")
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
#' @param model A list of the class \code{rtmpt_model}.
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
#' ## removing the process times (tau) for the failed (="minus") detection ("do" and "dn") 
#' new_model <- tau2zero(model = model, names = c("dn", "do"),
#'                       outcomes = c("minus", "minus"), values = 0)
#' new_model
#' 
#' @seealso \code{\link{delta2delta}}, \code{\link{theta2const}}, \code{\link{theta2theta}} and \code{\link{tau2tau}}
#' @author Raphael Hartmann
#' @export
tau2zero <- function(model, names, outcomes, values = 0) {
  
  if (!inherits(model, "rtmpt_model")) stop("model must be of class \"rtmpt_model\".")
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
#' @param model A list of the class \code{rtmpt_model}.
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
#' @return A list of the class \code{rtmpt_model}.
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
#' model <- to_rtmpt_model(mdl_file = mdl_2HTM)
#' 
#' ## make do = dn
#' new_model <- theta2theta(model = model, names = c("do", "dn"))
#' new_model
#' 
#' @seealso \code{\link{delta2delta}}, \code{\link{theta2const}}, \code{\link{tau2zero}} and \code{\link{tau2tau}}
#' @author Raphael Hartmann
#' @export
theta2theta <- function(model, names, keep_consts = FALSE) {
  
  if (!inherits(model, "rtmpt_model")) stop("model must be of class \"rtmpt_model\".")
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
#' @param model A list of the class \code{rtmpt_model}.
#' @param names Character vector giving the names of the processes for which the process completion times should be equal. If 
#'   \code{length(names) = 1} then the corresponding process completion times (for negative and positive outcomes) will be 
#'   estimates (i.e., they will be set to NA)
#' @param keep_zeros Can be one of the following
#'   \itemize{
#'    \item logical value: \code{FALSE} (default) means none of the zeros for \code{names} in the \code{model} will be kept; The times of
#'      the reference process (i.e., first of \code{names} in alphabetical order) will be set to \code{NA} (i.e., will be estimated) and the others
#'      will be set to the name of the reference process (i.e., will be set to equal the reference process times). \code{TRUE} means
#'      the zero(s) of the reference process times (if specified) is used for the same outcome of all other processes. 
#'    \item numeric value: index for \code{names}. If 1, the zero(s) of the first process in \code{names} (in original order defined by the user) is 
#'      used for the same outcome of all other processes in \code{names}. If 2, the zero(s) of the second process is used. And so on.
#'   }
#' @return A list of the class \code{rtmpt_model}.
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
#' model <- to_rtmpt_model(mdl_file = mdl_2HTM)
#' 
#' ## make do = dn
#' new_model <- tau2tau(model = model, names = c("do", "dn"))
#' new_model
#' 
#' @seealso \code{\link{delta2delta}}, \code{\link{theta2const}}, \code{\link{tau2zero}} and \code{\link{theta2theta}}
#' @author Raphael Hartmann
#' @export
tau2tau <- function(model, names, keep_zeros = FALSE) {
  
  if (!inherits(model, "rtmpt_model")) stop("model must be of class \"rtmpt_model\".")
  if (!("lines" %in% names(model)) || !("params" %in% names(model)) || !("responses" %in% names(model))) stop("No valid model file.")
  
  params_list <- model$params
  if (!is.list(params_list)) stop("params must be a list.")
  if ( !("probs" %in% names(params_list)) || 
       !("taus" %in% names(params_list)) ) stop("\"params\" must contain \"probs\" and \"taus\".")
  if (!is.vector(names)) stop("\"names\" must be a vector")
  if (!is.character(names)) stop("\"names\" must be characters.")
  if (any(!(names %in% names(params_list$probs)))) stop("\"names\" do not match with the names of the processes in the \"model\".")
  if (!is.logical(keep_zeros) & !is.numeric(keep_zeros)) stop("\"keep_zeros\" must either be a logical value or a numeric value smaller or equal the length of \"names\".")
  if (length(keep_zeros)!=1) stop("\"keep_zeros\" must be of length one")
  if (length(unique(names)) == 1) stop("\"names\" must be of length two or larger")
  names <- unique(names)
  
  sorted_names <- sort(names)
  if (keep_zeros == FALSE) {
    
    params_list[["taus"]][, names] <- sorted_names[1]
    params_list[["taus"]][, sorted_names[1]] <- NA
    
  } else if (keep_zeros == TRUE) {
    
    val <- params_list[["taus"]][, sorted_names[1]]

    if (any(is.na(val))) {
      if (all(is.na(val))) {
        params_list[["taus"]][, names] <- sorted_names[1]
        params_list[["taus"]][, sorted_names[1]] <- NA
      } else if (any(val==0)) {
        params_list[["taus"]][, names] <- sorted_names[1]
        params_list[["taus"]][, sorted_names[1]] <- NA
        ind_zero <- which(val == 0)
        params_list[["taus"]][ind_zero, names] <- val[ind_zero]
      }
    } else if (any(val %in% names(params_list[["taus"]]))) {
      if (val[1]==val[2]) {
        params_list[["taus"]][, names] <- sorted_names[1]
        params_list[["taus"]][, sorted_names[1]] <- NA
      } else if (any(val == 0)) {
        ind_zero <- which(val == 0)
        ind_proc <- which(as.character(val) %in% names(params_list[["taus"]]))
        params_list[["taus"]][, names] <- sorted_names[1]
        params_list[["taus"]][ind_proc, sorted_names[1]] <- NA
        params_list[["taus"]][ind_zero, sorted_names[1]] <- 0
        params_list[["taus"]][ind_zero, names] <- 0
      } else if (all(val %in% names(params_list[["taus"]])) & val[1]!=val[2]) {
        stop("Something is not right with the model. There cannot be two different process names in one column.")
      }
    } else if (all(val == 0)) {
      params_list[["taus"]][, names] <- 0
    }
    
  } else {
    
    val <- params_list[["taus"]][, names[keep_zeros]]
    
    if (any(is.na(val))) {
      if (all(is.na(val))) {
        params_list[["taus"]][, names] <- sorted_names[1]
        params_list[["taus"]][, sorted_names[1]] <- NA
      } else if (any(val==0)) {
        params_list[["taus"]][, names] <- sorted_names[1]
        params_list[["taus"]][, sorted_names[1]] <- NA
        ind_zero <- which(val == 0)
        params_list[["taus"]][ind_zero, names] <- val[ind_zero]
      }
    } else if (any(val %in% names(params_list[["taus"]]))) {
      if (val[1]==val[2]) {
        params_list[["taus"]][, names] <- sorted_names[1]
        params_list[["taus"]][, sorted_names[1]] <- NA
      } else if (any(val == 0)) {
        ind_zero <- which(val == 0)
        ind_proc <- which(as.character(val) %in% names(params_list[["taus"]]))
        params_list[["taus"]][, names] <- sorted_names[1]
        params_list[["taus"]][ind_proc, sorted_names[1]] <- NA
        params_list[["taus"]][ind_zero, sorted_names[1]] <- 0
        params_list[["taus"]][ind_zero, names] <- 0
      } else if (all(val %in% names(params_list[["taus"]])) & val[1]!=val[2]) {
        stop("Something is not right with the model. There cannot be two different process names in one column.")
      }
    } else if (all(val == 0)) {
      params_list[["taus"]][, names] <- 0
    }
    
  }
  
  for(i in 1:length(params_list[["taus"]][1,])) {
    if(all(!params_list[["taus"]][, i] %in% names(params_list[["taus"]]))) {
      params_list[["taus"]][, i] <- as.numeric(params_list[["taus"]][, i] )
    }
  }
  
  model$params <- params_list
  return(model)
}


#' @rdname tau2tau
#' @examples
#' 
#' ## make do = dn
#' new_model <- set_taus_equal(model = model, names = c("do", "dn"))
#' new_model
#' @export
set_taus_equal <- tau2tau








############ SUPPRESS TAUS AND SET PROBS TO CONSTANTS ############

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
#' @seealso \code{\link{delta2delta}}
#' @author Raphael Hartmann
#' @export
set_params <- function(model, parameter, names, values = NA) {
  
  warning("this function is deprecated and will be removed in a future version. Please use theta2const and tau2zero instead.")
  
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
