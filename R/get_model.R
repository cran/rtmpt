
#' Create a model list to fit an RT-MPT
#' 
#' Create a model list of the class \code{rtmpt_model} by providing either \code{eqn_file} or \code{mdl_file}.
#' If both are provided \code{mdl_file} will be used.
#' 
#' @param eqn_file Character string as shown in example 2 or path to the text file that specifies the 
#'   (RT-)MPT model with standard .eqn syntax (Heck et al., 2018; Hu, 1999). E.g. \code{studied ; hit ; (1-do)*g} for a correct guess 
#'   in the detect-guess 2HT model.
#' @param mdl_file Character string as shown in example 1 or path to the text file that specifies the 
#'   (RT-)MPT model and gives on each line the equation of one category using \code{+} to separate branches 
#'   and \code{*} to separate processes (Singmann and Kellen, 2013). E.g. \code{do+(1-do)*g} for the category "hit" in the detect-guess 
#'   2HT model.
#' @return A list of the class \code{rtmpt_model}.
#' @references
#' Heck, D. W., Arnold, N. R., & Arnold, D. (2018). TreeBUGS: An R package for hierarchical 
#'   multinomial-processing-tree modeling. \emph{Behavior Research Methods, 50(1)}, 264-284.
#'
#' Hu, X. (1999). Multinomial processing tree models: An implementation. 
#'   \emph{Behavior Research Methods, Instruments, & Computers, 31(4)}, 689-695.
#'
#' Singmann, H., & Kellen, D. (2013). MPTinR: Analysis of multinomial processing tree models in R. 
#'   \emph{Behavior Research Methods, 45(2)}, 560-575.
#' @examples 
#' ########################################################################################
#' # Detect-Guess variant of the Two-High Threshold model
#' #   with constant guessing and
#' #   suppressed process completion times for both failed detections.
#' # The encoding and motor execution times are assumed to be different for each response.
#' ########################################################################################
#' 
#' ## 1. using the mdl syntax
#' mdl_2HTM <- "
#' # targets
#' do+(1-do)*g     ; 0
#' (1-do)*(1-g)    ; 1
#'
#' # lures
#' (1-dn)*g        ; 0
#' dn+(1-dn)*(1-g) ; 1
#' 
#' # do: detect old; dn: detect new; g: guess
#' 
#' # OPTIONAL MPT CONSTRAINTS
#' #   set probabilities to constants:
#' const_prob: g=0.5
#' 
#' #   suppress process times:
#' suppress_process: dn-, do-
#' "
#' 
#' model <- to_rtmpt_model(mdl_file = mdl_2HTM)
#' model
#' 
#' ## 2. using the eqn syntax
#' eqn_2HTM <- "
#' # CORE MPT EQN
#' # tree ; cat ; mpt
#'      0 ;   0 ; do
#'      0 ;   0 ; (1-do)*g
#'      0 ;   1 ; (1-do)*(1-g)
#'        
#'      1 ;   2 ; (1-dn)*g
#'      1 ;   3 ; dn
#'      1 ;   3 ; (1-dn)*(1-g)
#' 
#' # OPTIONAL MPT CONSTRAINTS
#' #   set probabilities to constants:
#' const_prob: g=0.5
#' 
#' #   suppress process times:
#' suppress_process: dn-, do-
#' 
#' #     tree ; cat ;  MAP
#' resp:    0 ;   0 ;    0
#' resp:    0 ;   1 ;    1
#' resp:    1 ;   2 ;    0
#' resp:    1 ;   3 ;    1
#' # different motor execution times for old and new responses.
#' "
#' 
#' model <- to_rtmpt_model(eqn_file = eqn_2HTM)
#' model
#' 
#' @note Within a branch of a (RT-)MPT model it is not allowed to have the same process two or more times.
#' @seealso \code{\link{delta2delta}}, \code{\link{theta2const}}, \code{\link{tau2zero}}, \code{\link{theta2theta}}, and \code{\link{tau2tau}} for 
#'   functions to change the model
#' @author Raphael Hartmann
#' @export
to_rtmpt_model <- function(eqn_file = NULL, mdl_file = NULL) {
  
  if (is.null(eqn_file) && is.null(mdl_file)) stop("Neither an eqn- nor a mdl-file specified.")
  
  form <- ifelse(!is.null(mdl_file), 1, 2)
  
  if (!is.null(mdl_file)) {
    if (grepl(pattern = "\n", x = mdl_file)) {line_char <- strsplit(x = mdl_file, split = "\n")[[1]]
    } else line_char <- readLines(mdl_file)
  } else {
    if (grepl(pattern = "\n", x = eqn_file)) {line_char <- strsplit(x = eqn_file, split = "\n")[[1]]
    } else line_char <- readLines(eqn_file)
  }
  
  
  membership <- get_membership(line_char = line_char, form = form)
  
  
  
  RAW_MODEL <- make_raw_model(line_char = line_char, membership = membership, form = form)
  model_lines <- make_model(RAW_MODEL = RAW_MODEL, save_model = FALSE, form = form)
  
  
  ordered_probs <- get_ordered_probs(RAW_MODEL = RAW_MODEL, form = form)
  
  if(form == 1) {
    raw_model <- RAW_MODEL$resp
  } else if(form == 2) {
    raw_model <- RAW_MODEL$eqn
  }
  
  wrong_trees <- check_one(raw_model = raw_model, variables = ordered_probs)
  len_wt <- length(wrong_trees)
  if(len_wt != 1 | any(wrong_trees != 0)) {
    stop(paste0(ifelse(len_wt > 1, "Trees ", "Tree "), paste0(wrong_trees, collapse = " and "), " with ", ifelse(len_wt > 1, "labels ", "label ") , paste0("\'", unique(raw_model$TREE)[wrong_trees], "\'", collapse = " and "), ifelse(len_wt > 1, " do", " does"), " not sum to 1"))
  }
  
  
  probabilities <- get_probs(RAW_MODEL = RAW_MODEL, ordered_probs = ordered_probs)
  taus <- get_taus(RAW_MODEL = RAW_MODEL, ordered_probs = ordered_probs)
  
  params_list <- list()
  params_list$probs <- probabilities
  params_list$taus <- taus
  
  
  responses <- get_resp(RAW_MODEL = RAW_MODEL, form = form)
  
  
  model <- list()
  model$lines <- model_lines$lines
  model$params <- params_list
  model$responses <- responses$responses
  if(suppressWarnings(all(!is.na(as.numeric(model$responses$CAT))))) {
    model$responses$CAT <- as.numeric(model$responses$CAT)
  }
  
  class(model) <- "rtmpt_model"
  
  
  # test model for structure
  mdl_txt <- gsub("\\\\", "/", tempfile(pattern = "model", tmpdir = tempdir(), fileext = ".txt"))
  mdl_info <- gsub("\\\\", "/", tempfile(pattern = "model", tmpdir = tempdir(), fileext = ".info"))
  infofile <- get_infofile(model, mdl_txt = mdl_txt, mdl_info = mdl_info)
  
  
  # return
  return(model)
  
}

#' @export
print.rtmpt_model <- function(x, ...) {
  cat("\nMODEL OVERVIEW\n\n")
  
  cat("\nMDL syntax for the MPT part:\n")
  print(data.frame(MDL_lines=x$lines[1:(length(x$lines)-1)]))
  cat("\n* NOTE 1: Each line in the MDL syntax represents one response category.\n")
  cat("----------------------------\n")
  
  cat("\nProcess probability parameters (thetas):\n")
  print(x$params$probs)
  cat("\n* NOTE 1: NA means the parameter will be estimated.\n") 
  cat("* NOTE 2: A value larger than 0 and smaller than 1 means it will be held constant.\n")
  cat("* INFO:", length(which(!is.na(x$params$probs))), "of the process probabilities will be held constant.\n")
  cat("-------------------------------\n")
  
  cat("\nProcess completion time parameters (taus):\n")
  print(x$params$taus)
  cat("\n* NOTE 1: \"minus\" refers to the negative outcome (1-P) and \"plus\" to the positive.\n")
  cat("* NOTE 2: NA means the parameter will be estimated. 0 means it will be suppressed.\n")
  cat("* INFO:", length(which(x$params$taus==0)), "of the process completion times will be suppressed.\n")
  cat("-----------------------------------\n")
  
  cat("\nMapping of response categories and encoding plus motor execution times (deltas):\n")
  print(x$responses)
  cat("\n* NOTE 1: Unique representation of trees and categories.\n")
  cat("* NOTE 2: Each mapping number corresponds to a distinct encoding plus motor execution time.\n")
  cat("* INFO:", max(x$responses$MAP)+1,"distinct delta(s) assumed, namely with mapping [", paste(unique(x$responses$MAP), collapse = ", "),"].\n")
  cat("-----------------------------------\n\n")
}


