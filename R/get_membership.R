
get_membership <- function(line_char, form) {
  
  if (form == 1) {
    
    INorOUT <- function(row_char) {
      # prob_key <- c("prob", "porbabilities", "probability")
      # minus_key <- c("lambda_minus", "rate_minus", "process_minus", "tau_minus")
      # plus_key <- c("lambda_plus", "rate_plus", "process_minus", "tau_plus")
      suppr_key <- c("suppr","suppress_lambda", "suppress_rate", "suppress_process", "suppress_tau",
                     "suppr_lambda", "suppr_rate", "suppr_process", "suppr_tau", "suppress")
      const_key <- c("const_prob", "constant_prob", "const_probabilities", "constant_probabilities", "const", "constant")
      resp_key <- c("resp", "gamma", "mean")
      if ( (grepl(pattern = "#", x = row_char)) || row_char == "" ) {
        return(0)
      } else if ( any(sapply(X = const_key, FUN = function(y) {grepl(pattern = y, x = row_char)} )) ) {
        return(1)
      } else if ( any(sapply(X = suppr_key, FUN = function(y) {grepl(pattern = y, x = row_char)} )) ) {
        return(2)
      } else {
        return(3)
      }
    }
    
    membership <- as.numeric(sapply(X = line_char, FUN = INorOUT))
    
    # CONTROLS
    if (length(which(membership==1))>1) {stop("Too many lines for constant probabilities.")}
    if (length(which(membership==2))>1) {stop("Too many lines for suppressing process times.")}
    if (length(which(membership==3))<4) {warning("There seem to be too few equations.")}
    
  } else if (form == 2) {
    
    INorOUT <- function(row_char) {
      # prob_key <- c("prob", "porbabilities", "probability")
      # minus_key <- c("lambda_minus", "rate_minus", "process_minus", "tau_minus")
      # plus_key <- c("lambda_plus", "rate_plus", "process_minus", "tau_plus")
      suppr_key <- c("suppr", "suppress_lambda", "suppress_rate", "suppress_process", "suppress_tau",
                     "suppr_lambda", "suppr_rate", "suppr_process", "suppr_tau", "suppress")
      const_key <- c("const_prob", "constant_prob", "const_probabilities", "constant_probabilities", "const")
      resp_key <- c("resp", "gamma", "mean", "responses")
      if ( (grepl(pattern = "#", x = row_char)) || row_char == "" ) {
        return(0)
      } else if ( any(sapply(X = const_key, FUN = function(y) {grepl(pattern = y, x = row_char)} )) ) {
        return(1)
      } else if ( any(sapply(X = suppr_key, FUN = function(y) {grepl(pattern = y, x = row_char)} )) ) {
        return(2)
      } else if ( any(sapply(X = resp_key, FUN = function(y) {grepl(pattern = y, x = row_char)} )) ) {
        return(3)
      } else if ((3 == length(strsplit(gsub(" ", "", row_char), ";")[[1]])) | (3 == length(strsplit(gsub(" ", "", row_char), ",")[[1]]))) {
        return(4)
      } else {
        return(0)
      }
    }
    
    membership <- as.numeric(sapply(X = line_char, FUN = function(x) INorOUT(x)))
    
    # CONTROLS
    if (length(which(membership==1))>1) {stop("Too many lines for constant probabilities.")}
    if (length(which(membership==2))>1) {stop("Too many lines for suppressing process times.")}
    if (length(which(membership==4))<5) {warning("There seem to be too few equations.")}
    
  }
  
  
  return(membership)
  
}