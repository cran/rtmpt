
get_probs <- function(RAW_MODEL, ordered_probs) {
  estim <- NA
  probs_df <- data.frame(matrix(rep(estim, length(ordered_probs)), ncol = length(ordered_probs), nrow = 1))
  colnames(probs_df) <- ordered_probs
  #rownames(probs_df) <- "const"
  if ("const_probs" %in% names(RAW_MODEL)) {
    const_probs <- RAW_MODEL$const_probs
    const <- rep(estim, length(const_probs))
    for ( i in 1:length(const_probs)) {
      if (grepl(pattern = "=", x = const_probs[i])) {
        ending <- unlist(gregexpr(pattern = "=", const_probs[i]))-1
        const[i] <- as.numeric(substring(text = const_probs[i], first = (ending+2), last = 100))
        const_probs[i] <- substring(text = const_probs[i], first = 1, last = ending)
        if ((const[i] > 1) || (const[i] < 0)) {
          warning(const_probs[i], " = ", const[i], " is out of boundaries: Ignoring constant.")
          const[i] <- estim
        }
      }
    }
    
    
    for (i in 1:length(ordered_probs)) {
      if (ordered_probs[i] %in% const_probs) {
        ind <- which(const_probs == ordered_probs[i])[1]
        probs_df[1, i] <- const[ind]
      }
    }
    return(probs_df)
  } else return(probs_df)
  
}





get_taus <- function(RAW_MODEL, ordered_probs) {
  estim <- NA
  suppr <- 0
  rates_df <- data.frame(matrix(rep(estim, 2*length(ordered_probs)), ncol = length(ordered_probs), nrow = 2))
  colnames(rates_df) <- ordered_probs
  rownames(rates_df) <- c("minus", "plus")
  if ("suppr_taus" %in% names(RAW_MODEL)) {
    suppr_taus <- RAW_MODEL$suppr_taus
    suppr <- rep(estim, length(suppr_taus))
    for ( i in 1:length(suppr_taus)) {
      pattern = "[-]"
      if (grepl(pattern = pattern, x = suppr_taus[i])) {
        ending <- unlist(gregexpr(pattern = pattern, suppr_taus[i]))-1
        suppr[i] <- "minus"
        suppr_taus[i] <- substring(text = suppr_taus[i], first = 1, last = ending)
      }
      pattern = "[+]"
      if (grepl(pattern = pattern, x = suppr_taus[i])) {
        ending <- unlist(gregexpr(pattern = pattern, suppr_taus[i]))-1
        suppr[i] <- "plus"
        suppr_taus[i] <- substring(text = suppr_taus[i], first = 1, last = ending)
      }
    }
    
    
    for (i in 1:length(ordered_probs)) {
      if (ordered_probs[i] %in% suppr_taus) {
        ind <- which(suppr_taus == ordered_probs[i])
        if (length(ind) > 2) warning("\"", suppr_taus[ind[1]], "\" found to be suppressed more than twice. Only first two are used. Please check your equation file.")
        for (j in 1:length(ind)) {
          if (suppr[ind[j]] == "minus") rates_df[1, i] <- 0 else if (suppr[ind[j]] == "plus") rates_df[2, i] <- 0
        }
      }
    }
    return(rates_df)
  } else return(rates_df)
  
}
  
  
  
  

get_thresholds <- function(RAW_MODEL, ordered_probs) {
  thresh_df <- data.frame(matrix(rep(NA, length(ordered_probs)), ncol = length(ordered_probs), nrow = 1))
  colnames(thresh_df) <- ordered_probs
  return(thresh_df)
}

get_driftrate <- function(RAW_MODEL, ordered_probs) {
  drift_df <- data.frame(matrix(rep(NA, length(ordered_probs)), ncol = length(ordered_probs), nrow = 1))
  colnames(drift_df) <- ordered_probs
  return(drift_df)
}

get_startpoint <- function(RAW_MODEL, ordered_probs) {
  start_df <- data.frame(matrix(rep(NA, length(ordered_probs)), ncol = length(ordered_probs), nrow = 1))
  colnames(start_df) <- ordered_probs
  return(start_df)
}




