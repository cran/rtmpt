

get_resp <- function(RAW_MODEL, form) {
  
  resp_list <- list(responses=NULL)
  responses <- c()
  
  if (form == 1) {
    
    resp_values <- RAW_MODEL$resp$RESP
    min_resp <- min(unlist(resp_values))
    for (i in 1:length(resp_values)) responses[i] <- as.numeric(resp_values[i]) - min_resp
    if (max(responses) != (length(unique(responses))-1)) stop("Elements in \"resp\" must range from 0 to max(\"resp\").")
    resp_list$responses <- RAW_MODEL$resp

  } else if (form == 2) {
    
    uniq_resps <- unique(RAW_MODEL$eqn[, c(1,2)])
    if (suppressWarnings(all(!is.na(as.numeric(uniq_resps$CAT))))) {
      uniq_resps$CAT <- as.numeric(uniq_resps$CAT)
    }
    
    if ("resp" %in% names(RAW_MODEL)) {
      
      if (dim(RAW_MODEL$resp[, c(1,2)])[1] == dim(uniq_resps)[1]) {
        resp_values <- as.numeric(apply(X = uniq_resps, MARGIN = 1, FUN = function(x) {
          ind <- which(apply(X = RAW_MODEL$resp[, c(1,2)], MARGIN = 1, FUN = function(y) {all(x==y)}))
          return(as.numeric(RAW_MODEL$resp[ind, 3]))
        }))
        min_resp <- min(unlist(resp_values))
        for (i in 1:dim(uniq_resps)[1]) responses[i] <- as.numeric(resp_values[i]) - min_resp
        if (max(responses) != (length(unique(responses))-1)) stop("\"resp\" must be a sequence.")
        resp_list$responses <- uniq_resps
      } else warning("Incorrect number of rows in \"resp\". Must be ", dim(uniq_resps)[1], ". \"resp\" in equation file will be ignored.")
      
    } else if (suppressWarnings(all(!is.na(as.numeric(uniq_resps$CAT))))) {
      
      if(length(unique(uniq_resps$CAT)) < dim(uniq_resps)[1]) {
        resp_values <- as.numeric(uniq_resps$CAT)
        min_resp <- min(resp_values)
        for (i in 1:dim(uniq_resps)[1]) responses[i] <- resp_values[i] - min_resp
        resp_list$responses <- uniq_resps
      } else for (i in 1:dim(uniq_resps)[1]) {responses[i] <- 0; resp_list$responses <- uniq_resps}
      
    } else for (i in 1:dim(uniq_resps)[1]) {responses[i] <- 0; resp_list$responses <- uniq_resps}
    
  }

  resp_list$responses$RESP <- responses
  return(resp_list)
  
}

