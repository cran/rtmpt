
get_ordered_probs <- function(RAW_MODEL, form) {
  
  if (form == 1) {
    
    # Npaths <- nchar(RAW_MODEL$resp$MDL) - nchar(gsub("[+]", "", RAW_MODEL$resp$MDL)) + 1
    # max_paths <- max(Npaths)
    # path_df <- matrix(data = NA, nrow = length(Npaths), ncol = max_paths)
    # for (i in 1:length(Npaths)) {
    #   paths_string <- strsplit(RAW_MODEL$resp$MDL[i], "[+]")[[1]]
    #   if (Npaths[i] != max_paths) {paths_string <- c(paths_string, rep("", max_paths-Npaths[i]))}
    #   path_df[i,] <- paths_string
    # }
    # 
    # Nprobs <- apply(X = path_df, MARGIN = 2, FUN = function(x) {nchar(x) - nchar(gsub("[*]", "", x)) + 1})
    # max_probs <- max(Nprobs)
    # probs_df <- matrix(data = NA, nrow = length(Npaths), ncol = max_paths*max_probs)
    # for (i in 1:dim(path_df)[2]) {
    #   ind_start <- max_probs*(i-1)+1
    #   ind_end <- max_probs*i
    #   probs_string <- gsub(pattern = "[)]", replacement = "", x = gsub(pattern = "[(1-]", replacement = "", x = path_df[,i]))
    #   probs_temp <- sapply(X = probs_string, FUN = function(x) {strsplit(x, "[*]")[[1]]})
    #   probs_df[, ind_start:ind_end] <- t(sapply(X = probs_temp, FUN = function(x) {if (length(x) != max_probs) {c(x, rep("", max_probs-length(x)))} else x}))
    # }
    # 
    # 
    # temp <- as.vector(probs_df)
    # ind <- which(temp == "")
    # vec <- temp[-ind]
    
    vec <- as.character(unlist(sapply(X = RAW_MODEL$resp$MDL, FUN = function(x) {strsplit(x = gsub(pattern = "\\)", replacement = "", x = gsub(pattern = "\\(1-", replacement = "", x = x)), split = "[*]")[[1]]})))
    ordered_probs <- sort(unique(as.character(unlist(sapply(X = vec, FUN = function(x) {strsplit(x = x, split = "[+]")[[1]]})))))
    
  } else if (form == 2) {
    
    # unique_TREE_CAT <- unique(RAW_MODEL$eqn[, c(1,2)])
    # max_len <- max(sapply(X = RAW_MODEL$eqn$EQN, FUN = function(x) {length(strsplit(x, "[*]")[[1]])}))
    # probs_list_unordered <- lapply(X = RAW_MODEL$eqn$EQN, FUN = function(x) {
    #   vec <- rep("", max_len);
    #   y1 <- gsub("[(1-]", "", x);
    #   y2 <- gsub("[)]", "", y1); cc <- strsplit(y2, "[*]")[[1]];
    #   vec[1:length(cc)] <- cc; return(vec)
    # })
    # max_len_tot <- max_len*max(sapply(X = unique_TREE_CAT[,2], FUN = function(x) {length(which(RAW_MODEL$eqn$CAT == x))}))
    # probs_list <- lapply(X = seq(1,dim(unique_TREE_CAT)[1]), FUN = function(x) {
    #   vec <- c();
    #   ind1 <- which(RAW_MODEL$eqn[, 1]==unique_TREE_CAT[x, 1]);
    #   ind2 <- which(RAW_MODEL$eqn[, 2]==unique_TREE_CAT[x, 2]);
    #   ind <- intersect(ind1, ind2);
    #   vec <- unlist(lapply(X = ind, FUN = function(x) {probs_list_unordered[[x]]}));
    #   if (length(vec == max_len_tot)) vec <- c(vec, rep("", (max_len_tot-length(vec))));
    #   return(vec)
    # })
    # vec <- c()
    # for (j in 1:max_len_tot) {
    #   for (i in 1:length(probs_list)) {
    #     if (probs_list[[i]][j] != "") vec <- c(vec, probs_list[[i]][j])
    #   }
    # }
    
    ordered_probs <- sort(unique(as.character(unlist(sapply(X = RAW_MODEL$eqn$EQN, FUN = function(x) {strsplit(x = gsub(pattern = "\\)", replacement = "", x = gsub(pattern = "\\(1-", replacement = "", x = x)), split = "[*]")[[1]]})))))
    
  }
  
  
  
  # ordered_probs <- unique(vec)
  
  # CONTROLS
  # WARNINGS FOR NOT MATCHING PROBABILITY PARAMETERS
  if ("const_probs" %in% names(RAW_MODEL)) {
    const_probs <- RAW_MODEL$const_probs
    for ( i in 1:length(const_probs)) {
      if (grepl(pattern = "=", x = const_probs[i])) {
        ending <- unlist(gregexpr(pattern = "=", const_probs[i]))-1
        const_probs[i] <- substring(text = const_probs[i], first = 1, last = ending)
      }
      if (!(const_probs[i] %in% ordered_probs)) {
        warning("Definition of const_probs contain process name that is not in equations.\n Found \"", const_probs[i], "\" not to be used in equations. It will be ignored.")
      }
    }
    # for ( i in 1:length(ordered_probs)) {
    #   if (!(ordered_probs[i] %in% const_probs)) {
    #     warning("Equations contain more probability parameters than definition of const_probs.\n Found \"", unlist(probs_list)[i], "\" not in definition. It will be included with default value.")
    #   }
    # }
  }
  
  return(ordered_probs)
  
}