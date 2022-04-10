#' @importFrom stringr str_count
make_raw_model <- function(line_char, membership, form) {
  
  fnc1 <- function(row_char, start) {
    sub_char <- substring(text = row_char, first = start, last = 10000)
    vec_char <- strsplit(gsub("[[:space:]]", "", sub_char), ",")[[1]]
    return(vec_char)
  }
  fnc2 <- function(row_char, start) {
    sub_char <- substring(text = row_char, first = start, last = 10000)
    vec_char <- strsplit(gsub("[[:space:]]", "", sub_char), ",")[[1]]
    return(vec_char)
  }
  fnc3 <- function(row_char, start) {
    sub_char <- substring(text = row_char, first = start, last = 10000)
    if (2 == str_count(sub_char, ";")) {
      vec_char <- strsplit(gsub("[[:space:]]", "", sub_char), ";")[[1]]
    }
    if (2 == str_count(sub_char, ",")) {
      vec_char <- strsplit(gsub("[[:space:]]", "", sub_char), ",")[[1]]
    }
    return(vec_char)
  }
  fnc4 <- function(row_char) {
    if (2 == str_count(row_char, ";")) {
      vec_char <- strsplit(gsub("[[:space:]]", "", row_char), ";")[[1]]
    }
    if (2 == str_count(row_char, ",")) {
      vec_char <- strsplit(gsub("[[:space:]]", "", row_char), ",")[[1]]
    }
    return(vec_char)
  }
  fnc5 <- function(row_char) {
    output <- list(text = "", resp = NA)
    if (grepl(pattern = "[;]", x = row_char)) {
      vec_char <- strsplit(gsub("[[:space:]]", "", row_char), ";")[[1]]
      output$resp <- as.numeric(vec_char[2])
      output$text <- vec_char[1]
    } else if (grepl(pattern = "[,]", x = row_char)) {
      vec_char <- strsplit(gsub("[[:space:]]", "", row_char), ",")[[1]]
      output$resp <- as.numeric(vec_char[2])
      output$text <- vec_char[1]
    } else if (grepl(pattern = "#", x = row_char) || row_char == "") {
      output$text <- ""
    } else {
      vec_char <- gsub(" ", "", row_char)
      output$text <- vec_char
      output$resp <- 0
    }
    return(output)
  }
  
  model <- list()
  m <- list()
  if (form == 1) {
    for ( i in 1:3 ) {
      index <- which(membership == i)
      if (i == 1 && length(index) == 1) {
        for ( ind in 1:length(index) ) {
          row_char <- line_char[index[ind]]
          start <- max(unlist(gregexpr(pattern =":", row_char)))+1
          const_probs <- fnc1(row_char = row_char, start = start)
		  if(length(const_probs) > 0) {
		    model$const_probs <- const_probs
		  }
        }
      }
      if (i == 2 && length(index) == 1) {
        for ( ind in 1:length(index) ) {
          row_char <- line_char[index[ind]]
          start <- max(unlist(gregexpr(pattern =":", row_char)))+1
          suppr_taus <- fnc2(row_char = row_char, start = start)
		  if(length(suppr_taus) > 0) {
		    model$suppr_taus <- suppr_taus
		  }
        }
      }
      if (i == 3) {
        model$resp <- data.frame(TREE=rep(NA, length(index)), CAT=NA, MDL=NA, MAP=NA)
        len <- max(index) - min(index) + 1
        m$resp <- rep(0, len)
        m$mdl <- rep("", len)
        for ( ind in min(index):max(index) ) {
          row_char <- line_char[ind]
          if (ind %in% index) {
            output <- fnc5(row_char = row_char)
            m$mdl[ind-min(index)+1] <- output$text
            m$resp[ind-min(index)+1] <- output$resp
          } else {
            m$mdl[ind-min(index)+1] <- ""
            m$resp[ind-min(index)+1] <- NA
          }
        }
      }
    }
    ind <- 1; tree_val <- 0; cat_val <- 0; tree <- c()
    for (i in 1:length(m$mdl)) {
      if (m$mdl[i] != "") {
        model$resp$MAP[ind] <- m$resp[i]
        model$resp$MDL[ind] <- m$mdl[i]
        tree[ind] <- tree_val
        model$resp$CAT[ind] <- cat_val
        ind <- ind + 1
        cat_val <- cat_val + 1
      } else {
        tree_val <- tree_val + 1
      }
    }
    remove(m)
    uniq_tree <- unique(tree)
    uniq_new <- seq(1, length(uniq_tree)) - 1
    model$resp$TREE <- uniq_new[sapply(X = tree, FUN = function(x) {which(uniq_tree==x)})]


    
    return(model)
    
    
  } else if (form == 2) {
    for ( i in 1:4) {
      index <- which(membership == i)
      if (i == 1 && length(index) == 1) {
        for ( ind in 1:length(index) ) {
          row_char <- line_char[index[ind]]
          start <- max(unlist(gregexpr(pattern =":", row_char)))+1
          const_probs <- fnc1(row_char = row_char, start = start)
		  if(length(const_probs) > 0) {
		    model$const_probs <- const_probs
		  }
        }
      }
      if (i == 2 && length(index) == 1) {
        for ( ind in 1:length(index) ) {
          row_char <- line_char[index[ind]]
          start <- max(unlist(gregexpr(pattern =":", row_char)))+1
          suppr_taus <- fnc2(row_char = row_char, start = start)
		  if(length(suppr_taus) > 0) {
		    model$suppr_taus <- suppr_taus
		  }
        }
      }
      if (i == 3 && length(index) > 0) {
        model$resp <- data.frame(TREE=NA, CAT=NA, MAP=NA)
        for ( ind in 1:length(index) ) {
          row_char <- line_char[index[ind]]
          start <- max(unlist(gregexpr(pattern =":", row_char)))+1
          vec_char <- fnc3(row_char = row_char, start = start)
          model$resp[ ind, c(1,2)] <- vec_char[c(1,2)]
          model$resp[ ind, 3] <- as.numeric(vec_char[3])
        }
      }
      if (i == 4) {
        model$eqn <- data.frame(TREE=NA, CAT=NA, EQN=NA)
        for ( ind in 1:length(index) ) {
          row_char <- line_char[index[ind]]
          model$eqn[ ind, ] <- fnc4(row_char = row_char)
        }
      }
    }
    
    uniq_tree <- unique(model$eqn$TREE)
    uniq_cat <- unique(model$eqn$CAT)
    sorted_tree <- sapply(X = model$eqn$TREE, FUN = function(x) {which(uniq_tree==x)})
    sorted_cat <- sapply(X = model$eqn$CAT, FUN = function(x) {which(uniq_cat==x)})
    if (suppressWarnings(all(!is.na(as.numeric(model$eqn$CAT))))) {
      sorted_cat <- as.numeric(model$eqn$CAT)
    }
    model$eqn <- model$eqn[order(sorted_tree, sorted_cat),]
    
    if ("resp" %in% names(model)) {
      if (suppressWarnings(all(!is.na(as.numeric(model$resp$CAT))))) {
        sorted_tree <- sapply(X = model$resp$TREE, FUN = function(x) {which(uniq_tree==x)})
        model$resp$CAT <- as.numeric(model$resp$CAT)
        model$resp <- model$resp[order(sorted_tree, model$resp$CAT), ]
      }
    }
    
    
    # CONTROLS
    if ("resp" %in% names(model)) {
      if (!all(sort(unique(as.character(model$resp$TREE))) == sort(unique(as.character(model$eqn$TREE))))) {stop("TREE-Labels in eqn and resp do not match.")}
      if (!all(sort(unique(as.character(model$resp$CAT))) == sort(unique(as.character(model$eqn$CAT))))) {stop("CAT-Labels in eqn and resp do not match.")}
      if (!all(names(model$params$probs) == names(model$params$taus))) stop("Order of parameters for \"probs\" and \"taus\" in \"model\" not allowed.")
    }
    
    return(model)
    
  }
  
}




make_model <- function(RAW_MODEL, save_model = FALSE, form) {
  
  if (form == 1) {
    
    uniq_tree <- unique(RAW_MODEL$resp$TREE)
    ind <- lapply(X = uniq_tree, FUN = function(x) {which(RAW_MODEL$resp$TREE == x)})
    all_lines <- as.vector(sapply(X = ind, FUN = function(x) {c(paste0(RAW_MODEL$resp$MDL[x]), "")}))
    
  } else if (form == 2) {
    
    Rows <- dim(RAW_MODEL$eqn)[1]
    char_vec <- vector(mode = "list", length = Rows)
    for (i in 1:Rows) {
      char_vec[[i]] <- paste0(RAW_MODEL$eqn$TREE[i], " ", RAW_MODEL$eqn$CAT[i])
    }
    uniq_char_vec <- unique(char_vec)
    mapping_tree <- as.numeric(sapply(X = RAW_MODEL$eqn$TREE, FUN = function(x) {which(unique(RAW_MODEL$eqn$TREE) == x)}))
    mapping_cat <- sapply(X = char_vec, FUN = function(x) {which(uniq_char_vec == x)})
    model_lines <- sapply(X = seq(1, length(uniq_char_vec)), FUN = function(x) {ind <- which(mapping_cat==x); paste(RAW_MODEL$eqn$EQN[ind], collapse = "+")})
    comment_lines <- sapply(X = seq(1, length(unique(mapping_tree))), FUN = function(x) {
      ind <- which(mapping_tree==x); 
      paste("#", RAW_MODEL$eqn$TREE[ind[1]], "[", paste(unique(RAW_MODEL$eqn$CAT[ind]), collapse = " , "), "]")
    })
    
    all_lines <- c()
    index <- 1
    for (i in 1:length(comment_lines)) {
      ind <- which(mapping_tree==i)
      pat <- paste0(RAW_MODEL$eqn$TREE[ind[1]], " ")
      ind2 <- which(grepl(pattern = pat, x = substr(uniq_char_vec, 1, nchar(pat))))
      all_lines[index:(index+length(ind2)+1)] <- c(comment_lines[i], model_lines[ind2], "")
      index <- index+length(ind2)+2
    }
    
    # if (save_model) {
    #   last_dot <- max(unlist(gregexpr(pattern ="[.]", eqn_file)))
    #   filename <- paste0(substring(text = eqn_file, first = 1, last = (last_dot-1)), ".model")
    #   writeLines(all_lines, con = filename, sep = "\n")
    # }
    
  }
  
  model <- list()
  model$lines <- all_lines
  return(model)
  
}

#' @importFrom stringr str_replace_all
#' @importFrom Ryacas yac_str
check_one <- function(raw_model, variables) {
  
  new_vars <- paste0("var", 1:length(variables))
  
  unq_tree <- unique(raw_model[,1])
  
  nmbr <- vector(mode = "list", length = length(unq_tree))
  
  for(i in 1:length(nmbr)) {
    
    nmbr[[i]] <- raw_model[which(raw_model[,1] == unq_tree[i]), 3]
    for(j in 1:length(new_vars)) {
      nmbr[[i]] <- unlist(str_replace_all(nmbr[[i]], variables[j], new_vars[j]))
    }
    nmbr[[i]] <- paste(nmbr[[i]], collapse = "+")
    
    nmbr[[i]] <- yac_str(paste0("Simplify(", nmbr[[i]], ")"))
    
  }
  
  resp <- NULL
  if(all(unlist(nmbr)=="1")) {
    resp <- 0
  } else {
    resp <- which(unlist(nmbr)!="1")
  }
  return(resp)
  
}

