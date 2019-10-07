#################################
# Code carried together from    #
# MPTinR2 by Henrik Sinmann     #
# henrik.singmann@warwick.ac.uk #
# slightly modyfied by          #
# Raphael Hartmann              #
#################################



##############################
# function from
# mpt.classes.R
##
## S4 object definitions.
##
##############################
# restrictions
setClass("restriction", representation = representation(parameter = "character", "VIRTUAL"))
setClass("fixed.restriction", representation = representation(value = "numeric"), contains = "restriction")
# setClass("equality.restriction", representation = representation(value = "character"), contains = "restriction")
# setClass("inequality.restriction", representation = representation(exchange.parameter = "list", exchange.inverse = "list", compute.as = "expression"), contains = "restriction")
setClass("restrictions", representation = representation(fixed = "list", equality = "list", inequality = "list", raw = "list"))
setClassUnion("restrictionsOrNull", c("restrictions", "NULL"))

# MPT models
setClass("Rmpt.model", representation = representation(Rinitial.model = "list", Rcheck = "list", restrictions = "restrictionsOrNull", Rinitial.model.data.frame = "data.frame", model.data.frame = "data.frame", model.list = "list"))
setClass("Rbmpt.model", representation = representation(A = "array", B = "array", lbmpt = "character"), contains = "Rmpt.model")

# fitted MPT models
# setClass("mpt", representation(model = "Rmpt.model", observed.data = "array",  predicted.data = "array", C.matrix = "list", g2 = "numeric", log.likelihood = "numeric", estimates = "matrix", multifit = "logical", hessian = "list", default.ci = "numeric"))
# setClass("bmpt", representation(typeHessian = "character", hfail = "numeric", fia = "list", parametric.ci = "list", nonparametric.ci = "list"), contains = "mpt")

## make.mpt helper classes
setOldClass( c("file", "connection" ) )
setOldClass( c("url", "connection" ) )
setOldClass( c("textConnection", "connection" ) )
setClassUnion("characterOrConnection", c("character", "connection"))
# setClassUnion("characterOrNull", c("character", "NULL"))
#setClassUnion("listOrNull", c("list", "NULL"))


########################
# functions from
# methods.for.models.R
########################
if(!isGeneric("Rcheck")){
  if (is.function("Rcheck"))
    fun <- Rcheck
  else fun <- function(object) standardGeneric("Rcheck")
  setGeneric("Rcheck", fun)
}
setMethod("Rcheck", "Rmpt.model", function(object) object@Rcheck)
setGeneric("Rcheck<-", function(x, value) standardGeneric("Rcheck<-"))
setReplaceMethod("Rcheck", "Rmpt.model", function(x, value) {
  if (!is.list(value)) stop("the model Rcheck must be a list")
  x@Rcheck <- value
  x
})

if(!isGeneric("model.list")){
  if (is.function("model.list"))
    fun <- model.list
  else fun <- function(object) standardGeneric("model.list")
  setGeneric("model.list", fun)
}
setMethod("model.list", "Rmpt.model", function(object) object@model.list)
setGeneric("model.list<-", function(x, value) standardGeneric("model.list<-"))
setReplaceMethod("model.list", "Rmpt.model", function(x, value) {
  if (!is.list(value)) stop("the model.list must be a list")
  x@model.list <- value
  x
})

if(!isGeneric("model.data.frame")){
  if (is.function("model.data.frame"))
    fun <- model.data.frame
  else fun <- function(object) standardGeneric("model.data.frame")
  setGeneric("model.data.frame", fun)
}
setMethod("model.data.frame", "Rmpt.model", function(object) object@model.data.frame)
setGeneric("model.data.frame<-", function(x, value) standardGeneric("model.data.frame<-"))
setReplaceMethod("model.data.frame", "Rmpt.model", function(x, value) {
  if (!is.data.frame(value)) stop("the model.data.frame must be a data.frame")
  x@model.data.frame <- value
  x
})

if(!isGeneric("Rinitial.model.data.frame")){
  if (is.function("Rinitial.model.data.frame"))
    fun <- Rinitial.model.data.frame
  else fun <- function(object) standardGeneric("Rinitial.model.data.frame")
  setGeneric("Rinitial.model.data.frame", fun)
}
setMethod("Rinitial.model.data.frame", "Rmpt.model", function(object) object@Rinitial.model.data.frame)
setGeneric("Rinitial.model.data.frame<-", function(x, value) standardGeneric("Rinitial.model.data.frame<-"))
setReplaceMethod("Rinitial.model.data.frame", "Rmpt.model", function(x, value) {
  if (!is.data.frame(value)) stop("the Rinitial.model.data.frame must be a data.frame")
  x@Rinitial.model.data.frame <- value
  x
})

if(!isGeneric("Rinitial.model")){
  if (is.function("Rinitial.model")) fun <- Rinitial.model
  else fun <- function(object) standardGeneric("Rinitial.model")
  setGeneric("Rinitial.model", fun)
}
setMethod("Rinitial.model", "Rmpt.model", function(object) object@Rinitial.model)
setGeneric("Rinitial.model<-", function(x, value) standardGeneric("Rinitial.model<-"))
setReplaceMethod("Rinitial.model", "Rmpt.model", function(x, value) {
  if (!is.list(value)) stop("the initial model must be a list")
  x@Rinitial.model <- value
  x
})

if(!isGeneric("A")){
	if (is.function("A"))
		fun <- A
	else fun <- function(object) standardGeneric("A")
	setGeneric("A", fun)
}
setMethod("A", "Rbmpt.model", function(object) object@A)
setGeneric("A<-", function(x, value) standardGeneric("A<-"))
setReplaceMethod("A", "Rbmpt.model", function(x, value) {
	if (!is.array(value)) stop("A must be an array")
	x@A <- value
	x
})

if(!isGeneric("B")){
	if (is.function("B"))
		fun <- B
	else fun <- function(object) standardGeneric("B")
	setGeneric("B", fun)
}
setMethod("B", "Rbmpt.model", function(object) object@B)
setGeneric("B<-", function(x, value) standardGeneric("B<-"))
setReplaceMethod("B", "Rbmpt.model", function(x, value) {
	if (!is.array(value)) stop("B must be an array")
	x@B <- value
	x
})

if(!isGeneric("lbmpt")){
	if (is.function("lbmpt"))
		fun <- lbmpt
	else fun <- function(object) standardGeneric("lbmpt")
	setGeneric("lbmpt", fun)
}
setMethod("lbmpt", "Rbmpt.model", function(object) object@lbmpt)
setGeneric("lbmpt<-", function(x, value) standardGeneric("lbmpt<-"))
setReplaceMethod("lbmpt", "Rbmpt.model", function(x, value) {
	if (!is.character(value)) stop("lbmpt must be a character vector")
	x@lbmpt <- value
	x
})


##############################
# function from
# methods.for.restrictions.R
##############################
if(!isGeneric("restrictions")){
  if (is.function("restrictions"))
    fun <- restrictions
  else fun <- function(object) standardGeneric("restrictions")
  setGeneric("restrictions", fun)
}
setMethod("restrictions", "Rmpt.model", function(object) object@restrictions)
setGeneric("restrictions<-", function(x, value) standardGeneric("restrictions<-"))
setReplaceMethod("restrictions", "Rmpt.model", function(x, value) {
  if (class(value) != "restrictions" | is.null(value)) stop("the model restrictions must be of class restrictions or MULL")
  x@restrictions <- value
  x
})


######################
# functions from
# helper.functions.R
######################
.find.MPT.params <- function(model) {
  inobjects <- function(ex) {
    split.exp1 <- strsplit(as.character(ex),"[[:space:]]")
    split.exp2 <- sapply(split.exp1, strsplit, split = "[()+*-]")
    return(sort(unique(grep("[[:alpha:]]",unlist(split.exp2), value = TRUE))))
  }
  tmp <- sapply(model,inobjects)
  return(unique(sort(unique(unlist(tmp)))))
}

.count.branches <- function(model.df) {
  counters <- rep(1, max(model.df[,"category"]))
  for (branch in 1:dim(model.df)[1]) {
    tmp.cat <- model.df[branch,"category"]
    model.df[branch,"n.branch"] <- counters[tmp.cat]
    counters[tmp.cat] <- counters[tmp.cat] + 1
  }
  model.df
}


#####################
# functions from
# make.mpt.helper.R
#####################
.Rcheck.MPT.probabilities <- function(tree){
  tmp.env <- new.env()
  temp.param.names <- .find.MPT.params(tree)
  temp.param.val <- runif(length(temp.param.names))
  temp.branch <- sapply(tree,length)
  prob <- rep(NA,length(temp.branch))
  
  for (i in 1:length(temp.param.val)) {
    assign(temp.param.names[i],temp.param.val[i], envir = tmp.env)
  }
  temp.check <- sapply(unlist(tree),eval, envir = tmp.env)
  for (i in 1:length(temp.branch)){
    if (i==1) prob[1] <- sum(temp.check[1:temp.branch[1]])
    else prob[i] <- sum(temp.check[(sum(temp.branch[1:(i-1)])+1):sum(temp.branch[1:i])])
  }
  prob <- round(prob,digits=6)
  return(prob)
}

.make.model.list <- function(model.df) {
  parse.eqn <- function(x){
    branches <- unique(x[,2])
    l.tree <- length(branches)
    tree <- vector('expression', l.tree)
    for (branch in 1:l.tree) {
      tree[branch] <- parse(text = paste(x[x[,2] == branches[branch],"branches"], collapse = " + "))
    }
    tree
  }
  tmp.ordered <- model.df[order(model.df[,1]),]
  tmp.spl <- split(tmp.ordered, factor(tmp.ordered[,1]))
  tmp.spl <- lapply(tmp.spl, function(d.f) d.f[order(d.f[,2]),])
  model <- lapply(tmp.spl, parse.eqn)
  names(model) <- NULL
  model
}

.make.model.df <- function(model) {

  #require(stringr)
  oneLineDeparse <- function(expr){
    paste(deparse(expr), collapse="")
  }
  
  n.trees <- length(model)
  l.trees <- sapply(model, length)
  l.trees <- c(0, l.trees)
  
  fin.model <- vector("list", n.trees)
  
  for (tree in 1:n.trees) {
    utree <- unlist(model[[tree]])
    tree.df.unordered <- do.call("rbind",lapply(1:length(utree), function(t) data.frame(category = t, branches = oneLineDeparse(utree[[t]]), stringsAsFactors = FALSE)))
    
    tree.list <- vector("list", dim(tree.df.unordered)[1])
    for (c1 in 1:length(tree.list)) {
      category <- tree.df.unordered[c1,"category"]
      branch <- strsplit(tree.df.unordered[c1,"branches"], "\\+")
      branch <- gsub(" ", "", branch[[1]])
      tree.list[[c1]] <- data.frame(tree = tree, category = category, branches = branch, stringsAsFactors = FALSE)
    }
    tree.df <- do.call("rbind", tree.list)
    fin.model[[tree]] <- tree.df[rev(order(tree.df[["branches"]])),]
  }
  n.categories <- c(0,sapply(fin.model, function(x) max(x[["category"]])))
  n.cat.cumsum <- cumsum(n.categories)
  
  model.df <- do.call("rbind", fin.model)
  
  model.df[["category"]] <- model.df[,"category"] + n.cat.cumsum[model.df[,"tree"]]
  
  rownames(model.df) <- NULL
  
  .count.branches(model.df)
  
}

.make.Rmpt.cf <- function(model){
	
	bin.objects <- function(branch) {
		objects <- strsplit(branch, "\\*")[[1]]
		!(grepl("[()]", objects))
	}
	
	model.df.tmp <- model.data.frame(model)
	c.join <- 1
	
	while (length(unique(model.df.tmp[,"tree"])) > 1) {
		model.df.tmp[model.df.tmp$tree == 1, "branches"] <- paste("hank.join.", c.join, "*", model.df.tmp[model.df.tmp$tree == 1, "branches"], sep = "")
		model.df.tmp[model.df.tmp$tree == 2, "branches"] <- paste("(1-hank.join.", c.join, ")*", model.df.tmp[model.df.tmp$tree == 2, "branches"], sep = "")
		model.df.tmp[model.df.tmp$tree == 2, "tree"] <- rep(1, length(model.df.tmp[model.df.tmp$tree == 2, "tree"]))
		model.df.tmp[model.df.tmp$tree > 2, "tree"] <- model.df.tmp[model.df.tmp$tree > 2, "tree"] -1
		c.join <- c.join + 1
	}
	tree.ordered <- model.df.tmp
	tree.list <- lapply(1:dim(tree.ordered)[1], function(x) list(category = tree.ordered[x,"category"], branch = tree.ordered[x,"branches"], objects = strsplit(tree.ordered[x,"branches"], "\\*")[[1]], params = .find.MPT.params(tree.ordered[x,"branches"]), binary = bin.objects(tree.ordered[x,"branches"])))
	tmp.tree <- tree.list
	mpt.string <- c(tmp.tree[[1]][["objects"]], tmp.tree[[1]][["category"]])
	for (counter1 in 2:length(tmp.tree)) {
		if (length(tmp.tree[[counter1]][["binary"]]) == length(tmp.tree[[counter1-1]][["binary"]]) & tmp.tree[[counter1-1]][["binary"]][length(tmp.tree[[counter1-1]][["binary"]])] == TRUE & tmp.tree[[counter1]][["binary"]][length(tmp.tree[[counter1]][["binary"]])] == FALSE) {
			mpt.string <- c(mpt.string, tmp.tree[[counter1]][["category"]])
		} else {
		if (length(tmp.tree[[counter1]][["binary"]]) == length(tmp.tree[[counter1-1]][["binary"]]) & tmp.tree[[counter1-1]][["binary"]][length(tmp.tree[[counter1-1]][["binary"]])] == FALSE & tmp.tree[[counter1]][["binary"]][length(tmp.tree[[counter1]][["binary"]])] == TRUE) {
			change <- min(which((tmp.tree[[counter1]][["binary"]] == tmp.tree[[counter1-1]][["binary"]]) == FALSE))+1
			tmp.objects <- tmp.tree[[counter1]][["objects"]][change:(length(tmp.tree[[counter1]][["binary"]]))]
			mpt.string <- c(mpt.string, tmp.objects[tmp.tree[[counter1]][["binary"]][change:length(tmp.tree[[counter1]][["binary"]])]], tmp.tree[[counter1]][["category"]])
		} else {
		if (length(tmp.tree[[counter1]][["binary"]]) > length(tmp.tree[[counter1-1]][["binary"]])) {
			change <- min(which((tmp.tree[[counter1]][["binary"]] == tmp.tree[[counter1-1]][["binary"]][1:length(tmp.tree[[counter1]][["binary"]])]) == FALSE))+1
			if (change < (length(tmp.tree[[counter1-1]][["binary"]]))) {
				tmp.param <- tmp.tree[[counter1]][["objects"]][change:(length(tmp.tree[[counter1]][["binary"]]))]
			} else {
				tmp.new <- tmp.tree[[counter1]][["objects"]][(length(tmp.tree[[counter1-1]][["binary"]])):length(tmp.tree[[counter1]][["binary"]])]
				tmp.param <- tmp.new[tmp.tree[[counter1]][["binary"]][(length(tmp.tree[[counter1-1]][["binary"]])):length(tmp.tree[[counter1]][["binary"]])]]
			}
			mpt.string <- c(mpt.string, tmp.param, tmp.tree[[counter1]][["category"]])
		} else {
		if (length(tmp.tree[[counter1]][["binary"]]) < length(tmp.tree[[counter1-1]][["binary"]])) {
			change <- min(which((tmp.tree[[counter1]][["binary"]] == tmp.tree[[counter1-1]][["binary"]][1:length(tmp.tree[[counter1]][["binary"]])]) == FALSE))+1
			if (change <= length(tmp.tree[[counter1]][["binary"]])) {
			tmp.objects <- tmp.tree[[counter1]][["objects"]][change:(length(tmp.tree[[counter1]][["binary"]]))]
			} else tmp.objects <- NULL
			mpt.string <- c(mpt.string, tmp.objects[tmp.tree[[counter1]][["binary"]][change:length(tmp.tree[[counter1]][["binary"]])]], tmp.tree[[counter1]][["category"]])
		}
		}
		}
		}
	
	}
	mpt.string
}


##################
# functions from
# make.model.R
##################
.Rcheck.model <- function(model) {
  
  prob.tree.check <- .Rcheck.MPT.probabilities(model.list(model))
  if(all(prob.tree.check==1)) {
    prob.corr <- TRUE
  } else {
    prob.corr <- paste("Model not constructed well: Branch probabilities of tree(s) ", paste(which(prob.tree.check!=1), collapse= ", "), " do not sum to 1!", sep = "")
  }
  orig.params <- .find.MPT.params(model.list(model))	
  if (!is.null(restrictions(model))) {
	stop("restrictions of model should always be NULL\n")
    #for (restriction in fixed.restrictions(model)) orig.params <- orig.params[-which(parameter(restriction) == orig.params)]
    #fixed.parameters <- vapply(fixed.restrictions(model), parameter, "")
  } else {
    fixed.parameters <- NULL
  }
  l.orig.params <- length(orig.params)
  n.trees.orig <- length(model.list(model))
  n.categories <- length(unlist(model.list(model)))
  
  max.branches.per.category <- max(table(model.data.frame(model)[,2]))
  branches.per.category <- table(model.data.frame(model)[,2])
  
  original.parameters <- .find.MPT.params(Rinitial.model(model))
  
  suppressWarnings(lbmpt <- tryCatch(.make.Rmpt.cf(model), error = function(e) NULL))
  
  list(probabilities.eq.1 = prob.corr, n.trees = n.trees.orig, n.categories = n.categories, n.free.parameters = l.orig.params, free.parameters = orig.params, n.fixed.parameters = length(fixed.parameters), fixed.parameters = fixed.parameters, original.parameters = original.parameters, max.branches.per.category = max.branches.per.category, branches.per.category = branches.per.category, lbmpt = lbmpt)
}

setGeneric("make.Rmpt", function(model, restrictions = NULL, ...) standardGeneric("make.Rmpt"))
setMethod("make.Rmpt", signature(model = "characterOrConnection"), function(model, restrictions = NULL, model.type = c("easy", "eqn", "eqn2"), ...) {
  
  raw.model <- .read.mpt(model.filename = model, model.type = model.type)

  callGeneric(model = raw.model, restrictions = NULL, ...)
})
# setMethod("make.Rmpt", signature(model = "characterOrConnection", restrictions = "characterOrConnection"), function(model, restrictions = NULL, model.type = c("easy", "eqn", "eqn2"), ...) {
#   
#   raw.model <- .read.mpt(model.filename = model, model.type = model.type)
#   
#   restrictions <- .read.MPT.restrictions.file(restrictions)
# 
#   callGeneric(model = raw.model, restrictions = restrictions, ...)
# })
# setMethod("make.Rmpt", signature(model = "characterOrConnection", restrictions = "list"), function(model, restrictions, model.type = c("easy", "eqn", "eqn2"), ...) {
#   
#   raw.model <- .read.mpt(model.filename = model, model.type = model.type)
#   
#   restrictions <- .read.MPT.restrictions(restrictions)
# 
#   callGeneric(model = raw.model, restrictions = restrictions, ...)
# })
setMethod("make.Rmpt", signature(model = "list", restrictions = "restrictionsOrNull"), function(model, restrictions, ...) {
  
  Rinitial.model.list <- model
  
  model <- new("Rmpt.model", Rinitial.model = Rinitial.model.list, Rcheck = list(), restrictions = restrictions)
  
  Rinitial.model.data.frame(model) <- .make.model.df(Rinitial.model(model))
  
  if (is.null(restrictions(model))) {
	model.data.frame(model) <- Rinitial.model.data.frame(model)
  }
  else {
	stop("in RTMPT the restrictions of the internal model should always be NULL")
	# model.data.frame(model) <- .apply.restrictions(model)
  }
  
  model.list(model) <- .make.model.list(model.data.frame(model))
  
  first.checks <- .Rcheck.model(model)
  Rcheck(model) <- first.checks[-length(first.checks)]
  Rcheck(model)[["df"]] <- c(available = Rcheck(model)[["n.categories"]] - Rcheck(model)[["n.trees"]], used = Rcheck(model)[["n.free.parameters"]], model = (Rcheck(model)[["n.categories"]] - Rcheck(model)[["n.trees"]]) - Rcheck(model)[["n.free.parameters"]])
  if (!is.null(first.checks[["lbmpt"]])) {
	bmpt <- is.bmpt(first.checks[["lbmpt"]])
  }
  else {
	bmpt <- FALSE
  }
  if (bmpt) {
    #browser()
    model <- as(model, "Rbmpt.model")
    matrices <- .make.matrices(model)
    A(model) <- matrices[["A"]]
    storage.mode(A(model)) <- "integer"
    B(model) <- matrices[["B"]]
    storage.mode(B(model)) <- "integer"
    lbmpt(model) <- first.checks[["lbmpt"]]
    Rcheck(model)[["is.bmpt"]] <- TRUE
  }
  else {
	Rcheck(model)[["is.bmpt"]] <- FALSE
  }

  model
})

is.bmpt <- function(lbmpt) {
	
	is.category <- grepl("^[[:digit:]]+$", lbmpt)	
	type <- ifelse(is.category == 0, 1, 0)
	
	############################################################################################################
	## This code is adapted from Wu, Myung & Batchelder (2010a, 2010b) and based on Purdy & Batchelder (2009) ##
	############################################################################################################
	
	L <- length(lbmpt)
	code <- matrix(0, L, L)
	if (type[1] == 0 & L == 1) return(FALSE)	#This should return TRUE, but a model that small is uninteresting here (Henrik Singmann, 29-7-2011)
	if (type[1] == 0 & L != 1) return(FALSE)
	p <- 1
	u <- 1
	for (i in 2:L) {
		code[i,] <- code[p,]
		code[i,p] <- u
		if (type[i] == 1) {
			u <- 1
			p <- i
		} else {
			u <- -1
			ind <- i-1
			while (ind > 0) {
				if (ind <= 0 & i < L) return(FALSE)
				if (type[ind] == 1) {
					p <- ind
					break
				} else {
					if (type[ind] == 0) {
						if (type[ind-1] !=1) return(FALSE)
						type[c(ind-1,ind, ind+1)] <- -1
						ind <- ind-2
					} else {
						if (type[ind] == -1) {
							type[ind+1] <- -1
							while (type[ind] == -1) ind <- ind-1
							if (type[ind] != 1) return(FALSE)
							type[ind] <- -1
							ind <- ind-1
						}
					}
				}
			}
		}
	}
	if (ind > 0) return(FALSE)
	else return (TRUE)
}

.make.matrices <- function(model) {
	#browser()
	
	n.parameters <- Rcheck(model)[["n.free.parameters"]] + Rcheck(model)[["n.fixed.parameters"]]
	parameters <- c(Rcheck(model)[["free.parameters"]], Rcheck(model)[["fixed.parameters"]])
	
	A <- array(0, dim = c(Rcheck(model)[["n.categories"]], Rcheck(model)[["max.branches.per.category"]], n.parameters))
	dimnames(A)[[3]] <- parameters
	B <- A
	
	for (parameter in parameters) {
		for (branch in 1:dim(model.data.frame(model))[1]) {
			tmp.branch <- strsplit(model.data.frame(model)[branch,"branches"], split="\\*")[[1]]
			A[model.data.frame(model)[branch,"category"], model.data.frame(model)[branch,"n.branch"], parameter] <- sum(grepl(paste("^", parameter, "$", sep = ""), tmp.branch))
			B[model.data.frame(model)[branch,"category"], model.data.frame(model)[branch,"n.branch"], parameter] <- sum(grepl(paste("^\\(1-", parameter, "\\)$", sep = ""), tmp.branch))
		}
	}
	list(A = A, B = B)
}


##################
# functions from
# read.model.R
##################
.read.MPT.model <- function(model.filename) {
  whole <- readLines(model.filename)
  model <- vector("list", length(whole))
  c2 <- 1
  c3 <- 1
  s.flag <- FALSE
  for (c1 in 1:length(whole)) {
    if (!(grepl("^[[:space:]]*$", whole[c1]))) {
      if (grepl("^[[:space:]]*#", whole[c1])) next
      whole[c1] <- gsub("#.*", "", whole[c1])
      s.flag <- TRUE
      model[[c2]][c3] <- parse(text = whole[c1])[1]
      c3 <- c3 + 1
      fin <- c2
    }
    else {
      if (s.flag == TRUE) c2 <- c2 + 1
      c3 <- 1
      s.flag <- FALSE
    }
  }
  return (model[1:fin])
}

.read.mpt <- function(model.filename, model.type) {
  if (grepl("\\.eqn$", model.filename) || grepl("\\.EQN$", model.filename)) model.type <- "eqn"
  if (model.type[1] == "eqn") {
	stop("in RTMPT model.type cannot be eqn")
    #raw.model <- .read.EQN.model(model.filename)
  } else if (model.type[1] == "eqn2") {
	stop("in RTMPT model.type cannot be eqn2")
    #raw.model <- .read.EQN.model.2(model.filename)
  } else {
	raw.model <- .read.MPT.model(model.filename)
  }
  raw.model
}


####################
# recursive.tree.R
####################
setClassUnion("bmpt.leaf", "numeric")
setClassUnion("bmpt.tree", "bmpt.leaf")
setClass("strict.bmpt.tree", representation(top = "bmpt.tree", bottom = "bmpt.tree", parameter = "character", node = "numeric"))
setIs("strict.bmpt.tree", "bmpt.tree")

make.tree <- function(tree.df, envir) {
  if ((nrow(tree.df) == 1) & (tree.df[1, "branches"] == "")) return(tree.df[1, "category"])
  p <- str_split(tree.df[1,"branches"], "\\*")[[1]][1]
  counter <- 1
  #browser()
  for (c in 1:nrow(tree.df)) {
    if (str_detect(tree.df[c,"branches"], str_c("^", p, "\\b"))) counter <- counter + 1
    tree.df[c,"branches"] <- str_replace(tree.df[c,"branches"], str_c("^", p, "\\b\\*?|^\\(1-", p, "\\)\\**"), "")
  }
  #browser()   
  #node.counter <- node.counter + 1
  assign("node.counter", get("node.counter", envir = envir, inherits = TRUE) + 1, envir = envir, inherits = TRUE)
  #assign("node.counter", get("node.counter", envir = parent.frame(), inherits = TRUE) + 1, envir = parent.frame(), inherits = TRUE)
  #return(new("strict.bmpt.tree", node = node.counter, top = if (tree.df[1, "branches"] == "") tree.df[1, "category"] else make.tree(tree.df[seq_len(counter-1),], envir), bottom = if (tree.df[2, "branches"] == "") tree.df[2, "category"] else make.tree(tree.df[counter:nrow(tree.df),], envir), parameter = p))
  return(new("strict.bmpt.tree", node = get("node.counter", envir = envir, inherits = TRUE), 
		 top = if (tree.df[1, "branches"] == "") tree.df[1, "category"] else make.tree(tree.df[seq_len(counter-1),], envir), 
		 bottom = if (tree.df[2, "branches"] == "") tree.df[2, "category"] else make.tree(tree.df[counter:nrow(tree.df),], envir), parameter = p))
}

count.nodes <- function(tree) {
  count <- 1
  if (!is.numeric(tree@top)) count <- count + count.nodes(tree@top)
  if (!is.numeric(tree@bottom)) count <- count + count.nodes(tree@bottom)
  count
}

node.to.p <- function(tree) {
  return(c(tree@parameter, if (is.numeric(tree@top)) NULL else node.to.p(tree@top), if (is.numeric(tree@bottom)) NULL else node.to.p(tree@bottom)))
}

ar.node <- function(tree, x = list(NULL)) {
  if (is.numeric(tree)) {
    #browser()
    return(list(c(unlist(x), tree)))
  }
  xt <- lapply(x, function(y) c(y, -1))
  xb <- lapply(x, function(y) c(y, -3))
  #browser()
  for (c in 1:length(xt)) {
    names(xt[[c]])[length(xt[[c]])] <- tree@node
  }
  for (c in 1:length(xb)) {
    names(xb[[c]])[length(xb[[c]])] <- tree@node
  }
  return(c(ar.node(tree@top, xt), ar.node(tree@bottom, xb)))
}


#######################
# prep.hierarchical.R
#######################
prep.h <- function(mf, filename) {
  model <- make.Rmpt(mf)
  kernpar <- as.integer(Rcheck(model)[["n.free.parameters"]]+Rcheck(model)[["n.fixed.parameters"]])
  zweig <- as.integer(Rcheck(model)[["max.branches.per.category"]])
  # make recursive object:
  mo <- vector("list", Rcheck(model)[["n.trees"]])
  split.df <- split(model.data.frame(model), model.data.frame(model)$tree)
  #node.counter.env <- new.env(parent = .GlobalEnv)
  #browser()
  for (c in 1:Rcheck(model)[["n.trees"]]) {
    thisEnv <- environment()
	node.counter <- 0
    mo[[c]] <- make.tree(split.df[[c]], thisEnv)
  }
  node.per.tree <- vapply(mo, count.nodes, 0)
  branches.per.category <- Rcheck(model)[["branches.per.category"]]
  parameters <- c(Rcheck(model)[["free.parameters"]], Rcheck(model)[["fixed.parameters"]])
  l.nodes.to.p <- lapply(mo, node.to.p)
  tree.node.2.par <- matrix(-1, Rcheck(model)[["n.trees"]], max(node.per.tree))
  for (tree in 1:Rcheck(model)[["n.trees"]]) {
    for (node in 1:length(l.nodes.to.p[[tree]])) {
      tree.node.2.par[tree, node] <- which(parameters == l.nodes.to.p[[tree]][node])
    }    
  }
  AR <- array(0, dim = c(Rcheck(model)[["n.categories"]], Rcheck(model)[["max.branches.per.category"]], max(node.per.tree)))
  #browser()
  mo.nodes <- lapply(mo, ar.node)
  for (tree in mo.nodes) {
    tmp.cat <- 0
    tmp.tree <- tree[order(vapply(tree, function(x) x[length(x)], 0))]
    for (c.outer in 1:length(tmp.tree)) {
      if (tmp.cat != tmp.tree[[c.outer]][length(tmp.tree[[c.outer]])]) branch <- 1
      tmp.cat <- tmp.tree[[c.outer]][length(tmp.tree[[c.outer]])]
      tmp.num <- as.numeric(names(tmp.tree[[c.outer]]))
      for (c in 1:(length(tmp.tree[[c.outer]])-1)) {
        AR[tmp.cat, branch, tmp.num[c]] <- tmp.tree[[c.outer]][c] + 2
      }
      branch <- branch + 1
    }
  }
  #browser()
  cat(zweig, kernpar, max(node.per.tree), Rcheck(model)[["n.trees"]], Rcheck(model)[["n.categories"]], "\n", file = filename)
  cat(model.data.frame(model)[model.data.frame(model)$n.branch == 1,"tree"], "\n", append = TRUE, file = filename)
  cat(as.matrix(branches.per.category), "\n", append = TRUE, file = filename)
  write.table(tree.node.2.par, file = filename, col.names = FALSE, row.names = FALSE, append = TRUE)
  cat(node.per.tree, "\n", append = TRUE, file = filename)
  #cat(Rcheck(model)[["n.categories"]], Rcheck(model)[["max.branches.per.category"]], max(node.per.tree), "\n", file = "AR.txt")
  cat(AR, "\n", append = TRUE, file = filename)
  cat(Rcheck(model)[["original.parameters"]], "\n", append = TRUE, file = filename)
  #browser()
}


##################################
# use the functions from MPTinR2
##################################
get_infofile <- function(model, mdl_txt, mdl_info) {
  writeLines(text = model$lines, con = mdl_txt)
  
  prep.h(mdl_txt, mdl_info)
  # file.remove(mdl_txt)
  
  return(mdl_info)
}

