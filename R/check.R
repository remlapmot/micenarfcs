# internal function for checking input to main mice() function

check.visitSequence <- function(setup, where) {
  
  nwhere <- setup$nwhere
  nvar <- setup$nvar
  visitSequence <- setup$visitSequence

  # set default visit sequence, left to right
  if (is.null(visitSequence))
    visitSequence <- seq_len(ncol(where))[apply(where, 2, any)]
  
  if (!is.numeric(visitSequence)) {
    code <- match.arg(visitSequence, c("roman", "arabic", "monotone",
                                       "revmonotone"))
    visitSequence <- switch(
      code, 
      roman = seq_len(nvar)[nwhere > 0],
      arabic = rev(seq_len(nvar)[nwhere > 0]),
      monotone = order(nwhere)[(sum(nwhere == 0) + 1):length(nwhere)],
      revmonotone = rev(order(nwhere)[(sum(nwhere == 0) + 1):length(nwhere)]),
      seq_len(nvar)[nwhere > 0]
    )
  }
  
  # if (all(nwhere[visitSequence] == 0))
  #   stop("No locations to impute")
  flags <- nwhere == 0 & is.element(seq_len(nvar), visitSequence)
  if (any(flags))
    visitSequence <- visitSequence[!flags]
  visitSequence <- visitSequence[visitSequence <= nvar]
  visitSequence <- visitSequence[visitSequence >= 1]
  # if (length(visitSequence) == 0)
  #   stop("No locations to impute")
  setup$visitSequence <- visitSequence
  return(setup)
}

check.method <- function(setup, data) {
  # check method, set defaults if appropriate
  method <- setup$method
  defaultMethod <- setup$defaultMethod
  visitSequence <- setup$visitSequence
  nwhere <- setup$nwhere
  nvar <- setup$nvar
  
  # handle default
  if (all(method == "")) { 
    for (j in visitSequence) {
      y <- data[, j]
      if (is.numeric(y)) {
        method[j] <- defaultMethod[1]
      } else if (nlevels(y) == 2) {
        method[j] <- defaultMethod[2]
      } else if (is.ordered(y) && nlevels(y) > 2) {
        method[j] <- defaultMethod[4]
      } else if (nlevels(y) > 2) {
        method[j] <- defaultMethod[3]
      } else if (is.logical(y)) {
        method[j] <- defaultMethod[2]
      } else {
        method[j] <- defaultMethod[1]
      }
    }
  }
  
  # expand user's imputation method to all visited columns
  # single string supplied by user (implicit assumption of two columns)
  if (length(method) == 1) {
    if (is.passive(method))
      stop("Cannot have a passive imputation method for every column.")
    method <- rep(method, nvar)
  }
  
  # if user specifies multiple methods, check the length of the argument
  if (length(method) != nvar) {
    stop(paste0("The length of method (", length(method),
               ") does not match the number of columns in the data (", nvar,
               ")."))
  }
  
  # check whether the elementary imputation methods are on the search path
  active.check <- !is.passive(method) & nwhere > 0 & method != ""
  passive.check <- is.passive(method) & nwhere > 0 & method != ""
  check <- all(active.check) & any(passive.check)
  if (check) {
    fullNames <- rep.int("mice.impute.passive", length(method[passive.check]))
  } else {
    fullNames <- paste("mice.impute", method[active.check], sep = ".")
    if (length(method[active.check]) == 0) fullNames <- character(0)
  }
  notFound <- !vapply(fullNames, exists, logical(1), mode = "function", inherits = TRUE)
  if (any(notFound)) {
    stop(paste("The following functions were not found:",
               paste(fullNames[notFound], collapse = ", ")))
  }
  
  # type checks on built-in imputation methods
  for (j in visitSequence) {
    y <- data[, j]
    vname <- dimnames(data)[[2]][j]
    mj <- method[j]
    mlist <- list(m1 = c("logreg", "logreg.boot", "polyreg", "lda", "polr","logregSens"), 
                  m2 = c("norm", "norm.nob", "norm.predict", "norm.boot",
                         "mean", "2l.norm", "2l.pan","2l.lmm",
                         "2lonly.pan", "quadratic", "ri","normSens"), 
                  m3 = c("norm", "norm.nob", "norm.predict", "norm.boot",
                         "mean", "2l.norm", "2l.pan", "2l.lmm",
                         "2lonly.pan", "quadratic", "logreg", "logreg.boot",
                         "logregSens","normSens"))
    if (is.numeric(y) & (mj %in% mlist$m1)) {

      warning("Type mismatch for variable ", vname, "\nImputation method ",
              mj, " is for categorical data.",
              "\nIf you want that, turn variable ", vname, " into a factor,",
              "\nand store your data in a data frame.", call. = FALSE)
    } else if (is.factor(y) && nlevels(y) == 2 && (mj %in% mlist$m2)) {
      warning("Type mismatch for variable ", vname, "\nImputation method ",
              mj, " is not for factors.", call. = FALSE)
    } else if (is.factor(y) && nlevels(y) > 2 && (mj %in% mlist$m3)) {
      warning("Type mismatch for variable ", vname, "\nImputation method ",
              mj, " is not for factors with three or more levels.",
              call. = FALSE)
    }
  }
  setup$method <- method
  return(setup)
}


check.predictorMatrix <- function(setup) {
  ## checks and makes consistency edits of the predictormatrix
  pred <- setup$predictorMatrix
  varnames <- setup$varnames
  nwhere <- setup$nwhere
  nvar <- setup$nvar
  vis <- setup$visitSequence
  post <- setup$post
  
  if (!is.matrix(pred))
    stop("Argument 'predictorMatrix' not a matrix.")
  if (nvar != nrow(pred) || nvar != ncol(pred))
    stop(paste("The predictorMatrix has", nrow(pred), "rows and",
               ncol(pred), "columns. Both should be", nvar, "."))
  dimnames(pred) <- list(varnames, varnames)
  diag(pred) <- 0

  # inactivate predictors of complete variables
  for (j in seq_len(nvar)) {
    if (nwhere[j] == 0 && any(pred[j, ] != 0))
      pred[j, ] <- 0
  }
  
  setup$predictorMatrix <- pred
  setup$visitSequence <- vis
  setup$post <- post
  return(setup)
}


check.data <- function(setup, data, allow.na = FALSE, 
                       remove_collinear = TRUE, ...) {
  
  pred <- setup$predictorMatrix
  nvar <- setup$nvar
  varnames <- setup$varnames
  meth <- setup$method
  vis <- setup$visitSequence
  post <- setup$post
  
  # stop if the class variable is a factor
  isclassvar <- apply(pred == -2, 2, any)
  for (j in seq_len(nvar)) {
    if (isclassvar[j] && is.factor(data[,j])) 
      stop("Class variable (column ", j,
           ") cannot be factor. Convert to numeric by as.integer()")        
  }
  
  # remove constant variables but leave passive variables untouched
  for (j in seq_len(nvar)) {
    if (!is.passive(meth[j])) {
      d.j <- data[, j]
      v <- if (is.character(d.j)) NA else var(as.numeric(d.j), na.rm = TRUE)
      constant <- if (allow.na) {
                    if (is.na(v)) FALSE else v < 1000 * .Machine$double.eps
                  } else {
                    is.na(v) || v < 1000 * .Machine$double.eps
                  }
      didlog <- FALSE
      if (constant && any(pred[, j] != 0)) {
        out <- varnames[j]
        pred[, j] <- 0
        updateLog(out = out, meth = "constant")
        didlog <- TRUE
      }
      if (constant && meth[j] != "") {
        out <- varnames[j]
        pred[j, ] <- 0
        if (!didlog)
          updateLog(out = out, meth = "constant")
        meth[j] <- ""
        vis <- vis[vis != j]
        post[j] <- ""
      }
    }
  }
  
  ## remove collinear variables
  ispredictor <- apply(pred != 0, 2, any)
  if (any(ispredictor)) {
    droplist <- find.collinear(data[, ispredictor, drop = FALSE], ...)
  } else {
    droplist <- NULL
  }
  if (length(droplist) > 0 && remove_collinear) {
    for (k in seq_along(droplist)) {
      j <- which(varnames %in% droplist[k])
      didlog <- FALSE
      if (any(pred[, j] != 0)) {
        # remove as predictor
        out <- varnames[j]
        pred[, j] <- 0
        updateLog(out = out, meth = "collinear")
        didlog <- TRUE
      }
      if (meth[j] != "") {
        out <- varnames[j]
        pred[j, ] <- 0
        if (!didlog)
          updateLog(out = out, meth = "collinear")
        meth[j] <- ""
        vis <- vis[vis != j]
        post[j] <- ""
      }
    }
  }
  
  setup$predictorMatrix <- pred
  setup$visitSequence <- vis
  setup$post <- post
  setup$meth <- meth
  return(setup)
}




check.predictorSens <- function(setup) {
  pred <- setup$predictorSens
  varnames <- setup$varnames
  varnames2 <- paste(":",setup$varnames,sep="")
  
  nmis <- setup$nmis
  nvar <- setup$nvar
  vis <- setup$visitSequence
  method <- setup$method
  post <- setup$post
  if (!is.matrix(pred)) 
    stop("Argument 'predictorSens' not a matrix.")
  if (nvar != nrow(pred) | nvar != ncol(pred)) 
    stop(paste("The predictorSens has", nrow(pred), 
               "rows and", ncol(pred), "columns. Both should be", 
               nvar, "."))
  dimnames(pred) <- list(varnames, varnames2)
  diag(pred) <- 0
  isclassvar <- apply(pred == -2, 2, any)
  for (j in 1:nvar) {
    if (method[j] == "" & any(pred[, j] != 0) & nmis[j] > 
        0) {
      out <- varnames[j]
      updateLog(out = out)
      pred[, j] <- 0
      vis <- vis[vis != j]
      post[j] <- ""
      if (isclassvar[j]) 
        stop("Removed an incomplete class variable.")
    }
    if (nmis[j] == 0 & any(pred[j, ] != 0)) 
      pred[j, ] <- 0
  }
  setup$predictorSens <- pred
  setup$visitSequence <- vis
  setup$post <- post
  return(setup)
}


check.parmSens <- function(setup,data) {
  varnames <- setup$varnames
  nvar <- setup$nvar
  method <- setup$method
  predSens <- setup$predictorSens
  parmSens <- setup$parmSens
  formSens <- setup$formSens
  
  if (!is.list(parmSens)) 
    stop("Argument 'parmSens' not a list.")
  
  if (nvar != length(parmSens)) 
    stop(paste("The parmSens argument needs to be of length",nvar, ", i.e. with one component
               per variable in the order in they appear in the dataset."))
  
  for (j in 1:nvar) {
    
    if(any(unlist(parmSens[[j]])!="" ) && method[j]!="normSens" && method[j]!="logregSens")
      stop(paste("Entries in the parmSens argument other than the empty string '' can only be used 
                 with NARFCS methods, i.e. normSens or logregSens."))
    
    if(!any(unlist(parmSens[[j]])!="" ) && (method[j]=="normSens" || method[j]=="logregSens"))
      stop(paste("To use methods normSens or logregSens you need to provide sensitivity parameter values 
                 through the parmSens argument (see help). If no predictor variables
                 are included in the unidentifiable part of the model, you need to provide 
                 the sensitivity parameter value for the intercept term."))  
    
    if (any(unlist(parmSens[[j]])!="" ))  
    {
      if(! is.null(formSens) && nchar(formSens[j])>0)
      {
        tf<-terms(as.formula(formSens[j]),keep.order=T)
        npreds<-length(attr(tf,"term.labels"))+attr(tf,"intercept")
      }else
        npreds<-1+sum(predSens[j,]==1)
      
      if(length(parmSens[[j]])!= npreds)
        stop(paste("Component", j, "of parmSens needs to be a list of length equal,", npreds,
                   "i.e. 1 (for the intercept) + the number of variables specified as predictors 
                   in the unidentifiable part of the imputation model for variable", 
                   varnames[j],". In this list, the first element corresponds to the intercept shift 
                   and the following elements to the predictor variables, either in the order 
                   in which they appear in the dataset (if using the matrix syntax, 
                   i.e. predictorSens argument) or in the order in which they appear 
                   in the formula (if using formula syntax, i.e. formSens argument)."))
      
      if(length(parmSens[[j]][[1]])!=1)
        stop(paste("The sensitivity parameter corresponding to a shift in the
                   intercept must be a scalar. Thus the first entry number of component", j, 
                   "of parmSens needs to be of length 1."))
      
      if(npreds-1>0)
      {
        for(k in 1:(npreds-1))
        { 
          if(! is.null(formSens) && nchar(formSens[j])>0)
            fvar<-data[,attr(tf,"term.labels")[k]]
          else fvar<- data[,which(predSens[j,]!=0)[k]]
          
          if(is.factor(fvar))
          {
            npars<-length(levels(fvar))-1
            if(npars!=length(parmSens[[j]][[k+1]]))
              stop(paste("For factor predictors, a sensitivity parameter
                         needs to be supplied for each non-reference category.
                         Thus, entry number", k+1,"of component", j, "of parmSens needs to be 
                         of length", npars, "."))
          } else
            if(length(parmSens[[j]][[k+1]])!=1)
              stop(paste("For non-factor predictors only one sensitivity parameter
                         should be supplied. Thus, entry number", k+1,"of component", j, 
                         "of parmSens needs to be of length 1."))
          }      
        }
      
        }
    }
  
  return(setup)
}
