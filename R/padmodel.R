padModel <- function(data, method, predictorMatrix, visitSequence, 
                     form, post, nvar,predictorSens, formSens,
                     parmSens) {
  # Called by mice().  Augments the imputation model by including 
  # dummy variables. Adapts data, predictorMatrix, method and
  # visitSequence.  Returns a list whose components make up the padded model.
  
  categories <- data.frame(
    is.factor = factor(rep.int(FALSE, nvar), levels = c("TRUE", "FALSE")), 
    n.dummy   = rep.int(0, nvar), 
    is.dummy  = factor(rep.int(FALSE, nvar), levels = c("TRUE", "FALSE")), 
    father    = rep.int(0, nvar))
  
  
  categoriesSens<-categories 
  
  # make explicit local copy
  data <- data
  dataSens<-data
  
  for (j in seq_len(nvar)) {
    if (is.factor(data[, j]) && any(predictorMatrix[, j] != 0)) {
      
      categories[j, 1] <- TRUE
      data[, j] <- C(data[, j], contr.treatment)
      
      n.dummy <- length(levels(data[, j])) - 1
      categories[j, 2] <- n.dummy
      
      # add n.dummy rows and columns to curreny predictorMatrix
      predictorMatrix <- rbind(predictorMatrix,
                               matrix(0, 
                                      ncol = ncol(predictorMatrix), 
                                      nrow = n.dummy))
      predictorMatrix <- cbind(predictorMatrix, 
                               matrix(rep(predictorMatrix[, j], times = n.dummy), 
                                      ncol = n.dummy))
      
      # remove j as predictor
      predictorMatrix[seq_len(nvar), j] <- rep.int(0, times = nvar)
      form <- c(form, rep.int("", n.dummy))
      
      if (any(visitSequence == j)) {
        # make column j a predictor for its dummies
        idx <- (ncol(predictorMatrix) - n.dummy + 1):ncol(predictorMatrix)
        predictorMatrix[idx, j] <- rep.int(1, times = n.dummy)
        
        # extend visitSequence
        newcol <- ncol(predictorMatrix) - n.dummy + 1
        nloops <- sum(visitSequence == j)
        for (ii in seq_len(nloops)) {
          idx2 <- seq_along(visitSequence)[visitSequence == j][ii]
          visitSequence <- append(visitSequence, newcol, idx2)
        }
      }
      
      # add n.dummy columns to data
      data <- cbind(data, matrix(0, ncol = n.dummy, nrow = nrow(data)))
      idx <- (ncol(predictorMatrix) - n.dummy + 1):ncol(predictorMatrix)

      # set missing entries in new columns
      data[is.na(data[, j]), idx] <- NA
      
      # and merge model expansion by model.matrix() for observed values
      cat.column <- data[!is.na(data[, j]), j]
      data[!is.na(data[, j]), idx] <- model.matrix(~ cat.column - 1)[, -1]
      
      # name new columns
      names(data)[idx] <- paste(attr(data, "names")[j], seq_len(n.dummy), sep = ".")
      
      # extend method and post arrays
      method <- c(method, rep.int("dummy", n.dummy))
      post <- c(post, rep.int("", n.dummy))
      
      # append administrative info
      categories <- rbind(
        categories, 
        data.frame(is.factor = rep.int(FALSE, n.dummy), 
                   n.dummy   = rep.int(0, n.dummy), 
                   is.dummy  = rep.int(TRUE, n.dummy), 
                   father    = rep.int(j, n.dummy)))
  
    # perform the same manipulations for predictorSens, formSens, dataSens and categoriesSens
    # this is done for all the same variables as if a variable is included in the 
    # unidentifiable part it must have been included in the identifiable part
    # N.B. no further modifications are done to method, visitSequence or post
    categoriesSens[j, 1] <- TRUE
    dataSens[, j] <- C(dataSens[, j], contr.treatment)
    n.dummy <- length(levels(dataSens[, j])) - 1
    categoriesSens[j, 2] <- n.dummy
    predictorSens <- rbind(predictorSens, matrix(0, ncol = ncol(predictorSens), nrow = n.dummy))
    predictorSens <- cbind(predictorSens, matrix(rep(predictorSens[, j], times = n.dummy), ncol = n.dummy))
    predictorSens[1:nvar, j] <- rep(0, times = nvar)
    formSens <- c(formSens, rep("", n.dummy))
    idxSens<-(ncol(predictorSens) - n.dummy + 1):ncol(predictorSens)
    if (any(visitSequence == j)) {
      predictorSens[idxSens, j] <- rep(1, times = n.dummy)

    }
    dataSens <- (cbind(dataSens, matrix(0, ncol = n.dummy, nrow = nrow(dataSens))))
    dataSens[is.na(dataSens[, j]), idxSens] <- NA
    cat.column <- dataSens[!is.na(dataSens[, j]), j]
    dataSens[!is.na(dataSens[, j]), idxSens] <- model.matrix(~cat.column - 1)[, -1]
    names(dataSens)[idxSens] <- paste(attr(dataSens,"names")[j], (1:n.dummy), sep = ".")
    categoriesSens <- rbind(
      categoriesSens, 
      data.frame(is.factor = rep(FALSE, n.dummy), 
                 n.dummy = rep(0, n.dummy), 
                 is.dummy = rep(TRUE, n.dummy), 
                 father = rep(j, n.dummy)))                                                                                                       
                                                                                                                           
    }}
  
  varnames <- dimnames(data)[[2]]
  dimnames(predictorMatrix) <- list(varnames, varnames)
  varnamesSens <- dimnames(dataSens)[[2]]
  dimnames(predictorSens) <- list(varnamesSens, varnamesSens)
  names(method) <- varnames
  names(form) <- varnames
  names(formSens) <- varnamesSens
  names(post) <- varnames
  names(visitSequence) <- varnames[visitSequence]
  dimnames(categories)[[1]] <- dimnames(data)[[2]]
  dimnames(categoriesSens)[[1]] <- dimnames(dataSens)[[2]]

  if (anyDuplicated(names(data)))
    stop("Column names of padded data not unique")
  
  return(list(data = as.data.frame(data), 
              predictorMatrix = predictorMatrix, 
              method = method, 
              visitSequence = visitSequence,
              form = form, post = post, categories = categories,
              dataSens = as.data.frame(dataSens),
              categoriesSens = categoriesSens,
              predictorSens = predictorSens,formSens = formSens,parmSens=parmSens))
  
}
