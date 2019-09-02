#'Imputation by a linear mixed model
#'
#'Imputes missing data from a two-level normal model using \code{lme4::lmer()}.
#'
#'Draws are performed as follows. Using standard asymptotic results, the 
#'conditional (posterior) distribution of the parameters is approximated by a normal
#'distribution with mean and variance equal to the estimates of fixed and random effects
#'and of their variance--covariance, respectively, given restricted maximum likelihood 
#'(REML) estimates of the variance components.
#'@name mice.impute.2l.lmm
#'@inheritParams mice.impute.pmm
#'@param type Vector of length \code{ncol(x)} identifying random and class
#'variables.  Random variables are identified by a '2'. The class variable
#'(only one is allowed) is coded as '-2'. Fixed effects are indicated by 
#'a '1'.
#'@return Vector with imputed data, same type as \code{y}, and of length 
#'\code{sum(wy)}
#'@author Margarita Moreno-Betancur (2013)
#'@references
#'Moreno-Betancur M., Chavance M. (2016) Sensitivity analysis of incomplete 
#'longitudinal data departing from the missing at random assumption: Methodology 
#'and application in a clinical trial with drop-outs.
#'\emph{Statistical methods in medical research}, 25 (4), 1471-1489
#'[Epub ahead of print, May 22 2013]
#'@family univariate \code{2l} functions
#'@keywords datagen
#'@export

mice.impute.2l.lmm <- function(y, ry, x, type, wy = NULL,...) {
    
    if (!requireNamespace("lme4", quietly = TRUE))
      stop("Please install package 'lme4'", call. = FALSE)
    
    if (is.null(wy)) wy <- !ry
    
    # define cluster, random and fixed effects
    nameG <- names(x)[type == -2]
    namesZ <- names(x)[type == 2]
    namesX <- names(x)[type > 0]
    nam<-"y"
    data<-data.frame(cbind(y=y,x))
    
    ### CREATE MODEL FORMULA ###
    if(length(nameG) == 0)
      stop(paste("No grouping factor to include random effects in model."))
    
    if (length(namesZ) == 0) {
      formulaGLMER = as.formula(paste(nam,"~1+",paste(namesX,collapse="+"),
                                      paste("+(1|",nameG,")")))
    } else if (length(namesZ) != 0)
      formulaGLMER = as.formula(paste(nam,"~1+",paste(namesX,collapse="+"),
                                      paste("+(1+",paste(namesZ,collapse="+"),"|",nameG,")")))
    
    ### FIT LMM FROM AVAILABLE CASES ###
    fit <- lme4::lmer(formulaGLMER, data[ry==1,])
    ids <- sort(unique(unlist(data[ry==1, nameG], FALSE, FALSE)))
    
    ### RETRIEVE FIXED EFFECTS ESTIMATES ###
    Beta    <- lme4::fixef(fit)
    VarBeta <- t(Matrix::chol(as.matrix(vcov(fit))))
  
    ### RETRIEVE RANDOM EFFECTS PREDICTIONS AND CONDITIONAL VAR-COV MATRIX ###
    Bmat    <- lme4::ranef(fit, condVar = TRUE)[[1]]  # nLevels x nRanef matrix
    nLevels <- nrow(Bmat)
    nRanef  <- ncol(Bmat)
    B       <- unlist(Bmat, FALSE, FALSE)
    VarB    <- attr(Bmat, "postVar", exact = TRUE)
    
    ### DRAW FIXED AND RANDOM EFFECTS FROM THEIR ASYMPTOTIC DISTRIBUTION ###
    beta.star <- Beta + VarBeta %*% rnorm(ncol(VarBeta))
    B.disturbance <- as.vector(t(
      do.call(cbind, lapply(seq(dim(VarB)[3]),
                            function(x) {
                              L <- t(Matrix::chol(VarB[ , , x]))
                              B.disturbance <- L %*% rnorm(ncol(L))
                            }))))
    b.star    <- B + B.disturbance
    
    ### GET TRANSPOSE AND DIAGONAL OPERATORS FOR SPARSE MATRICES ###
    ### FROM MATRIX PACKAGE ###
    tr<-get("t",env=environment(rep2abI),mode="function")
    diago<-get("diag",env=environment(rep2abI),mode="function")
    
    ### CALCULATE DEGREES OF FREEDOM OF ESTIMATOR OF RESIDUAL VARIANCE ###
    A<-getME(fit,"A")
    X1<-getME(fit,"X")
    M<-cBind(tr(A),X1)
    I <- matrix(0,nrow=nrow(as(A,"matrix")),ncol=nrow(as(A,"matrix")))
    I[row(I)==col(I)] <- 1
    Zero<-matrix(0,nrow=nrow(I),ncol=ncol(X1))
    N<-cBind(rBind(tr(A), I),rBind(X1,Zero))
    ddl<-nrow(X1)-sum(diago(M%*%solve(tr(N)%*%N)%*%tr(M)))
    # = number of observations used to fit the model - trace of the hat matrix
    
    ### DRAW THE RESIDUAL VARIANCE AND A GAUSSIAN ERROR ##
    sigma<-sqrt((attr(VarCorr(fit),"sc",exact=TRUE)^2)*ddl/rchisq(1,ddl))
    e<-rnorm(sum(!ry),0,sigma)
    
    ### OBTAIN MODEL MATRICES FOR IMPUTATION ###
    data2 <- data
    data2[!ry, nam] <- 0
    ids2  <- sort(unique(unlist(data2[!ry, nameG], FALSE, FALSE)))
    if (any(!ids2%in%ids)) #Note: different from survtd implementation because there everyone has missing data
      stop(paste0("Trying to impute for individuals who did not have ",
                  "random effect estimates from available case data."))
    formulaGLMER2 = as.formula(paste(nam,"~1+",paste(namesX,collapse="+")))
    X <- model.matrix(formulaGLMER2, data2)
    
    formulaGLMER3 = as.formula(paste(nam,"~-1+as.factor(",nameG,")+as.factor(",nameG,"):",paste(namesZ,collapse="+")))
    Z <- Matrix::sparse.model.matrix(formulaGLMER3, data2)
    
    ### VALUES TO IMPUTE MISSING OUTCOMES = LINEAR PREDICTOR  ###
    fixedpart <- X[!ry, ] %*% beta.star
    vec <- fixedpart + Z[!ry, ] %*% b.star +e
    return(vec@x)
}
  