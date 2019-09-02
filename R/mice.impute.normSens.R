# --------------------MICE.IMPUTE.NORMSENS----------------------------------

#'For NARFCS procedure only: Not-at-Random Imputation by Bayesian linear regression
#'
#'Imputes univariate missing data using Bayesian linear regression analysis in 
#'exactly the same way as in \code{mice.impute.norm} except the drawn imputations are 
#'then modified by adding a quantity as specified for the 
#'unidentified part of the model through \code{predictorSens} or \code{formSens}, and 
#'corresponding sensitivity parameter values through \code{parmSens}.
#'This results in "not-at-random" imputations.
#'
#'@aliases mice.impute.normSens normSens
#'@param y Vector to be imputed
#'@param ry Vector of missing data pattern (\code{FALSE}=missing,
#'\code{TRUE}=observed). This function fits the imputation model to all the data 
#'(observed and imputed).
#'(\code{TRUE}) and missing values (\code{FALSE}) in \code{y}.
#'@param x Numeric design matrix with \code{length(y)} rows with predictors for 
#'\code{y}. Matrix \code{x} may have no missing values.
#'@param wy Logical vector of length \code{length(y)}. A \code{TRUE} value 
#'indicates locations in \code{y} for which imputations are created.
#'@param xSens Matrix of covariates for unidentified part of the model.
#'@param parmS Vector of sensitivity parameters for the unidentified part of the model.
#'@param ... Other named arguments.
#'@return Vector with imputed data, same type as \code{y}, and of length 
#'\code{sum(wy)}
#'@author M. Moreno-Betancur, F. Leacy, D. Tompsett, I. White, 2017
#'@seealso \code{link{mice.impute.norm}}
#'@references 
#'
#'Moreno-Betancur M, Leacy FP, Tompsett D, White I. "mice: The NARFCS procedure for sensitivity analyses" 
#'(available at: \url{https://rawgit.com/moreno-betancur/NARFCS/master/README.html})
#'
#'Leacy FP. Multiple imputation under missing not at random assumptions via fully conditional 
#'specification 2016 (PhD thesis).
#'
#'Tompsett D, Leacy FP, Moreno-Betancur M, White I. On the use of the not at random fully conditional
#'specification procedure (NARFCS) in practice (submitted).
#'@family univariate imputation functions
#'@export


mice.impute.normSens<-
  function (y, ry, x, wy = NULL, xSens, parmS, ...) 
  { 
    if (is.null(wy)) wy <- !ry
    x <- cbind(1, as.matrix(x))
    xSens <- cbind(1, as.matrix(xSens))
    parm <- .norm.draw(y, ry, x, ...)
    
    return(x[wy, ] %*% parm$beta + rnorm(sum(wy)) * parm$sigma+ as.matrix(xSens[wy, ])%*%parmS)
    
  }

