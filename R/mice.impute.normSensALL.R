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
#'@aliases mice.impute.normSensALL normSensALL
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


mice.impute.normSensALL<-
  function (y, ry, x, wy = NULL, xSens, parmS, ...) 
  { 
    if (is.null(wy)) wy <- !ry
    x <- cbind(1, as.matrix(x))
    xSens <- cbind(1, as.matrix(xSens))
    parm <- .norm.draw2(y, ry, x, xSens, parmS,...)
    
    return(x[wy, ] %*% parm$beta + rnorm(sum(wy)) * parm$sigma+ as.matrix(xSens[wy, ])%*%parmS)
    
  }

#' Draws values of beta and sigma by Bayesian linear regression
#' 
#' This function draws random values of beta and sigma under the Bayesian 
#' linear regression model as described in Rubin (1987, p. 167). This function
#' can be called by user-specified imputation functions.
#' 
#'@aliases norm.draw2 .norm.draw2
#'@param y Incomplete data vector of length \code{n}
#'@param ry Vector of missing data pattern (\code{FALSE}=missing,
#'\code{TRUE}=observed)
#'@param x Matrix (\code{n} x \code{p}) of complete covariates.
#'@param ridge A small numerical value specifying the size of the ridge used. 
#' The default value \code{ridge = 1e-05} represents a compromise between stability
#' and unbiasedness. Decrease \code{ridge} if the data contain many junk variables.
#' Increase \code{ridge} for highly collinear data. 
#'@param ... Other named arguments.
#'@return A \code{list} containing components \code{coef} (least squares estimate),
#'\code{beta} (drawn regression weights) and \code{sigma} (drawn value of the 
#'residual standard deviation).
#'@references
#'Rubin, D.B. (1987). \emph{Multiple imputation for nonresponse in surveys}. New York: Wiley.
#'@author Margarita Moreno-Betancur
#'@export
norm.draw2 <- function(y, ry, x,xSens, parmS,  ridge = 1e-05, ...) 
  return(.norm.draw2(y, ry, x,xSens, parmS,  ridge = 1e-05, ...))

###'@rdname norm.draw2  
###'@export
.norm.draw2 <- function(y, ry, x, xSens, parmS, ridge = 1e-05, ...) {
  xobs <- x
  yobs <- y-(rowSums(as.matrix(xSens)%*%parmS)*(!ry))
  xtx <- crossprod(xobs)  # SvB 21/01/2014
  pen <- ridge * diag(xtx)
  if (length(pen) == 1)
    pen <- matrix(pen)
  v <- solve(xtx + diag(pen))
  coef <- t(yobs %*% xobs %*% v)
  residuals <- yobs - xobs %*% coef
  df <- max(length(ry) - ncol(x), 1)  # SvB 31/10/2012
  sigma.star <- sqrt(sum((residuals)^2)/rchisq(1, df))  # SvB 01/02/2011
  beta.star <- coef + (t(chol(sym(v))) %*% rnorm(ncol(x))) * sigma.star
  parm <- list(coef, beta.star, sigma.star)  # SvB 10/2/2010
  names(parm) <- c("coef", "beta", "sigma")  # SvB 10/2/2010
  return(parm)
}
