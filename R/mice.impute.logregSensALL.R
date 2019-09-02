# -------------------------MICE.IMPUTE.LOGREGSENS-------------------------

#'For NARFCS procedure only: Not-at-Random Imputation by logistic regression.
#'
#'Imputes univariate missing data using using Bayesian logistic regression in 
#'exactly the same way as in \code{mice.impute.logreg} except the linear predictor
#'of the logistic model is modified before drawing the binary imputation by adding 
#'a quantity as specified for the unidentified part of the model through 
#'\code{predictorSens} or \code{formSens}, and corresponding sensitivity parameter 
#'values through \code{parmSens}. This results in "not-at-random" imputations.
#'
#'As for \code{logreg}, perfect prediction is handled by the data augmentation method.
#'
#'@aliases mice.impute.logregSens
#'@param y Vector to be imputed
#'@param ry Vector of missing data pattern (\code{FALSE}=missing,
#'\code{TRUE}=observed). This function fits the imputation model to all the data 
#'(observed and imputed).
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
#'@seealso \code{link{mice.impute.logreg}}
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
#'@keywords datagen
#'@export

mice.impute.logregSensALL<-
  function (y, ry, x, wy = NULL, xSens,parmS,...) 
  {
    if (is.null(wy)) wy <- !ry
    
    rySens <- ry
    aug <- augment(y, ry, x, wy)
    x <- aug$x
    y <- aug$y
    ry <- aug$ry
    wy <- aug$wy
    w <- aug$w
    
    x <- cbind(1, as.matrix(x))
    xSens <- cbind(1, as.matrix(xSens))
    expr <- expression(glm.fit(x, y, family = binomial(link = logit), 
                               weights = w,offset=rowSums(as.matrix(xSens)%*%parmS)*(!ry)))
fit <- suppressWarnings(eval(expr))
fit.sum <- summary.glm(fit)
beta <- coef(fit)
rv <- t(chol(sym(fit.sum$cov.unscaled)))
beta.star <- beta + rv %*% rnorm(ncol(rv))

p <- 1/(1 + exp(-(x[wy, ] %*% beta.star+as.matrix(xSens[!rySens,])%*%parmS)))

vec <- (runif(nrow(p)) <= p)
vec[vec] <- 1
if (is.factor(y)) {
  vec <- factor(vec, c(0, 1), levels(y))
}
return(vec)
  }

