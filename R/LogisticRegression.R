#' Sigmoid function
#' @description Sigmoid function, otherwise known as a logistic function
#' @param z input
#'
#' @details
#' The sigmoid (or logistic) function,
#' \deqn{
#' \sigma(z) = 1 / (1 + e^{-z})
#' }
#' is used in logistic regression to model probabilities, commonly the probability of predicting the positive class,
#' i.e. \eqn{p(y=1)}.
#'
#' @return output value of the function
#' @export
#' @examples
#' theta = c(1, 0.1)
#' X = matrix(rnorm(4), 2, 2)
#' sigmoid(X %*% theta)
sigmoid = function(z) { 1/( 1 + exp(-z) ) }

#' Bernoulli distribution log-likelihood
#'
#' @description Log-likelihood of the Bernoulli distribution, commonly used in optimisation procedures for maximisation.
#'
#' @param theta parameters relating to \code{x}
#' @param x explanatory variables in model matrix
#' @param y binary response variable
#'
#' @details The log-likelihood is written
#' \deqn{
#' \sum_{i=1}^n y \log(\sigma(X\beta)) + (1-y) \log(1-\sigma(X\beta))
#' }
#' @return single value, the log-likelihood of the Bernoulli distribution with fixed \code{x} and \code{y}
#' @export
#'
#' @examples
#' theta = c(1, 0.1)
#' X = matrix(rnorm(4), 50, 2)
#' y = sample(0:1, 50, replace=TRUE)
#' loglik_lr(theta, X, y)
loglik_lr = function(theta, x, y) {sum(
  y*as.vector(log(sigmoid(x %*% theta)))+log(1-as.vector(sigmoid(x %*% theta)))*(1-y)
)}

#' Logistic Regression
#'
#' @description \code{lr} is used to fit a logistic regression model for a binary response variable.
#'
#' @param formula an object of class \code{\link[stats]{formula}}, a symbolic description of the model to be fitted.
#' The specified names need to also be in \code{data}.
#' @param data a required data frame containing the variables in the model
#' @param init optional initial conditions to be passed to optimisation of the log-likelihood
#'
#' @details The form of the \code{formula} argument will be of the form \code{response ~ predictor1 + predictor2 + ...},
#' with \code{predictor1} and \code{predictor2} being named columns of the data frame in \code{data}.
#'
#' The log-likelihood (from \code{loglik_lr}) is maximised using \code{\link{irls.lr}} with initial estimates given by \code{init}.
#' If no initial values are supplied, this uses a vector of zeros instead.
#'
#' @return An S3 object of class 'lr', which is a list containing
#' \item{\code{coefficients}}{a vector of coefficients corresponding to covariates specified in \code{formula}}
#' \item{\code{data}}{the \code{data} input to the function}
#' \item{\code{formula}}{the \code{formula} input to the function}
#' \item{\code{X}}{the model matrix \code{X}}
#' \item{\code{val}}{value of the final log-likelihood at the values of \code{coefficients}, given by \code{loglik_lr}}
#' \item{\code{its}}{number of iterations performed to retrieve the maximised log-likelihood}
#' @export
#' @importFrom stats model.matrix
#' @importFrom Matrix Matrix
#' @examples
#' y = sample(0:1, 50, replace=TRUE)
#' d = data.frame(y = y, x = rnorm(10*y + 15))
#' fit = lr(y ~ x, data = d)
lr <- function(formula, data, init=NULL){
  y = data[,all.vars(formula)[1]]
  X = model.matrix(formula, data)
  if((sum(X==0)/length(X))>=0.5) X = Matrix::Matrix(X, sparse=TRUE)
  est = irls.lr(y, X, init=init)
  out = list("coefficients" = est$par, "data" = data, "formula" = formula, "X" = X,
             "val" = est$val, "iters" = est$iters)
  class(out) = "lr"
  return(out)
}

#' Iteratively Re-weighted Least Squares
#'
#' @description Optimisation procedure based on iteratively re-weighted least squares, called by \code{\link{lr}}.
#'
#' @param y response vector
#' @param x model matrix of covariates
#' @param init initial estimate of \code{theta}
#' @param tol tolerance parameter, default \code{1e-6}
#' @param maxiter optional number of iterations
#'
#' @details
#' Iterative method used to find parameter estimates of \eqn{\beta} from the least squares problem
#' \deqn{
#' \beta = argmin (z - X\beta)^T W (z- X\beta),
#' }
#' where \eqn{W} is a diagonal matrix of weights with \eqn{i}-th diagonal element being
#' \deqn{
#' \sigma(x_i ; \beta) (1 - \sigma(x_i ; \beta)),
#' }
#' and \eqn{z} is the vector
#' \deqn{
#' z = X \beta + W^{-1} (y - \sigma(x_i ; \beta)).
#' }
#' The parameter vector \eqn{\beta} is updated iteratively with a Newton update of the form
#' \deqn{
#' \beta = (X^T W X)^{-1} X^T W z.
#' }
#'
#' @return a list containing
#' \item{\code{par}}{a vector of estimates of \code{theta}, the parameters being optimised}
#' \item{\code{val}}{the value of the log-likelihood at the final \code{theta} estimate}
#' \item{\code{iters}}{the number of iterations needed to converge}
#'
#' @importFrom Matrix solve
#' @importFrom Matrix t
#' @export
irls.lr = function(y, x, init=NULL, tol=1e-6, maxiter=100){
  theta = if(is.null(init)) rep(0, ncol(x)) else init
  W = Matrix::Diagonal(nrow(x))
  diag(W) = sigmoid(x %*% theta) * (1-sigmoid(x %*% theta))
  z = x %*% theta + Matrix::solve(W) %*% as.matrix(y - as.vector(sigmoid(x %*% theta)))
  for(i in 1:maxiter){
    thetanew = Matrix::solve(t(x)%*%(W%*%x)) %*% (t(x)%*%(W%*%z))
    if( all(abs(theta-thetanew) < tol) ) break
    theta = thetanew
    diag(W) = sigmoid(x %*% theta) *(1-sigmoid(x %*% theta))
    z = x %*% theta + Matrix::solve(W) %*% as.matrix(y-as.vector(sigmoid(x %*% theta)))
  }
  return(list("par"=as.matrix(theta), "val"=loglik_lr(theta,x,y), "iters"=i))
}

#' Print Logistic Regression Model
#'
#' @description Displays important output from a \code{\link{lr}} object; the parameter estimates and the maximised log-likelihood value.
#'
#' @param x an object of class "\code{lr}", the output from \code{\link{lr}}
#' @param ... further arguments
#'
#' @return a printed, named vector of coefficients and log-likelihood value at these estimates
#' @export
#' @examples
#' y = sample(0:1, 50, replace=TRUE)
#' d = data.frame(y = y, x = rnorm(10*y + 15))
#' fit = lr(y ~ x, data = d)
#' print(fit)
print.lr <- function(x,...){
  cat(" \n")
  cat("Logistic Regression model \n")
  cat("_________________________ \n")
  cat("Coefficients: \n")
  coef = t(as.matrix(x$coefficients))
  colnames(coef) = c("Intercept", colnames(x$X)[2:ncol(x$X)])
  print(coef)
  cat("Log-likelihood value: ", x$val, "\n")
}



#' Predict from a Logistic Regression model
#'
#' @description Predicted values of a certain type based on a \code{\link{lr}} object.
#'
#' @param object an object of class "\code{lr}", the output from \code{\link{lr}}
#' @param newdata an optional data frame to predict from. If omitted, predictions will be from the original data frame used to fit the model.
#' @param thresh an optional threshold parameter, see 'Details'.
#' @param type type of predictions, can be one of "preds", "probs" or "vals", see 'Details'.
#' @param ... further arguments
#'
#' @details
#' Predictions from a logistic regression model are formed by first multiplying the model matrix by the parameter vector, i.e.
#' \deqn{
#' f(x; \beta) := X \beta.
#' }
#' These values are outputted when the \code{type} argument is "vals". Predictions for the positive class can be obtained when these values
#' are above a certain threshold value, the argument \code{thresh}, i.e.
#' \deqn{
#' f(x; \beta) \ge t,
#' }
#' where the threshold \eqn{t} is usually equal to zero. This is how the predictions are made when the \code{type} argument is "preds".
#' To get probabilities, the \code{\link{sigmoid}} function is used on the values, i.e.
#' \deqn{
#' p(y=1) = 1/(1 + e^{-X \beta}).
#' }
#' @return a vector of predictions, probabilities or values, depending on the input to \code{type}.
#' @importFrom stats model.matrix
#' @export
#' @examples
#' y = sample(0:1, 50, replace=TRUE)
#' d = data.frame(y = y, x = rnorm(10*y + 15))
#' fit = lr(y ~ x, data = d)
#' predict(fit, type="preds")
predict.lr <- function(object, newdata = NULL, thresh = 0, type = "preds", ...){
  if(is.null(newdata)) {
    X = model.matrix(object$formula, object$data, ...)
  } else{
    X = model.matrix(object$formula, newdata, ...)
  }
  ind = grep(type, c("preds","probs","vals"), ignore.case=TRUE)
  if(length(ind)>1 | length(ind)==0) stop("Supply either 'preds', 'probs' or 'vals' to argument 'type'.")

  p = X %*% object$coefficients
  if(ind==3) return(p)
  if(ind==2) return(sigmoid(p))
  if(ind==1){
    pred = numeric(length(p))
    pred[p > thresh] = 1
    return(pred)
  }
}

#' AUC Calculator
#'
#' @param p a vector of predictions
#' @param y a vector of the true observed response variable
#'
#' @description
#' See \code{\link{roc.lr}} for details. This function is used to only output the AUC value for use in other functions, such as \code{\link{cv.lr}}.
#' @return the AUC value
#' @export
auc_in <- function(p, y){
  thresh = seq(min(p), max(p), length=100)
  TPR = FPR = rep(NA, length(thresh))
  for(i in 1:length(thresh)){
    f_pred = numeric(length(p))
    f_pred[p >= thresh[i]] = 1

    TPR[i] = sum(f_pred[y==1])/length(y[y==1])
    FPR[i] = sum(f_pred[y!=1])/length(y[y==0])
  }
  height = (TPR[-1]+TPR[-length(TPR)])/2
  width = -diff(FPR)
  return(sum(height*width))
}

#' Match Error Function
#'
#' @description Used to get error function in \code{\link{cv.lr}}.
#' @param type type of error function to be output, either "mse", "auc" or "log"
#'
#' @return function of either mean squared error or \code{\link{auc_in}}
#' @export
get_error <- function(type="mse"){
  if(type=="mse") return( function(f, y) mean((f-y)^2, na.rm=TRUE))
  if(type=="auc") return( auc_in )
  if(type=="log") return( log_score )
}

#' Cross-Validation of Logistic Regression Model
#'
#' @description Implementation of cross-validation for a \code{\link{lr}} object, calculation of error across
#' a number of subsets of the inputted data set.
#'
#' @param lrfit an object of class "\code{lr}", the output to \code{\link{lr}}
#' @param metric which metric to calculate, one of "mse", "auc" or "both". See 'Details'.
#' @param leave_out number of points to leave out for cross-validation.
#' @param verbose logical; whether to print information about number of iterations completed.
#' @param seed optional; number to be passed to \code{\link[base]{set.seed}} before shuffling the data set
#'
#' @details
#' \eqn{k}-fold cross-validation, where \eqn{k} is the input to the \code{leave_out} argument.
#' This can be used to judge the out-of-sample predictive power of the model by subsetting the original
#' data set into two partitions; fitting the model for the (usually larger) one, and testing the
#' predictions of that model on the (usually smaller) partition. The position of the \eqn{k}
#' points separated from the data set are selected uniformly at random.
#'
#' The error metrics available are that of mean squared error, AUC, or log score; selected by the \code{metric}
#' argument being one of "mse", "auc", "log" or "all". See \code{\link{roc.lr}} for details on AUC.
#' If \code{metric} is "all", then a vector will be output containing all three metrics.
#'
#' Note that the output from \code{metric = "auc"} has non-deterministic elements due to the shuffling
#' of the data set. To mitigate this, include a number to the \code{seed} argument.
#' @export
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom stats predict
#' @return error value or vector consisting of the average of the chosen \code{metric}
#'
cv.lr <- function(lrfit, metric = "mse", leave_out = nrow(lrfit$data)/10, verbose=TRUE, seed=1){

  ind = grep(metric, c("mse", "auc", "log", "all"), ignore.case=TRUE)
  if(length(ind)>1 | length(ind)==0) stop("Supply either 'mse', 'auc', 'log', or 'all' to argument 'type'.")
  if(ind==1) method = get_error("mse")
  if(ind==2) method = get_error("auc")
  if(ind==3) method = get_error("log")
  if(ind==4) method = c(get_error("mse"), get_error("auc"), get_error("log"))

  set.seed(seed)
  data = lrfit$data[sample(1:nrow(lrfit$data), nrow(lrfit$data)),]

  iters = nrow(data) %/% leave_out
  final_iter = nrow(data) %% leave_out
  fin = if(round(final_iter)!=0) 1 else 0
  if(verbose) cat("   --- CROSS VALIDATION PROCESS ---  \n")
  if(ind!=4) error = rep(NA, iters+fin)
  if(ind==4) {
    error = array(NA, c(iters+fin, 3))
    colnames(error) = c("mse", "auc", "log")
  }
  start_points = c(0, (1:iters)*leave_out)
  if(verbose) progress = txtProgressBar(1, iters+fin, style=3)
  for(i in 1:iters){
    cross_out = data[(start_points[i]+1):(start_points[i+1]),]
    cross_in = data[-((start_points[i]+1):(start_points[i+1])),]
    cross_fit = lr(lrfit$formula, data = cross_in)

    y_out = cross_out[,all.vars(lrfit$formula)[1]]
    if(ind==1) p = predict(cross_fit, newdata=cross_out, type="preds")
    if(ind==2) p = predict(cross_fit, newdata=cross_out, type="vals")
    if(ind==3) p = predict(cross_fit, newdata=cross_out, type="probs")
    if(ind==4) {
      p1 = predict(cross_fit, newdata=cross_out, type="preds")
      p2 = predict(cross_fit, newdata=cross_out, type="vals")
      p3 = predict(cross_fit, newdata=cross_out, type="probs")
      error[i, 1] = method[[1]](p1, y_out)
      error[i, 2] = method[[2]](p2, y_out)
      error[i, 3] = method[[3]](p3, y_out)
    } else error[i] = method(p, y_out)
    if(verbose) setTxtProgressBar(progress, i)
  }
  if(fin) {
    cross_out = data[(start_points[iters+1]+1):(start_points[iters+1]+final_iter),]
    cross_in = data[-((start_points[iters+1]+1):(start_points[iters+1]+final_iter)),]
    cross_fit = lr(lrfit$formula, data = cross_in)

    y_out = cross_out[,all.vars(lrfit$formula)[1]]
    if(ind==1) p = predict(cross_fit, newdata=cross_out, type="preds")
    if(ind==2) p = predict(cross_fit, newdata=cross_out, type="vals")
    if(ind==3) p = predict(cross_fit, newdata=cross_out, type="probs")
    if(ind==4) {
      p1 = predict(cross_fit, newdata=cross_out, type="preds")
      p2 = predict(cross_fit, newdata=cross_out, type="vals")
      p3 = predict(cross_fit, newdata=cross_out, type="probs")
      error[iters+1, 1] = method[[1]](p1, y_out)
      error[iters+1, 2] = method[[2]](p2, y_out)
      error[iters+1, 3] = method[[3]](p3, y_out)
    } else error[iters+1] = method(p, y_out)
    if(verbose) setTxtProgressBar(progress, iters+1)
  }
  if(verbose) cat("\n")
  if(ind == 4) return(colMeans(error)) else return(mean(error))

}

#' Receiver Operating Characteristic (ROC) curve
#'
#' @description A method of judging the predictive performance of a model by plotting and/or averaging the
#' probability of predicting the positive class correctly, over multiple thresholds. See 'Details'.
#'
#' @param lrfit an object of class "\code{lr}", the output from \code{\link{lr}}
#' @param newdata an optional data frame to predict from. If ignored, the default data frame is that used to fit the original model.
#' @param plot logical; if \code{TRUE}, then the ROC curve is plotted
#' @param len optional; number of different thresholds to use.
#' @export
#' @details
#' A positive prediction from a logistic regression model is made when
#' \deqn{
#' f(x; \beta) := X \beta \ge t.
#' }
#'
#' where \eqn{t} is some threshold. See \code{\link{predict.lr}} for details. A different threshold
#' \eqn{t_0} will yield a different set of predictions. For a given sequence \eqn{t_j in [min(t), max(t)]},
#'  for \eqn{j=1,\dots,J}, the True Positive Rate (TPR) and False Positive Rate (FPR) can be calculated as
#' \deqn{
#'  TPR(j) = \sum I(f(x_i;\beta) \geq t_j)/\sum I(y_i = 1),
#' }
#' \deqn{
#'  FPR(j) = \frac{\sum I(f(x_i;\beta) < t_j)}{\sum I(y_i = 1)}.
#' }
#' The ROC curve is plotted from the pairs \eqn{(FPR(j), TPR(j))}, and the AUC is calculated as the area under this curve, i.e.
#' \deqn{
#'  AUC = \int_{j=1}^J TPR(FPR(j))dj.
#' }
#' @return the AUC value, and a plot of the ROC curve if \code{plot=TRUE}
#' @importFrom graphics abline grid points
#' @importFrom stats predict
roc.lr <- function(lrfit, newdata=NULL, plot=TRUE, len=50){
  f = predict(lrfit, newdata = newdata, type="vals")
  if(is.null(newdata)) insample = lrfit$data[,all.vars(lrfit$formula)[1]]
  if(!is.null(newdata)) insample = newdata[,all.vars(lrfit$formula)[1]]

  thresh = seq(min(f), max(f), length=len)
  TPR = FPR = rep(NA, length(thresh))
  for(i in 1:length(thresh)){
    f_pred = numeric(length(f))
    f_pred[f >= thresh[i]] = 1

    insample1s = insample[insample==1]
    insample0s = insample[insample==0]

    TPR[i] = sum(f_pred[insample==1])/length(insample[insample==1])
    FPR[i] = sum(f_pred[insample!=1])/length(insample[insample==0])
  }
  if(plot){
    plot(FPR, TPR, type="l", lwd=2)
    grid()
    abline(0, 1, col="red", lty=2)
    if(len<=50) points(FPR, TPR, pch=20)
  }
  height = (TPR[-1]+TPR[-length(TPR)])/2
  width = -diff(FPR)
  sum(height*width)
}


#' Log Score
#'
#' @description Log score diagnostic based on probabilities of predicting the correct class
#'
#' @param p a vector of probabilities
#' @param y a vector of the true observed response variable
#'
#' @details
#' The log score is defined as
#' \deqn{
#' LS = 1/n \sum - log p (z)
#' }
#' where \eqn{p(z)} is the probability of predicting the correct value \eqn{z}. This is averaged over all data points.
#'
#' The log score penalises probabilities that the model assigns to the correct class that are low, and rewards
#' those that are high in the correct place.
#'
#' @return a single value, the mean of the negative log of the probabilities for predicting the correct class
#' @export
log_score <- function(p, y){

  predicted_probs = cbind(1-p, p)
  correct_ind = y+1

  correct_probs = rep(NA, length(p))
  for(i in 1:length(p)) correct_probs[i] = predicted_probs[i,correct_ind[i]]

  mean(-log(correct_probs),na.rm=TRUE)
}




