#' Print GEV model
#'
#' @description Print elements of the GEV model.
#'
#' @param x object of class "\code{gev}" fit by \code{\link{gev_fit}}
#' @param ... further arguments
#'
#' @details This prints the maximum likelihood estimates of parameters for each parameter in the GEV distribution,
#' and the value of the log-likelihood here.
#'
#' @return printed elements
#' @export
print.gev <- function(x, ...){
  cat("-- GEV Model -- \n \n")

  mu = t(as.matrix(x$par$mu))
  colnames(mu) = colnames(x$modX$mu)
  sig = t(as.matrix(x$par$sig))
  colnames(sig) = colnames(x$modX$sig)
  xi = t(as.matrix(x$par$xi))
  colnames(xi) = colnames(x$modX$xi)

  cat("Mu (location) parameters: \n")
  print(mu)
  cat("\n Sigma (scale) parameters: \n")
  print(sig)
  cat("\n Xi (shape) parameters: \n \n")
  print(xi)
  cat("Log-likelihood:", x$loglik, "\n")
}

#' Plot GEV model diagnostics
#'
#' @description Plot the probability and the quantile-quantile plot for the GEV distribution, for a model
#' fit by \code{\link{gev_fit}}.
#'
#' @param x object of class "\code{gev}" fit by \code{\link{gev_fit}}
#' @param ... further arguments for plotting (will apply to both plots)
#'
#' @details The probability plot plots the model calculated
#' \deqn{
#' G(z_i) = \exp \{ - ( 1 + \xi (z_i - \mu)/\sigma)^{-1/\xi} \}
#' }
#' against the empirical
#' \deqn{
#' H(z_i) = i / (m + 1).
#' }
#' If the model is a good fit, these should be approximately equal, and so the points plotted should lie on a
#' diagonal line.
#'
#' The quantile-quantile plot plots
#' \deqn{
#' G^{-1} = \mu - (\sigma / \xi) (1 - (- log(i/(m+1)))^{-\xi} )
#' }
#' against the observed maxima \eqn{z_i}. Similar to the probability plot, this should also lie on the diagonal line
#' for a good fit, as these points should be approximately equal.
#' @return printed elements
#' @importFrom graphics par
#' @export
plot.gev = function(x,  ...){
  p = gev_parsum(x, x$data)
  old_par = list(mfrow=par()$mfrow)
  par(mfrow=c(1,2))
  gev_prob_plot(x$response, p, ...)
  gev_quantile_plot(x$response, p, ...)
  par(old_par)
}

#' Ghat function
#' @noRd
Ghat = function(z, mu, sig, xi){
  xizeros = round(xi, 6) == 0
  out = rep(NA, length(z))
  out[!xizeros] = exp(- (1 + xi*(z-mu)/sig)^(-1/xi) )
  out[xizeros] = exp(-exp(- (z-mu)/sig))
  out
}

#' Probability plot
#' @importFrom graphics abline plot
#' @noRd
gev_prob_plot = function(z, pars, ...){
  o = order(z)
  Gh = Ghat(z[o], pars$musum[o], pars$sigsum[o], pars$xisum[o])
  plot((1:length(z)/length(z)), sort(Gh), main = "Probability Plot", xlab ="Empirical", ylab = "Model", ...)
  abline(0, 1, col="red", lty=2, lwd=2)
}

#' Quantile Quantile plot
#' @importFrom graphics abline plot
#' @noRd
gev_quantile_plot = function(z, pars, ...){
  o = order(z)
  Ghinv = pars$musum - (pars$sigsum/pars$xisum) * (1 - (-log((1:length(z))/(length(z)+1)))^(-pars$xisum))
  plot(sort(Ghinv), z[o], main = "Quantile-Quantile Plot", xlab = "Empirical", ylab = "Model", ...)
  abline(0, 1, col="red", lty=2, lwd=2)
}










