#' GEV gradient
#' @noRd
GEVingrad <- function(mu, sig, xi,z){
  .e1 <- z - mu
  .e3 <- xi * .e1/sig
  .e4 <- 1 + .e3
  .e5 <- 1/xi
  .e6 <- 1 + .e5
  .e8 <- 1/.e4^.e5
  .e9 <- sig * .e4
  .e10 <- .e8 - xi * .e6
  list(mu = -(.e10/.e9), sig = -((.e10 * .e1/.e9 + 1)/sig),
       xi = -(((.e8 - 1) * log1p(.e3)/xi - .e1/(sig * .e4^.e6))/xi +
                .e6 * .e1/.e9))
}

#' GEV gradient with log(sigma)
#' @noRd
GEVingradLog <- function(mu, logsig, xi, z){
  # Taken from Deriv(GEVinlikLog,c("mu","logsig","xi"), combine="list")
  .e1 <- exp(logsig)
  .e2 <- z - mu
  .e4 <- xi * .e2/.e1
  .e5 <- 1 + .e4
  .e6 <- 1/xi
  .e7 <- 1 + .e6
  .e8 <- .e5 * .e1
  .e10 <- 1/.e5^.e6
  .e11 <- .e10 - xi * .e7
  list(mu = -(.e11/.e8), sig = -(.e11 * .e2/.e8 + 1),
       xi = -(((.e10 - 1) * log1p(.e4)/xi - .e2/(.e5^.e7 *
                                                   .e1))/xi + .e7 * .e2/.e8))
}

#' GEV gradient with log(sigma), when shape xi=0
#' @noRd
GEVingradLogxi0 = function(mu, logsig, xi, z){
  # from Deriv(GEVinlikLogxi0, c("mu","logsig","xi"), combine = "list")
  .e1 <- exp(logsig)
  .e2 <- z - mu
  .e5 <- exp(-(.e2/.e1)) - 1
  list(mu = -(.e5/.e1), sig = -(.e5 * .e2/.e1 + 1), xi = rep(0,length(mu)) )
}

#' GEV gradient when shape xi=0
#' @noRd
GEVingradxi0 = function(mu, sig, xi, z) {
  # From Deriv(GEVinlikxi0, c("mu","sig","xi"), combine = "list")
  .e1 <- z - mu
  .e4 <- exp(-(.e1/sig)) - 1
  list(mu = -(.e4/sig), sig = -((.e4 * .e1/sig + 1)/sig), xi = rep(0,length(mu)))
}

#' GEV gradient function
#'
#' @description Function to calculate gradient of parameters, accounting for non-stationarity
#'
#' @param theta vector of parameters, for each of location, scale and shape (in that order)
#' @param muX model matrix for the location parameter
#' @param sigX model matrix for the scale parameter
#' @param xiX model matrix for the shape parameter
#' @param z vector of observed maxima
#' @param mulink string describing what link function to use for the location parameter
#' @param siglink string describing what link function to use for the scale parameter
#' @param xilink string describing what link function to use for the shape parameter
#'
#' @details For the purposes of optimisation, this function takes the same arguments as \code{\link{gev_likfit}}
#' so it can be passed into the \code{\link[stats]{optim}} function.
#'
#' @return vector of gradients for each parameter
#' @export
gev_grad <- function(theta, muX, sigX, xiX, z,
                           mulink="identity", siglink="exponential", xilink="identity"){
  # Set up parameters
  mup = dim(muX)[2]
  sigp = dim(sigX)[2]
  xip = dim(xiX)[2]

  bindex <- c(rep(1,mup), rep(2,sigp), rep(3,xip))

  # Set parameters
  mubeta  <- theta[bindex==1]
  sigbeta <- theta[bindex==2]
  xibeta  <- theta[bindex==3]

  # Set up link function
  mulinkF = getLinkFunction(mulink)
  siglinkF = getLinkFunction(siglink)
  xilinkF = getLinkFunction(xilink)

  # get matrices for parameters with link functions
  mumat <- as.vector(mulinkF(mubeta, muX))
  sigmat <- as.vector(siglinkF(sigbeta, sigX))
  ximat <- as.vector(xilinkF(xibeta, xiX))
  gradients = list("mu"=mumat, "sig"=sigmat, "xi"=ximat)

  # Criteria for shape parameter to be zero (or close enough)
  xizeros = round(ximat,6)==0

  # hard code if siglink is exponential (which is likely)
  # and separating likelihood functions for when shape==0
  if(siglink!="exponential") {
    gradientsno0s <- GEVingrad(mumat[!xizeros],sigmat[!xizeros],ximat[!xizeros],z[!xizeros])
  } else if(siglink=="exponential") {
    gradientsno0s <- GEVingradLog(mumat[!xizeros],log(sigmat[!xizeros]),ximat[!xizeros],z[!xizeros])
  }

  if(siglink!="exponential") {
    gradients0s <- GEVingradxi0(mumat[xizeros],sigmat[xizeros],ximat[xizeros],z[xizeros])
  } else if(siglink=="exponential") {
    gradients0s <- GEVingradLogxi0(mumat[xizeros],log(sigmat[xizeros]),ximat[xizeros],z[xizeros])
  }

  # Separate gradients for different parameters
  gradients$mu[xizeros] = gradients0s$mu; gradients$mu[!xizeros] = gradientsno0s$mu
  gradients$sig[xizeros] = gradients0s$sig; gradients$sig[!xizeros] = gradientsno0s$sig
  gradients$xi[xizeros] = gradients0s$xi; gradients$xi[!xizeros] = gradientsno0s$xi

  # Multiply by model matrices
  gradients$mu <- colSums(as.matrix(muX * gradients$mu))
  gradients$sig <- colSums(as.matrix(sigX * gradients$sig))
  gradients$xi <- colSums(as.matrix(xiX * gradients$xi))

  # Return gradients as a vector
  return(unlist(gradients))
}

#' GEV likelihood function
#'
#' @description Function to calculate the log-likelihood of the GEV distribution for a given set of parameters,
#' accounting for non-stationarity
#'
#' @param theta vector of parameters, for each of location, scale and shape (in that order)
#' @param muX model matrix for the location parameter
#' @param sigX model matrix for the scale parameter
#' @param xiX model matrix for the shape parameter
#' @param z vector of observed maxima
#' @param mulink string describing what link function to use for the location parameter
#' @param siglink string describing what link function to use for the scale parameter
#' @param xilink string describing what link function to use for the shape parameter
#'
#' @details For the purposes of optimisation, this function takes the same arguments as \code{\link{gev_grad}}
#' so it can be passed into the \code{\link[stats]{optim}} function.
#'
#' @return vector of gradients for each parameter
#' @export
gev_likfit = function(theta, z, muX, sigX, xiX,
                         mulink="identity", siglink="exponential", xilink="identity"){
  mup = dim(muX)[2]
  sigp = dim(sigX)[2]
  xip = dim(xiX)[2]
  bindex <- c(rep(1,mup), rep(2,sigp), rep(3,xip))

  mubeta  <- theta[bindex==1]
  sigbeta <- theta[bindex==2]
  xibeta  <- theta[bindex==3]

  mulinkF = getLinkFunction(mulink)
  siglinkF = getLinkFunction(siglink)
  xilinkF = getLinkFunction(xilink)

  mumat <- mulinkF(mubeta, muX)
  sigmat <- siglinkF(sigbeta, sigX)
  ximat <- xilinkF(xibeta, xiX)

  z0 <- 1 + ximat*((z-mumat)/sigmat)

  if(any(z0<0)) return(-1e20)
  if(any(sigmat<0)) return(-1e20)
  if(any(ximat<0)) return(-1e20)

  loglike = -sum(log(sigmat) + (1 + 1/ximat)*log(z0) + (z0^(-1/ximat) ) )

  return(loglike)
}

#' From a string, output a link function as a function
#' @noRd
getLinkFunction = function(link="identity"){
  if(link == "identity") outf = function(B,X) X %*% B
  if(link == "exponential") outf = function(B,X) exp(X%*%B)
  return(outf)
}

#' Get initial conditions for GEV optimisation
#' @import stats
#' @noRd
gev_init = function(z,mup,sigp,xip,siglink="exponential"){

  init = numeric(sum(mup,sigp,xip))
  init[mup+1] = sqrt(6*var(z))/pi

  init[mup+sigp+1] = 0.1
  init[1] = mean(z) - 0.57722*init[mup+1]
  if(mup>1) if(siglink=="exponential" && init[mup+1]!=0) init[mup+1] = log(init[mup+1])

  return(init)
}

#' GEV optimisation function
#'
#' @description Optimise a non-stationary GEV distribution
#'
#' @param z vector of observed maxima
#' @param muX model matrix for the location parameter
#' @param sigX model matrix for the scale parameter
#' @param xiX model matrix for the shape parameter
#' @param gr gradient function to be used in \code{\link[stats]{optim}}, if \code{NULL}, a numerical gradient approximation is used
#' @param se logical; if \code{TRUE}, standard error of parameters will be output too
#' @param mulink string describing what link function to use for the location parameter
#' @param siglink string describing what link function to use for the scale parameter
#' @param xilink string describing what link function to use for the shape parameter
#' @param init initial conditions for optimisation
#
#'
#' @details This function uses the \code{\link[stats]{optim}} function with the BFGS method to
#' optimise a parameter vector corresponding to a given set of covariates for each GEV parameter.
#' @import stats
#' @return vector of gradients for each parameter
#' @export
gev_optim = function(z, muX, sigX, xiX, gr = gev_grad, se = TRUE,
                         mulink="identity", siglink="exponential", xilink="identity",
                         init=NULL){
  # Fit the GEV distribution with allowances of different model matrices, these are specified in advance
  # fit_GEV will create model matrices based on link functions and the data frame

  # Get initial conditions
  mup = dim(muX)[2]; sigp = dim(sigX)[2]; xip = dim(xiX)[2]
  if(is.null(init)) init = gev_init(z,mup,sigp,xip)
  bindex = c(rep(1, mup), rep(2, sigp), rep(3, xip))

  # Use optim with fnscale = -1 to maximise the likelihood from GEVlikX_diffX
  maxl = optim(par = init, fn = gev_likfit, gr=gr, hessian = TRUE, method = "BFGS",
               control=list(fnscale=-1, maxit=1000), z=z, muX=muX, sigX = sigX, xiX=xiX,
               mulink=mulink, siglink=siglink, xilink=xilink)
  val = maxl$value

  if(!se) se = NULL
  if(se){
    cov = solve(maxl$hessian)
    se = sqrt(diag(abs(cov)))
  }

  # Output normal optim output, but also details of the model, so that when the model is put
  # into other functions, they can read the link functions and model matrices etc
  return(list("value"=maxl$value, "par"=maxl$par, "se"=se, "mulink"=mulink,
              "par_list" = list("mu" = maxl$par[bindex==1], "sig" = maxl$par[bindex==2], "xi" = maxl$par[bindex==3]),
              "se_list" = list("mu" = se[bindex==1], "sig" = se[bindex==2], "xi" = se[bindex==3]),
              "siglink"=siglink, "xilink"=xilink, "modX" = list("mu"=muX, "sig"=sigX, "xi"=xiX)))

}

#' Fit a GEV model
#'
#' @description \code{gev_fit} is used to fit a Generalised Extreme Value (GEV) linear model. This is
#' a wrapper function, using \code{\link{gev_likfit}} and \code{\link{gev_grad}} to obtain parameter
#' estimates via maximum likelihood.
#'
#' @param formulas a list with three elements, all objects of class "\code{\link[stats]{formula}}",
#' for the location, scale and shape parameter in that order. A symbolic description of the model
#' @param data data frame containing all variables listed in \code{formulas}
#' @param response the observed maxima corresponding to those variables in \code{data}
#' @param rl_n the value of n of which to calculate the return level of period n, if \code{NULL}, the return
#' level is not calculated
#' @param links a list with three elements, all strings detailing which link function to use for the
#' location, scale and shape parameter in that order
#
#'
#' @details The \code{formulas} object have formulas specified symbolically, of the form \code{~ predictor1 + predictor2 + ...}
#' where both \code{predictor1} and \code{predictor2} are variables in \code{data}.
#'
#' The GEV model is fit by using \code{\link[stats]{optim}} with the BFGS method.
#'
#' @return an S3 object of class "\code{gev}", containing
#' \item{\code{loglik}}{value of the final log-likelihood at the values of \code{par}}
#' \item{\code{par}}{list of maximum likelihood parameter estimates for each GEV distribution parameter}
#' \item{\code{par_se}}{standard errors of the parameters in \code{par}}
#' \item{\code{modX}}{list of model matrices used for each parameter}
#' \item{\code{link}}{list of link functions inputted}
#' \item{\code{formula}}{list of formulae inputted}
#' \item{\code{big_formula}}{combination of unique elements in all formulae used}
#' \item{\code{data}}{data frame inputted}
#' \item{\code{response}}{response vector (observed maxima) inputted}
#' @import stats
#' @export
gev_fit = function(formulas, data, response, rl_n = 100,
                   links = list("identity","exponential","identity")){

    muX = model.matrix(formulas[[1]], data = data)
  sigX = model.matrix(formulas[[2]], data = data)
  xiX = model.matrix(formulas[[3]], data = data)
  mul = links[[1]]; sigl = links[[2]]; xil = links[[3]]

  GEVmod = gev_optim(response, muX, sigX, xiX,
                         mulink = mul, siglink=sigl, xilink=xil)

  f1 = labels(terms(formulas[[1]]))
  f2 = labels(terms(formulas[[2]]))
  f3 = labels(terms(formulas[[3]]))

  fall = unique(c(f1,f2,f3))

  if(length(fall)>0)  big.formula = reformulate(fall)
  if(length(fall)==0) big.formula = ~1

  out = list("loglik" = GEVmod$value, "par" = GEVmod$par_list, "par_se" = GEVmod$se_list,
             "modX" = GEVmod$modX, "link" = links, "formula" = formulas,
             "big_formula"=big.formula, "data" = data, "response" = response)

  if(!is.null(rl_n)){
    unique.data = model.matrix(big.formula, data = data)
    unique.data = unique.data[,colnames(unique.data)[colnames(unique.data)!="(Intercept)"], drop=FALSE]
    unique.data = data.frame(unique(unique.data))
    rl = gev_rl(rl_n, out, unique.data)
    if(ncol(rl) > 2) colnames(rl)[1:(ncol(rl)-1)] = fall
    out = c(out, "rl"=rl)
  } else out = c(out, "rl"=NULL)
  class(out) <- "gev"
  return(out)
}

#' Calculate Return Values
#'
#' @description Calculates return values for each set of covariates specified in both \code{data} and the model.
#'
#' @param n value for which the return value is calculated for
#' @param gevmodel object of class "\code{gev}", fit by \code{\link{gev_fit}}
#' @param data data frame to find return values for, using parameter estimates from \code{gevmodel}
#
#'
#' @details The \eqn{n}-level return value is defined as the largest expected return over the next \eqn{n} time period.
#' For example, if the maxima is measured yearly, this would be the \eqn{n}-year return value.
#'
#' @return a new data frame with an additional column \code{return}, containing the return values for the corresponding
#' other entries in the data frame
#' @export
gev_rl = function(n, gevmodel, data){
  parsums = gev_parsum(gevmodel, data)
  xizeros = round(parsums$xisum, 6) == 0

  if(length(parsums$xisum)!=1){
    z_n = numeric(length(parsums$xisum))
    z_n[!xizeros] = parsums$musum - (parsums$sigsum/parsums$xisum)*(1-(-log(1-1/n))^(-parsums$xisum))
    z_n[xizeros] = parsums$musum - parsums$sigsum * log(-log(1-(1/n)))
  } else {
    z_n = rep(NA, length(parsums$musum))
    z_n[!xizeros] = parsums$musum - (parsums$sigsum/parsums$xisum)*(1-(-log(1-1/n))^(-parsums$xisum))
    z_n[xizeros] = parsums$musum - parsums$sigsum * log(-log(1-(1/n)))
  }
  retlevels = data.frame(data, return=z_n)
  return(retlevels)
}

#' GEV parameter sums
#' @import stats
#' @noRd
gev_parsum = function(gevmodel, data){


  data = data.frame(data)

  mu  = gevmodel$par$mu
  sig = gevmodel$par$sig
  xi  = gevmodel$par$xi

  munames  = colnames(gevmodel$modX[[1]])[colnames(gevmodel$modX[[1]])!="(Intercept)"]
  signames = colnames(gevmodel$modX[[2]])[colnames(gevmodel$modX[[2]])!="(Intercept)"]
  xinames  = colnames(gevmodel$modX[[3]])[colnames(gevmodel$modX[[3]])!="(Intercept)"]

  muformula = if(length(munames)>0) reformulate(munames) else ~1
  sigformula = if(length(signames)>0) reformulate(signames) else ~1
  xiformula = if(length(xinames)>0) reformulate(xinames) else ~1

  muX = model.matrix(muformula, data = data)
  sigX = model.matrix(sigformula, data = data)
  xiX = model.matrix(xiformula, data = data)

  fmu = fsig = fxi = identity
  fmu = getLinkFunction(gevmodel$link[[1]])
  fsig = getLinkFunction(gevmodel$link[[2]])
  fxi = getLinkFunction(gevmodel$link[[3]])

  musum = fmu(as.matrix(mu), muX)
  sigsum = fsig(as.matrix(sig), sigX)
  xisum = fxi(as.matrix(xi), xiX)

  return(list("musum"=musum,"sigsum"=sigsum,"xisum"=xisum))
}
