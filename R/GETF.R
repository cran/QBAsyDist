#' @title Generalized tick-exponential family
#'
#' @description Density, cumulative distribution function, quantile function and random sample generation
#' from the generalized tick-exponential family (GTEF) of densities discusse in Gijbels et al. (2019b).
#'
##' @param y,q These are each a vector of quantiles.
#' @param eta This is the location parameter \eqn{\eta}.
#' @param phi This is the scale parameter  \eqn{\phi}.
#' @param alpha This is the index parameter  \eqn{\alpha}.
#' @param p This is the shape parameter, which must be positive.
#' @param g This is the "link" function. The function \eqn{g} is to be differentiated. Therefore, \eqn{g} must be written as a function. For example, {g<-function(y)\{log(y)\}} for log link function.
#' @return \code{\link{dGTEF}} provides the density, \code{\link{pGTEF}} provides the cumulative distribution function, \code{\link{qGTEF}} provides the quantile function, and \code{\link{rGTEF}} generates a random sample from the generalized tick-exponential family of densities.
#' The length of the result is determined by \eqn{n} for \code{\link{rGTEF}}, and is the maximum of the lengths of the numerical arguments for the other functions.
#'
#'
#' @references{
#'  Gijbels, I., Karim, R. and Verhasselt, A. (2019b). Quantile estimation in a generalized  asymmetric distributional setting. To appear in \emph{Springer Proceedings in Mathematics & Statistics, Proceedings of `SMSA 2019', the 14th Workshop on Stochastic Models, Statistics and their Application}, Dresden, Germany, in March 6--8, 2019. Editors: Ansgar Steland, Ewaryst Rafajlowicz, Ostap Okhrin.
#' }
#'
#'
#' @name GTEF
NULL
#' @import Deriv
#' @rdname GTEF
#'
#' @export
dGTEF<-function(y,eta,phi,alpha,p,g){
  g.prime<- Deriv::Deriv(g)
  return((alpha*(1-alpha)*g.prime(y)/(2*phi*gamma(1+1/p)))*ifelse(y>eta,exp(-alpha^p*((g(y)-g(eta))/phi)^p),
                                                                  exp(-(1-alpha)^p*((g(eta)-g(y))/phi)^p)))}


#' @import zipfR
#' @rdname GTEF
#'
#' @export
pGTEF<-function(q,eta,phi,alpha,p,g){
  mu=g(eta)
  z<-g(q)
  return(pAEPD(z,mu,phi,alpha,p))
  }

#' @param lower	This is the lower limit of the domain (support of the random variable) \eqn{f_{\alpha}^g(y;\eta,\phi)}, default {-Inf}.
#' @param upper	This is the upper limit of the domain (support of the random variable) \eqn{f_{\alpha}^g(y;\eta,\phi)}, default {Inf}.
#' @param beta This is a vector of probabilities.
#' @import zipfR
#' @import GoFKernel
#' @rdname GTEF
#'
#' @export
qGTEF<-function(beta,eta,phi,alpha,p,g,lower = -Inf, upper = Inf){
  qGTEF.inv<-function(beta,eta,phi,alpha,p,g,lower = -Inf, upper = Inf){
    mu=g(eta)
       qf_z<-ifelse(beta<alpha,mu-(phi/(1-alpha))*(Igamma.inv(1/p, gamma(1/p)*(alpha-beta)/alpha,lower=TRUE))^(1/p),
                      mu+(phi/alpha)*(Igamma.inv(1/p, gamma(1/p)*(beta-alpha)/(1-alpha),lower=TRUE))^(1/p))
    g.inv<-GoFKernel::inverse(g,lower,upper)
    return(g.inv(qf_z))
  }
  qf.y<-NA
  for (j in 1:length(beta)){
    qf.y[j]<-qGTEF.inv(beta=beta[j],eta=eta,phi=phi,alpha=alpha,p=p,g=g,lower = lower, upper = upper )
  }
  return(qf.y)
}



#' @param n This is the number of observations, which must be a positive integer that has length 1.
#' @rdname GTEF
#'
#'
#' @examples
#' # For identiy link function
#' y=rnorm(100)
#' g_id<-function(y){y}
#' dGTEF(y,eta=0,phi=1,alpha=0.5,p=2,g=g_id)
#'
#' # cumulative distribution function
#' pGTEF(q=y,eta=10,phi=1,alpha=0.5,p=2,g=g_id)
#'
#' # Quantile function
#' beta=c(0.25,0.5,0.75)
#' qGTEF(beta=beta,eta=10,phi=1,alpha=0.5,p=2,g=g_id)
#'
#' # random sample generation
#' rGTEF(n=100,eta=10,phi=1,alpha=.5,p=2,g=g_id,lower = -Inf, upper = Inf)
#'
#' # For log link function
#' y=rexp(100)
#' g_log<-function(y){log(y)}
#' dGTEF(y,eta=10,phi=1,alpha=0.5,p=2,g=g_log)
#'
#' # cumulative distribution function
#' pGTEF(q=y,eta=10,phi=1,alpha=0.5,p=2,g=g_log)
#'
#' # Quantile function
#' g_log<-function(y){log(y)}
#' #' beta=c(0.25,0.5,0.75)
#' qGTEF(beta=beta,eta=10,phi=1,alpha=0.5,p=2,g=g_log,lower = 0, upper = Inf)
#'
#' # random sample generation
#' rGTEF(n=100,eta=10,phi=1,alpha=.5,p=2,g=g_log,lower = 0, upper = Inf)
#'
#'
#'
#'
#' @export
rGTEF<-function(n,eta,phi,alpha,p,g,lower = -Inf, upper = Inf){
  mu=g(eta)
  g.inv<-GoFKernel::inverse(g,lower,upper)
  u<-runif(n, min = 0, max = 1)
  qf_z<-NA
  rnum<-NA
  for (i in 1:length(u)) {
    qf_z[i]<-ifelse(u[i]<alpha,mu-(phi/(1-alpha))*(Igamma.inv(1/p, gamma(1/p)*(alpha-u[i])/alpha,lower=TRUE))^(1/p),
                    mu+(phi/alpha)*(Igamma.inv(1/p, gamma(1/p)*(u[i]-alpha)/(1-alpha),lower=TRUE))^(1/p))
    rnum[i]<-g.inv( qf_z[i])
  }
  return(rnum)}

#' @title Log-likelihood function for the generalized tick-exponential family (GTEF) of distributions.
#' @description Log-Likelihood function \eqn{\ell_n(\eta,\phi,\alpha,p)=\ln[L_n(\eta,\phi,\alpha,p)]}
#' in the generalized tick-exponential family of densities discussed in Gijbels et al. (2019b).
#' @param y This is a vector of quantiles.
#' @param eta This is the location parameter \eqn{\eta}.
#' @param phi This is the scale parameter  \eqn{\phi}.
#' @param alpha This is the index parameter  \eqn{\alpha}.
#' @param g This is the "link" function. The function \eqn{g} is to be differentiated. Therefore, \eqn{g} must be written as a function. For example, {g<-function(y)\{log(y)\}} for log link function.
#' @param p This is the shape parameter, which must be positive.
#' @return \code{\link{LogLikGAD}} provides the realized value of the Log-likelihood function of the generalized quantile-based asymmetric family of distributions.
#'
#' @references{
#'  Gijbels, I., Karim, R. and Verhasselt, A. (2019b). Quantile estimation in a generalized  asymmetric distributional setting. To appear in \emph{Springer Proceedings in Mathematics & Statistics, Proceedings of `SMSA 2019', the 14th Workshop on Stochastic Models, Statistics and their Application}, Dresden, Germany, in March 6--8, 2019. Editors: Ansgar Steland, Ewaryst Rafajlowicz, Ostap Okhrin.
#' }
#'
#'
#' @name LogLikGTEF
NULL
#' @rdname LogLikGTEF
#' @examples
#' # Examples
#' y<-rnorm(100)
#' g_id<-function(y){y}
#' g_log<-function(y){log(y)}
#' LogLikGTEF(y,eta=0,phi=1,alpha=0.5,p=2,g=g_id) # For identity-link
#' LogLikGTEF(rexp(100,0.1),eta=10,phi=1,alpha=0.5,p=2,g=g_log) # For log-link
#'
#' @export
LogLikGTEF<- function(y,eta,phi,alpha,p,g){
  LL<-log(dGTEF(y,eta,phi,alpha,p,g))
  return(sum(LL[!is.infinite(LL)]))
}


#' @title Maximum likelihood estimation (MLE) for the generalized tick-exponential family (GTEF) of distributions.
#' @description The log-likelihood function \eqn{\ell_n(\eta,\phi,\alpha,p)=\ln[L_n(\eta,\phi,\alpha,p)]}
#' and parameter estimation of \eqn{ \theta=(\eta,\phi,\alpha,p)} in the generalized tick-exponential family of distributions
#' by using the maximum likelihood estimation are discussed in Gijbels et al. (2019b).
#' @param y This is a vector of quantiles.
#' @param g This is the "link" function. The function \eqn{g} is to be differentiated. Therefore, \eqn{g} must be written as a function. For example, {g<-function(y)\{log(y)\}} for log link function.
#' @param lower	This is the lower limit of the domain (support of the random variable) \eqn{f_{\alpha}^g(y;\eta,\phi)}, default {-Inf}.
#' @param upper	This is the upper limit of the domain (support of the random variable) \eqn{f_{\alpha}^g(y;\eta,\phi)}, default {Inf}.

#'
#' @return The maximum likelihood estimate of parameter \eqn{\theta=(\eta,\phi,\alpha,p)} of the generalized tick-exponential family of distributions.
#'
#' @import ald
#' @import GoFKernel
#'
#'
#' @references{
#'  Gijbels, I., Karim, R. and Verhasselt, A. (2019b). Quantile estimation in a generalized  asymmetric distributional setting. To appear in \emph{Springer Proceedings in Mathematics & Statistics, Proceedings of `SMSA 2019', the 14th Workshop on Stochastic Models, Statistics and their Application}, Dresden, Germany, in March 6--8, 2019. Editors: Ansgar Steland, Ewaryst Rafajlowicz, Ostap Okhrin.
#' }
#'
#'
#' @name mleGTEF
NULL
#' @rdname mleGTEF
#' @examples
#' # Example
#' rnum=rnorm(100)
#' g_id<-function(y){y}
#' g_log<-function(y){log(y)}
#' mleGTEF(rnum,g_id) # For identity-link
#' mleGTEF(rexp(100),g_log,lower = 0, upper = Inf) # For log-link
#' @export
mleGTEF<-function(y,g,lower = -Inf, upper = Inf){
  z=g(y)
  g.inv<-GoFKernel::inverse(g,lower,upper)
  mle_Z<-mleAEPD(z)
  eta.MLE<-g.inv(mle_Z$mu.MLE)
  LL=LogLikGTEF(y,eta=eta.MLE,phi=mle_Z$phi.MLE,alpha=mle_Z$alpha.MLE,p=mle_Z$p.MLE,g)
  list(eta.MLE=eta.MLE,phi.MLE=mle_Z$phi.MLE,alpha.MLE=mle_Z$alpha.MLE,p.MLE=mle_Z$p.MLE,LogLikelihood=LL)
}





