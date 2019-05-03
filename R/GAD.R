#' @title Generalized quantile-based asymmetric family
#'
#' @description Density, cumulative distribution function, quantile function and random sample generation
#' from the generalized quantile-based asymmetric family of densities defined in Gijbels et al. (2019b).
#'
##' @param y,q These are each a vector of quantiles.
#' @param eta This is the location parameter \eqn{\eta}.
#' @param phi This is the scale parameter  \eqn{\phi}.
#' @param alpha This is the index parameter  \eqn{\alpha}.
#' @param g This is the "link" function. The function \eqn{g} is to be differentiated. Therefore, \eqn{g} must be written as a function. For example, {g<-function(y)\{log(y)\}} for log link function.
#' @param f This is the reference density function \eqn{f} which is a standard version of a unimodal and symmetric around 0 density.
#'
#' @references{
#'  Gijbels, I., Karim, R. and Verhasselt, A. (2019b). Quantile estimation in a generalized  asymmetric distributional setting. To appear in \emph{Springer Proceedings in Mathematics & Statistics, Proceedings of `SMSA 2019', the 14th Workshop on Stochastic Models, Statistics and their Application}, Dresden, Germany, in March 6--8, 2019. Editors: Ansgar Steland, Ewaryst Rafajlowicz, Ostap Okhrin.
#' }
#'
#'
#' @name GAD
NULL
#' @import Deriv
#' @rdname GAD
#'
#' @export
dGAD<-function(y,eta,phi,alpha,f,g){
  g.prime<- Deriv::Deriv(g)
  return(2*(alpha*(1-alpha)*g.prime(y)/phi)*ifelse(y>eta,f(alpha*(g(y)-g(eta))/phi),f((1-alpha)*(g(eta)-g(y))/phi)))}




#' @param F This is the cumulative distribution function \eqn{F} of the unimodal and symmetric around 0 reference density function \eqn{f}.
#' @import zipfR
#' @rdname GAD
#'
#' @export
pGAD<-function(q,eta,phi,alpha,F,g){pQBAD(g(q),g(eta),phi,alpha,F)}

#' @param QF This is the quantile function of the reference density \eqn{f}.
#' @param lower	This is the lower limit of the domain (support of the random variable) \eqn{f_{\alpha}^g(y;\eta,\phi)}, default {-Inf}.
#' @param upper	This is the upper limit of the domain (support of the random variable) \eqn{f_{\alpha}^g(y;\eta,\phi)}, default {Inf}.
#' @param beta This is a vector of probabilities.
#' @import zipfR
#' @import GoFKernel
#' @rdname GAD
#'
#' @export
qGAD<-function(beta,eta,phi,alpha,F,g,QF=NULL,lower = -Inf, upper = Inf){
  qGAD.inv<-function(beta,eta,phi,alpha,F,g,QF=NULL,lower = lower, upper = upper){
    qf_z<-qQBAD(beta=beta,mu=g(eta),phi=phi,alpha=alpha,F=F,QF=QF)
    g.inv<-GoFKernel::inverse(g,lower,upper)
    return(g.inv(qf_z))
  }

  qf.y<-NA
  for (j in 1:length(beta)){
    qf.y[j]<-qGAD.inv(beta=beta[j],eta=eta,phi=phi,alpha=alpha,F=F,g=g,lower = lower, upper = upper)
  }
  return(qf.y)
}


#' @param n This is the number of observations, which must be a positive integer that has length 1.
#' @rdname GAD
#'
#'
#' @examples

#' # Example 1: Let F be a standard normal cumulative distribution function then
#' f_N<-function(s){dnorm(s, mean = 0,sd = 1)} # density function of N(0,1)
#' F_N<-function(s){pnorm(s, mean = 0,sd = 1)} # distribution function of N(0,1)
#' QF_N<-function(beta){qnorm(beta, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)}
#'
#' # For identiy link function
#' g_id<-function(y){y}
#' # For log-link function
#' g_log<-function(y){log(y)}
#'
#' rnum<-rnorm(100)
#' beta=c(0.25,0.50,0.75)
#'
#' # Density
#' dGAD(y=rnorm(100),eta=10,phi=1,alpha=0.5,f=f_N,g=g_id) # For identity link
#' dGAD(y=rexp(100,0.1),eta=10,phi=1,alpha=0.5,f=f_N,g=g_log) # For log-link
#'
#' # Distribution function
#' pGAD(q=rnorm(100),eta=0,phi=1,alpha=.5,F=F_N,g=g_id) # For identity link
#' pGAD(q=rexp(100,0.1),eta=10,phi=1,alpha=.5,F=F_N,g=g_log) # For log-link
#'
#' # Quantile function
#' qGAD(beta=beta,eta=0,phi=1,alpha=0.5,F=F_N,g=g_id) # For identity link
#' qGAD(beta=beta,eta=10,phi=1,alpha=0.5,F=F_N,g=g_log,lower = 0, upper = Inf)  # For log-link
#'
#' # random sample generation
#' rGAD(n=100,eta=0,phi=1,alpha=.5,F=F_N,g=g_id ,lower = -Inf, upper = Inf,QF=NULL) # For identity link
#' rGAD(n=100,eta=10,phi=1,alpha=.5,F=F_N,g=g_log ,lower =0, upper = Inf,QF=NULL)   # For log-link
#'
#'
#' # Example 2: Let F be a standard Laplace cumulative distribution function then
#' f_La<-function(s){0.5*exp(-abs(s))} # density function of Laplace(0,1)
#' F_La<-function(s){0.5+0.5*sign(s)*(1-exp(-abs(s)))} # distribution function of Laplace(0,1)
#' QF_La<-function(beta){-sign(beta-0.5)*log(1-2*abs(beta-0.5))}
#'
#' # For identiy link function
#' g_log<-function(y){log(y)}
#' beta=c(0.25,0.50,0.75)
#'
#' # Density
#' dGAD(y=rnorm(100),eta=10,phi=1,alpha=0.5,f=f_La,g=g_id) # For identity-link
#' dGAD(y=rexp(100,0.1),eta=10,phi=1,alpha=0.5,f=f_La,g=g_log) # For log-link
#'
#' # Distribution function
#' pGAD(q=rnum,eta=0,phi=1,alpha=.5,F=F_La,g=g_id) # For identity-link
#' pGAD(q=rexp(100,0.1),eta=10,phi=1,alpha=.5,F=F_La,g=g_log) # For log-link
#'
#' # Quantile function
#' qGAD(beta=beta,eta=0,phi=1,alpha=0.5,F=F_La,g=g_id,lower = -Inf, upper = Inf) # For identity link
#' qGAD(beta=beta,eta=10,phi=1,alpha=0.5,F=F_La,g=g_log,lower = 0, upper = Inf) # For log-link
#'
#' # random sample generation
#' rGAD(n=100,eta=0,phi=1,alpha=.5,F=F_La,g=g_id) # For identity link
#' rGAD(n=100,eta=10,phi=1,alpha=.5,F=F_La,g=g_log ,lower =0, upper = Inf,QF=NULL)   # For log-link
#'
#'
#'
#'
#' @export
rGAD<-function(n,eta,phi,alpha,F,g,lower = -Inf, upper = Inf,QF=NULL){
  g.inv<-GoFKernel::inverse(g,lower,upper)
  qf_z<-rQBAD(n=n,mu=g(eta),phi=phi,alpha=alpha,F=F,QF=QF)
  rnum<-NA
  for (i in 1:length(qf_z)){
    rnum[i]<-g.inv(qf_z[i])}
  return(rnum)
}

#' @title Log-likelihood function for the generalized quantile-based asymmetric family of distributions.
#' @description Log-Likelihood function \eqn{\ell_n(\eta,\phi,\alpha)=\ln[L_n(\eta,\phi,\alpha)]}
#' in the three parameter generalized quantile-based asymmetric family of densities defined in Gijbels et al. (2019b).
#' @param y This is a vector of quantiles.
#' @param eta This is the location parameter \eqn{\eta}.
#' @param phi This is the scale parameter  \eqn{\phi}.
#' @param alpha This is the index parameter  \eqn{\alpha}.
#' @param g This is the "link" function. The function \eqn{g} is to be differentiated. Therefore, \eqn{g} must be written as a function. For example, {g<-function(y)\{log(y)\}} for log link function.
#'
#' @param f This is the reference density function \eqn{f} which is a standard version of a unimodal and symmetric around 0 density.
#' @return \code{\link{LogLikGAD}} provides the realized value of the Log-likelihood function of the generalized quantile-based asymmetric family of distributions.
#'
#' @references{
#'
#'
#'  Gijbels, I., Karim, R. and Verhasselt, A. (2019b). Quantile estimation in a generalized  asymmetric distributional setting. To appear in \emph{Springer Proceedings in Mathematics & Statistics, Proceedings of `SMSA 2019', the 14th Workshop on Stochastic Models, Statistics and their Application}, Dresden, Germany, in March 6--8, 2019. Editors: Ansgar Steland, Ewaryst Rafajlowicz, Ostap Okhrin.
#' }
#'
#'
#' @name LogLikGAD
NULL
#' @rdname LogLikGAD
#' @examples
#' # Example 1: Let F be a standard normal cumulative distribution function then
#' f_N<-function(s){dnorm(s, mean = 0,sd = 1)} # density function of N(0,1)
#' y<-rnorm(100)
#' g_id<-function(y){y}
#' g_log<-function(y){log(y)}
#' LogLikGAD(y,eta=0,phi=1,alpha=0.5,f=f_N,g=g_id) # For identity-link
#' LogLikGAD(rexp(100,0.1),eta=10,phi=1,alpha=0.5,f=f_N,g=g_log) # For log-link
#'
#'
#' # Example 2: Let F be a standard Laplace cumulative distribution function then
#' f_La<-function(s){0.5*exp(-abs(s))} # density function of Laplace(0,1)
#' LogLikGAD(y,eta=0,phi=1,alpha=0.5,f=f_La,g=g_id) # For identity-link
#' LogLikGAD(rexp(100,0.1),eta=10,phi=1,alpha=0.5,f=f_La,g=g_log) # For log-link
#' @export
LogLikGAD<- function(y,eta,phi,alpha,f,g){
  LL<-log(dGAD(y,eta,phi,alpha,f,g))
    return(sum(LL[!is.infinite(LL)]))
}



#' @title Maximum likelihood estimation (MLE) for the generalized quantile-based asymmetric family of distributions (GAD).
#' @description The log-likelihood function \eqn{\ell_n(\eta,\phi,\alpha)=\ln[L_n(\eta,\phi,\alpha)]}
#' and parameter estimation of \eqn{ \theta=(\eta,\phi,\alpha)} in the three parameter generalized quantile-based asymmetric family of densities
#' by using the maximum likelihood estimation are discussed in Gijbels et al. (2019b).
#' @param y This is a vector of quantiles.
#' @param f This is the reference density function \eqn{f} which is a standard version of a unimodal and symmetric around 0 density.
#' @param g This is the "link" function. The function \eqn{g} is to be differentiated. Therefore, \eqn{g} must be written as a function. For example, {g<-function(y)\{log(y)\}} for log link function.
#' @param lower	This is the lower limit of the domain (support of the random variable) \eqn{f_{\alpha}^g(y;\eta,\phi)}, default {-Inf}.
#' @param upper	This is the upper limit of the domain (support of the random variable) \eqn{f_{\alpha}^g(y;\eta,\phi)}, default {Inf}.
#'
#' @return The maximum likelihood estimate of parameter \eqn{\theta=(\eta,\phi,\alpha)} of the generalized quantile-based asymmetric family of densities
#'
#'
#'
#' @references{
#'
#'
#'  Gijbels, I., Karim, R. and Verhasselt, A. (2019b). Quantile estimation in a generalized  asymmetric distributional setting. To appear in \emph{Springer Proceedings in Mathematics & Statistics, Proceedings of `SMSA 2019', the 14th Workshop on Stochastic Models, Statistics and their Application}, Dresden, Germany, in March 6--8, 2019. Editors: Ansgar Steland, Ewaryst Rafajlowicz, Ostap Okhrin.
#' }
#'
#'
#' @import GoFKernel
#'
#' @name mleGAD
NULL
#' @rdname mleGAD
#' @examples
#' # Example 1: Let F be a standard normal cumulative distribution function then
#' f_N<-function(s){dnorm(s, mean = 0,sd = 1)} # density function of N(0,1)
#' y<-rnorm(100)
#' g_id<-function(y){y}
#' g_log<-function(y){log(y)}
#' mleGAD(y,f=f_N,g=g_id) # For identity-link
#' mleGAD(rexp(100,0.1),f=f_N,g=g_log,lower = 0, upper = Inf) # For log-link
#'
#'
#' # Example 2: Let F be a standard Laplace cumulative distribution function then
#' f_La<-function(s){0.5*exp(-abs(s))} # density function of Laplace(0,1)
#' mleGAD(y,f=f_La,g=g_id) # For identity-link
#' mleGAD(rexp(100,0.1),f=f_La,g=g_log,lower = 0, upper = Inf) # For log-link
#' @export
mleGAD<-function(y,f,g,lower = -Inf, upper = Inf){
  z=g(y)
  g.inv<-GoFKernel::inverse(g,lower,upper)
  mle_Z<-mleQBAD(z,f)
  eta.MLE<-g.inv(mle_Z$mu.MLE)
  LL=LogLikGAD(y,eta=eta.MLE,phi=mle_Z$phi.MLE,alpha=mle_Z$alpha.MLE,f,g)
  list(eta.MLE=eta.MLE,phi.MLE=mle_Z$phi.MLE,alpha.MLE=mle_Z$alpha.MLE,LogLikelihood=LL)
}



