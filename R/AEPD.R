#' @title Quantile-based asymmetric exponential power distribution
#'
#'
#' @description Density, cumulative distribution function, quantile function and random sample generation
#' from the quantile-based asymmetric exponential power distribution (AEPD) studied in Gijbels et al. (2019b). An alternative form of the density AEPD is also studied in Komunjer (2007).
#' @param y,q These are each a vector of quantiles.
#' @param mu This is the location parameter \eqn{\mu}.
#' @param phi This is the scale parameter  \eqn{\phi}.
#' @param alpha This is the index parameter  \eqn{\alpha}.
#' @param p This is the shape parameter, which must be positive.
#' @param n This is the number of observations, which must be a positive integer that has length 1.
#' @param beta This is a vector of probabilities.
#'
#' @return \code{\link{dAEPD}} provides the density, \code{\link{pAEPD}} provides the cumulative distribution function, \code{\link{qAEPD}} provides the quantile function, and \code{\link{rAEPD}} generates a random sample from the quantile-based asymmetric exponential power distribution.
#'
#' @references{
#'  Gijbels, I., Karim, R. and Verhasselt, A. (2019b). Quantile estimation in a generalized  asymmetric distributional setting. To appear in \emph{Springer Proceedings in Mathematics & Statistics, Proceedings of `SMSA 2019', the 14th Workshop on Stochastic Models, Statistics and their Application}, Dresden, Germany, in March 6--8, 2019. Editors: Ansgar Steland, Ewaryst Rafajlowicz, Ostap Okhrin.
#'
#'
#'  Komunjer, I., (2007). Asymmetric power distribution: theory and applications to risk measurement. \emph{Journal of Applied Econometrics}, \bold{22}(5), 891-921.
#'
#' }
#'
#' @name AEPD
NULL
#' @rdname AEPD
#' @export
dAEPD<-function(y,mu,phi,alpha,p){(alpha*(1-alpha)/(2*phi*gamma(1+1/p)))*ifelse(y>mu,exp(-alpha^p*((y-mu)/phi)^p), exp(-(1-alpha)^p*((mu-y)/phi)^p))}

#' @import zipfR
#' @rdname AEPD
#'
#' @export
pAEPD<-function(q,mu,phi,alpha,p){ifelse(q>mu,
      alpha+((1-alpha)/gamma(1/p))*zipfR::Igamma(1/p,alpha^p*(q-mu)^p/phi^p, lower=TRUE),
      alpha-(alpha/gamma(1/p))*zipfR::Igamma(1/p,(1-alpha)^p*((mu-q)^p/phi^p), lower=TRUE))}


#' @import zipfR
#' @rdname AEPD
#' @export
qAEPD<-function(beta,mu,phi,alpha,p){
  qf<-NA;
  for (i in 1:length(beta)) {
    qf[i]<-ifelse(beta[i]<alpha,mu-(phi/(1-alpha))*(zipfR::Igamma.inv(1/p, gamma(1/p)*(alpha-beta[i])/alpha,lower=TRUE))^(1/p),
                  mu+(phi/alpha)*(zipfR::Igamma.inv(1/p, gamma(1/p)*(beta[i]-alpha)/(1-alpha),lower=TRUE))^(1/p))}
  return(qf)
}




#' @import  zipfR
#' @rdname AEPD
#' @examples
#' # Quantile-based asymmetric exponential power distribution
#' # Density
#' rnum<-rnorm(100)
#' dAEPD(y=rnum,mu=0,phi=1,alpha=.5,p=2)
#'
#' # Distribution function
#' pAEPD(q=rnum,mu=0,phi=1,alpha=.5,p=2)
#'
#' # Quantile function
#' beta<-c(0.25,0.5,0.75)
#' qAEPD(beta=beta,mu=0,phi=1,alpha=.5,p=2)
#'
#' # random sample generation
#' rAEPD(n=100,mu=0,phi=1,alpha=.5,p=2)
#'
#' @export
rAEPD<-function(n,mu,phi,alpha,p){
  u<-runif(n, min = 0, max = 1)
  y<-NA;
  for (i in 1:length(u)) {
    y[i]<-ifelse(u[i]>alpha,mu+(phi/alpha)*(zipfR::Igamma.inv(1/p, gamma(1/p)*(u[i]-alpha)/(1-alpha),lower=TRUE))^(1/p),
                 mu-(phi/(1-alpha))*(zipfR::Igamma.inv(1/p, gamma(1/p)*(alpha-u[i])/alpha,lower=TRUE))^(1/p))}
  return(y)
}


#' @title Log-likelihood function for the quantile-based asymmetric exponential power distribution (AEPD) of distributions.
#' @description Log-Likelihood function \eqn{\ell_n(\mu,\phi,\alpha,p)=\ln[L_n(\mu,\phi,\alpha,p)]}
#' in the quantile-based asymmetric exponential power distribution (AEPD) of densities defined in Gijbels et al. (2019b).
#' @param y This is a vector of quantiles.
#' @param mu This is the location parameter \eqn{\mu}.
#' @param phi This is the scale parameter  \eqn{\phi}.
#' @param alpha This is the index parameter  \eqn{\alpha}.
#' @param p This is the shape parameter, which must be positive.
#' @return \code{\link{LogLikAEPD}} provides the realized value of the Log-likelihood function of the quantile-based asymmetric exponential power distribution.
#' @references{
#'  Gijbels, I., Karim, R. and Verhasselt, A. (2019b). Quantile estimation in a generalized  asymmetric distributional setting. To appear in \emph{Springer Proceedings in Mathematics & Statistics, Proceedings of `SMSA 2019', the 14th Workshop on Stochastic Models, Statistics and their Application}, Dresden, Germany, in March 6--8, 2019. Editors: Ansgar Steland, Ewaryst Rafajlowicz, Ostap Okhrin.
#'
#' }
#'
#' @name LogLikAEPD
NULL
#' @rdname LogLikAEPD
#'
#'
#' @examples
#' # Example
#' y<-rnorm(100)
#' LogLikAEPD(rexp(100,0.1),mu=10,phi=1,alpha=0.5,p=2)
#'
#'
#' @export
LogLikAEPD<- function(y,mu,phi,alpha,p){
  LL<-sum(log(dAEPD(y,mu,phi,alpha,p)))
  return(sum(LL[!is.infinite(LL)]))
}



#' @title Maximum likelihood estimation (MLE) for the quantile-based asymmetric exponential power distribution.
#' @description The log-likelihood function \eqn{\ell_n(\mu,\phi,\alpha,p)=\ln[L_n(\mu,\phi,\alpha,p)]}
#' and parameter estimation of \eqn{ \theta=(\mu,\phi,\alpha,p)} in the three parameter quantile-based asymmetric exponential power distribution
#' by using the maximum likelihood estimation are discussed in Gijbels et al. (2019b).
#' @param y This is a vector of quantiles.
#' @return The maximum likelihood estimate of parameter \eqn{\theta=(\mu,\phi,\alpha,p)} of the quantile-based asymmetric exponential power distribution.
#'
#' @import ald
#' @import nloptr
#'
#' @references{
#'  Gijbels, I., Karim, R. and Verhasselt, A. (2019b). Quantile estimation in a generalized  asymmetric distributional setting. To appear in \emph{Springer Proceedings in Mathematics & Statistics, Proceedings of `SMSA 2019', the 14th Workshop on Stochastic Models, Statistics and their Application}, Dresden, Germany, in March 6--8, 2019. Editors: Ansgar Steland, Ewaryst Rafajlowicz, Ostap Okhrin.
#' }
#' @name mleAEPD
NULL
#' @rdname mleAEPD
#' @examples
#' # Example
#' rnum=rnorm(100)
#' mleAEPD(rnum)
#' @export
mleAEPD<-function(y){
  theta0<-ald::mleALD(y, initial = NA)$par
  fn <- function(theta){
    mu=theta[1]
    phi=theta[2]
    alpha=theta[3]
    p=theta[4]
    LL<-log(dAEPD(y,mu,phi,alpha,p))
    return(-sum(LL[!is.infinite(LL)]))}
  Est.par=nloptr::bobyqa(x0=c(theta0,2),fn = fn ,lower = c(min(y),00.001,0.001,0.001), upper = c(max(y),Inf,1,Inf))$par
  estimated.LL<-LogLikAEPD(y=y,mu=Est.par[1],phi=Est.par[2],alpha=Est.par[3],p=Est.par[4])
  list(mu.MLE=Est.par[1], phi.MLE=Est.par[2],alpha.MLE=Est.par[3],p.MLE=Est.par[4],LogLikelihood=estimated.LL)
}



