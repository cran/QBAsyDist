#' @title Quantile-based asymmetric Laplace distribution
#' @description Density, cumulative distribution function, quantile function and random sample generation
#' for the quantile-based asymmetric Laplace distribution (ALaD) discussed in Yu and Zhang (2005) and Gijbels et al. (2019a).
#' @param y,q These are each a vector of quantiles.
#' @param mu This is the location parameter \eqn{\mu}.
#' @param phi This is the scale parameter  \eqn{\phi}.
#' @param alpha This is the index parameter  \eqn{\alpha}.
#' @param n This is the number of observations, which must be a positive integer that has length 1.
#' @param beta This is a vector of probabilities.
#'
#' @return \code{\link{dALaD}} provides the density, \code{\link{pALaD}} provides the cumulative distribution function, \code{\link{qALaD}} provides the quantile function, and \code{\link{rALaD}} generates a random sample from the quantile-based asymmetric Laplace distribution.
#' The length of the result is determined by \eqn{n} for \code{\link{rALaD}}, and is the maximum of the lengths of the numerical arguments for the other functions.
#'
#' @references{
#'  Gijbels, I., Karim, R. and Verhasselt, A. (2019a). On quantile-based asymmetric family of distributions: properties and inference. \emph{International Statistical Review}, \url{https://doi.org/10.1111/insr.12324}.
#'
#'
#'  Yu., K, and Zhang, J. (2005). A three-parameter asymmetric Laplace distribution and its extension. \emph{Communications in Statistics–Theory and Methods}, \bold{34}(9-10), 1867–1879.
#' }
#'
#'
#' @seealso \code{\link{dQBAD}},  \code{\link{pQBAD}},  \code{\link{qQBAD}},  \code{\link{rQBAD}}
#'
#' @name ALaD
NULL
#' @rdname ALaD
#'
#' @export

dALaD<-function(y,mu,phi,alpha){alpha*(1-alpha)*(1/phi)*ifelse(y>mu,exp(-alpha*(y-mu)/phi),exp(-(1-alpha)*(mu-y)/phi))}


#' @rdname ALaD
#' @export
pALaD<-function(q,mu,phi,alpha){ifelse(q>mu,(1-(1-alpha)*exp(-alpha*(q-mu)/phi)),(alpha*exp((1-alpha)*(q-mu)/phi)))}


#' @examples
#' # Density
#' rnum<-rnorm(100)
#' dALaD(y=rnum,mu=0,phi=1,alpha=.5)
#'
#' # Distribution function
#'  pALaD(q=rnum,mu=0,phi=1,alpha=.5)
#'
#' # Quantile function
#' beta<-c(0.25,0.5,0.75)
#' qALaD(beta=beta,mu=0,phi=1,alpha=.5)
#'
#' # random sample generation
#' rALaD(n=100,mu=0,phi=1,alpha=.5)
#'
#'
#' @import stats
#' @rdname ALaD
#' @export
qALaD<-function(beta,mu,phi,alpha){ifelse(beta<alpha,mu+(phi/(1-alpha))*log(beta/alpha),(mu-(phi/alpha)*log((1-beta)/(1-alpha))))}


#' @import stats
#' @rdname ALaD
#' @export
rALaD<-function(n,mu,phi,alpha){
  u=runif(n, min = 0, max = 1)
  y=ifelse(u>alpha,mu-(phi/alpha)*log((1-u)/(1-alpha)),mu+phi/(1-alpha)*log(u/alpha));
  return(y)
}


#' @title Moments estimation for the quantile-based asymmetric Laplace distribution.
#' @description Mean, variance, skewness, kurtosis and moments about the location parameter (i.e., \eqn{\alpha}th quantile) of the quantile-based asymmetric Laplace distribution studied in Gijbels et al. (2019a) useful for quantile regression with location parameter equal to \eqn{\mu}, scale parameter \eqn{\phi} and index parameter \eqn{\alpha}.
#' @param mu This is the location parameter \eqn{\mu}.
#' @param phi This is the scale parameter  \eqn{\phi}.
#' @param alpha This is the index parameter  \eqn{\alpha}.
#' @return \code{\link{meanALaD}} provides the mean, \code{\link{varALaD}} provides the variance, \code{\link{skewALaD}} provides the skewness, \code{\link{kurtALaD}} provides the kurtosis, and  \code{\link{momentALaD}} provides the \eqn{r}th moment about the location parameter \eqn{\mu} of the quantile-based asymmetric Laplace distribution.
#'
#'
#' @references{
#'  Gijbels, I., Karim, R. and Verhasselt, A. (2019a). On quantile-based asymmetric family of distributions: properties and inference. \emph{International Statistical Review}, \url{https://doi.org/10.1111/insr.12324}.
#' }
#'
#' @name momentALaD
NULL
#' @rdname momentALaD
#' @export
meanALaD<- function(mu,phi,alpha){mu+phi*(1-2*alpha)/(alpha*(1-alpha))}


#' @rdname momentALaD
#' @export
varALaD<- function(mu,phi,alpha){(phi^2/(alpha^2*(1-alpha)^2))*((1-2*alpha)^2+2*alpha*(1-alpha))}

#' @rdname momentALaD
#' @export
skewALaD<- function(alpha){
  nu<-(1-2*alpha)*(1-2*alpha)^2*2
  deno<-((1-2*alpha)^2+alpha*(1-alpha)*2)^(3/2)
  skew<-(nu/deno)
  return(skew)
}



#' @rdname momentALaD
#'
#' @export
kurtALaD<- function(alpha){
  nu<-((1-alpha)^5+alpha^5)*24-(1-2*alpha)^2*(4*(1-2*alpha+2*alpha^2)*6-6*(1-3*alpha+3*alpha^2)*2+3*(1-2*alpha)^2)
  deno<-((1-2*alpha)^2+alpha*(1-alpha)*2)^(2)
  kurto<-(nu/deno)
  return(kurto)
}

#' @param r This is a value which is used to calculate \eqn{r}th moment about \eqn{\mu}.
#' @import stats
#' @examples
#' # Example
#' meanALaD(mu=0,phi=1,alpha=0.5)
#' varALaD(mu=0,phi=1,alpha=0.5)
#' skewALaD(alpha=0.5)
#' kurtALaD(alpha=0.5)
#' momentALaD(phi=1,alpha=0.5,r=1)
#' @rdname momentALaD
#' @export
momentALaD<- function(phi,alpha,r){
  integrand <- function(x) {x^r*0.5*exp(-abs(x))}
  mu.r<-2*stats::integrate(integrand, lower = 0, upper = Inf)$ value
  return(phi^r*((1-alpha)^(r+1)+(-1)^r*alpha^(r+1))*mu.r/(alpha^r*(1-alpha)^r))
}


#' @title Method of moments (MoM) estimation for the quantile-based asymmetric Laplace distribution.
#' @description Parameter estimation in the quantile-based asymmetric Laplace distribution by using method of moments is studied in Gijbels et al. (2019a).
#' @param y This is a vector of quantiles.
#' @param alpha This is the index parameter  \eqn{\alpha}.  If \eqn{\alpha} is unknown, the it should be NULL which is default option. In this case, the sample skewness will be used to estimate \eqn{\alpha}. If \eqn{\alpha} is known, then the value of \eqn{\alpha} has to be specified in the function.
#'
#'
#'
#' @import stats
#' @return \code{\link{momALaD}} provides the method of moments estimates of the unknown parameters of the distribution.
#'
#'
#' @references{
#'  Gijbels, I., Karim, R. and Verhasselt, A. (2019a). On quantile-based asymmetric family of distributions: properties and inference. \emph{International Statistical Review}, \url{https://doi.org/10.1111/insr.12324}.
#' }
#'
#' @name momALaD
NULL
#' @rdname momALaD
#' @examples
#' # Example
#' y=rnorm(100)
#' momALaD(y=y,alpha=0.5) # If alpha is known with alpha=0.5
#' momALaD(y=y) # If alpha is unknown
#'
#' @export
momALaD<-function(y,alpha=NULL){
  muMoM<-function(y,f,alpha=NULL){
    if (is.null(alpha)){
      m1<-mean(y)
      m2<-sum(y^2)/length(y)
      cm2<-sum((y-m1)^2)/length(y)
      cm3<-sum((y-m1)^3)/length(y)

      mu1<-1
      mu2<-2
      mu3<-6
      diffskew<-matrix(NA,99,2)
      for (al.i in 1:99){
        al=al.i/100
        nu<-(1-2*al)*((1-2*al)^2*(mu3-3*mu1*mu2+2*mu1^3)+al*(1-al)*(2*mu3-3*mu1*mu2))
        deno<-((1-2*al)^2*(mu2-mu1^2)+al*(1-al)*mu2)^(3/2)
        skewPop<-(nu/deno)
        diffskew[al.i,1]<-abs(skewPop- cm3/cm2^(3/2))
        diffskew[al.i,2]<-al}
      alpha<-diffskew[diffskew[,1]== min( diffskew[,1]), ][2]
      k1<-(1-2*alpha)*mu1/(alpha*(1-alpha))
      k2<-((1-alpha)^3+alpha^3)*mu2/(alpha^2*(1-alpha)^2)
      mu<-m1-(k1/sqrt(k2-k1^2))*sqrt(m2-m1^2)


    } else {
      m1<-mean(y)
      m2<-sum(y^2)/length(y)
      mu1<-1
      mu2<-2
      k1<-(1-2*alpha)*mu1/(alpha*(1-alpha))
      k2<-((1-alpha)^3+alpha^3)*mu2/(alpha^2*(1-alpha)^2)
      mu<-m1-(k1/sqrt(k2-k1^2))*sqrt(m2-m1^2)
    }
    return(mu)
  }
  phiMoM<-function(y,alpha=NULL){
    if (is.null(alpha)){
      m1<-mean(y)
      m2<-sum(y^2)/length(y)
      cm2<-sum((y-m1)^2)/length(y)
      cm3<-sum((y-m1)^3)/length(y)
      mu1<-1
      mu2<-2
      mu3<-6

      diffskew<-matrix(NA,99,2)
      for (al.i in 1:99){
        al=al.i/100
        nu<-(1-2*al)*((1-2*al)^2*(mu3-3*mu1*mu2+2*mu1^3)+al*(1-al)*(2*mu3-3*mu1*mu2))
        deno<-((1-2*al)^2*(mu2-mu1^2)+al*(1-al)*mu2)^(3/2)
        skewPop<-(nu/deno)
        diffskew[al.i,1]<-abs(skewPop- cm3/cm2^(3/2))
        diffskew[al.i,2]<-al}
      alpha<-diffskew[diffskew[,1]== min( diffskew[,1]), ][2]
      k1<-(1-2*alpha)*mu1/(alpha*(1-alpha))
      k2<-((1-alpha)^3+alpha^3)*mu2/(alpha^2*(1-alpha)^2)
      phi<-(1/sqrt(k2-k1^2))*sqrt(m2-m1^2)


    } else {
      m1<-mean(y)
      m2<-sum(y^2)/length(y)
      mu1<-1
      mu2<-2
      k1<-(1-2*alpha)*mu1/(alpha*(1-alpha))
      k2<-((1-alpha)^3+alpha^3)*mu2/(alpha^2*(1-alpha)^2)
      phi<-(1/sqrt(k2-k1^2))*sqrt(m2-m1^2)
    }
    return(phi)
  }

  alphaMoM<-function(y){
    m1<-mean(y)
    cm2<-sum((y-m1)^2)/length(y)
    cm3<-sum((y-m1)^3)/length(y)
    mu1<-1
    mu2<-2
    mu3<-6
    diffskew<-matrix(NA,99,2)
    for (al.i in 1:99){
      al=al.i/100
      nu<-(1-2*al)*((1-2*al)^2*(mu3-3*mu1*mu2+2*mu1^3)+al*(1-al)*(2*mu3-3*mu1*mu2))
      deno<-((1-2*al)^2*(mu2-mu1^2)+al*(1-al)*mu2)^(3/2)
      skewPop<-(nu/deno)
      diffskew[al.i,1]<-abs(skewPop- cm3/cm2^(3/2))
      diffskew[al.i,2]<-al
    }
    return(diffskew[diffskew[,1]== min( diffskew[,1]), ][2])
  }

  if (is.null(alpha)){
    list(mu.MoM=muMoM(y),phi.MoM=phiMoM(y),alpha.MoM=alphaMoM(y))
  } else {
    list(mu.MoM=muMoM(y),phi.MoM=phiMoM(y))}
}






#' @title Log-likelihood function for the quantile-based asymmetric Laplace distribution.
#' @description Log-Likelihood function \eqn{\ell_n(\mu,\phi,\alpha)=\ln[L_n(\mu,\phi,\alpha)]}
#' of the quantile-based asymmetric Laplace distribution discussed in Gijbels et al. (2019a).
#' @param y This is a vector of quantiles.
#' @param mu This is the location parameter \eqn{\mu}.
#' @param phi This is the scale parameter  \eqn{\phi}.
#' @param alpha This is the index parameter  \eqn{\alpha}.
#' @return \code{\link{LogLikALaD}} provides the value of the Log-likelihood function of the quantile-based asymmetric Laplace distribution.
#' @references{
#'  Gijbels, I., Karim, R. and Verhasselt, A. (2019a). On quantile-based asymmetric family of distributions: properties and inference. \emph{International Statistical Review}, \url{https://doi.org/10.1111/insr.12324}.
#' }
#'
#' @examples
#' # Example
#' y<-rnorm(100)
#' LogLikALaD(y,mu=0,phi=1,alpha=0.5)
#'
#' @export
LogLikALaD<- function(y,mu,phi,alpha){
  LL<-log(dALaD(y,mu,phi,alpha))
  return(sum(LL[!is.infinite(LL)]))
}



#' @title Maximum likelihood estimation (MLE) for the quantile-based asymmetric Laplace distribution.
#' @description The log-likelihood function \eqn{\ell_n(\mu,\phi,\alpha)=\ln[L_n(\mu,\phi,\alpha)]}
#' and parameter estimation of \eqn{ \theta=(\mu,\phi,\alpha)} in the quantile-based asymmetric Laplace distribution
#' by using the maximum likelihood estimation are discussed in Gijbels et al. (2019a).
#'  See also in Yu and Zhang (2005). The linear programing (LP) algorithm is used to obtain a solution to the maximization problem.
#'  The LP algorithm can be found in  Koenker (2005). See also \code{\link[ald:mleALD]{mleALD}} in the Package \pkg{ald}.
#' @param y This is a vector of quantiles.
#' @return The maximum likelihood estimate of parameter \eqn{ \theta=(\mu,\phi,\alpha)} of the quantile-based asymmetric Laplace distribution.
#' @import zipfR
#' @import ald
#'
#' @references{
#'  Gijbels, I., Karim, R. and Verhasselt, A. (2019a). On quantile-based asymmetric family of distributions: properties and inference. \emph{International Statistical Review}, \url{https://doi.org/10.1111/insr.12324}.
#'
#'
#'  Koenker, R. (2005). \emph{Quantile Regression}. Cambridge University Press.
#'
#'
#'  Yu., K, and Zhang, J. (2005). A three-parameter asymmetric Laplace distribution and its extension. \emph{Communications in Statistics–Theory and Methods}, \bold{34}(9-10), 1867–1879.
#'
#' }
#'
#' @name mleALaD
NULL
#' @rdname mleALaD
#' @examples
#' ## Example:
#' y=rnorm(100)
#' mleALaD(y)
#'
#' @export
mleALaD<-function(y){
  y<-sort(y)
  fit.ALD<-ald::mleALD(y, initial = NA)
  LogLik.ALaD<-sum(log(dALaD(y, mu = fit.ALD$par[1], phi = fit.ALD$par[2], alpha =fit.ALD$par[3])))
  return(list(LogLikelihood.ALaD=round(LogLik.ALaD,4),mu.ALaD=round(fit.ALD$par[1],4),phi.ALaD=round(fit.ALD$par[2],4),alpha.ALaD=round(fit.ALD$par[3],4)))
}
