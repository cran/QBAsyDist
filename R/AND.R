#' @title Quantile-based asymmetric normal distribution
#'
#' @description Density, cumulative distribution function, quantile function and random sample generation
#' from the quantile-based asymmetric normal distribution (AND) introduced in Gijbels et al. (2019a).
#' @param y,q These are each a vector of quantiles.
#' @param mu This is the location parameter \eqn{\mu}.
#' @param phi This is the scale parameter  \eqn{\phi}.
#' @param alpha This is the index parameter  \eqn{\alpha}.
#' @param n This is the number of observations, which must be a positive integer that has length 1.
#' @param beta This is a vector of probabilities.
#'
#' @return \code{\link{dAND}} provides the density, \code{\link{pAND}} provides the cumulative distribution function, \code{\link{qAND}} provides the quantile function, and \code{\link{rAND}} generates a random sample from the quantile-based asymmetric normal distribution.
#'
#' @references{
#'  Gijbels, I., Karim, R. and Verhasselt, A. (2019a). On quantile-based asymmetric family of distributions: properties and inference. \emph{International Statistical Review}, \url{https://doi.org/10.1111/insr.12324}.
#' }
#'
#'
#' @name AND
NULL

#' @seealso \code{\link{dQBAD}},  \code{\link{pQBAD}},  \code{\link{qQBAD}},  \code{\link{rQBAD}}
#' @rdname AND
#' @export

dAND<-function(y,mu,phi,alpha){2*alpha*(1-alpha)*(1/phi)*ifelse(y>mu,dnorm(alpha*(y-mu)/phi,0,1),dnorm((1-alpha)*(y-mu)/phi,0,1))}



#' @import zipfR
#' @rdname AND
#'
#' @export
pAND<-function(q,mu,phi,alpha){ifelse(q>mu,alpha+((1-alpha)/sqrt(pi))*zipfR::Igamma(0.5,(alpha^2/2)*(q-mu)^2/phi^2,lower=TRUE),(alpha/sqrt(pi))*zipfR::Igamma(0.5,((1-alpha)^2/2)*(q-mu)^2/phi^2,lower=FALSE))}

#' @import zipfR
#' @rdname AND
#' @export
qAND<-function(beta,mu,phi,alpha){
  qf<-NA;
  for (i in 1:length(beta)) {
    qf[i]<-ifelse(beta[i]<alpha,mu-sqrt((2*phi^2/(1-alpha)^2)*zipfR::Igamma.inv(0.5, beta[i]*sqrt(pi)/alpha,lower=FALSE)),
                  mu+sqrt((2*phi^2/alpha^2)*zipfR::Igamma.inv(0.5, sqrt(pi)*(beta[i]-alpha)/(1-alpha))))}
  return(qf)
}




#' @import zipfR
#' @rdname AND
#' @examples
#' # Quantile-based asymmetric normal distribution (AND)
#' # Density
#' rnum<-rnorm(100)
#' dAND(y=rnum,mu=0,phi=1,alpha=.5)
#'
#' # Distribution function
#' pAND(q=rnum,mu=0,phi=1,alpha=.5)
#'
#' # Quantile function
#' beta<-c(0.25,0.5,0.75)
#' qAND(beta=beta,mu=0,phi=1,alpha=.5)
#'
#' # random sample generation
#' rAND(n=100,mu=0,phi=1,alpha=.5)
#'
#' @export
rAND<-function(n,mu,phi,alpha){
  u<-runif(n, min = 0, max = 1)
  y<-NA;
  for (i in 1:length(u)) {
    y[i]<-ifelse(u[i]>alpha,mu+sqrt((2*phi^2/alpha^2)*zipfR::Igamma.inv(0.5, sqrt(pi)*(u[i]-alpha)/(1-alpha))),
                 mu-sqrt((2*phi^2/(1-alpha)^2)*zipfR::Igamma.inv(0.5, u[i]*sqrt(pi)/alpha,lower=FALSE)))}
  return(y)
}



#' @title Moments estimation for the quantile-based asymmetric normal distribution.
#' @description Mean, variance, skewness, kurtosis and moments about the location parameter (i.e., \eqn{\alpha}th quantile) of the quantile-based asymmetric normal distribution introduced in Gijbels et al. (2019a) useful for quantile regression with location parameter equal to \eqn{\mu}, scale parameter \eqn{\phi} and index parameter \eqn{\alpha}.
#' @param mu This is the location parameter \eqn{\mu}.
#' @param phi This is the scale parameter  \eqn{\phi}.
#' @param alpha This is the index parameter  \eqn{\alpha}.
#' @return \code{\link{meanAND}} provides the mean, \code{\link{varAND}} provides the variance, \code{\link{skewAND}} provides the skewness, \code{\link{kurtAND}} provides the kurtosis, and  \code{\link{momentAND}} provides the \eqn{r}th moment about the location parameter \eqn{\mu} of the quantile-based asymmetric normal distribution.
#'
#'
#' @references{
#'  Gijbels, I., Karim, R. and Verhasselt, A. (2019a). On quantile-based asymmetric family of distributions: properties and inference. \emph{International Statistical Review}, \url{https://doi.org/10.1111/insr.12324}.
#' }
#'
#' @name momentAND
NULL
#' @rdname momentAND
#' @export
meanAND<- function(mu,phi,alpha){mu+phi*(1-2*alpha)/(alpha*(1-alpha))*sqrt(2/pi)}

#' @rdname momentAND
#' @export
varAND<- function(mu,phi,alpha){(phi^2/(alpha^2*(1-alpha)^2))*((1-2*alpha)^2*(1- (sqrt(2/pi))^2)+alpha*(1-alpha))}


#' @rdname momentAND
#' @export
skewAND<- function(alpha){
  nu<-(1-2*alpha)*((1-2*alpha)^2*((2*sqrt(2/pi))-3*(sqrt(2/pi))+2*(sqrt(2/pi))^3)+alpha*(1-alpha)*(2*(2*sqrt(2/pi))-3*(sqrt(2/pi))))
  deno<-((1-2*alpha)^2*(1-(sqrt(2/pi))^2)+alpha*(1-alpha))^(3/2)
  skew<-(nu/deno)
  return(skew)
}



#' @rdname momentAND
#'
#' @export
kurtAND<- function(alpha){
  nu<-((1-alpha)^5+alpha^5)*3-(1-2*alpha)^2*(4*(1-2*alpha+2*alpha^2)*(sqrt(2/pi))*(2*sqrt(2/pi))-6*(1-3*alpha+3*alpha^2)*(sqrt(2/pi))^2+3*(1-2*alpha)^2*(sqrt(2/pi))^4)
  deno<-((1-2*alpha)^2*(1-(sqrt(2/pi))^2)+alpha*(1-alpha))^(2)
  kurto<-(nu/deno)
  return(kurto)
}

#' @param r This is a value which is used to calculate \eqn{r}th moment about \eqn{\mu}.
#' @import stats
#' @examples
#' # Example
#' meanAND(mu=0,phi=1,alpha=0.5)
#' varAND(mu=0,phi=1,alpha=0.5)
#' skewAND(alpha=0.5)
#' kurtAND(alpha=0.5)
#' momentAND(phi=1,alpha=0.5,r=1)
#'
#'
#' @rdname momentAND
#' @export
momentAND<-function(phi,alpha,r){
  integrand <- function(x) {x^r*(exp(-x^2/2))/sqrt(2*pi)}
  mu.r<-2*stats::integrate(integrand, lower = 0, upper = Inf)$ value
  return(phi^r*((1-alpha)^(r+1)+(-1)^r*alpha^(r+1))*mu.r/(alpha^r*(1-alpha)^r))
}


#' @title Method of moments (MoM) estimation for the quantile-based asymmetric normal distribution.
#' @description Parameter estimation in the quantile-based asymmetric normal distribution by using method of moments are discussed in Gijbels et al. (2019a).
#' @param y This is a vector of quantiles.
#' @param alpha This is the index parameter  \eqn{\alpha}.  If \eqn{\alpha} is unknown, indicate NULL which is the default option. In this case, the sample skewness will be used to estimate \eqn{\alpha}. If \eqn{\alpha} is known, then the value of \eqn{\alpha} has to be specified in the function.
#'
#'
#'
#' @import stats
#' @return \code{\link{momAND}} provides the method of moments estimates of the unknown parameters of the distribution.
#'
#'
#' @references{
#'  Gijbels, I., Karim, R. and Verhasselt, A. (2019a). On quantile-based asymmetric family of distributions: properties and inference. \emph{International Statistical Review}, \url{https://doi.org/10.1111/insr.12324}.
#' }
#'
#'
#' @name momAND
NULL
#' @rdname momAND
#' @examples
#' # Example
#' y=rnorm(100)
#' momAND(y=y,alpha=0.5) # If alpha is known wtih alpha=0.5
#' momAND(y=y) # If alpha is unknown
#'
#' @export
momAND<-function(y,alpha=NULL){
  muMoM<-function(y,f,alpha=NULL){
    if (is.null(alpha)){
      m1<-mean(y)
      m2<-sum(y^2)/length(y)
      cm2<-sum((y-m1)^2)/length(y)
      cm3<-sum((y-m1)^3)/length(y)
      mu1<-sqrt(2/pi)
      mu2<-1
      mu3<-2*sqrt(2/pi)
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
      mu1<-sqrt(2/pi)
      mu2<-1
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
      mu1<-sqrt(2/pi)
      mu2<-1
      mu3<-2*sqrt(2/pi)

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
      mu1<-sqrt(2/pi)
      mu2<-1
      k1<-(1-2*alpha)*mu1/(alpha*(1-alpha))
      k2<-((1-alpha)^3+alpha^3)*mu2/(alpha^2*(1-alpha)^2)
      phi<-(1/sqrt(k2-k1^2))*sqrt(m2-m1^2)
    }
    return(phi)
  }

  alphaMoM<-function(y){
    m1<-mean(y)
    m2<-sum(y^2)/length(y)
    cm2<-sum((y-m1)^2)/length(y)
    cm3<-sum((y-m1)^3)/length(y)
    mu1<-sqrt(2/pi)
    mu2<-1
    mu3<-2*sqrt(2/pi)
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



#' @title Log-likelihood function for the quantile-based asymmetric normal distribution.
#' @description The log-likelihood function \eqn{\ell_n(\mu,\phi,\alpha)=\ln[L_n(\mu,\phi,\alpha)]} in the quantile-based asymmetric normal distribution
#' is presented in Gijbels et al. (2019a).
#' @param y This is a vector of quantiles.
#' @param mu This is the location parameter \eqn{\mu}.
#' @param phi This is the scale parameter  \eqn{\phi}.
#' @param alpha This is the index parameter  \eqn{\alpha}.
#' @return \code{\link{LogLikAND}} provides the value of the Log-likelihood function of the quantile-based asymmetric normal distribution.
#'
#' @references{
#'  Gijbels, I., Karim, R. and Verhasselt, A. (2019a). On quantile-based asymmetric family of distributions: properties and inference. \emph{International Statistical Review}, \url{https://doi.org/10.1111/insr.12324}.
#' }
#'
#' @examples
#' # Example
#' y<-rnorm(100)
#' LogLikAND(y,mu=0,phi=1,alpha=0.5)
#'
#' @export
LogLikAND<- function(y,mu,phi,alpha){
  LL<-log(dAND(y,mu,phi,alpha))
  return(sum(LL[!is.infinite(LL)]))
}





#' @title Maximum likelihood estimation (MLE) for the quantile-based asymmetric normal distribution.
#' @description The log-likelihood function \eqn{\ell_n(\mu,\phi,\alpha)=\ln[L_n(\mu,\phi,\alpha)]}
#' and parameter estimation of \eqn{ \theta=(\mu,\phi,\alpha)} in the asymmetric normal distribution
#' by using the maximum likelihood estimation are discussed in Gijbels et al. (2019a).
#' @param y This is a vector of quantiles.
#' @param alpha This is the index parameter  \eqn{\alpha}.
#' @return The maximum likelihood estimate of parameter \eqn{\theta=(\mu,\phi,\alpha)} of the quantile-based asymmetric normal distribution.
#' @import zipfR
#' @import ald
#' @references{
#'  Gijbels, I., Karim, R. and Verhasselt, A. (2019a). On quantile-based asymmetric family of distributions: properties and inference. \emph{International Statistical Review}, \url{https://doi.org/10.1111/insr.12324}.
#' }
#'
#'
#'
#' @name mleAND
NULL
#' @rdname mleAND
#' @examples
#' \donttest{
#' # Maximum likelihood estimation
#' y=rnorm(100)
#' mleAND(y)
#' mleAND(y,alpha=0.5)
#' }
#' @export
mleAND<-function(y,alpha=NULL){
  f_N<-function(s){dnorm(s, mean = 0,sd = 1)} # density function of N(0,1)
 return(mleQBAD(y,f=f_N,alpha))}


