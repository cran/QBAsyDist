#' @title Quantile-based asymmetric logistic distribution
#'
#' @description Density, cumulative distribution function, quantile function and random sample generation
#' from the quantile-based asymmetric logistic distribution (ALoD) proposed in Gijbels et al. (2019a).
#' @param y,q These are each a vector of quantiles.
#' @param mu This is the location parameter \eqn{\mu}.
#' @param phi This is the scale parameter  \eqn{\phi}.
#' @param alpha This is the index parameter  \eqn{\alpha}.
#' @param n This is the number of observations, which must be a positive integer that has length 1.
#' @param beta This is a vector of probabilities.
#'
#' @return \code{\link{dALoD}} provides the density, \code{\link{pALoD}} provides the cumulative distribution function, \code{\link{qALoD}} provides the quantile function, and \code{\link{rALoD}} generates a random sample from the quantile-based asymmetric logistic distribution.
#' The length of the result is determined by \eqn{n} for \code{\link{rALoD}}, and is the maximum of the lengths of the numerical arguments for the other functions.
#'
#' @references{
#'  Gijbels, I., Karim, R. and Verhasselt, A. (2019a). On quantile-based asymmetric family of distributions: properties and inference. \emph{International Statistical Review}, to appear.
#' }
#'
#' @name ALoD
NULL

#' @seealso \code{\link{dQBAD}},  \code{\link{pQBAD}},  \code{\link{qQBAD}},  \code{\link{rQBAD}}
#' @rdname ALoD
#' @export
dALoD<-function(y,mu,phi,alpha){2*alpha*(1-alpha)/phi*ifelse(y>mu,exp(-alpha*(y-mu)/phi)/(1+exp(-alpha*(y-mu)/phi))^2,exp(-(1-alpha)*(mu-y)/phi)/(1+exp(-(1-alpha)*(mu-y)/phi))^2)}



#' @rdname ALoD
#'
#' @export
pALoD<-function(q,mu,phi,alpha){ifelse(q>mu,2*alpha-1+2*(1-alpha)/(1+exp(-alpha*(q-mu)/phi)),
                                       (2*alpha)/(1+exp(-(1-alpha)*(q-mu)/phi)))}

#' @rdname ALoD
#' @export
qALoD<-function(beta,mu,phi,alpha){
  qf<-NA;
  for (i in 1:length(beta)) {
    qf[i]<-ifelse(beta[i]<alpha,mu-(phi/(1-alpha))*log(2*alpha/beta[i]-1),
                  mu-(phi/alpha)*log((1-beta[i])/(beta[i]-2*alpha+1)))}
  return(qf)
}




#' @rdname ALoD
#' @examples
#' # Quantile-based asymmetric logistic distribution (ALoD)
#' # Density
#' rnum<-rnorm(100)
#' dALoD(y=rnum,mu=0,phi=1,alpha=.5)
#'
#' # Distribution function
#' pALoD(q=rnum,mu=0,phi=1,alpha=.5)
#'
#' # Quantile function
#' beta<-c(0.25,0.5,0.75)
#' qALoD(beta=beta,mu=0,phi=1,alpha=.5)
#'
#' # random sample generation
#' rALoD(n=100,mu=0,phi=1,alpha=.5)
#'
#' @export
rALoD<-function(n,mu,phi,alpha){
  u<-runif(n, min = 0, max = 1)
  y<-NA;
  for (i in 1:length(u)) {
    y[i]<-ifelse(u[i]>alpha,mu-(phi/alpha)*log((1-u[i])/(u[i]-2*alpha+1)),
                 mu-(phi/(1-alpha))*log(2*alpha/u[i]-1))}
  return(y)
}


#' @title Moments estimation for the quantile-based asymmetric logistic distribution.
#' @description Mean, variance, skewness, kurtosis and moments about the location parameter (i.e., \eqn{\alpha}th quantile) of the quantile-based asymmetric logistic distribution defined in Gijbels et al. (2019a) useful for quantile regression with location parameter equal to \eqn{\mu}, scale parameter \eqn{\phi} and index parameter \eqn{\alpha}.
#' @param mu This is the location parameter \eqn{\mu}.
#' @param phi This is the scale parameter  \eqn{\phi}.
#' @param alpha This is the index parameter  \eqn{\alpha}.
#' @return \code{\link{meanALoD}} provides the mean, \code{\link{varALoD}} provides the variance, \code{\link{skewALoD}} provides the skewness, \code{\link{kurtALoD}} provides the kurtosis, and  \code{\link{momentALoD}} provides the \eqn{r}th moment about the location parameter \eqn{\mu} of the quantile-based asymmetric logistic distribution.
#'
#'
#' @references{
#'  Gijbels, I., Karim, R. and Verhasselt, A. (2019a). On quantile-based asymmetric family of distributions: properties and inference. \emph{International Statistical Review}, to appear.
#' }
#'
#'
#' @name momentALoD
NULL
#' @rdname momentALoD
#' @export
meanALoD<- function(mu,phi,alpha){mu+phi*(1-2*alpha)/(alpha*(1-alpha))*2*log(2)}

#' @rdname momentALoD
#' @export
varALoD<- function(mu,phi,alpha){(phi^2/(alpha^2*(1-alpha)^2))*((1-2*alpha)^2*(pi^2/3-(2*log(2))^2)+alpha*(1-alpha)*pi^2/3)}


#' @rdname momentALoD
#' @export
skewALoD<- function(alpha){
  f_Lo<-function(s){exp(-s)/(1+exp(-s))^2} # density function of LD(0,1)
  mu_1=mu_k(f=f_Lo,k=1)
  mu_2=mu_k(f=f_Lo,k=2)
  mu_3=mu_k(f=f_Lo,k=3)
  skew<-skewQBAD(alpha,mu_1,mu_2,mu_3)
  return(skew)
}



#' @rdname momentALoD
#'
#' @export
kurtALoD<- function(alpha){
  f_Lo<-function(s){exp(-s)/(1+exp(-s))^2} # density function of LD(0,1)
  mu_1=mu_k(f=f_Lo,k=1)
  mu_2=mu_k(f=f_Lo,k=2)
  mu_3=mu_k(f=f_Lo,k=3)
  mu_4=mu_k(f=f_Lo,k=4)
  kurto<-kurtQBAD(alpha,mu_1,mu_2,mu_3,mu_4)
  return(kurto)
}

#' @param r This is a value which is used to calculate the \eqn{r}th moment about \eqn{\mu}.
#' @import stats
#' @examples
#' # Example
#' meanALoD(mu=0,phi=1,alpha=0.5)
#' varALoD(mu=0,phi=1,alpha=0.5)
#' skewALoD(alpha=0.5)
#' kurtALoD(alpha=0.5)
#' momentALoD(phi=1,alpha=0.5,r=1)
#'
#'
#' @rdname momentALoD
#' @export
momentALoD<-function(phi,alpha,r){
  f_Lo<-function(s){exp(-s)/(1+exp(-s))^2} # density function of LD(0,1)
  return(momentQBAD(phi,alpha,f=f_Lo,r))
}



#' @title Method of moments (MoM) estimation for the quantile-based asymmetric logistic distribution.
#' @description Parameter estimation in the quantile-based asymmetric logistic distribution by using method of moments are studied in Gijbels et al. (2019a).
#' @param y This is a vector of quantiles.
#' @param alpha This is the index parameter  \eqn{\alpha}.  If \eqn{\alpha} is unknown, indicate NULL which is the default option. In this case, the sample skewness will be used to estimate \eqn{\alpha}. If \eqn{\alpha} is known, then the value of \eqn{\alpha} has to be specified in the function.
#'
#'
#'
#' @import stats
#' @return \code{\link{momALoD}} provides the method of moments estimates of the unknown parameters of the distribution.
#'
#'
#' @references{
#'  Gijbels, I., Karim, R. and Verhasselt, A. (2019a). On quantile-based asymmetric family of distributions: properties and inference. \emph{International Statistical Review}, to appear.
#' }
#'
#'
#' @name momALoD
NULL
#' @rdname momALoD
#' @examples
#' # Example
#' y=rnorm(100)
#' momALoD(y=y,alpha=0.5) # If alpha is known with alpha=0.5
#' momALoD(y=y) # If alpha is unknown
#'
#' @export
momALoD<-function(y,alpha=NULL){
  muMoM<-function(y,f,alpha=NULL){
    if (is.null(alpha)){
      m1<-mean(y)
      m2<-sum(y^2)/length(y)
      cm2<-sum((y-m1)^2)/length(y)
      cm3<-sum((y-m1)^3)/length(y)
      f_Lo<-function(s){exp(-s)/(1+exp(-s))^2} # density function of LD(0,1)
      mu1=mu_k(f=f_Lo,k=1)
      mu2=mu_k(f=f_Lo,k=2)
      mu3=mu_k(f=f_Lo,k=3)
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
      f_Lo<-function(s){exp(-s)/(1+exp(-s))^2} # density function of LD(0,1)
      mu1=mu_k(f=f_Lo,k=1)
      mu2=mu_k(f=f_Lo,k=2)

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
      f_Lo<-function(s){exp(-s)/(1+exp(-s))^2} # density function of LD(0,1)
      mu1=mu_k(f=f_Lo,k=1)
      mu2=mu_k(f=f_Lo,k=2)
      mu3=mu_k(f=f_Lo,k=3)
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
      f_Lo<-function(s){exp(-s)/(1+exp(-s))^2} # density function of LD(0,1)
      mu1=mu_k(f=f_Lo,k=1)
      mu2=mu_k(f=f_Lo,k=2)
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
    f_Lo<-function(s){exp(-s)/(1+exp(-s))^2} # density function of LD(0,1)
    mu1=mu_k(f=f_Lo,k=1)
    mu2=mu_k(f=f_Lo,k=2)
    mu3=mu_k(f=f_Lo,k=3)
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






#' @title Log-likelihood function for the quantile-based asymmetric logistic distribution.
#' @description The log-likelihood function \eqn{\ell_n(\mu,\phi,\alpha)=\ln[L_n(\mu,\phi,\alpha)]}
#' in the quantile-based asymmetric logistic distribution is presented in Gijbels et al. (2019a).
#'
#' @param y This is a vector of quantiles.
#' @param mu This is the location parameter \eqn{\mu}.
#' @param phi This is the scale parameter  \eqn{\phi}.
#' @param alpha This is the index parameter  \eqn{\alpha}.
#'
#' @return \code{\link{LogLikALoD}} provides the value of the Log-likelihood function of the quantile-based asymmetric logistic distribution.
#'
#'
#' @references{
#'  Gijbels, I., Karim, R. and Verhasselt, A. (2019a). On quantile-based asymmetric family of distributions: properties and inference. \emph{International Statistical Review}, to appear.
#' }
#'
#' @examples
#' # Example
#' y<-rnorm(100)
#' LogLikALoD(y,mu=0,phi=1,alpha=0.5)
#'
#' @export
LogLikALoD<- function(y,mu,phi,alpha){
  LL<-log(dALoD(y,mu,phi,alpha))
  return(sum(LL[!is.infinite(LL)]))
}


#' @title Maximum likelihood estimation (MLE) for the quantile-based asymmetric logistic distribution.
#' @description The log-likelihood function \eqn{\ell_n(\mu,\phi,\alpha)=\ln[L_n(\mu,\phi,\alpha)]}
#' and parameter estimation of \eqn{ \theta=(\mu,\phi,\alpha)} in the quantile-based asymmetric logistic distribution
#' by using the maximum likelihood estimation are discussed in Gijbels et al. (2019a).
#'
#' @param y This is a vector of quantiles.
#' @return The maximum likelihood estimate of parameter \eqn{\theta=(\mu,\phi,\alpha)} of the quantile-based asymmetric family of distributions.
#'
#'
#' @references{
#'  Gijbels, I., Karim, R. and Verhasselt, A. (2019a). On quantile-based asymmetric family of distributions: properties and inference. \emph{International Statistical Review}, to appear.
#' }
#'
#'
#' @name mleALoD
NULL
#' @rdname mleALoD
#' @examples
#' # Example
#' rnum=rnorm(100)
#' mleALoD(rnum)
#'
#' @export
mleALoD<-function(y){
  f_Lo<-function(s){exp(-s)/(1+exp(-s))^2} # density function of LD(0,1)
  return(mleQBAD(y,f=f_Lo))
}
