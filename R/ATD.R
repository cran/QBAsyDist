#' @title Quantile-based asymmetric Student's-t distribution
#'
#' @description Density, cumulative distribution function, quantile function and random sample generation
#' from the quantile-based asymmetric Student's-\eqn{t} distribution (ATD) proposed in Gijbels et al. (2019a).
#' @param y,q These are each a vector of quantiles.
#' @param mu This is the location parameter \eqn{\mu}.
#' @param phi This is the scale parameter  \eqn{\phi}.
#' @param alpha This is the index parameter  \eqn{\alpha}.
#' @param nu This is the degrees of freedom parameter \eqn{\nu}, which must be positive.
#' @param n This is the number of observations, which must be a positive integer that has length 1.
#' @param beta This is a vector of probabilities.
#'
#' @return \code{\link{dATD}} provides the density, \code{\link{pATD}} provides the cumulative distribution function, \code{\link{qATD}} provides the quantile function, and \code{\link{rATD}} generates a random sample from the quantile-based asymmetric Student's-\eqn{t} distribution.
#' The length of the result is determined by \eqn{n} for \code{\link{rATD}}, and is the maximum of the lengths of the numerical arguments for the other functions.
#'
#'
#'
#' @references{
#'  Gijbels, I., Karim, R. and Verhasselt, A. (2019a). On quantile-based asymmetric family of distributions: properties and inference. \emph{International Statistical Review}, \url{https://doi.org/10.1111/insr.12324}.
#' }
#'
#'
#'
#' @name ATD
NULL

#' @seealso \code{\link{dQBAD}},  \code{\link{pQBAD}},  \code{\link{qQBAD}},  \code{\link{rQBAD}}
#' @rdname ATD
#' @export
dATD<-function(y,mu,phi,alpha,nu){(2*alpha*(1-alpha)/(sqrt(nu)*phi*beta(.5,nu/2)))*ifelse(y>mu,(1+(alpha*(y-mu)/phi)^2/nu)^(-(nu+1)/2),(1+((1-alpha)*(y-mu)/phi)^2/nu)^(-(nu+1)/2))}



#' @import zipfR
#' @rdname ATD
#'
#' @export
pATD<-function(q,mu,phi,alpha,nu){ifelse(q>mu,
                                         alpha+(1-alpha)*zipfR::Rbeta(((alpha^2/nu)*((q-mu)/phi)^2)/(1+((alpha^2/nu)*((q-mu)/phi)^2)), 0.5, nu/2, lower=TRUE),
                                         alpha-alpha*zipfR::Rbeta((((1-alpha)^2/nu)*((q-mu)/phi)^2)/(1+(((1-alpha)^2/nu)*((q-mu)/phi)^2)), 0.5, nu/2, lower=TRUE))}

#' @import zipfR
#' @rdname ATD
#' @export
qATD<-function(beta,mu,phi,alpha,nu){
  qf<-NA;
  for (i in 1:length(beta)) {
    qf[i]<-ifelse(beta[i]<alpha,mu-(phi/(1-alpha))*sqrt(nu*zipfR::Rbeta.inv((alpha-beta[i])/alpha, 0.5, nu/2, lower=TRUE)/(1-zipfR::Rbeta.inv((alpha-beta[i])/alpha, 0.5, nu/2, lower=TRUE))),
                  mu+(phi/alpha)*sqrt(nu*zipfR::Rbeta.inv((beta[i]-alpha)/(1-alpha), 0.5, nu/2, lower=TRUE)/(1-zipfR::Rbeta.inv((beta[i]-alpha)/(1-alpha), 0.5, nu/2, lower=TRUE))))}
  return(qf)
}



#' @import zipfR
#' @rdname ATD
#' @examples
#' # Quantile-based asymmetric Student's-\eqn{t} distribution (ATD)
#' # Density
#' rnum<-rnorm(100)
#' dATD(rnum,mu=0,phi=1,alpha=0.5,nu=10)
#'
#' # Distribution function
#' pATD(rnum,mu=0,phi=1,alpha=0.5,nu=10)
#'
#' # Quantile function
#' beta<-c(0.25,0.5,0.75)
#' qATD(beta=beta,mu=0,phi=1,alpha=.5,nu=10)
#'
#' # random sample generation
#' rATD(n=100,mu=0,phi=1,alpha=.5,nu=10)
#'
#' @export
rATD<-function(n,mu,phi,alpha,nu){
  u<-runif(n, min = 0, max = 1)
  y<-NA;
  for (i in 1:length(u)) {
    y[i]<-ifelse(u[i]<alpha,mu-(phi/(1-alpha))*sqrt(nu*zipfR::Rbeta.inv((alpha-u[i])/alpha, 0.5, nu/2, lower=TRUE)/(1-zipfR::Rbeta.inv((alpha-u[i])/alpha, 0.5, nu/2, lower=TRUE))),
                 mu+(phi/alpha)*sqrt(nu*zipfR::Rbeta.inv((u[i]-alpha)/(1-alpha), 0.5, nu/2, lower=TRUE)/(1-zipfR::Rbeta.inv((u[i]-alpha)/(1-alpha), 0.5, nu/2, lower=TRUE))))}
  return(y)
}


#' @title Moments estimation for the quantile-based asymmetric Student's-\eqn{t} distribution.
#' @description Mean, variance, skewness, kurtosis and moments about the location parameter (i.e., \eqn{\alpha}th quantile) of the quantile-based asymmetric Student's-\eqn{t} distribution defined in Gijbels et al. (2019a) useful for quantile regression with location parameter equal to \eqn{\mu}, scale parameter \eqn{\phi} and index parameter \eqn{\alpha}.
#' @param mu This is the location parameter \eqn{\mu}.
#' @param phi This is the scale parameter  \eqn{\phi}.
#' @param alpha This is the index parameter  \eqn{\alpha}.
#' @param nu This is the degrees of freedom parameter \eqn{\nu}, which must be positive.
#'
#' @return \code{\link{meanATD}} provides the mean, \code{\link{varATD}} provides the variance, \code{\link{skewATD}} provides the skewness, \code{\link{kurtATD}} provides the kurtosis, and  \code{\link{momentATD}} provides the \eqn{r}th moment about the location parameter \eqn{\mu} of the quantile-based asymmetric Student's-\eqn{t} distribution.
#'
#'
#'
#'
#' @references{
#'  Gijbels, I., Karim, R. and Verhasselt, A. (2019a). On quantile-based asymmetric family of distributions: properties and inference. \emph{International Statistical Review}, \url{https://doi.org/10.1111/insr.12324}.
#' }
#'
#'
#' @name momentATD
NULL
#' @rdname momentATD
#' @export
meanATD<- function(mu,phi,alpha,nu){
  mu_1<-sqrt(nu)*gamma((nu-1)/2)/(sqrt(pi)*gamma(nu/2))
  return(meanQBAD(mu,phi,alpha,mu_1=mu_1))
}

#' @rdname momentATD
#' @export
varATD<- function(mu,phi,alpha,nu){
  mu_1<-sqrt(nu)*gamma((nu-1)/2)/(sqrt(pi)*gamma(nu/2))
  mu_2<-nu/(nu-2)
  return(varQBAD(mu,phi,alpha,mu_1=mu_1,mu_2=mu_2))
  }


#' @rdname momentATD
#' @export
skewATD<- function(alpha,nu){
  mu_1<-sqrt(nu)*gamma((nu-1)/2)/(sqrt(pi)*gamma(nu/2))
  mu_2<-nu/(nu-2)
  mu_3=nu*sqrt(nu)*gamma((nu-1)/2)/(sqrt(pi)*gamma(nu/2))
  skew<-skewQBAD(alpha,mu_1=mu_1,mu_2=mu_2,mu_3=mu_3)
  return(skew)
}



#' @rdname momentATD
#'
#' @export
kurtATD<- function(alpha,nu){
  mu_1<-sqrt(nu)*gamma((nu-1)/2)/(sqrt(pi)*gamma(nu/2))
  mu_2<-nu/(nu-2)
  mu_3=nu*sqrt(nu)*gamma((nu-1)/2)/(sqrt(pi)*gamma(nu/2))
  mu_4=3*nu^2/((nu-2)*(nu-4))
  kurto<-kurtQBAD(alpha,mu_1,mu_2,mu_3,mu_4)
  return(kurto)
}

#' @param r This is a value which is used to calculate the \eqn{r}th moment \eqn{(r\in\{1,2,3,4\})} about \eqn{\mu}.
#' @import stats
#' @examples
#' # Example
#' meanATD(mu=0,phi=1,alpha=0.5,nu=10)
#' varATD(mu=0,phi=1,alpha=0.5,nu=10)
#' skewATD(alpha=0.5,nu=10)
#' kurtATD(alpha=0.5,nu=10)
#' momentATD(phi=1,alpha=0.5,nu=10,r=1)
#'
#'
#' @rdname momentATD
#' @export
momentATD<-function(phi,alpha,nu,r){
  quan<-phi^2*((1-alpha)^(r+1)+(-1)^r*alpha^(r+1))/((1-alpha)^r*alpha^r)
  if (r == 1) {
    mu_1<-sqrt(nu)*gamma((nu-1)/2)/(sqrt(pi)*gamma(nu/2))
    return(quan*mu_1)
  }
  else if (r == 2) {
    mu_2<-nu/(nu-2)
    return(quan*mu_2)
  }
  else if (r == 3){
    mu_3=nu*sqrt(nu)*gamma((nu-1)/2)/(sqrt(pi)*gamma(nu/2))
    return(quan*mu_3)
  }
  else if (r == 4){
    mu_4=3*nu^2/((nu-2)*(nu-4))
    return(quan*mu_4)
  }
  else if (r >4){
     return("This function is useful for r={1,2,3,4}")
  }
}



#' @title Method of moments (MoM) estimation for the quantile-based asymmetric Student's-\eqn{t} distribution.
#' @description Parameter estimation in the quantile-based asymmetric Student's-\eqn{t} distribution by using method of moments are discussed in Gijbels et al. (2019a). We here used the first four sample moments to estimate parameter \eqn{\theta=(\mu,\phi,\alpha,\nu)} under the assumption that the first four population moments exist, which needs to assume \eqn{\nu>4}.
#' @param y This is a vector of quantiles.
#' @param alpha This is the index parameter  \eqn{\alpha}.  If \eqn{\alpha} is unknown, indicate NULL which is the default option. In this case, the sample skewness will be used to estimate \eqn{\alpha}. If \eqn{\alpha} is known, then the value of \eqn{\alpha} has to be specified in the function.
#'
#'
#'
#' @import stats
#' @return \code{\link{momATD}} provides the method of moments estimates of the unknown parameters of the distribution.
#'
#'
#'
#'
#' @references{
#'  Gijbels, I., Karim, R. and Verhasselt, A. (2019a). On quantile-based asymmetric family of distributions: properties and inference. \emph{International Statistical Review}, \url{https://doi.org/10.1111/insr.12324}.
#' }
#'
#'
#'
#' @name momATD
NULL
#' @rdname momATD
#' @examples
#' # Example
#' y=rnorm(100)
#' momATD(y=y,alpha=0.5) # If alpha is known with alpha=0.5
#' momATD(y=y) # If alpha is unknown
#'
#' @export
momATD<-function(y,alpha=NULL){
  muMoM<-function(y,alpha,nu){
      m1<-mean(y)
      m2<-sum(y^2)/length(y)
      mu1<-sqrt(nu)*gamma((nu-1)/2)/(sqrt(pi)*gamma(nu/2))
      mu2<-nu/(nu-2)
      k1<-(1-2*alpha)*mu1/(alpha*(1-alpha))
      k2<-((1-alpha)^3+alpha^3)*mu2/(alpha^2*(1-alpha)^2)
      mu<-m1-(k1/sqrt(k2-k1^2))*sqrt(m2-m1^2)
    return(mu)
  }
  phiMoM<-function(y,alpha,nu){
      m1<-mean(y)
      m2<-sum(y^2)/length(y)
      mu1<-sqrt(nu)*gamma((nu-1)/2)/(sqrt(pi)*gamma(nu/2))
      mu2<-nu/(nu-2)
      k1<-(1-2*alpha)*mu1/(alpha*(1-alpha))
      k2<-((1-alpha)^3+alpha^3)*mu2/(alpha^2*(1-alpha)^2)
      phi<-(1/sqrt(k2-k1^2))*sqrt(m2-m1^2)
    return(phi)
  }

  if (is.null(alpha)){
    alpha_nu<-function(theta){
      m1<-mean(y)
      cm2<-sum((y-m1)^2)/length(y)
      cm3<-sum((y-m1)^3)/length(y)
      cm4<-sum((y-m1)^4)/length(y)
      sample.skewness<-cm3/cm2^(3/2)
      sample.kurtosis<-cm4/cm2^2
      sum((sample.skewness-skewATD(theta[1],theta[2]))^2+(sample.kurtosis-kurtATD(theta[1],theta[2]))^2)
      }
    theta0<-ald::mleALD(y, initial = NA)$par
    alpha_nu<-optim(c(theta0[3],5), alpha_nu)$par
    list(mu.MoM=muMoM(y,alpha=alpha_nu[1],nu=alpha_nu[2]),phi.MoM=phiMoM(y,alpha=alpha_nu[1],nu=alpha_nu[2]),alpha.MoM=alpha_nu[1],nu.MoM=alpha_nu[2])
  } else {
      nu_fn<-function(theta,alpha){
        m1<-mean(y)
        cm2<-sum((y-m1)^2)/length(y)
        cm3<-sum((y-m1)^3)/length(y)
        cm4<-sum((y-m1)^4)/length(y)
        sample.skewness<-cm3/cm2^(3/2)
        sample.kurtosis<-cm4/cm2^2
        sum((sample.skewness-skewATD(alpha,theta))^2+(sample.kurtosis-kurtATD(alpha,theta))^2)
      }
      nu.est<-optimize(nu_fn,c(4.1,500),alpha=alpha)$minimum
    list(mu.MoM=muMoM(y,alpha=alpha,nu=nu.est),phi.MoM=phiMoM(y,alpha=alpha,nu=nu.est),nu.MoM=nu.est)}
}







#' @title Log-likelihood function for the quantile-based asymmetric Student's-\eqn{t} distribution.
#' @description The log-likelihood function \eqn{\ell_n(\mu,\phi,\alpha,\nu)=\ln[L_n(\mu,\phi,\alpha,\nu)]}
#' and parameter estimation of \eqn{ \theta=(\mu,\phi,\alpha,\nu)} in
#' the quantile-based asymmetric Student's-\eqn{t} distribution by using the maximum likelihood estimation
#' are discussed in Gijbels et al. (2019a).
#'
#' @param y This is a vector of quantiles.
#' @param mu This is the location parameter \eqn{\mu}.
#' @param phi This is the scale parameter  \eqn{\phi}.
#' @param alpha This is the index parameter  \eqn{\alpha}.
#' @param nu This is the degrees of freedom parameter \eqn{\nu}, which must be positive.
#'
#' @return \code{\link{LogLikATD}} provides the value of the Log-likelihood function of the quantile-based asymmetric Student's-\eqn{t} distribution.
#'
#'
#' @references{
#'  Gijbels, I., Karim, R. and Verhasselt, A. (2019a). On quantile-based asymmetric family of distributions: properties and inference. \emph{International Statistical Review}, \url{https://doi.org/10.1111/insr.12324}.
#' }
#'
#'
#' @examples
#' # Example
#' y<-rnorm(100)
#' LogLikATD(y,mu=0,phi=1,alpha=0.5,nu=10)
#'
#' @export
LogLikATD<- function(y,mu,phi,alpha,nu){
  LL<-log(dATD(y,mu,phi,alpha,nu))
  return(sum(LL[!is.infinite(LL)]))
}





#' @title Maximum likelihood estimation (MLE) for the quantile-based asymmetric Student's-\eqn{t} distribution.
#' @description The log-likelihood function \eqn{\ell_n(\mu,\phi,\alpha,\nu)=\ln[L_n(\mu,\phi,\alpha,\nu)]}
#' and parameter estimation of \eqn{ \theta=(\mu,\phi,\alpha,\nu)} in the quantile-based asymmetric Student's-\eqn{t} distribution.
#' by using the maximum likelihood estimation are discussed in Gijbels et al. (2019a).
#' @param y This is a vector of quantiles.
#' @return The maximum likelihood estimate of parameter \eqn{\theta=(\mu,\phi,\alpha,\nu)} of the quantile-based asymmetric Student's-\eqn{t} distribution.
#'
#' @import ald
#' @import nloptr
#'
#'
#' @references{
#'  Gijbels, I., Karim, R. and Verhasselt, A. (2019a). On quantile-based asymmetric family of distributions: properties and inference. \emph{International Statistical Review}, \url{https://doi.org/10.1111/insr.12324}.
#' }
#'
#'
#'
#' @name mleATD
NULL
#' @rdname mleATD
#' @examples
#' \donttest{
#' # Example
#' y=rnorm(20)
#' mleATD(y)
#'
#'}
#' @export
mleATD<-function(y){
  theta0<-ald::mleALD(y, initial = NA)$par
  fn <- function(theta){
    mu=theta[1]
    phi=theta[2]
    alpha=theta[3]
    nu=theta[4]
    LL<-log(dATD(y,mu,phi,alpha,nu))
    return(-sum(LL[!is.infinite(LL)]))}
  Est.par=nloptr::bobyqa(x0=c(theta0,2),fn = fn ,lower = c(min(y),00.001,0.001,0.001), upper = c(max(y),Inf,1,Inf))$par
  estimated.LL<-LogLikATD(y=y,mu=Est.par[1],phi=Est.par[2],alpha=Est.par[3],nu=Est.par[4])
  list(mu.MLE=Est.par[1], phi.MLE=Est.par[2],alpha.MLE=Est.par[3],nu.MLE=Est.par[4],LogLikelihood=estimated.LL)
}



