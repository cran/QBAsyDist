#' @title Quantile-based asymmetric family of distributions
#'
#' @description Density, cumulative distribution function, quantile function and random sample generation
#' from the quantile-based asymmetric family of densities defined in Gijbels et al. (2019a).
#'
#' @param y,q These are each a vector of quantiles.
#' @param mu This is the location parameter \eqn{\mu}.
#' @param phi This is the scale parameter  \eqn{\phi}.
#' @param alpha This is the index parameter  \eqn{\alpha}.
#'
#' @param f This is the reference density function \eqn{f} which is a standard version of a unimodal and symmetric around 0 density.
#'
#'
#'
#' @references{
#'  Gijbels, I., Karim, R. and Verhasselt, A. (2019a). On quantile-based asymmetric family of distributions: properties and inference. \emph{International Statistical Review}, \url{https://doi.org/10.1111/insr.12324}.
#' }
#'
#'
#' @name QBAD
NULL
#' @rdname QBAD
#'
#'
#' @export
dQBAD<-function(y,mu,phi,alpha,f){ifelse(y>mu,((2*alpha*(1-alpha)/phi)*f(alpha*(y-mu)/phi)), ((2*alpha*(1-alpha)/phi)*f((1-alpha)*(mu-y)/phi)))}


#' @param F This is the cumulative distribution function \eqn{F} of a unimodal and symmetric around 0 reference density function \eqn{f}.
#'
#' @rdname QBAD
#' @export
pQBAD<-function(q,mu,phi,alpha,F){ifelse(q>mu,(2*alpha-1+2*(1-alpha)*F(alpha*(q-mu)/phi)), (2*alpha*F((1-alpha)*(q-mu)/phi)))}

#' @param beta This is a vector of probabilities.
#' @param QF This is the quantile function of the reference density \eqn{f}.
#' @import GoFKernel
#' @rdname QBAD
#' @export
qQBAD<-function(beta,mu,phi,alpha,F,QF=NULL){
  if (is.null(QF)){
    QF<-GoFKernel::inverse(F, lower = -Inf, upper = Inf)
    qf<-NA;
    for (i in 1:length(beta)) {
      qf[i]<-ifelse(beta[i]<alpha,(mu+(phi/(1-alpha))*QF(beta[i]/(2*alpha))),
                    mu+(phi/alpha)*QF((1+beta[i]-2*alpha)/(2*(1-alpha))))}
  } else {

    qf<-ifelse(beta<alpha,(mu+(phi/(1-alpha))*QF(beta/(2*alpha))),
               mu+(phi/alpha)*QF((1+beta-2*alpha)/(2*(1-alpha))))
  }
  return(qf)
}






#' @param n This is the number of observations, which must be a positive integer that has length 1.
#' @import stats
#' @return \code{\link{dQBAD}} provides the density, \code{\link{pQBAD}} provides the cumulative distribution function, \code{\link{qQBAD}} provides the quantile function, and \code{\link{rQBAD}} generates a random sample from the quantile-based asymmetric family of distributions.
#' The length of the result is determined by \eqn{n} for \code{\link{rQBAD}}, and is the maximum of the lengths of the numerical arguments for the other functions.
#'
#' @rdname QBAD
#'
#' @examples
#' # Example 1: Let F be a standard normal cumulative distribution function then
#' f_N<-function(s){dnorm(s, mean = 0,sd = 1)} # density function of N(0,1)
#' F_N<-function(s){pnorm(s, mean = 0,sd = 1)} # distribution function of N(0,1)
#' QF_N<-function(beta){qnorm(beta, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)}
#' rnum<-rnorm(100)
#' beta=c(0.25,0.50,0.75)
#'
#' # Density
#' dQBAD(y=rnum,mu=0,phi=1,alpha=.5,f=f_N)
#'
#' # Distribution function
#' pQBAD(q=rnum,mu=0,phi=1,alpha=.5,F=F_N)
#'
#' # Quantile function
#' qQBAD(beta=beta,mu=0,phi=1,alpha=.5,F=F_N,QF=QF_N)
#' qQBAD(beta=beta,mu=0,phi=1,alpha=.5,F=F_N)
#'
#' # random sample generation
#' rQBAD(n=100,mu=0,phi=1,alpha=.5,QF=QF_N)
#' rQBAD(n=100,mu=0,phi=1,alpha=.5,F=F_N)
#'
#'
#'
#' # Example 2: Let F be a standard Laplace cumulative distribution function then
#' f_La<-function(s){0.5*exp(-abs(s))} # density function of Laplace(0,1)
#' F_La<-function(s){0.5+0.5*sign(s)*(1-exp(-abs(s)))} # distribution function of Laplace(0,1)
#' QF_La<-function(beta){-sign(beta-0.5)*log(1-2*abs(beta-0.5))}
#' rnum<-rnorm(100)
#' beta=c(0.25,0.50,0.75)
#'
#' # Density
#' dQBAD(y=rnum,mu=0,phi=1,alpha=.5,f=f_La)
#'
#' # Distribution function
#' pQBAD(q=rnum,mu=0,phi=1,alpha=.5,F=F_La)
#'
#' # Quantile function
#' qQBAD(beta=c(0.25,0.50,0.75),mu=0,phi=1,alpha=.5,F=F_La,QF=QF_La)
#' qQBAD(beta=c(0.25,0.50,0.75),mu=0,phi=1,alpha=.5,F=F_La)
#'
#' # random sample generation
#' rQBAD(n=100,mu=0,phi=1,alpha=.5,QF=QF_La)
#' rQBAD(n=100,mu=0,phi=1,alpha=.5,F=F_La)
#'
#' @export
rQBAD<-function(n,mu,phi,alpha,F,QF=NULL){
  u<-runif(n, min = 0, max = 1)
  if (is.null(QF)){
    QF<-GoFKernel::inverse(F, lower = -Inf, upper = Inf)
    r<-NA;
    for (i in 1:length(u)) {
      r[i]<-ifelse(u[i]<alpha,(mu+(phi/(1-alpha))*QF(u[i]/(2*alpha))),
                   (mu+(phi/alpha)*QF((1+u[i]-2*alpha)/(2*(1-alpha)))))}
  } else {
    r<-ifelse(u<alpha,(mu+(phi/(1-alpha))*QF(u/(2*alpha))),
              (mu+(phi/alpha)*QF((1+u-2*alpha)/(2*(1-alpha)))))
  }
  return(r)
}

#' @title Log-likelihood function for the quantile-based asymmetric family of distributions.
#' @description Log-Likelihood function \eqn{\ell_n(\mu,\phi,\alpha)=\ln[L_n(\mu,\phi,\alpha)]}
#' in the three parameter quantile-based asymmetric family of densities defined in Section 3.2 of Gijbels et al. (2019a).
#' @param y This is a vector of quantiles.
#' @param mu This is the location parameter \eqn{\mu}.
#' @param phi This is the scale parameter  \eqn{\phi}.
#' @param alpha This is the index parameter  \eqn{\alpha}.
#'
#' @param f This is the reference density function \eqn{f} which is a standard version of a unimodal and symmetric around 0 density.
#' @return \code{\link{LogLikQBAD}} provides the realized value of the Log-likelihood function of quantile-based asymmetric family of distributions.
#'
#'
#' @references{
#'  Gijbels, I., Karim, R. and Verhasselt, A. (2019a). On quantile-based asymmetric family of distributions: properties and inference. \emph{International Statistical Review}, \url{https://doi.org/10.1111/insr.12324}.
#' }
#'
#'
#' @examples
#' # Example 1: Let F be a standard normal cumulative distribution function then
#' f_N<-function(s){dnorm(s, mean = 0,sd = 1)} # density function of N(0,1)
#' y<-rnorm(100)
#' LogLikQBAD(y,mu=0,phi=1,alpha=0.5,f=f_N)
#'
#' # Example 2: Let F be a standard Laplace cumulative distribution function then
#' f_La<-function(s){0.5*exp(-abs(s))} # density function of Laplace(0,1)
#' LogLikQBAD(y,mu=0,phi=1,alpha=0.5,f=f_La)
#' @export
LogLikQBAD<- function(y,mu,phi,alpha,f){
  LL<-(log(ifelse(y>mu,((2*alpha*(1-alpha)/phi)*f(alpha*(y-mu)/phi)), ((2*alpha*(1-alpha)/phi)*f((1-alpha)*(mu-y)/phi)))))
  return(sum(LL[!is.infinite(LL)]))
}

#' @title Moment estimation for the quantile-based asymmetric family of distributions.
#' @description Mean, variance, skewness, kurtosis and moments about the location parameter (i.e., \eqn{\alpha}th quantile) of the quantile-based asymmetric family of densities defined in Gijbels et al. (2019a) useful for quantile regression with location parameter equal to \eqn{\mu}, scale parameter \eqn{\phi} and index parameter \eqn{\alpha}.
#' @param f This is the reference density function \eqn{f} which is a standard version of a unimodal and symmetric around 0 density.
#' @param r This is a value which is used to calculate the \eqn{r}th moment about \eqn{\mu}.
#' @param k This is an integer value (\eqn{k=1,2,3,\ldots}) for calculating \eqn{\mu_k=\int_{0}^{\infty} 2s^k f(s) ds} and \eqn{\gamma_k=\int_{0}^{\infty}s^{k-1}\frac{[f'(s)]^2}{f(s)}ds}.
#' @import stats
#' @return \code{\link{mu_k}} provides the quantity \eqn{\int_{0}^{\infty} 2s^k f(s) ds}, \code{\link{gamma_k}} provides the quantity \eqn{\int_{0}^{\infty}s^{k-1}\frac{[f'(s)]^2}{f(s)}ds},  \code{\link{meanQBAD}} provides the mean, \code{\link{varQBAD}} provides the variance, \code{\link{skewQBAD}} provides the skewness, \code{\link{kurtQBAD}} provides the kurtosis, and  \code{\link{momentQBAD}} provides the \eqn{r}th moment about the location parameter \eqn{\mu} of the asymmetric family of distributions.
#'
#'
#'
#' @references{
#'  Gijbels, I., Karim, R. and Verhasselt, A. (2019a). On quantile-based asymmetric family of distributions: properties and inference. \emph{International Statistical Review}, \url{https://doi.org/10.1111/insr.12324}.
#' }
#'
#'
#'
#' @name momentQBAD
NULL
#' @rdname momentQBAD
#' @export
mu_k<-function(f,k){
  integrand <- function(x) {x^k*f(x)}
  return(2*stats::integrate(integrand, lower = 0, upper = Inf)$ value)}


#' @import Deriv
#' @rdname momentQBAD
#' @export
gamma_k<-function(f,k){
  fn <- function(x) {1/f(x)}
  df<- Deriv::Deriv(f)
  gamma_fn<-function(x){x^(k-1)*(df(x))^2*(fn(x))}
  return(stats::integrate(gamma_fn, lower = 0, upper = 30, stop.on.error=FALSE)$ value)
}

#' @param mu This is the location parameter \eqn{\mu}.
#' @param phi This is the scale parameter  \eqn{\phi}.
#' @param alpha This is the index parameter  \eqn{\alpha}.
#'
#' @param mu_1 This is the quantity \eqn{\int_{0}^{\infty} 2 s f(s) ds}.
#' @rdname momentQBAD
#' @export
meanQBAD<- function(mu,phi,alpha,mu_1){mu+phi*(1-2*alpha)/(alpha*(1-alpha))*mu_1}

#' @param mu_2 This is the quantity \eqn{\int_{0}^{\infty} 2 s^2 f(s) ds}.
#' @rdname momentQBAD
#' @export
varQBAD<- function(mu,phi,alpha,mu_1,mu_2){(phi^2/(alpha^2*(1-alpha)^2))*((1-2*alpha)^2*(mu_2-mu_1^2)+alpha*(1-alpha)*mu_2)}

#' @param mu_3 This is the quantity \eqn{\int_{0}^{\infty} 2 s^3 f(s) ds}.
#' @rdname momentQBAD
#' @export
skewQBAD<- function(alpha,mu_1,mu_2,mu_3){
  nu<-(1-2*alpha)*((1-2*alpha)^2*(mu_3-3*mu_1*mu_2+2*mu_1^3)+alpha*(1-alpha)*(2*mu_3-3*mu_1*mu_2))
  deno<-((1-2*alpha)^2*(mu_2-mu_1^2)+alpha*(1-alpha)*mu_2)^(3/2)
  skew<-(nu/deno)
  return(skew)
}


#' @param mu_4 This is the quantity \eqn{\int_{0}^{\infty} 2 s^4 f(s) ds}.
#' @rdname momentQBAD
#'
#' @export
kurtQBAD<- function(alpha,mu_1,mu_2,mu_3,mu_4){
  nu<-((1-alpha)^5+alpha^5)*mu_4-(1-2*alpha)^2*(4*(1-2*alpha+2*alpha^2)*mu_1*mu_3-6*(1-3*alpha+3*alpha^2)*mu_1^2*mu_2+3*(1-2*alpha)^2*mu_1^4)
  deno<-((1-2*alpha)^2*(mu_2-mu_1^2)+alpha*(1-alpha)*mu_2)^(2)
  kurto<-(nu/deno)
  return(kurto)
}


#' @import stats
#' @examples
#' # Example 1: Let F be a standard normal cumulative distribution function then
#' f_N<-function(s){dnorm(s, mean = 0,sd = 1)} # density function of N(0,1)
#' mu_k(f=f_N,k=1)
#' gamma_k(f=f_N,k=1)
#' mu.1_N=sqrt(2/pi)
#' mu.2_N=1
#' mu.3_N=2*sqrt(2/pi)
#' mu.4_N=4
#' meanQBAD(mu=0,phi=1,alpha=0.5,mu_1=mu.1_N)
#' varQBAD(mu=0,phi=1,alpha=0.5,mu_1=mu.1_N,mu_2=mu.2_N)
#' skewQBAD(alpha=0.5,mu_1=mu.1_N,mu_2=mu.2_N,mu_3=mu.3_N)
#' kurtQBAD(alpha=0.5,mu_1=mu.1_N,mu_2=mu.2_N,mu_3=mu.3_N,mu_4=mu.4_N)
#' momentQBAD(phi=1,alpha=0.5,f=f_N,r=1)
#'
#'
#' # Example 2: Let F be a standard Laplace cumulative distribution function then
#' f_La<-function(s){0.5*exp(-abs(s))} # density function of Laplace(0,1)
#' mu_k(f=f_La,k=1)
#' gamma_k(f=f_La,k=1)
#' mu.1_La=1
#' mu.2_La=2
#' mu.3_La=6
#' mu.4_La=24
#' meanQBAD(mu=0,phi=1,alpha=0.5,mu_1=mu.1_La)
#' varQBAD(mu=0,phi=1,alpha=0.5,mu_1=mu.1_La,mu_2=mu.2_La)
#' skewQBAD(alpha=0.5,mu_1=mu.1_La,mu_2=mu.2_La,mu_3=mu.3_La)
#' kurtQBAD(alpha=0.5,mu_1=mu.1_La,mu_2=mu.2_La,mu_3=mu.3_La,mu_4=mu.4_La)
#' momentQBAD(phi=1,alpha=0.5,f=f_La,r=1)
#' @rdname momentQBAD
#' @export
momentQBAD<- function(phi,alpha,f,r){
  integrand <- function(x) {x^r*f(x)}
  mu.r<-2*stats::integrate(integrand, lower = 0, upper = Inf)$ value
  return(phi^r*((1-alpha)^(r+1)+(-1)^r*alpha^(r+1))*mu.r/(alpha^r*(1-alpha)^r))
}




#' @title Method of moments (MoM) estimation for the quantile-based asymmetric family of distributions.
#' @description Parameter estimation in the quantile-based asymmetric family of densities by using method of moments are discussed in Section 3.1 of Gijbels et al. (2019a).
#' @param y This is a vector of quantiles.
#' @param alpha This is the index parameter  \eqn{\alpha}.  If \eqn{\alpha} is unknown, indicate NULL which is default option. In this case, the sample skewness will be used to estimate \eqn{\alpha}. If \eqn{\alpha} is known, then the value of \eqn{\alpha} has to be specified in the function.
#'
#'
#' @param f This is the reference density function \eqn{f} which is a standard version of a unimodal and symmetric around 0 density.
#'
#' @import stats
#' @return \code{\link{momQBAD}} provides the method of moments estimates of the unknown parameters of the distribution.
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
#' @name momQBAD
NULL
#' @rdname momQBAD
#' @examples
#' # Example 1: Let F be a standard normal cumulative distribution function then
#' f_N<-function(s){dnorm(s, mean = 0,sd = 1)} # density function of N(0,1)
#' y=rnorm(100)
#' momQBAD(y=y,f=f_N,alpha=0.5) # If alpha is known with alpha=0.5
#' momQBAD(y=y,f=f_N) # If alpha is unknown
#'
#'
#' # Example 2: Let F be a standard Laplace cumulative distribution function then
#' f_La<-function(s){0.5*exp(-abs(s))} # density function of Laplace(0,1)
#' momQBAD(y=y,f=f_La,alpha=0.5) # If alpha is known with alpha=0.5
#' momQBAD(y=y,f=f_La) # If alpha is unknown
#' @export
momQBAD<-function(y,f,alpha=NULL){
  muMoM<-function(y,f,alpha=NULL){
    if (is.null(alpha)){
      m1<-mean(y)
      m2<-sum(y^2)/length(y)
      cm2<-sum((y-m1)^2)/length(y)
      cm3<-sum((y-m1)^3)/length(y)
      cm4<-sum((y-m1)^4)/length(y)
      integrand1 <- function(x) {x*f(x)}
      integrand2 <- function(x) {x^2*f(x)}
      integrand3 <- function(x) {x^3*f(x)}
      mu1<-2*stats::integrate(integrand1, lower = 0, upper = Inf)$ value
      mu2<-2*stats::integrate(integrand2, lower = 0, upper = Inf)$ value
      mu3<-2*stats::integrate(integrand3, lower = 0, upper = Inf)$ value
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
      integrand1 <- function(x) {x*f(x)}
      integrand2 <- function(x) {x^2*f(x)}
      mu1<-2*stats::integrate(integrand1, lower = 0, upper = Inf)$ value
      mu2<-2*stats::integrate(integrand2, lower = 0, upper = Inf)$ value

      k1<-(1-2*alpha)*mu1/(alpha*(1-alpha))
      k2<-((1-alpha)^3+alpha^3)*mu2/(alpha^2*(1-alpha)^2)
      mu<-m1-(k1/sqrt(k2-k1^2))*sqrt(m2-m1^2)
    }
    return(mu)
  }
  phiMoM<-function(y,f,alpha=NULL){
    if (is.null(alpha)){
      m1<-mean(y)
      m2<-sum(y^2)/length(y)
      cm2<-sum((y-m1)^2)/length(y)
      cm3<-sum((y-m1)^3)/length(y)
      integrand1 <- function(x) {x*f(x)}
      integrand2 <- function(x) {x^2*f(x)}
      integrand3 <- function(x) {x^3*f(x)}
      mu1<-2*stats::integrate(integrand1, lower = 0, upper = Inf)$ value
      mu2<-2*stats::integrate(integrand2, lower = 0, upper = Inf)$ value
      mu3<-2*stats::integrate(integrand3, lower = 0, upper = Inf)$ value
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
      integrand1 <- function(x) {x*f(x)}
      integrand2 <- function(x) {x^2*f(x)}
      mu1<-2*stats::integrate(integrand1, lower = 0, upper = Inf)$ value
      mu2<-2*stats::integrate(integrand2, lower = 0, upper = Inf)$ value
      k1<-(1-2*alpha)*mu1/(alpha*(1-alpha))
      k2<-((1-alpha)^3+alpha^3)*mu2/(alpha^2*(1-alpha)^2)
      phi<-(1/sqrt(k2-k1^2))*sqrt(m2-m1^2)
    }
    return(phi)
  }

  alphaMoM<-function(y,f){
    m1<-mean(y)
    cm2<-sum((y-m1)^2)/length(y)
    cm3<-sum((y-m1)^3)/length(y)
    integrand1 <- function(x) {x*f(x)}
    integrand2 <- function(x) {x^2*f(x)}
    integrand3 <- function(x) {x^3*f(x)}
    mu1<-2*stats::integrate(integrand1, lower = 0, upper = Inf)$ value
    mu2<-2*stats::integrate(integrand2, lower = 0, upper = Inf)$ value
    mu3<-2*stats::integrate(integrand3, lower = 0, upper = Inf)$ value
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
    list(mu.MoM=muMoM(y,f),phi.MoM=phiMoM(y,f),alpha.MoM=alphaMoM(y,f))
  } else {
    list(mu.MoM=muMoM(y,f),phi.MoM=phiMoM(y,f))}
}



#' @title Maximum likelihood estimation (MLE) for the quantile-based asymmetric family of distributions.
#' @description The log-likelihood function \eqn{\ell_n(\mu,\phi,\alpha)=\ln[L_n(\mu,\phi,\alpha)]}
#' and parameter estimation of \eqn{ \theta=(\mu,\phi,\alpha)} in the three parameter quantile-based asymmetric family of densities
#' by using the maximum likelihood estimation are discussed in Section 3.2 of Gijbels et al. (2019a).
#' @param y This is a vector of quantiles.
#' @param f This is the reference density function \eqn{f} which is a standard version of a unimodal and symmetric around 0 density.
#' @param alpha This is the index parameter  \eqn{\alpha}.
#'
#' @return The maximum likehood estimate of paramter \eqn{\theta=(\mu,\phi,\alpha)} of the quantile-based asymmetric family of densities
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
#' @name mleQBAD
NULL
#' @rdname mleQBAD
#' @examples
#' # Example 1: Let F be a standard normal cumulative distribution function then
#' f_N<-function(s){dnorm(s, mean = 0,sd = 1)} # density function of N(0,1)
#' rnum=rnorm(100)
#' mleQBAD(rnum,f=f_N)
#' mleQBAD(rnum,f=f_N,alpha=.5)
#'
#' # Example 2: Let F be a standard Laplace cumulative distribution function then
#' f_La<-function(s){0.5*exp(-abs(s))} # density function of Laplace(0,1)
#' mleQBAD(rnum,f=f_La)
#' mleQBAD(rnum,f=f_La,alpha=.5)
#' @export
 mleQBAD<-function(y,f,alpha=NULL){
   if (is.null(alpha)){
theta0<-ald::mleALD(y, initial = NA)$par
fn <- function(theta){
  mu=theta[1]
  phi=theta[2]
  alpha=theta[3]
  LL<-log(dQBAD(y,mu,phi,alpha,f))
  return(-sum(LL[!is.infinite(LL)]))}
Est.par=nloptr::bobyqa(x0=theta0,fn = fn ,lower = c(min(y),00.001,0.001), upper = c(max(y),Inf,1))$par
estimated.LL<-LogLikQBAD(y=y,mu=Est.par[1],phi=Est.par[2],alpha=Est.par[3],f)
list(mu.MLE=Est.par[1], phi.MLE=Est.par[2],alpha.MLE=Est.par[3],LogLikelihood=estimated.LL)
 } else {
   theta0<-c(median(y),sd(y),alpha)
   fn <- function(theta){
     mu=theta[1]
     phi=theta[2]
     alpha=alpha
     LL<-log(dQBAD(y,mu,phi,alpha,f))
     return(-sum(LL[!is.infinite(LL)]))}
   Est.par=nloptr::bobyqa(x0=theta0,fn = fn ,lower = c(min(y),00.001,alpha), upper = c(max(y),Inf,alpha))$par
   estimated.LL<-LogLikQBAD(y=y,mu=Est.par[1],phi=Est.par[2],alpha=alpha,f)
   list(mu.MLE=Est.par[1], phi.MLE=Est.par[2],alpha.MLE=Est.par[3],LogLikelihood=estimated.LL)}
 }


