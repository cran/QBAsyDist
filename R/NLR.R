#' @title Nonparametric likelihood ratio test to test for symmetry.
#' @description The nonparametric likelihood ratio test to test for symmetry is discussed in Gijbels et al. (2019d).
#' @param y This is a vector of quantiles.
#' @param f This is the reference density function \eqn{f} which is a standard version of a unimodal and symmetric around 0 density.
#' @param F This is the cumulative distribution function \eqn{F} of a unimodal and symmetric around 0 reference density function \eqn{f}.
#' @param method The method to be used for drawing bootstrap samples. The default method is a smooth bootstrap procedure. The density function \eqn{f}, the cumulative distribution function \eqn{F} and the quantile function \eqn{QF} are required for parametric bootstrap procedure. If \eqn{QF} is not given, then the numerical \eqn{QF} will be used.
#' @param nboot The number of bootstrap samples desired. The default number is 500.
#' @param QF This is the quantile function of the reference density \eqn{f}.
#' @return The nonparametric likelihood ratio test statistic with bootstrap \eqn{P}-value.
#'
#'
#' @references{
#'  Gijbels, I., Karim, R. and Verhasselt, A. (2019d). Test of symmetry and a quantile-based asymmetric family of densities, \emph{Manuscript}.
#' }
#'
#' @import scdensity
#'
#' @name NLRTest
NULL
#' @rdname NLRTest
#'
#' @examples
#' \donttest{
#' # Example 1: Let F be a standard normal cumulative distribution function then
#' f_N<-function(s){dnorm(s, mean = 0,sd = 1)} # density function of N(0,1)
#' F_N<-function(s){pnorm(s, mean = 0,sd = 1)} # distribution function of N(0,1)
#' QF_N<-function(beta){qnorm(beta, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)}
#' # Example: STRength dataset
#' my.sample<-c(1.901,2.132,2.203,2.228,2.257,2.350,2.361,2.396,2.397,
#'              2.445,2.454,2.474,2.518,2.522,2.525,2.532,2.575,2.614,2.616,
#'              2.618,2.624,2.659,2.675,2.738,2.740,2.856,2.917,2.928,2.937,
#'              2.937,2.977,2.996,3.030,3.125,3.139,3.145,3.220,3.223,3.235,
#'             3.243,3.264,3.272,3.294,3.332,3.346,3.377,3.408,3.435,3.493,
#'              3.501,3.537,3.554,3.562,3.628,3.852,3.871,3.886,3.971,4.024,
#'              4.027,4.225,4.395,5.020)
#'
#'
#' NLRTest(my.sample,f=NULL,F=NULL,QF=NULL,method=c("smooth"),nboot=500)
#' NLRTest(my.sample,f=f_N,F=F_N,QF=QF_N,method=c("parametric"),nboot=500)
#' NLRTest(my.sample,f=f_N,F=F_N,QF=NULL,method=c("parametric"),nboot=500)
#' }
#' @export
NLRTest<-function(y,f=NULL,F=NULL,QF=NULL,method=c("smooth","parametric"),nboot=500){
method <- match.arg(method)
LLRatio_SemiPar<-function(y){
  f.hat_H1<-scdensity(y,constraint = c("unimodal"))
  LLik_H1<-sum(log(f.hat_H1$fhat(y)))

  # Under H0
  f.hat_H0<-scdensity(y,constraint = c("unimodal", "symmetric"),opts = list(pointOfSymmetry = median(y)))
  LLik_H0<-sum(log(f.hat_H0$fhat(y)))
  LLRatio=-2*(LLik_H0-LLik_H1)
  return(LLRatio=LLRatio)}

################### Parametric bootstrap p-value

Pvalue.boot<-function(y,Tstatistic,nboot,f,F, QF = NULL){
  n=length(y)
  fit<-mleQBAD(y,f,alpha=0.5)
  Tstat<-function(y){Tstatistic(y)}
  boot.stat_SPar<-NA
  for (i in 1:nboot){
    #  print(i)
    tryCatch({
      boot.sample<-rQBAD(n, mu=fit$mu.MLE,phi=fit$phi.MLE,alpha=0.5, F, QF)
      boot.stat_SPar[i]<-Tstat(boot.sample)
    }, error=function(e){})
  }
  boot.pvalue=2*min(length(which((boot.stat_SPar-Tstat(y))>(Tstat(y))))/nboot,length(which((boot.stat_SPar-Tstat(y))<(Tstat(y))))/nboot)
list(t1=boot.stat_SPar,t0=Tstat(y),boot.pvalue=boot.pvalue)}


################### Smooth bootstrap p-value
### bootstrapping P.value

rf.hat <- function(n, x0,fhat)
{
  y0=fhat
  diffs <- diff(x0)
  # Making sure we have equal increments
  stopifnot(all(abs(diff(x0) - mean(diff(x0))) < 1e-9))
  total <- sum(y0)
  y0_pro <- y0 / total
  ydistr <- cumsum(y0_pro)
  yunif <- runif(n)
  indices <- sapply(yunif, function(y) min(which(ydistr > y)))
  x <- x0[indices]

  return(x)
}

Pvalue.Smooth.boot<-function(y,Tstatistic,nboot){

  f.hat_H0<-scdensity(y,constraint = c("unimodal", "symmetric"),opts = list(pointOfSymmetry = median(y)))
  fhat_H0<-f.hat_H0$fhat(y)
  x0=f.hat_H0$extra$conCheckGrid
  y0=f.hat_H0$fhat(x0)
  Tstat<-function(y){Tstatistic(y)}
  boot.stat_SPar<-NA
  for (i in 1:nboot){
    #  print(i)
    tryCatch({
      boot.sample<-rf.hat(n=length(y), x0,y0)
      boot.stat_SPar[i]<-Tstat(boot.sample)
    }, error=function(e){})
  }
  boot.pvalue=2*min(length(which((boot.stat_SPar-Tstat(y))>(Tstat(y))))/nboot,length(which((boot.stat_SPar-Tstat(y))<(Tstat(y))))/nboot)
list(t1=boot.stat_SPar,t0=Tstat(y),boot.pvalue=boot.pvalue)}


if (method=="parametric"){ return(Pvalue.boot(y,Tstatistic=LLRatio_SemiPar,nboot=nboot,f,F, QF = NULL))}
if (method=="smooth"){ return(Pvalue.Smooth.boot(y,Tstatistic=LLRatio_SemiPar,nboot=nboot))}
}




