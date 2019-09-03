#' @title Likelihood ratio test to test for symmetry.
#' @description The likelihood ratio test to test for symmetry, in the context of a framework of quantile-based asymmetric family of densities is discussed in Gijbels et al. (2019d).
#' @param y This is a vector of quantiles.
#' @param f This is the reference density function \eqn{f} which is a standard version of a unimodal and symmetric around 0 density.
#'
#' @return The likelihood ratio test statistic with \eqn{P}-value.
#'
#' @references{
#'  Gijbels, I., Karim, R. and Verhasselt, A. (2019d). Test of symmetry and a quantile-based asymmetric family of densities, \emph{Manuscript}.
#' }
#' @name LRTest
NULL
#' @examples
#' # Example: Let F be a standard normal cumulative distribution function then
#' f_N<-function(s){dnorm(s, mean = 0,sd = 1)} # density function of N(0,1)
#' rnum=rnorm(100)
#' LRTest(rnum,f=f_N)
#' @export
LRTest<-function(y,f){
  LL_H1<-mleQBAD(y,f=f)$LogLikelihood
  LL_H0<-mleQBAD(y,f=f,alpha=0.5)$LogLikelihood


  LogLRatio=-2*(LL_H0-LL_H1)
  pvalue=1-pchisq(LogLRatio,df=1)
  list(LRTStat=LogLRatio, P.value=pvalue)
}


