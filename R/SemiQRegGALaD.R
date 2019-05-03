#' @title Semiparametric quantile regression in generalized Laplace distributional settings.
#' @description The local polynomial technique is used to estimate location and scale functions of the quantile-based asymmetric Laplace distribution as discussed in Gijbels et al. (2019c). Using these estimates, the quantile function of the generalized asymmetric Laplace distribution will be estimated. A detailed study can be found in Gijbels et al. (2019b).
#' @param y The is a response variable.
#' @param x This is a conditioning covariate.
#' @param p1 This is the order of the Taylor expansion for the location function (i.e.,\eqn{\mu(X)}) in local polynomial fitting technique. The default value is 1.
#' @param p2 This is the order of the Taylor expansion for the log of scale function (i.e., \eqn{\ln[\phi(X)]}) in local polynomial fitting technique. The default value is 1.
#' @param h This is the bandwidth parameter \eqn{h}.
#' @param alpha This is the index parameter  \eqn{\alpha} of the generalized asymmetric Laplace density. The default value of \eqn{\alpha} is NULL in the code \code{\link{SemiQRegGALaD}}. In this case, the \eqn{\alpha} will be estimated based on the residuals form local linear mean regression.
#' @param g This is the "link" function. The function \eqn{g} is to be differentiated. Therefore, \eqn{g} must be written as a function. For example, {g<-function(y)\{log(y)\}} for log link function.
#' @param lower	This is the lower limit of the domain (support of the random variable) \eqn{f_{\alpha}^g(y;\eta,\phi)}, default {-Inf}.
#' @param upper	This is the upper limit of the domain (support of the random variable) \eqn{f_{\alpha}^g(y;\eta,\phi)}, default {Inf}.
#' @param beta This is a specific probability for estimating \eqn{\beta}th quantile function.
#' @param m This is the number of grid points at which the functions are to be evaluated. The default value is 101.
#' @return The code \code{\link{SemiQRegGALaD}} provides the realized value of the \eqn{\beta}th conditional quantile estimator by using semiparametric quantile regression technique discussed in Gijbels et al. (2019b) and Gijbels et al. (2019c).
#' @import locpol
#' @rdname SemiQRegGALaD
#'
#'
#'
#' @references{
#'  Gijbels, I., Karim, R. and Verhasselt, A. (2019b). Quantile estimation in a generalized  asymmetric distributional setting. To appear in \emph{Springer Proceedings in Mathematics & Statistics, Proceedings of `SMSA 2019', the 14th Workshop on Stochastic Models, Statistics and their Application}, Dresden, Germany, in March 6--8, 2019. Editors: Ansgar Steland, Ewaryst Rafajlowicz, Ostap Okhrin.
#'
#' Gijbels, I., Karim, R. and Verhasselt, A. (2019c).  Semiparametric quantile regression using quantile-based asymmetric family of densities. Manuscript.
#'
#'
#' }
#'

#'
#' @import GoFKernel
#' @examples
#' \donttest{
#'
#' data(LocomotorPerfor)
#' x=log(LocomotorPerfor$Body_Mass)
#' y=LocomotorPerfor$MRRS
#'
#' # For log-link function
#' g_log<-function(y){log(y)}
#' h_ROT =  0.9030372
#' fit<-SemiQRegGALaD(beta=0.90,x,y,p1=1,p2=1,h=h_ROT,g=g_log,lower=0)
#' plot(x,y)
#' lines(fit$x0,fit$qf_g)
#'
#' }
#' @name SemiQRegGALaD
NULL
#' @rdname SemiQRegGALaD
#' @export
SemiQRegGALaD<-function(beta,x, y,p1=1,p2=1,  h,alpha=NULL,g,lower = -Inf, upper = Inf, m = 101){
  z=g(y)
  g.inv<-GoFKernel::inverse(g,lower,upper)
  SemiQReg_Z<-SemiQRegALaD(beta,x, z,  p1,p2,  h,alpha, m)
  qf_g<-NA
for(i in 1:length(SemiQReg_Z$fit_beta_ALaD)){
  qf_g[i]<-g.inv(SemiQReg_Z$fit_beta_ALaD[i])}
  list(x0=SemiQReg_Z$x0,qf_g=qf_g)
}
