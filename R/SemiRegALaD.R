#' @title Semiparametric quantile regression in quantile-based asymmetric Laplace distributional settings.
#' @description The local polynomial technique is used to estimate location and scale function of the quantile-based asymmetric Laplace distribution discussed in Gijbels et al. (2019c).
#' The semiparametric quantile estimation technique is used to estimate \eqn{\beta}th conditional quantile function in quantile-based asymmetric Laplace distributional setting discussed in Gijbels et al. (2019b) and Gijbels et al. (2019c).
#' @param y The is a response variable.
#' @param x This a conditioning covariate.
#' @param p1 This is the order of the Taylor expansion for the location function (i.e.,\eqn{\mu(X)}) in local polynomial fitting technique. The default value is 1.
#' @param p2 This is the order of the Taylor expansion for the log of scale function (i.e., \eqn{\ln[\phi(X)]}) in local polynomial fitting technique. The default value is 1.
#' @param h This is the bandwidth parameter \eqn{h}.
#' @param alpha This is the index parameter  \eqn{\alpha} of the quantile-based asymmetric Laplace density. The default value is 0.5 in the codes code \code{\link{locpolALaD_x0}} and code \code{\link{locpolALaD}}. The default value of \eqn{\alpha} is NULL in the code \code{\link{SemiQRegALaD}}. In this case, the \eqn{\alpha} will be estimated based on the residuals of local linear mean regression.
#' @param x0 This is a grid-point \eqn{x_0} at which the function is to be estimated.
#' @param tol the desired accuracy. See details in \code{\link[stats]{optimize}}.
#' @return The code \code{\link{locpolALaD_x0}} provides the realized value of the local maximum likelihood estimator of \eqn{\widehat{\theta}_{rj}(x_0)} for \eqn{(r\in \{1,2\}; j=1,2,...,p_r)} with the estimated approximate asymptotic bias and variance at the grind point \eqn{x_0} discussed in Gijbels et al. (2019c).
#' @import quantreg
#'
#'
#' @references{
#'  Gijbels, I., Karim, R. and Verhasselt, A. (2019b). Quantile estimation in a generalized  asymmetric distributional setting. To appear in \emph{Springer Proceedings in Mathematics & Statistics, Proceedings of `SMSA 2019', the 14th Workshop on Stochastic Models, Statistics and their Application}, Dresden, Germany, in March 6--8, 2019. Editors: Ansgar Steland, Ewaryst Rafajlowicz, Ostap Okhrin.
#'
#' Gijbels, I., Karim, R. and Verhasselt, A. (2019c).  Semiparametric quantile regression using quantile-based asymmetric family of densities. Manuscript.
#'
#'
#' }
#' @name SemiQRegALaD
NULL
#' @rdname SemiQRegALaD
#'
#' @examples
#' \donttest{
#' data(Hurricane)
#' locpolALaD_x0(Hurricane$Year, Hurricane$WmaxST, p1=1,p2=1,h=2.18,
#'  alpha=0.16,x0=median(Hurricane$Year))
#'  }
#' @export
locpolALaD_x0<-function(x,y,p1=1,p2=1,h, alpha=0.5,x0,tol=1e-08)
{
  n<-length(y)
  x<-as.matrix(x)
  # r<-matrix(NA,n,1)
  u1<-matrix(0,n,p1)
  for (j in 1:p1)
  {
    for (i in 1:n)
    {
      u1[i,j]<-(x[i,1]-x0)^j
    }
  }

  X1<-cbind(1, u1)
  z <- x - x0
  K_h <- dnorm(z/h)

  if (p2==0){
    X2<-matrix(1, length(x),1)

    b1=quantreg::rq(y ~ X1[,-1], weights = K_h , tau = alpha, ci = FALSE)$coefficients
    r<-{ifelse(y-X1 %*% b1<0,alpha-1 ,alpha)}
    theta2<-log(sum(((y-X1 %*% b1)*r)*K_h)/sum(K_h))
  }  else  if (p2==1){
    u2<-matrix(0,n,p2)
    for (j in 1:p2)
    {
      for (i in 1:n)
      {
        u2[i,j]<-(x[i,1]-x0)^j
      }
    }

    X2<-cbind(1, u2)


    fit <- quantreg::rq(y ~ X1[,-1], weights = K_h , tau = alpha, ci = FALSE)
    theta1IniValue<-as.vector(fit$coefficients)
    b1=theta1IniValue
    fitted.values = X1 %*% b1


    r1<-function(b1){ifelse(y-X1 %*% b1<0,alpha-1 ,alpha)}
    r=r1(b1)


    rho_hat<-((y-fitted.values)*r)
    #	theta2IniValue<-exp(lm(formula = log(rho_hat+.1)  ~ u2)$ coefficients)
    #      theta2IniValue<-log(sum(rho_hat*K_h)/sum(K_h))
    theta2_initial<-c(log(sum(rho_hat*K_h)/sum(K_h)),rep(0.005,p2))
    for(j in 1:100) {
      ### Likelihood function for Theta2
      rho_hat<-((y-X1 %*% b1)*r1(b1))


      ### theta2 estimation
      ## Newton-Raphson Method
      NRM<-function(theta21_old,tol=1e-07,N=100){
        i=1;thetha_new=theta21_old
        u=X2[,-1]
        pp=numeric(N)
        while (i<=N){
          c0<-sum(exp(-thetha_new*u)*rho_hat*K_h);
          c1<-sum(exp(-thetha_new*u)*u*rho_hat*K_h);
          c2<-sum(exp(-thetha_new*u)*u^2*rho_hat*K_h);
          f_theta21<-sum(K_h)*c1/c0-sum(u*K_h);
          df_theta21<-(sum(K_h)*c1^2-sum(K_h)*c0*c2)/c0^2;
          thetha_new<-theta21_old-f_theta21/df_theta21
          pp[i]=thetha_new
          i=i+1
          if (abs(thetha_new-theta21_old)<tol) break
          theta21_old=thetha_new
        }
        return(pp[(i-1)])
      }

      theta21<-NRM(0)
      theta20<--log(sum(K_h)/sum(exp(-theta21*X2[,-1])*rho_hat*K_h))

      theta2<-c(theta20,theta21)


      b1_old   = b1
      theta2_old   = theta2
      quantreg::rq(y ~ X1[,-1], weights = K_h/exp(theta2[1]+theta2[2]*(x-x0)) , tau = alpha, ci = FALSE)

      b1<-as.vector(quantreg::rq(y ~ X1[,-1], weights = K_h/exp(theta2[1]+theta2[2]*(x-x0)) , tau = alpha, ci = FALSE)$coefficients)
      if(sqrt(crossprod(b1-b1_old)+crossprod(theta2-theta2_old)) < tol) break
    }
  } else {
    u2<-matrix(0,n,p2)
    for (j in 1:p2)
    {
      for (i in 1:n)
      {
        u2[i,j]<-(x[i,1]-x0)^j
      }
    }

    X2<-cbind(1, u2)
    fit <- quantreg::rq(y ~ X1[,-1], weights = K_h , tau = alpha, ci = FALSE)
    b1<-as.vector(fit$coefficients)
    fitted_p1 = X1 %*% b1
    r_p1<-matrix(NA,n,1)
    for (i in 1:length(y))
    {
      if(y[i]-fitted_p1[i]<0) r_p1[i]=alpha-1
      else r_p1[i]=alpha

    }
    rho_hat_p1<-((y-fitted_p1)*r_p1)
    LogLikeTheta2_p2 <- function(theta2_p2) sum(X2%*%theta2_p2+rho_hat_p1/exp(X2%*%theta2_p2))
    y_log_rho<-log(rho_hat_p1)
    y_log_rho[!is.finite(y_log_rho)]<- NA
    theta2IValue<-as.vector(lm(y_log_rho~X2[,-1],na.action=na.omit))$coefficients

    theta20<-optim(theta2IValue, LogLikeTheta2_p2)$par


    fn<-function(theta){
      mu<- as.matrix(X1)%*%as.vector(theta[1:(p1+1)])
      phi<-exp(as.matrix(X2)%*%as.vector(theta[(p1+2):(p1+p2+2)]))
      LL<-log(dALaD(y,mu,phi,alpha))
      return(-sum(LL[!is.infinite(LL)]))}

    theta0<-c(b1,theta20)
    Est.par=optim(par=theta0,fn = fn)$par
    b1<-Est.par[1:(p1+1)]
    theta2<-Est.par[(p1+2):(p1+p2+2)]
  }


  ### Bias and Variance Estiamtion
  r_p1_2<-matrix(NA,n,1)
  u_p1_2<-matrix(0,n,p1+2)
  u_p2_2<-matrix(0,n,p2+2)
  for (j in 1:(p1+2))
  {
    for (i in 1:n)
    {
      u_p1_2[i,j]<-(x[i]-x0)^j
    }
  }

  for (j in 1:(p2+2))
  {
    for (i in 1:n)
    {
      u_p2_2[i,j]<-(x[i]-x0)^j
    }
  }


  X_p1_2<-cbind(1, u_p1_2)
  X_p2_2<-cbind(1, u_p2_2)
  K_h <- dnorm((x-x0)/h)
  fit_p1_2 <- quantreg::rq(y ~ X_p1_2[,-1], weights = K_h , tau = alpha, ci = FALSE)
  b_p1_2<-as.vector(fit_p1_2$coefficients)
  fitted_p1_2 = X_p1_2 %*% b_p1_2
  for (i in 1:length(y))
  {
    if(y[i]-fitted_p1_2[i]<0) r_p1_2[i]=alpha-1
    else r_p1_2[i]=alpha

  }
  rho_hat_p1_2<-((y-fitted_p1_2)*r_p1_2)
  LogLikeTheta2 <- function(theta2_p2_2) sum(X_p2_2%*%theta2_p2_2+rho_hat_p1_2/exp(X_p2_2%*%theta2_p2_2))
  y_log_rho<-log(rho_hat_p1_2)
  y_log_rho[!is.finite(y_log_rho)]<- NA
  theta2IValue<-as.vector(lm(y_log_rho~X_p2_2[,-1],na.action=na.omit))$coefficients

  theta2_p2_2<-optim(theta2IValue, LogLikeTheta2)$par


  v11_x0<-sum((ifelse(fit_p1_2$residuals>0,alpha,alpha-1))^2*K_h/exp(2*X_p2_2%*%theta2_p2_2))/sum(K_h)
  v22_x0<-sum(((-1+rho_hat_p1_2/exp(2*X_p2_2%*%theta2_p2_2))^2)*K_h)/sum(K_h)

  delta_1<-(ifelse(fit_p1_2$residuals>0,alpha,alpha-1))/exp(X_p2_2%*%theta2_p2_2)
  delta_2<--1+rho_hat_p1_2/exp(X_p2_2%*%theta2_p2_2)
  W<-diag(c(K_h))
  S_11<-t(X1)%*%W%*%X1
  S_22<-t(X2)%*%W%*%X2
  W_K_2<-diag(c(K_h^2))
  S_11_bar_n<-t(X1)%*%W_K_2%*%X1
  S_22_bar_n<-t(X2)%*%W_K_2%*%X2
  Bias_1_x0<-(solve(S_11))%*%(t(X1)%*%W%*%delta_1)/v11_x0
  Bias_2_x0<-(solve(S_22))%*%(t(X2)%*%W%*%delta_2)/v22_x0
  Var_1_x0<-solve(S_11)%*%S_11_bar_n%*%solve(S_11)/v11_x0
  Var_2_x0<-solve(S_22)%*%S_22_bar_n%*%solve(S_22)/v22_x0

  list(theta1_x0=b1,theta2_x0=theta2,Bias_1_x0=Bias_1_x0,Bias_2_x0=Bias_2_x0,Var_1_x0=Var_1_x0,Var_2_x0=Var_2_x0)
}

#' @param m This is the number of grid points at which the functions are to be evaluated. The default value is 101.
#' @return The code \code{\link{locpolALaD}} provides the realized value of the local maximum likelihood estimator of \eqn{\widehat{\theta}_{r0}(x_0)} for \eqn{(r\in \{1,2\})} with the estimated approximate asymptotic bias and variance at all \eqn{m} grind points \eqn{x_0} discussed in Gijbels et al. (2019c).
#'
#' @rdname SemiQRegALaD
#'
#' @examples
#' \donttest{
#' data(Hurricane)
#' locpolALaD(Hurricane$Year, Hurricane$WmaxST, p1=1,p2=1,h=2.18, alpha=0.16)
#' }
#' @export
locpolALaD<-function(x, y,p1=1,p2=1, h, alpha=0.5,m = 101)
{

  xx <- seq(min(x)+.01*h, max(x)-.01*h, length = m)
  theta_10<-matrix(NA,length(xx),1)

  theta_20<-matrix(NA,length(xx),1)
  Bias_10<-matrix(NA,length(xx),1)
  Bias_20<-matrix(NA,length(xx),1)
  Var_10<-matrix(NA,length(xx),1)
  Var_20<-matrix(NA,length(xx),1)

  for (ii in 1:length(xx)) {

    fit_x0<-locpolALaD_x0(x,y,p1,p2,h,alpha,x0=xx[ii])
    theta_10[ii] <- fit_x0$ theta1_x0[[1]]
    theta_20[ii] <-fit_x0$theta2_x0[[1]]
    Bias_10[ii]<- fit_x0$Bias_1_x0[1,1]
    Bias_20[ii]<- fit_x0$Bias_2_x0[1,1]
    Var_10[ii]<- fit_x0$Var_1_x0[1,1]
    Var_20[ii]<- fit_x0$Var_2_x0[1,1]
  }
  list(x0 = xx, theta_10= theta_10, theta_20= theta_20,Bias_10=Bias_10,Bias_20=Bias_20,Var_10=Var_10,Var_20=Var_20)
}



#' @param beta This is a specific probability for estimating \eqn{\beta}th quantile function.
#' @return The code \code{\link{SemiQRegALaD}} provides the realized value of the \eqn{\beta}th conditional quantile estimator by using semiparametric quantile regression technique discussed in Gijbels et al. (2019b) and Gijbels et al. (2019c).
#' @import locpol
#' @rdname SemiQRegALaD
#' @examples
#' \donttest{
#' ## For Hurricane Data
#' data(Hurricane)
#' Hurricane<-Hurricane[which(Hurricane$Year>1970),]
#'
#' plot(Hurricane$Year,Hurricane$WmaxST)
#'
#' h=2.181082
#' alpha=0.1649765
#' gridPoints=101
#' fit_ALaD <-locpolALaD(Hurricane$Year, Hurricane$WmaxST, p1=1,p2=1,h=h, alpha=alpha, m = gridPoints)
#' str(fit_ALaD)
#' par(mgp=c(2,.4,0),mar=c(5,4,4,1)+0.01)
#'
#' # For phi plot
#' plot(fit_ALaD$x0,exp(fit_ALaD$theta_20),ylab=expression(widehat(phi)(x[0])),xlab="Year",
#' type="l",font.lab=2,cex.lab=1.5,bty="l",cex.axis=1.5,lwd =3)
#'
#' ## For theta2 plot
#' plot(fit_ALaD$x0,fit_ALaD$theta_20,ylab=expression(bold(widehat(theta[2]))(x[0])),
#' xlab="Year",type="l",col=c(1), lty=1, font.lab=1,cex.lab=1.5,bty="l",cex.axis=1.3,lwd =3)
#'
#'
#'
#' #### Estimated Quantile lines by ALaD
#' par(mgp=c(2.5, 1, 0),mar=c(5,4,4,1)+0.01)
#' # X11()
#' plot(Hurricane$Year, Hurricane$WmaxST, xlab = "Year",ylim=c(20,210),
#' ylab = "Maximum Wind Spreed",font.lab=1,cex.lab=1.3,bty="l",pch=20,cex.axis=1.3)
#'
#' lines(fit_ALaD$x0,fit_ALaD$theta_10, type='l',col=c(4),lty=1,lwd =3)
#'
#' #####  Conditioanl Quantile line for ALaD
#'
#' lines(fit_ALaD$x0,SemiQRegALaD(beta=0.50,Hurricane$Year, Hurricane$WmaxST,
#' p1=1,p2=1, h=h,alpha=alpha,m=gridPoints)$fit_beta_ALaD,type='l',col=c(1),lty=1,lwd =3)
#'
#'
#' lines(fit_ALaD$x0,SemiQRegALaD(beta=0.90,Hurricane$Year, Hurricane$WmaxST,
#' p1=1,p2=1, h=h,alpha=alpha,m=gridPoints)$fit_beta_ALaD,type='l',col=c(14),lty=1,lwd =3)
#' lines(fit_ALaD$x0,SemiQRegALaD(beta=0.95,Hurricane$Year, Hurricane$WmaxST,
#' p1=1,p2=1, h=h,alpha=alpha,m=gridPoints)$fit_beta_ALaD,type='l',col=c(19),lty=1,lwd =3)
#'
#' # Add local linear mean regression line
#' library(locpol)
#' fit_mean<-locpol(WmaxST~Year, data=Hurricane,kernel=gaussK,deg=1,
#' xeval=NULL,xevalLen=101)
#'
#' lines(fit_mean$lpFit[,1], fit_mean$lpFit[,2],type='l',col=c(2),lty=1,lwd =3)
#' axis(1, at = c(1975, 1985, 1995,2005,2015),cex.axis=1.3)
#' axis(2, at = c(25, 75, 125,175),cex.axis=1.3)
#'
#' legend("topright", legend = c(expression(beta==0.1650), expression(beta==0.50),
#' "Mean line",expression(beta==0.90), expression(beta==0.95)), col = c(4,1,2,14,19),
#'  lty=c(1,1,1,1,1), inset = 0, lwd = 3,cex=1.2)
#' }
#' @export
SemiQRegALaD<-function(beta,x, y,  p1=1,p2=1,  h,alpha=NULL, m = 101){
  if (is.null(alpha)){
    alpha_fn<-function(x,y){
      Newdata<-data.frame(x,y)
      fit_mean<-locpol(y~x, data=Newdata,weig=rep(1,nrow(Newdata)) ,kernel=gaussK,deg=1,xeval=NULL,xevalLen=101)
      MeanRegResidual<-fit_mean$residuals
      MeanRegResidual<-sort(MeanRegResidual)
      alpha<-mleALaD(MeanRegResidual)$alpha
      list(alpha=alpha)
    }

    alpha<-alpha_fn(x,y)$alpha}
  fit_theta <-locpolALaD(x, y,p1,p2, h, alpha,m)
  theta1=fit_theta$theta_10
  theta2=fit_theta$theta_20
  quanty<-matrix(NA,length(theta1),1);
  for (i in 1:length(theta1))
  {
    if (beta <alpha)
    {
      quanty[i]<-theta1[i]+(exp(theta2[i])/(1-alpha))*log(beta/alpha);

    }
    else
    {
      quanty[i]<-theta1[i]-(exp(theta2[i])/alpha)*log((1-beta)/(1-alpha));
    }
  }
  C_alpha_beta<-qALaD(beta,mu=0,phi=1,alpha)
  Bias_q_x0<-fit_theta$Bias_10+C_alpha_beta*fit_theta$Bias_20*exp(fit_theta$theta_20)
#  Var_q_x0<-fit_theta$Var_10+(C_alpha_beta)^2*fit_theta$Var_20*exp(2*fit_theta$theta_20)*exp(fit_theta$Bias_20)
  Var_q_x0<-fit_theta$Var_10+(C_alpha_beta)^2*fit_theta$Var_20*exp(2*fit_theta$theta_20)
    list(x0=fit_theta$x0,fit_beta_ALaD=quanty,Bias_q_x0=Bias_q_x0,Var_q_x0=Var_q_x0,fit_alpha_ALaD=theta1,alpha=alpha,beta=beta)
}

