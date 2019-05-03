#' @title Semiparametric quantile regression in quantile-based asymmetric normal distributional settings.
#' @description The local polynomial technique is used to estimate location and scale function of the quantile-based asymmetric normal distribution discussed in Gijbels et al. (2019c).
#' The semiparametric quantile estimation technique is used to estimate \eqn{\beta}th conditional quantile function in quantile-based asymmetric normal distributional setting discussed in Gijbels et al. (2019b) and Gijbels et al. (2019c).
#' @param y The is a response variable.
#' @param x This a conditioning covariate.
#' @param p1 This is the order of the Taylor expansion for the location function (i.e.,\eqn{\mu(X)}) in local polynomial fitting technique. The default value is 1.
#' @param p2 This is the order of the Taylor expansion for the log of scale function (i.e., \eqn{\ln[\phi(X)]}) in local polynomial fitting technique. The default value is 1.
#' @param h This is the bandwidth parameter \eqn{h}.
#' @param alpha This is the index parameter  \eqn{\alpha} of the quantile-based asymmetric normal density. The default value is 0.5 in the codes code \code{\link{locpolAND_x0}} and code \code{\link{locpolAND}}. The default value of \eqn{\alpha} is NULL in the code \code{\link{SemiQRegAND}}. In this case, \eqn{\alpha} will be estimated based on the residuals from local linear mean regression.
#' @param x0 This is a grid-point \eqn{x_0} at which the function is to be estimated.
#' @param tol the desired accuracy. See details in \code{\link[stats]{optimize}}.
#' @return The code \code{\link{locpolAND_x0}} provides the realized value of the local maximum likelihood estimator of \eqn{\widehat{\theta}_{rj}(x_0)} for \eqn{(r\in \{1,2\}; j=1,2,...,p_r)} with the estimated approximate asymptotic bias and variance at the grind point \eqn{x_0} discussed in Gijbels et al. (2019c).
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
#' @name SemiQRegAND
NULL
#' @rdname SemiQRegAND
#'
#' @examples
#' \donttest{
#' data(LocomotorPerfor)
#' x=log(LocomotorPerfor$Body_Mass)
#' y=log(LocomotorPerfor$MRRS)
#' h_ROT =  0.9030372
#' locpolAND_x0(x, y, p1=1,p2=1,h=h_ROT,alpha=0.50,x0=median(x))
#'}
#' @export
locpolAND_x0<-function(x,y,p1=1,p2=1,h, alpha=0.5,x0,tol=1e-08)
{
  n<-length(y)
  x<-as.matrix(x)
  u1<-matrix(0,n,p1)
  for (j in 1:p1)
  {
    for (i in 1:n)
    {
      u1[i,j]<-(x[i,1]-x0)^j
    }
  }
  X1<-cbind(1, u1)
  irls <-function(x,y,p,h, alpha,x0,maxit=25, tol=1e-08)
  {
    n<-length(y)
    x<-as.matrix(x)
    u<-matrix(0,n,p)
    for (j in 1:p)
    {
      for (i in 1:n)
      {
        u[i,j]<-(x[i,1]-x0)^j
      }
    }

    X<-cbind(1, u)
    z <- x - x0
    K_h <- dnorm(z/h)

    ##### Weight Matrix Estimation
    W <-as.vector(array(NA,length(y)))
    U <-as.vector(array(NA,length(y)))
    b = rep(0,ncol(X))
    for(j in 1:maxit)
    {
      mu = X %*% b
      for (i in 1:length(y))
      {
        if(y[i]-mu[i]>0) W[i]=alpha^2*K_h[i] else W[i]=(1-alpha)^2*K_h[i]
      }
      b_old   = b
      b      = solve(crossprod(X,W*X), crossprod(X,W*y), tol=2*.Machine$double.eps)
      if(sqrt(crossprod(b-b_old)) < tol) break
    }
    fitted.values = X %*% b
    for (i in 1:length(y))
    {
      if(y[i]-fitted.values[i]>0) W[i]=alpha^2*K_h[i]
      else W[i]=(1-alpha)^2*K_h[i]

    }
    for (i in 1:length(y))
    {
      if(y[i]-fitted.values[i]>0) U[i]=alpha^2
      else U[i]=(1-alpha)^2

    }
    SSE<-crossprod((y-fitted.values)^2,W)
    #	phi2=phi^2
    phi2<-sum((y-fitted.values)^2*U*K_h)/(sum(K_h))
    phi<-sqrt(phi2)
    list(b1=b,phi=phi,phi2=phi2,iterations=j)
  }

  if (p2==0) {
    b1=irls(x,y,p=p1,h,alpha,x0)$b
    theta2<-(1/2)*log(irls(x,y,p=p1,h,alpha,x0)$phi2)
    X2<-rep(1, length(x))
  }  else if(p2==1){
    u2<-matrix(0,n,p2)
    for (j in 1:p2)
    {
      for (i in 1:n)
      {
        u2[i,j]<-(x[i,1]-x0)^j
      }
    }

    X2<-cbind(1, u2)
    z <- x - x0
    K_h <- dnorm(z/h)


    b1=irls(x,y,p=p1,h,alpha,x0)$b

    theta2_initial=c((1/2)*log(irls(x,y,p=p1,h,alpha,x0)$phi2),rep(0.005,p2))

    for(j in 1:100) {
      W<-ifelse(y>X1%*% as.vector(b1),alpha^2*K_h,(1-alpha)^2*K_h)

      SqError<-((y-X1%*%as.vector(b1))^2*W)
      #	theta2IniValue<-lm(formula = log(SqError)  ~ u2)$ coefficients

      fr<- function(theta2) {   ## (-ve) likelihood function
        theta20 <- theta2[1]
        theta21 <- theta2[2]
        sum(theta20+theta21*(x-x0)*K_h)+sum(SqError/(2*exp(2*theta20+2*theta21*(x-x0))))
      }


      gr<- function(theta2) {   ## gradian
        theta20 <- theta2[1]
        theta21 <- theta2[2]
        c(-sum(SqError/exp(2*theta20+2*theta21*(x-x0)))+sum(K_h),sum(SqError*(x-x0)/(exp(2*theta20+2*theta21*(x-x0))))-sum((x-x0)*K_h))
      }



      # theta2<-optim(c(theta2_initial,0.05), fr, gr, method = "Nelder-Mead")$par
      theta2<-optim(theta2_initial, fr, gr, method = "BFGS")$par

      #theta2<-(1/2)*log(local_irls(x,y,p=p1,h,alpha,x0)$phi2)


      W1<-ifelse(y>X1%*% as.vector(b1),alpha^2*K_h/(2*exp(2*X2%*% as.vector(theta2))),(1-alpha)^2*K_h/(2*exp(2*X2%*% as.vector(theta2))))

      b1_old   = b1
      theta2_old   = theta2
      b1=solve(crossprod(X1,diag(c(W1))%*%X1), crossprod(X1,diag(c(W1))%*%y))
      #  print(j);print(b1);print(theta2)
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


    theta1_initial=irls(x,y,p=p1,h,alpha,x0)$b

    theta2_initial=c((1/2)*log(irls(x,y,p=p1,h,alpha,x0)$phi2),rep(0.005,p2))

    fn<-function(theta){
      mu<- as.matrix(X1)%*%as.vector(theta[1:(p1+1)])
      phi<-exp(as.matrix(X2)%*%as.vector(theta[(p1+2):(p1+p2+2)]))
      LL<-log(dAND(y,mu,phi,alpha))
      return(-sum(LL[!is.infinite(LL)]))}

    theta0<-c(theta1_initial,theta2_initial)
    Est.par=optim(par=theta0,fn = fn)$par

    b1<-Est.par[1:(p1+1)]
    theta2<-Est.par[(p1+2):(p1+p2+2)]
  }

  ### Bias_x0 and Variance_x0
  ### Bias and Variance Estiamtion
  a=2;
  u_p1_2<-matrix(0,n,p1+a)
  u_p2_2<-matrix(0,n,p2+a)
  for (j in 1:(p1+a))
  {
    for (i in 1:n)
    {
      u_p1_2[i,j]<-(x[i]-x0)^j
    }
  }

  for (j in 1:(p2+a))
  {
    for (i in 1:n)
    {
      u_p2_2[i,j]<-(x[i]-x0)^j
    }
  }


  X1_p1_2<-cbind(1, u_p1_2)
  X2_p2_2<-cbind(1, u_p2_2)
  K_h <- dnorm((x-x0)/h)
  theta1_p1_2=irls(x,y,p=(p1+a),h,alpha,x0)$b

  theta2_p2_2=c((1/2)*log(irls(x,y,p=p1,h,alpha,x0)$phi2),rep(0.005,(p2+a)))

  fn<-function(theta){
    mu<- as.matrix(X1_p1_2)%*%as.vector(theta[1:(p1+3)])
    phi<-exp(as.matrix(X2_p2_2)%*%as.vector(theta[(p1+4):(p1+p2+6)]))
    LL<-log(dAND(y,mu,phi,alpha))
    return(-sum(LL[!is.infinite(LL)]))}

  theta0_2<-c(theta1_p1_2, theta2_p2_2)
  Est.par_2=optim(par=theta0_2,fn = fn)$par

  b1_p1_2<- Est.par_2[1:(p1+a+1)]
  b2_p2_2<- Est.par_2[(p1+a+2):(p1+p2+2+2*a)]
  e_1_hat<-  X1_p1_2[,(p1+1):(p1+a)]%*% b1_p1_2[(p1+1):(p1+a)]
  e_2_hat<-  X2_p2_2[,(p2+1):(p2+a)]%*% b2_p2_2[(p2+1):(p2+a)]
  yhat_p1_2<-as.matrix(X1)%*%as.vector(b1)+e_1_hat
  yhat_p2_2<-as.matrix(X2)%*%as.vector(theta2)+e_2_hat

  delta_1<-((y-yhat_p1_2)*ifelse((y-yhat_p1_2)>0,alpha^2,(1-alpha)^2))/exp(2*yhat_p2_2)
  delta_2<--1+((y-yhat_p1_2)^2*ifelse((y-yhat_p1_2)>0,alpha^2,(1-alpha)^2))/exp(2*yhat_p2_2)
  v11_x0<-sum(((y-X1_p1_2%*% b1_p1_2)*ifelse((y-X1_p1_2%*% b1_p1_2)>0,alpha^2,(1-alpha)^2))^2*K_h/exp(2*X2_p2_2%*% b2_p2_2))/sum(K_h)

  v22_x0<-sum((-1+(y-X1_p1_2%*% b1_p1_2)^2*ifelse((y-X1_p1_2%*% b1_p1_2)>0,alpha^2,(1-alpha)^2)/exp(2*X2_p2_2%*% b2_p2_2))^2*K_h)/sum(K_h)


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
#' @return The code \code{\link{locpolAND}} provides the realized value of the local maximum likelihood estimator of \eqn{\widehat{\theta}_{r0}(x_0)} for \eqn{(r\in \{1,2\})} with the estimated approximate asymptotic bias and variance at all \eqn{m} grind points \eqn{x_0} discussed in Gijbels et al. (2019c).
#'
#' @rdname SemiQRegAND
#'
#' @examples
#' \donttest{
#' data(LocomotorPerfor)
#' x=log(LocomotorPerfor$Body_Mass)
#' y=log(LocomotorPerfor$MRRS)
#' h_ROT =  0.9030372
#' locpolAND(x, y, p1=1,p2=1,h=h_ROT, alpha=0.50)
#' }
#' @export
locpolAND<-function(x, y,p1,p2, h, alpha,m = 101)
{

  xx <- seq(min(x)+.01*h, max(x)-.01*h, length = m)
  theta_10<-matrix(NA,length(xx),1)

  theta_20<-matrix(NA,length(xx),1)
  Bias_10<-matrix(NA,length(xx),1)
  Bias_20<-matrix(NA,length(xx),1)
  Var_10<-matrix(NA,length(xx),1)
  Var_20<-matrix(NA,length(xx),1)

  for (ii in 1:length(xx)) {

    fit_x0<-locpolAND_x0(x,y,p1,p2,h,alpha,x0=xx[ii])
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
#' @return The code \code{\link{SemiQRegAND}} provides the realized value of the \eqn{\beta}th conditional quantile estimator by using semiparametric quantile regression technique discussed in Gijbels et al. (2019b) and Gijbels et al. (2019c).
#' @import zipfR
#' @import locpol
#' @rdname SemiQRegAND
#' @examples
#' \donttest{
#' # Data
#' data(LocomotorPerfor)
#' x=log(LocomotorPerfor$Body_Mass)
#' y=log(LocomotorPerfor$MRRS)
#' h_ROT =  0.9030372
#' gridPoints=101
#' alpha= 0.5937
#' plot(x,y)
#' # location and scale functions estimation at the grid point x0
#' gridPoints=101
#' fit_AND <-locpolAND(x, y, p1=1,p2=1,h=h_ROT, alpha=alpha, m = gridPoints)
#' par(mgp=c(2,.4,0),mar=c(5,4,4,1)+0.01)
#'
#' # For phi plot
#' plot(fit_AND$x0,exp(fit_AND$theta_20),ylab=expression(widehat(phi)(x[0])),
#' xlab="log(Body mass)",type="l",font.lab=2,cex.lab=1.5,
#' bty="l",cex.axis=1.5,lwd =3)

#'
#' ## For theta2 plot
#' plot(fit_AND$x0,fit_AND$theta_20,ylab=expression(bold(widehat(theta[2]))(x[0])),
#' xlab="log(Body mass)",type="l",col=c(1), lty=1, font.lab=1,cex.lab=1.5,
#' bty="l",cex.axis=1.3,lwd =3)
#'
#'
#'
#### Estimated Quantile lines by ALaD
#' par(mgp=c(2.5, 1, 0),mar=c(5,4,4,1)+0.01)
#' # X11(width=7, height=7)
#' plot(x,y, ylim=c(0,4.5),xlab = "log(Body mass (kg))",
#' ylab = "log(Maximum relative running speed)",font.lab=1.5,
#' cex.lab=1.5,bty="l",pch=20,cex.axis=1.5)
#'
# alpha quantle line
#' lines(fit_AND$x0,fit_AND$theta_10, type='l',col=c(4),lty=6,lwd =3)
#' lines(fit_AND$x0,SemiQRegAND(beta=0.50,x, y,
#' p1=1,p2=1, h=h_ROT,alpha=alpha,m=gridPoints)$fit_beta_AND,
#' type='l',col=c(1),lty=5,lwd =3)
#' lines(fit_AND$x0,SemiQRegAND(beta=0.90,x, y,
#' p1=1,p2=1, h=h_ROT,alpha=alpha,m=gridPoints)$fit_beta_AND,type='l',col=c(14),lty=4,lwd =3)
#' lines(fit_AND$x0,SemiQRegAND(beta=0.10,x, y,
#' p1=1,p2=1, h=h_ROT,alpha=alpha,m=gridPoints)$fit_beta_AND,type='l',
#' col=c(19),lty=2,lwd =3)
#'
#' legend("topright", legend = c(expression(beta==0.10),
#'                              expression(beta==0.50), expression(beta==0.5937),
#'                              expression(beta==0.90)), col = c(19,1,4,14), lty=c(2,5,6,4),
#'       adj = c(.07, 0.5),, inset = c(0.05, +0.01), lwd = 3,cex=1.2)
#'
#' }
#' @export
SemiQRegAND<-function(beta,x, y, p1=1,p2=1,h,alpha=NULL, m = 101){
  if (is.null(alpha)){
    alpha_fn<-function(x,y){
      Newdata<-data.frame(x,y)
      fit_mean<-locpol::locpol(y~x, data=Newdata,weig=rep(1,nrow(Newdata)) ,kernel=gaussK,deg=1,xeval=NULL,xevalLen=101)
      MeanRegResidual<-fit_mean$residuals
      MeanRegResidual<-sort(MeanRegResidual)
      alpha<-mleAND(MeanRegResidual)$alpha
      list(alpha=alpha)
    }

    alpha<-alpha_fn(x,y)$alpha}
  fit_theta <-locpolAND(x, y,p1,p2, h, alpha,m)
  mu=fit_theta$theta_10
  phi=exp(fit_theta$theta_20)
  quanty<-matrix(NA,length(mu),1);
  for (i in 1:length(mu))
  {
    if (beta <alpha)
    {
      quanty[i]<-mu[i]-sqrt((2*phi[i]^2/(1-alpha)^2)*Igamma.inv(0.5, beta*sqrt(pi)/alpha,lower=FALSE));

    }
    else
    {
      quanty[i]<-mu[i]+sqrt((2*phi[i]^2/alpha^2)*Igamma.inv(0.5, sqrt(pi)*(beta-alpha)/(1-alpha)));
    }
  }
  C_alpha_beta<-qAND(beta,mu=0,phi=1,alpha)
  Bias_q_x0<-fit_theta$Bias_10+C_alpha_beta*fit_theta$Bias_20*exp(fit_theta$theta_20)
  Var_q_x0<-fit_theta$Var_10+(C_alpha_beta)^2*fit_theta$Var_20*exp(2*fit_theta$theta_20)*exp(fit_theta$Bias_20)
  list(x0=fit_theta$x0,fit_beta_AND=quanty,Bias_q_x0=Bias_q_x0,Var_q_x0=Var_q_x0,fit_alpha_AND=mu,alpha=alpha,beta=beta)
}

