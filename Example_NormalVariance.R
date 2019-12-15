# Confidence limits for the variance (theta) of a normal 
# distribution with known mean based on n iid observations.
# Example 2 in
# Higher-order approximate confidence intervals
# Pinheiro, E.C, Ferrari, S.L.P., Medeiros, F.M.C.
# 2019
# Date: 04/12/2019

n=20 #sample size
thetahat=1 #sum(x^2)/n
alpha=0.05 #To construct gamma = 1-alpha confidence limits

k<-function(alpha)
{
  return( sqrt(2)*( qnorm(alpha)+(sqrt(2)*(qnorm(alpha)^2-1))/(3*sqrt(n))+(qnorm(alpha)^3-7*qnorm(alpha))/(18*n) ) )
}
U<-function(theta)
{
  return( (n/(2*theta)*(thetahat/theta - 1)) )
}
Ulower<-function(theta)
{
  return( (n/(2*theta))*( thetahat/theta - 1 - k(1-alpha/2)/sqrt(n) ) )
}
Uupper<-function(theta)
{
  return( (n/(2*theta))*( thetahat/theta - 1 - k(alpha/2)/sqrt(n) ) )
}
UMedian<-function(theta)
{
  return( (n/(2*theta))*( thetahat/theta - 1 - k(1/2)/sqrt(n) ) )
}

theta<-seq(0.3,5,0.01)
plot(theta,U(theta),type="l",col="blue",lwd=3,lty=1,xlim=c(0.3,5),ylim=c(-8.5,6),xlab=expression(theta),ylab="score")
abline(h=0)
lines(theta,Ulower(theta),type="l",col="red" ,lwd=3,lty=2)
lines(theta,Uupper(theta),type="l",col="red" ,lwd=3,lty=3)
lines(theta,UMedian(theta),type="l",col="blue" ,lwd=3,lty=4)
legend(bty="n",x=3.2,y=7,legend=c(expression(U(theta),U[0.5](theta),U[alpha/2](theta),U[1-alpha/2](theta))),col=c("blue","blue","red","red"),lwd=c(3,3,3,3),lty=c(1,4,2,3),cex=1.3)

# MLE
MLE <- thetahat
MLE
# Median unbiased estimate
MUE <- MLE/(1+k(0.5)/sqrt(n))
MUE
MUE <- uniroot(UMedian, c(0.1,8), tol = 1e-15)
MUE[1]

#1 - alpha confidence limits
#
#MLE + asymp normality
#
LowerMLE <- MLE - qnorm(1-alpha/2)*MLE*sqrt(2)/sqrt(n)
LowerMLE
UpperMLE <- MLE + qnorm(1-alpha/2)*MLE*sqrt(2)/sqrt(n)
UpperMLE
#
#MUE with asymp normality
#
LowerMUE <- as.numeric(MUE[1]) - qnorm(1-alpha/2)*as.numeric(MUE[1])/sqrt(n)
LowerMUE
UpperMUE <- as.numeric(MUE[1]) + qnorm(1-alpha/2)*as.numeric(MUE[1])/sqrt(n)
UpperMUE
#
#HOACI 
#
LowerHOACI <- uniroot(Ulower, c(0.1,20), tol = 1e-15)
LowerHOACI[1]
UpperHOACI <- uniroot(Uupper, c(0.1,8), tol = 1e-15)
UpperHOACI[1]
#
#Exact - uses the exact distribution of n*thetahat/theta, a chi-squared dist with n df.
#
LowerExact<- n*thetahat/qchisq(1-alpha/2,n)
LowerExact
UpperExact<- n*thetahat/qchisq(alpha/2,n)
UpperExact

