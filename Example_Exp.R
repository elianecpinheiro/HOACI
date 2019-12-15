# Confidence limits for the parameter theta of an exponential 
# distribution with mean 1/theta based on n iid observations.
# Example 1 in
# Higher-order approximate confidence intervals
# Pinheiro, E.C, Ferrari, S.L.P., Medeiros, F.M.C.
# 2019
# Date: 04/12/2019

setEPS()
n=5 #sample size
xbar=1 #sample mean
alpha=0.01 #To construct gamma = 1-alpha confidence limits
k<-function(alpha)
{
  return( 1-(1/sqrt(n))*(qnorm(alpha)-(qnorm(alpha)^2-1)/(3*sqrt(n))+(qnorm(alpha)^3-7*qnorm(alpha))/(36*n)))
}
U<-function(theta)
{
  return(n*(1/theta-xbar))
}
Ulower<-function(theta)
{
  return(n*(1/theta*k(1-alpha/2)-xbar))
}
Uupper<-function(theta)
{
  return(n*(1/theta*k(alpha/2)-xbar))
}
UMedian<-function(theta)
{
  return(n*(1/theta*k(0.5)-xbar))
}

theta<-seq(0.01,4.5,0.01)
plot(theta,U(theta),type="l",col="blue",lwd=3,lty=1,xlim=c(0.1,4),ylim=c(-6,10),xlab=expression(theta),ylab="score")
abline(h=0)
lines(theta,Ulower(theta),type="l",col="red" ,lwd=3,lty=2)
lines(theta,Uupper(theta),type="l",col="red" ,lwd=3,lty=3)
lines(theta,UMedian(theta),type="l",col="blue" ,lwd=3,lty=4)
legend(bty="n",x=2.5,y=10.2,legend=c(expression(U(theta),U[0.5](theta),U[alpha/2](theta),U[1-alpha/2](theta))),col=c("blue","blue","red","red"),lwd=c(3,3,3,3),lty=c(1,4,2,3),cex=1.3)

# MLE
MLE <- 1/xbar
MLE
# Median unbiased estimate
MUE <- uniroot(UMedian, c(0.1,3), tol = 1e-15)
MUE[1]

#1 - alpha confidence limits
#
#MLE with asymp normality
#
LowerMLE <- MLE - qnorm(1-alpha/2)*MLE/sqrt(n)
UpperMLE <- MLE + qnorm(1-alpha/2)*MLE/sqrt(n)
#
#MUE with asymp normality
#
LowerMUE <- as.numeric(MUE[1]) - qnorm(1-alpha/2)*as.numeric(MUE[1])/sqrt(n)
UpperMUE <- as.numeric(MUE[1]) + qnorm(1-alpha/2)*as.numeric(MUE[1])/sqrt(n)
#
#HOACI - solve modified score equation 
#
LowerHOACISo <- uniroot(Ulower, c(0.1,4), tol = 1e-15)
LowerHOACISo <- as.numeric(LowerHOACISo[1])
UpperHOACISo <- uniroot(Uupper, c(0.1,4), tol = 1e-15)
UpperHOACISo <- as.numeric(UpperHOACISo[1])
#
#HOACI - analytic formula 
#
LowerHOACIAn <- MLE*k(1-alpha/2)
UpperHOACIAn <- MLE*k(alpha/2)
#
#Exact - uses the exact distribution of 2*n*theta*xbar, a chi-squared dist with 2n df.
#
LowerExact<-qchisq(alpha/2,2*n)/(2*n*xbar)
UpperExact<-qchisq(1-alpha/2,2*n)/(2*n*xbar)

print(c("n=",n,"alpha=",alpha))
print(c("LowerMLE=",LowerMLE,"UpperMLE=",UpperMLE))
print(c("LowerMUE=",LowerMUE,"UpperMUE=",UpperMUE))
print(c("LowerHOACISo=",LowerHOACISo,"UpperHOACISo=",UpperHOACISo))
print(c("LowerHOACIAn=",LowerHOACIAn,"UpperHOACIAn=",UpperHOACIAn))
print(c("LowerExact=",LowerExact,"UpperExact=",UpperExact))

