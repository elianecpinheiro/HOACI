# Confidence limits for the parameter theta of an exponential 
# distribution with mean 1/theta based on n iid observations.
# Example 1 in
# Higher-order approximate confidence intervals
# Pinheiro, E.C, Ferrari, S.L.P., Medeiros, F.M.C.
# 2020
# Date: 07/07/2020


#-----------------------------------------------------------------------------------------
# Plot 
#-----------------------------------------------------------------------------------------

setEPS()
n=5 #sample size
xbar=1 #sample mean
alpha=0.01 #To construct gamma = 1-alpha confidence limits
k<-function(alpha){
  return( 1-(1/sqrt(n))*(qnorm(alpha)-(qnorm(alpha)^2-1)/(3*sqrt(n))+(qnorm(alpha)^3-7*qnorm(alpha))/(36*n)))
}
U<-function(theta){
  return(n*(1/theta-xbar))
}
Ulower<-function(theta){
  return(n*(1/theta*k(1-alpha/2)-xbar))
}
Uupper<-function(theta){
  return(n*(1/theta*k(alpha/2)-xbar))
}
UMedian<-function(theta){
  return(n*(1/theta*k(0.5)-xbar))
}

theta<-seq(0.01,4.5,0.01)

alpha=0.01 #To construct gamma = 1-alpha confidence limits
par(cex=1.5)
plot(theta,U(theta),type="l",col="blue",lwd=3,lty=1,xlim=c(0.1,4),ylim=c(-6,10),xlab=expression(theta),ylab="score")
abline(h=0)
lines(theta,Ulower(theta),type="l",col="red" ,lwd=3,lty=2)
lines(theta,Uupper(theta),type="l",col="red" ,lwd=3,lty=3)
lines(theta,UMedian(theta),type="l",col="blue" ,lwd=3,lty=4)
legend(bty="n",x=2,y=10.2,
       legend=c(expression(U(theta),tilde(U)[0.5](theta),tilde(U)[alpha/2](theta),tilde(U)[1-alpha/2](theta))),
       col=c("blue","blue","red","red"),lwd=c(3,3,3,3),lty=c(1,4,2,3),cex=0.8,y.intersp=1.5)

#-----------------------------------------------------------------------------------------
# Table 
#-----------------------------------------------------------------------------------------

# MLE
MLE <- 1/xbar
MLE

# 1 - alpha confidence limits

liminf<-0.01  
limsup<-3.5
step<-0.001
theta <- seq(liminf,limsup,step) 

# thetahat <- 1
ALPHA<-c(0.10,0.05,0.01) 
N<-c(3,5,7)
table<-matrix(,nrow=21,ncol=6)
for(i in 1:length(N)){
  for(j in 1:length(ALPHA)){
    n<-N[i]
    alpha<-ALPHA[j] #To construct gamma = 1-alpha confidence limits

    #MLE with asymp normality
    #
    LowerMLE <- MLE - qnorm(1-alpha/2)*MLE/sqrt(n)
    UpperMLE <- MLE + qnorm(1-alpha/2)*MLE/sqrt(n)
    #
    #MBR with asymp normality
    #
    # Median bias reduction
    MBR <- as.numeric(uniroot(UMedian, c(0.1,4), tol = 1e-15)[1])
    LowerMBR <- MBR - qnorm(1-alpha/2)*MBR/sqrt(n)
    UpperMBR <- MBR + qnorm(1-alpha/2)*MBR/sqrt(n)
    #
    #HOACI - solve modified score equation 
    #
    LowerHOACISo <- as.numeric(uniroot(Ulower, c(0.1,4), tol = 1e-15)[1])
    UpperHOACISo <- as.numeric(uniroot(Uupper, c(0.1,4), tol = 1e-15)[1])
    #
    #HOACI - analytic formula 
    #
    LowerHOACIAn <- MLE*k(1-alpha/2)
    UpperHOACIAn <- MLE*k(alpha/2)
    #
    #Exact - uses the exact distribution of 2*n*theta*xbar, a chi-squared dist with 2n df.
    #
    LowerExact<-qchisq(  alpha/2,2*n)/(2*n*xbar)
    UpperExact<-qchisq(1-alpha/2,2*n)/(2*n*xbar)
    #
    #r and r*
    #
    r <- (MLE-theta)/abs(MLE-theta)*sqrt(2*n*(log(MLE/theta) - 1 + theta/MLE ))
    Ustar <- sqrt(n)/MLE*abs((MLE-theta))
    rstar <- r + log(Ustar/abs(r))/r
    res <- cbind(theta,r,rstar)
    
    u <- qnorm(1-alpha/2)
    Lowerr<-NULL
    l<-1
    while(res[l,2]>u){
      Lowerr<-res[l,1]
      l<-l+1
    } 
    Lowerrstar<-NULL
    l<-1
    while(res[l,3]>u){
      Lowerrstar<-res[l,1]
      l<-l+1
    } 
    
    u <- -qnorm(1-alpha/2)
    Upperr<-NULL
    l<-length(theta)
    while(res[l,2]<u){
      Upperr<-res[l,1]
      l<-l-1
    } 
    Upperrstar<-NULL
    l<-length(theta)
    while(res[l,3]<u){
      Upperrstar<-res[l,1]
      l<-l-1
    }

    table[7*i-6,2*j-1]  <-LowerMLE
    table[7*i-6,2*j]    <-UpperMLE
    table[7*i-5,2*j-1]  <-LowerMBR
    table[7*i-5,2*j]    <-UpperMBR
    table[7*i-4,2*j-1]  <-LowerHOACISo
    table[7*i-4,2*j]    <-UpperHOACISo
    table[7*i-3,2*j-1]  <-LowerMLE
    table[7*i-3,2*j]    <-UpperMLE
    table[7*i-2,2*j-1]  <-Lowerr
    table[7*i-2,2*j]    <-Upperr
    table[7*i-1,2*j-1]  <-Lowerrstar
    table[7*i-1,2*j]    <-Upperrstar
    table[7*i  ,2*j-1]  <-LowerExact
    table[7*i  ,2*j]    <-UpperExact
  }
}
rownames(table)<-c("MLE     n=3","MBR     n=3","QBR     n=3","score   n=3","r       n=3","rstar   n=3","Exact   n=3",
                   "MLE     n=5","MBR     n=5","QBR     n=5","score   n=5","r       n=5","rstar   n=5","Exact   n=5",
                   "MLE     n=7","MBR     n=7","QBR     n=7","score   n=7","r       n=7","rstar   n=7","Exact   n=7")
colnames(table)<-c("90%","90%","95%","95%","99%","99%")
# 1 - alpha confidence limits
round(table,2)


