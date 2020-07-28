# Confidence limits for the variance (theta) of a normal 
# distribution with known mean based on n iid observations.
# Example 2 in
# Higher-order approximate confidence intervals
# Pinheiro, E.C, Ferrari, S.L.P., Medeiros, F.M.C.
# 2020
# Date: 15/07/2020

#-----------------------------------------------------------------------------------------
# Plot 

k<-function(alpha){
  return( sqrt(2)*( qnorm(alpha)+(sqrt(2)*(qnorm(alpha)^2-1))/(3*sqrt(n))+(qnorm(alpha)^3-7*qnorm(alpha))/(18*n) ) )
}
U<-function(theta){
  return( (n/(2*theta)*(thetahat/theta - 1)) )
}
Ulower<-function(theta){
  return( (n/(2*theta))*( thetahat/theta - 1 - k(1-alpha/2)/sqrt(n) ) )
}
Uupper<-function(theta){
  return( (n/(2*theta))*( thetahat/theta - 1 - k(alpha/2)/sqrt(n) ) )
}
UMedian<-function(theta){
  return( (n/(2*theta))*( thetahat/theta - 1 - k(1/2)/sqrt(n) ) )
}
Urlower<-function(theta){
  r <- (1-theta)/abs(1-theta)*sqrt(n*(-log(MLE/theta)-1+MLE/theta))
  return( r + qnorm(  alpha/2) )
}
Urupper<-function(theta){
  r <- (1-theta)/abs(1-theta)*sqrt(n*(-log(MLE/theta)-1+MLE/theta))
  return(  r + qnorm(1-alpha/2) )
}
Urstarlower<-function(theta){
  r <- (1-theta)/abs(1-theta)*sqrt(n*(-log(MLE/theta)-1+MLE/theta))
  Ustar <- sqrt(n/2)*abs(-1+MLE/theta)
  rstar <- r + log(Ustar/abs(r))/r
  return( rstar + qnorm(  alpha/2) )
}
Urstarupper<-function(theta){
  r <- (1-theta)/abs(1-theta)*sqrt(n*(-log(MLE/theta)-1+MLE/theta))
  Ustar <- sqrt(n/2)*abs(-1+MLE/theta)
  rstar <- r + log(Ustar/abs(r))/r
  return(  rstar + qnorm(1-alpha/2) )
}
UScorelower<-function(theta){
  return( sqrt(n/2)*(-1+MLE/theta) + qnorm(  alpha/2) )
}
UScoreupper<-function(theta){
  return(  sqrt(n/2)*(-1+MLE/theta) + qnorm(1-alpha/2) )
}

n<-15 #sample size
thetahat<-1 #sum(x^2)/n

theta<-seq(0.3,5,0.01)

alpha<-0.01 #To construct gamma = 1-alpha confidence limits
par(cex=1.5)
plot(theta,U(theta),type="l",col="blue",lwd=3,lty=1,xlim=c(0.3,5),ylim=c(-10,11),xlab=expression(theta),ylab="score")
abline(h=0)
lines(theta,Ulower(theta),type="l",col="red" ,lwd=3,lty=2)
lines(theta,Uupper(theta),type="l",col="red" ,lwd=3,lty=3)
lines(theta,UMedian(theta),type="l",col="blue" ,lwd=3,lty=4)
legend(bty="n",x=2.5,y=13,
       legend=c(expression(U(theta),tilde(U)[0.5](theta),tilde(U)[alpha/2](theta),tilde(U)[1-alpha/2](theta))),
       col=c("blue","blue","red","red"),lwd=c(3,3,3,3),lty=c(1,4,2,3),cex=0.8,y.intersp=1.5)

#-----------------------------------------------------------------------------------------
# Table 

# MLE
MLE <- thetahat
MLE

# 1 - alpha confidence limits

liminf<-0.1  
limsup<-20
step<-0.001
theta <- seq(liminf,limsup,step) 

ALPHA<-c(0.10,0.05,0.01)
qnorm(1-ALPHA/2)
N<-c(10,15,20)

table<-matrix(,nrow=7*length(N),ncol=6)
for(i in 1:length(N)){
  for(j in 1:length(ALPHA)){
    n<-N[i]
    alpha<-ALPHA[j] #To construct gamma = 1-alpha confidence limits

    #
    #MLE with asymp normality
    #
    LowerMLE <- MLE - qnorm(1-alpha/2)*MLE*sqrt(2)/sqrt(n)
    UpperMLE <- MLE + qnorm(1-alpha/2)*MLE*sqrt(2)/sqrt(n)
    #
    #MBR with asymp normality
    #
    # Median bias reduction
    MBR <- MLE/(1+k(0.5)/sqrt(n))
    LowerMBR <- MBR - qnorm(1-alpha/2)*MBR*sqrt(2)/sqrt(n)
    UpperMBR <- MBR + qnorm(1-alpha/2)*MBR*sqrt(2)/sqrt(n)
    #
    #Exact - uses the exact distribution of n*thetahat/theta, a chi-squared dist with n df.
    #
    LowerExact<- n*MLE/qchisq(1-alpha/2,n)
    UpperExact<- n*MLE/qchisq(alpha/2,n)
    #
    #HOACI 
    #
    options(show.error.messages = FALSE)
    
    ans<-try(uniroot(Ulower, c(0.1,1e+15), tol = 1e-15))
    if(inherits(ans, "try-error")){
      #error handling code, maybe just skip this iteration using
      LowerHOACI<-Inf
    }
    else{
      LowerHOACI <- as.numeric(ans[1])
    }
    ans<-try(uniroot(Uupper, c(0.1,1e+15), tol = 1e-15))
    if(inherits(ans, "try-error")){
      #error handling code, maybe just skip this iteration using
      UpperHOACI<-Inf
    }
    else{
      UpperHOACI <- as.numeric(ans[1])
    }
    #
    # r
    #
    ans<-try(uniroot(Urlower, c(0.1,1e+15), tol = 1e-15))
    if(inherits(ans, "try-error")){
      #error handling code, maybe just skip this iteration using
      Lowerr<-Inf
    }
    else{
      Lowerr <- as.numeric(ans[1])
    }
    ans<-try(uniroot(Urupper, c(0.1,1e+15), tol = 1e-15))
    if(inherits(ans, "try-error")){
      #error handling code, maybe just skip this iteration using
      Upperr<-Inf
    }
    else{
      Upperr <- as.numeric(ans[1])
    }
    #
    # r*
    #
    ans<-try(uniroot(Urstarlower, c(0.1,1e+15), tol = 1e-15))
    if(inherits(ans, "try-error")){
      #error handling code, maybe just skip this iteration using
      Lowerrstar<-Inf
    }
    else{
      Lowerrstar <- as.numeric(ans[1])
    }
    ans<-try(uniroot(Urstarupper, c(0.1,1e+15), tol = 1e-15))
    if(inherits(ans, "try-error")){
      #error handling code, maybe just skip this iteration using
      Upperrstar<-Inf
    }
    else{
      Upperrstar <- as.numeric(ans[1])
    }
    #
    #score
    #
    ans<-try(uniroot(UScorelower, c(0.1,1e+15), tol = 1e-15))
    if(inherits(ans, "try-error")){
      #error handling code, maybe just skip this iteration using
      LowerScore<-Inf
    }
    else{
      LowerScore <- as.numeric(ans[1])
    }
    ans<-try(uniroot(UScoreupper, c(0.1,1e+15), tol = 1e-15))
    if(inherits(ans, "try-error")){
      #error handling code, maybe just skip this iteration using
      UpperScore<-Inf
    }
    else{
      UpperScore <- as.numeric(ans[1])
    }
    
    
    table[7*i-6,2*j-1]  <-LowerMLE
    table[7*i-6,2*j]    <-UpperMLE
    table[7*i-5,2*j-1]  <-LowerMBR
    table[7*i-5,2*j]    <-UpperMBR
    table[7*i-4,2*j-1]  <-LowerHOACI
    table[7*i-4,2*j]    <-UpperHOACI
    table[7*i-3,2*j-1]  <-LowerScore
    table[7*i-3,2*j]    <-UpperScore
    table[7*i-2,2*j-1]  <-Lowerr
    table[7*i-2,2*j]    <-Upperr
    table[7*i-1,2*j-1]  <-Lowerrstar
    table[7*i-1,2*j]    <-Upperrstar
    table[7*i  ,2*j-1]  <-LowerExact
    table[7*i  ,2*j]    <-UpperExact
    
  }
}
rownames<-c()
for(i in 1:length(N)){
  rownames<-c(rownames,
              paste("MLE       n=",N[i], sep = "", collapse = NULL),
              paste("MBR       n=",N[i], sep = "", collapse = NULL),
              paste("QBR       n=",N[i], sep = "", collapse = NULL),
              paste("score     n=",N[i], sep = "", collapse = NULL),
              paste("r         n=",N[i], sep = "", collapse = NULL),
              paste("rstar     n=",N[i], sep = "", collapse = NULL),
              paste("Exact     n=",N[i], sep = "", collapse = NULL)
  )
}
rownames(table)<-rownames
colnames(table)<-c("90%","90%","95%","95%","99%","99%")
round(table,2)
