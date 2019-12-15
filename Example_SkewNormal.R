# Confidence limits for the skewness parameter theta of a 
# standard skew normal distribution based on n iid observations.
# Example 3 in
# Higher-order approximate confidence intervals
# Pinheiro, E.C, Ferrari, S.L.P., Medeiros, F.M.C.
# 2019
# Data from Example 1 of Sartori (Journal of Statistical Planning and Inference, 2006),
# and simulated data with skewness parameter = 1
# Based on the code for evaluating the median bias-corrected estimate
# by Euloge Clovis Kenne Pagui.
# Date: 04/12/2019

library(sn)
library(nleqslv)
library(e1071)      ## to also calculate the skewness

#####################################################################################
#                   negative log-likelihood                  
#####################################################################################

nlogL = function(theta,data)
  -sum(dsn(data,alpha=theta,log=TRUE))


#####################################################################################
#                  numerical expectation of  y^kz(,)^h                   
#####################################################################################

a.kh<-function(theta,k,h)
{
  integrate(function(x) 2*dnorm(x)*pnorm(theta*x)*x^k*zeta(1,theta*x)^h,-Inf,Inf)$value
}

#####################################################################################
#                      median unbiased modication 
#####################################################################################

modification.mc <- function(para)
{
  a.kh(para,3,3)/(a.kh(para,2,2)*6)
}


#####################################################################################
#                    score function                        
#####################################################################################

score = function(theta,data){
  n = length(data)
  res <- sum(sapply(1:n, function(x) data[x]*zeta(1,theta*data[x]) ) ) 
  return(res)
}

#####################################################################################
#                    score median corrected                  
#####################################################################################

score.mc = function(theta,data){
  n = length(data)
  res <- score(theta,data) + modification.mc(theta) 
  return(res)
}

#####################################################################################
#                    score qbr corrected                  
#####################################################################################

score.qbr = function(theta,data, up){
  n = length(data)
  a22 = a.kh(theta,2,2)
  a33 = a.kh(theta,3,3)
  a44 = a.kh(theta,4,4)
  k2 = n*a22
  k3 = n*a33
  k4 = n*a44 - 3*n*a22^2
  
  M  = - sqrt(k2)*up - ( k3/(6*k2) )*(up^2-1) - 
    ( k4/(24*k2^(3/2)) )*( up^3-3*up ) +
    ( k3^2/(36*k2^(5/2)) )*(2*up^3 - 5*up)
  res <- score(theta,data) + M 
  return(res)
}

#####################################################################################
#                    Initial value                    
#####################################################################################

theta.start = function(y)
{
  n = length(y)
  gamma     = sum((y-mean(y))/n*var(y)^3)
  abs.gamma = min(0.9952717, abs(gamma))
  r = (2*abs.gamma/(4-pi))^(1/3)
  if(gamma < 0)
  {
    thetahat = - r/sqrt(2/pi - (1-2/pi)*r^2)
  }
  else
  {
    thetahat = r/sqrt(2/pi - (1-2/pi)*r^2)
  }   
  thetahat
}   

#####################################################################################
#                 Data 
#####################################################################################

#Sartori's data

#Original data
data <- c(0.195, 0.847, 0.726, -0.139, 1.788, 0.570, 2.069, 0.452, 0.868, 
          1.199, 0.894, 0.887, 1.258, 0.918, 0.496, 0.183, 0.119, 1.207,
          0.446, 2.579)

#Modified data (note the change of sign of obs 4)
data <- c(0.195, 0.847, 0.726, 0.139, 1.788, 0.570, 2.069, 0.452, 0.868, 
          1.199, 0.894, 0.887, 1.258, 0.918, 0.496, 0.183, 0.119, 1.207,
          0.446, 2.579)

#Simulated data
n     = 20
set.seed(2019, kind = NULL, normal.kind = NULL)
data = rsn(n,alpha=1)

#####################################################################################
#                 Main Function 
#####################################################################################

skew.sim=function(data, alpha)
{
  est=matrix(0, nrow=1, ncol=4)
  n = length(data)  
  xstart <- theta.start(data)
  mle=nlminb(xstart,nlogL,data=data)
  conv=mle$convergence
  up = qnorm(alpha/2)
  mle.upper = nleqslv(xstart, score.qbr, data=data, up=up)
  convupper=mle.upper$termcd
  up = qnorm(1-alpha/2)
  mle.lower = nleqslv(xstart, score.qbr, data=data, up=up)
  convlower=mle.lower$termcd 
  if(-mle$objective < -nlogL(Inf,data)) 
  {
    mle$objective= nlogL(Inf,data=data)
    mle$par= Inf
    mle.upper$x = Inf
  }    
  if(-mle$objective < -nlogL(-Inf,data=data)) 
  {
    mle$objective= nlogL(-Inf,data=data)
    mle$par= -Inf
    mle.lower$x = -Inf
  }
  mle.mc = nleqslv(xstart,score.mc,data=data)
  convmc=mle.mc$termcd
  
  est[1,1]=mle$par       ## maximum likelihood
  est[1,2]=mle.mc$x      ## Median bias reduction
  est[1,3]=mle.lower$x   ## qbr lower
  est[1,4]=mle.upper$x   ## qbr lower
  
  colnames(est)=c("mle","mle.mc", "qbr.lower", "qbr.upper")
  list(est=est,conv=conv,convmc=convmc,convlower=convlower,convupper=convupper,n=n,data=data)
}


#####################################################################################
#                         Results
#####################################################################################

alpha = 0.1
res.sn = skew.sim(data, alpha)

mle_notinf = 100 - length(res.sn$count)/100
mle_notinf

convergenceQBR <- data.frame(cbind(res.sn$convlower,res.sn$convlower))
convQBR         <- (convergenceQBR[][1]==1&convergenceQBR[][2]==1)
conv_QBR        <- length(convQBR[convQBR==TRUE])
conv_QBR


#----------------------------------------------------------------------------
#                         confidence limits
#----------------------------------------------------------------------------

res = data.frame(res.sn$est)
res
################ confidence limits - MLE ###################

if(res$mle==Inf){
  li.mv = -Inf
  ls.mv = Inf
  int.mv = data.frame(li.mv, ls.mv)
}else{
  a22_mv = a.kh(res$mle,2,2)
  li.mv = res$mle + qnorm(alpha/2)*sqrt(1/(n*a22_mv))
  ls.mv = res$mle + qnorm(1-alpha/2)*sqrt(1/(n*a22_mv))
  int.mv = data.frame(li.mv, ls.mv)
}
int.mv

############# confidence limits - median ###################

a22_med = a.kh(res$mle.mc,2,2)
li.med = res$mle.mc + qnorm(alpha/2)*sqrt(1/(n*a22_med))
ls.med = res$mle.mc + qnorm(1-alpha/2)*sqrt(1/(n*a22_med))
int.med = data.frame(li.med, ls.med)
int.med


############# confidence limits - qbr ###################

liqbr = res$qbr.lower
lsqbr = res$qbr.upper
int.qbr = c(liqbr, lsqbr)

#############  summary of results ###################

ic      = c(int.mv,int.med,int.qbr)
int.res = matrix(ic, ncol=6, byrow=T)
colnames(int.res) = c("Lower.ML ", "Upper.ML","Lower.M", "Upper.M", "Lower.QBR", "Upper.QBR")
int.res

#__________________________________________________________________
#      SCORE FUNCTION PLOTS: Sartori's data 
#__________________________________________________________________
x <- seq(0,15,0.01)
u <- umed <- ulower <- uupper <- NULL     
for(i in 1:length(x)){
  u[i]      = score(x[i],  res.sn$data)
  umed[i]   = score.mc(x[i], res.sn$data)
  ulower[i] = score.qbr(x[i], res.sn$data, qnorm(1-alpha/2))
  uupper[i] = score.qbr(x[i], res.sn$data, qnorm(alpha/2))
} 

plot(x,u,type="l",col="blue",lwd=3,lty=1,xlim=c(0,15),ylim=c(-0.5,5),xlab=expression(theta),ylab="score")
abline(h=0)
lines(x,umed,type="l",col="blue" ,lwd=3,lty=4)
lines(x,ulower,type="l",col="red" ,lwd=3,lty=2)
lines(x,uupper,type="l",col="red" ,lwd=3,lty=3)
legend("topright", c(expression(U(theta)), expression(U[0.5](theta)), 
                     expression(U[alpha/2](theta)), expression(U[1-alpha/2](theta))), 
       bty="n", lty=c(1,4,2,3), col=c("blue","blue","red","red") 
)

#__________________________________________________________________
#      SCORE FUNCTION PLOTS: Simulated data
#__________________________________________________________________
x <- seq(0,5,0.001)
u <- umed <- ulower <- uupper <- NULL     
for(i in 1:length(x)){
  u[i]      = score(x[i],  res.sn$data)
  umed[i]   = score.mc(x[i], res.sn$data)
  ulower[i] = score.qbr(x[i], res.sn$data, qnorm(1-alpha/2))
  uupper[i] = score.qbr(x[i], res.sn$data, qnorm(alpha/2))
} 

plot(x,u,type="l",col="blue",lwd=3,lty=1,xlim=c(0,3),ylim=c(-4,4),xlab=expression(theta),ylab="score")
abline(h=0)
lines(x,umed,type="l",col="blue" ,lwd=3,lty=4)
lines(x,ulower,type="l",col="red" ,lwd=3,lty=2)
lines(x,uupper,type="l",col="red" ,lwd=3,lty=3)
legend("topright", c(expression(U(theta)), expression(U[0.5](theta)), 
                     expression(U[alpha/2](theta)), expression(U[1-alpha/2](theta))), 
       bty="n", lty=c(1,4,2,3), col=c("blue","blue","red","red") 
)

