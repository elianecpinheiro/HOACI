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
# Date: 28/07/2020

library(sn)
library(nleqslv)
library(e1071)      ## to also calculate the skewness
library(stringr)

#####################################################################################
#                   negative log-likelihood                  
#####################################################################################

nlogL = function(theta,data)
  -sum(dsn(data,alpha=theta,log=TRUE))


#####################################################################################
#                  numerical expectation of  y^kz(,)^h                   
#####################################################################################

a.kh<-function(theta,k,h){
  integrate(function(x) 2*dnorm(x)*pnorm(theta*x)*x^k*zeta(1,theta*x)^h,-Inf,Inf)$value
}

#####################################################################################
#                      median unbiased modication 
#####################################################################################

modification.mc <- function(para){
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

theta.start = function(y){
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
n = length(data)

#Modified data (note the change of sign of obs 4)
#data <- c(0.195, 0.847, 0.726, 0.139, 1.788, 0.570, 2.069, 0.452, 0.868, 
#          1.199, 0.894, 0.887, 1.258, 0.918, 0.496, 0.183, 0.119, 1.207,
#          0.446, 2.579)
#n = length(data)


#####################################################################################
#                 Main Function 
#####################################################################################
skew.sim=function(data, alpha){
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

alpha = 0.05
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
################ confidence limits - MLE ###################

if(res$mle==Inf){
  li.mv = str_c("NaN")
  ls.mv = str_c("NaN")
  int.mv = data.frame(li.mv, ls.mv)
}else{
  a22_mv = a.kh(res$mle,2,2)
  li.mv = res$mle + qnorm(alpha/2)*sqrt(1/(n*a22_mv))
  ls.mv = res$mle + qnorm(1-alpha/2)*sqrt(1/(n*a22_mv))
  int.mv = data.frame(li.mv, ls.mv)
}


############# confidence limits - MBR ###################

a22_med = a.kh(res$mle.mc,2,2)
li.med = res$mle.mc + qnorm(alpha/2)*sqrt(1/(n*a22_med))
ls.med = res$mle.mc + qnorm(1-alpha/2)*sqrt(1/(n*a22_med))
int.med = data.frame(li.med, ls.med)


############# confidence limits - QBR ###################

liqbr = res$qbr.lower
lsqbr = res$qbr.upper
int.qbr = c(liqbr, lsqbr)

#############  summary of results ###################
ic      = c(int.mv,int.med,int.qbr)
int.res = data.frame(ic)
if(res$mle==Inf|res$mle==-Inf){
colnames(int.res) = c("Lower.ML ", "Upper.ML","Lower.MBR", "Upper.MBR", "Lower.QBR", "Upper.QBR")
rownames(int.res) = expression(theta)
cat("\n=======================================================================")
cat(paste("\n                          ",(1-alpha)*100,"% CI - Two-sided\n"))
cat("=======================================================================")
cat("\n")
print(int.res)
}else
{
  colnames(int.res) = c("Lower.ML ", "Upper.ML","Lower.MBR", "Upper.MBR", "Lower.QBR", "Upper.QBR")
  rownames(int.res) = expression(theta)
  cat("\n=======================================================================")
  cat(paste("\n                          ",(1-alpha)*100,"% CI - Two-sided\n"))
  cat("=======================================================================")
  cat("\n")
  print(int.res)
}

#__________________________________________________________________
#      SCORE FUNCTION PLOTS: Sartori's data 
#__________________________________________________________________
x <- seq(0,20,0.01)
u <- umed <- ulower <- uupper <- NULL     
for(i in 1:length(x)){
  u[i]      = score(x[i],  res.sn$data)
  umed[i]   = score.mc(x[i], res.sn$data)
  ulower[i] = score.qbr(x[i], res.sn$data, qnorm(1-alpha/2))
  uupper[i] = score.qbr(x[i], res.sn$data, qnorm(alpha/2))
} 

par(cex=1.5)
plot(NA,NA,type="l",col="blue",lwd=3,lty=1,xlim=c(0,20),ylim=c(-1,2),
     xlab=expression(theta),ylab="score",axes=F, frame.plot = TRUE)
abline(h=0)
lines(x,u,     type="l",col="blue",lwd=3,lty=1)
lines(x,umed,  type="l",col="blue",lwd=3,lty=4)
lines(x,ulower,type="l",col="red" ,lwd=3,lty=2)
lines(x,uupper,type="l",col="red" ,lwd=3,lty=3)
legend(bty="n",x=10,y=2.2, 
       legend=c(expression(U(theta),tilde(U)[0.5](theta),tilde(U)[alpha/2](theta),tilde(U)[1-alpha/2](theta))),
       col=c("blue","blue","red","red"), lwd=c(3,3,3,3), lty=c(1,4,2,3),cex=0.8,y.intersp=1.5)
y=seq(-1, 2, by = 1)
my.aty <- round(y,1)
axis(2, at = my.aty, labels = my.aty)
x=seq(0, 20, by = 5)
my.atx <- round(x,1)
axis(1, at = my.atx, labels = my.atx)


