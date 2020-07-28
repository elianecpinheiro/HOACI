# Confidence limits for the skewness parameter theta of a 
# standard skew normal distribution based on n iid observations.
# Example 3 in
# Higher-order approximate confidence intervals
# Pinheiro, E.C, Ferrari, S.L.P., Medeiros, F.M.C.
# 2020
# Example 3: Simulation
# Date: 28/07/2020


library(sn)
library(nleqslv)
library(e1071)      ##to also calculate the skewness

#####################################################################################
#                   negative log-likelihood                  
#####################################################################################

nlogL = function(theta,data)
  -sum(dsn(data,alpha=theta,log=TRUE))


#####################################################################################
#                  numerical expectation of  y^kz(,)^h                   
#####################################################################################

a.kh<-function(theta,k,h){
  integrate(function(x) 2*dnorm(x)*pnorm(theta*x)*x^k*zeta(1,theta*x)^h,-Inf,Inf,  
            rel.tol = .Machine$double.eps^0.8, abs.tol = .Machine$double.eps^0.8)$value
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
  gamma3 = mean(((y-mean(y))/sqrt(var(y)))^3)
  abs.gamma3 = min(0.99, abs(mean(((y-mean(y))/sqrt(var(y)))^3)))
  abs.delta = sqrt((pi/2)*(abs.gamma3^(2/3)/(abs.gamma3^(2/3)+((4-pi)/2)^(2/3))))
  if(gamma3 < 0)
  {
    thetahat = - abs.delta/sqrt(1+ abs.delta^2)
  }
  else
  {
    thetahat =  abs.delta/sqrt(1+ abs.delta^2)
  }   
  thetahat
}   


#####################################################################################
#                 Main Function 
#####################################################################################
  
 skew.sim = function(alpha){
  est=matrix(0, nrow=Nsim, ncol=4)
  
  conv <-  convmc <- convlower <- convupper <- rep(NA,Nsim)
  count = NULL
  
  for (i in 1:Nsim)
  {
    data <- datasim[i,]
    xstart <- theta.start(data)
    mle = nlminb(xstart,nlogL,data=data)
    conv[i]=mle$convergence
    if (-mle$objective < -nlogL(Inf,data))
    {
      # print(c(-mle$objective,-nlogL(Inf,data),-nlogL(theta,data)))
      mle$objective= nlogL(Inf,data=data)
      mle$par= Inf
    }
    if (-mle$objective < -nlogL(-Inf,data=data))
    {
      mle$objective= nlogL(-Inf,data=data)
      mle$par= -Inf
    }
    if(mle$par== Inf|mle$par == -Inf)
    {
      count = c(count,i)
    }
    
    up = 0
    mle.mc = nleqslv(xstart, score.qbr, data=data, up=up)
    convmc[i]=mle.mc$termcd
    
    if(mle$par== Inf | mle$par == -Inf){
      up = qnorm(1-alpha/2)
      mle.lower = nleqslv(0.5, score.qbr, data=data, up=up)
      convlower[i]=mle.lower$termcd    
    
      up = qnorm(alpha/2)
      mle.upper = nleqslv(0.5, score.qbr, data=data, up=up)
      convupper[i]=mle.upper$termcd     
   }else{
    
    up = qnorm(1-alpha/2)
    mle.lower = nleqslv(xstart, score.qbr, data=data, up=up)
    convlower[i]=mle.lower$termcd    
    
    up = qnorm(alpha/2)
    mle.upper = nleqslv(xstart, score.qbr, data=data, up=up)
    convupper[i]=mle.upper$termcd       
   }
  
    est[i,1]=mle$par       ## maximum likelihood
    est[i,2]=mle.mc$x      ## Median bias reduction
    est[i,3]=mle.lower$x   ## qbr lower
    est[i,4]=mle.upper$x   ## qbr upper
    
  }
  colnames(est)=c("mle", "mle.mc", "qbr.lower", "qbr.upper")
  list(datasim=datasim, est=est,conv=conv,convmc=convmc,convlower=convlower,convupper=convupper,count=count,theta=theta,n=n)
}

#####################################################################################
#                                  Simulation
#####################################################################################
 
 n     = 20
 theta = 10
 Nsim  = 10000
 
 set.seed(2020, kind = NULL, normal.kind = NULL)
 datasim <- matrix(rsn(n*Nsim,alpha=theta),nrow=Nsim,ncol=n, byrow=T)
 
 t0<-Sys.time()
 
 alpha = 0.05
 simul.sn = skew.sim(alpha)
 
 t1<-Sys.time()
 t1-t0
 
 results <- data.frame(simul.sn$est, simul.sn$conv, simul.sn$convmc, simul.sn$convlower, simul.sn$convupper)
 colnames(results) <- c("mle", "mle.mc", "qbr.lower", "qbr.upper",
                        "coverage.mle", "coverage.mbr", "coverage.lower", "coverage.upper")
 
 ###################################################################
 
 #------------------------------------------------------------------------------------------
 #                   Infinite MLE
 #------------------------------------------------------------------------------------------
 
 results.mle.inf = results[results$mle==Inf|results$mle==-Inf,]
 qbr.comp1 = results.mle.inf[results.mle.inf$qbr.lower < results.mle.inf$qbr.upper,]
 res.pro1 = round(rbind(nrow(results.mle.inf), nrow(qbr.comp1)),4)
 
 
 #------------------------------------------------------------------------------------------
 #                   finite MLE
 #------------------------------------------------------------------------------------------
 results.mle.fin = results[results$mle!=Inf & results$mle!=-Inf,]
 qbr.comp2 = results.mle.fin[results.mle.fin$qbr.lower<results.mle.fin$qbr.upper,]
 res.pro2 = round(rbind(nrow(results.mle.fin), nrow(qbr.comp2)),4)
 
 
 res.sum = cbind(res.pro1, res.pro2, rbind(res.pro1[1,]+res.pro2[1,],res.pro1[2,]+res.pro2[2,]))
 colnames(res.sum) = c("%MLE = Inf", " %MLE < Inf", "Total")
 rownames(res.sum) = c("MLE", "QBR computed")
 
 #####################################################################################
 #                         Results - coverage probability
 #####################################################################################
 
 # ----------------- coverage probability - MLE ------------------------
 
 res.mle = data.frame(results$mle[!is.infinite(results$mle)])
 
 ic.mv <- function(Dresult){
   a22_mv = a.kh(Dresult,2,2)
   
   li = Dresult + qnorm(alpha/2)*sqrt(1/(n*a22_mv))
   ls = Dresult + qnorm(1-alpha/2)*sqrt(1/(n*a22_mv))
   return(c(li, ls))
 }
 
 int.mv = data.frame(t(apply(res.mle, 1, ic.mv)))
 names(int.mv) = c("li", "ls")
 
 cp_bi_mv  = int.mv$li<=theta & theta<=int.mv$ls
 cp_bi_mv = length(cp_bi_mv[cp_bi_mv==TRUE])/nrow(res.mle)
 
 #----------------- coverage probability - MBR ------------------------
 
 res.med = data.frame(results$mle.mc)
 
 ic.med <- function(Dresult){
   a22_med = a.kh(Dresult,2,2)
   
   li = Dresult + qnorm(alpha/2)*sqrt(1/(n*a22_med))
   ls = Dresult + qnorm(1-alpha/2)*sqrt(1/(n*a22_med))
   return(c(li, ls))
 } 
 
 int.med = data.frame(t(apply(res.med, 1, ic.med)))
 names(int.med) = c("li", "ls")
 
 cp_bi_med  = int.med$li<=theta & theta<=int.med$ls
 cp_bi_med = length(cp_bi_med[cp_bi_med==TRUE])/nrow(res.med)
 
 #----------------- coverage probability - QBR ------------------------
 
 res.qbr    = results[results$qbr.lower < results$qbr.upper,]
 cp_bi_qbr  = res.qbr[res.qbr$qbr.lower<=theta & theta<=res.qbr$qbr.upper,]
 cp_bi_qbr  = nrow(cp_bi_qbr)/nrow(res.qbr)
 
 
 #-----------------  summary of results -------------------------------
 print.result<-function(){
   cat("\n-------------------------------------------\n")
   cat(paste("         n =", n, "and theta =", theta))
   cat("\n-------------------------------------------\n")
   print(res.sum)
   cp     = c(cp_bi_mv,cp_bi_med, cp_bi_qbr)
   cp_res = matrix(cp, ncol=3, byrow=T)
   colnames(cp_res) = c("ML ","MBR", "QBR")
   rownames(cp_res) = c("cp %")
   cat("\n-------------------------------------------\n")
   cat("        coverage probability                         ")
   cat(paste("\n       ",(1-alpha)*100,"% CI - Two-sided\n"))
   cat("-------------------------------------------")
   cat("\n")
   print(round(100*cp_res,2))
 }
 
 print.result()

 