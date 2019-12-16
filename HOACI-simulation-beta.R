# Pinheiro, Ferrari and Medeiros (2019) - Higher-order approximate confidence intervals
#
# Monte Carlo simulation experiment to evaluate the finite-sample performance of confidence 
# intervals based on the alpha-quantile modified score (modified CIs - QBR), Wald-type CIs, 
# i.e. CIs that use the asymptotic normality of MLE (usual CIs - MLE), and those that are 
# constructed similarly to Wald-type CIs with median bias reduced estimates of all the 
# parameters in place of MLEs (adjusted CIs - MBR), using the beta regression model with
# parameters beta and gamma (Example 7).
#
# Define:
# prob  - nominal confidence level
# psi   - index of the interest parameter theta=(beta,gamma) 
# Rep   - number of replicates
# nobs  - sample size
# beta  - parameter vector of the mean (mu) regression model
# gamma - parameter vector of the precision (phi) regression model
prob <- 0.995 # CL=99.5% for one-sided interval and CL=99% for two-sided interval
psi <- 4      # psi=1 for beta_0, psi=2 for beta_1, psi=3 for gamma_0, psi=4 for gamma_1
Rep <- 100000
nobs <- 25
beta <- c(1,1)
gamma <- c(1,2)

# See the working directory
getwd() 
# Set the working directory on the terminal where the HOACI repository is located
setwd("put here the path where the HOACI repository is located") 

# Install packages :
# install.packages("nleqslv", repos="https://cran.rstudio.com")
# library(nleqslv)
# install.packages("betareg", repos="https://cran.rstudio.com")
# library(betareg)
# install.packages("ssym", repos="https://cran.rstudio.com")
# library(ssym)
# install.packages("doParallel", repos="https://cran.rstudio.com")
# library(doParallel)
# install.packages("tensorA", repos="https://cran.rstudio.com")
# library(tensorA)

require(doParallel)
require(doRNG)
require(setRNG)

## ----Defining functions-------------------------------------------------------
source("HOACI-functions.R",local=TRUE)
simula_hoaci <- function(s_vY, s_mX, s_mZ, parameter, psi, prob, r){
  require(nleqslv)
  require(betareg)
  require(tensorA)
  
  filenamefail <-paste("fail_psi",psi,"_prob",prob*1000,".txt",sep="") 
  
  FAIL=FALSE
  
  y<-s_vY
  X<-s_mX
  Z<-s_mZ
  
  nobs = nrow(X)
  iota =  matrix(rep(1,nobs),nobs,1);
  q<-ncol(X)
  m<-ncol(Z)
  numpar <- q+m
  beta<-parameter[1:q]
  gamma<-parameter[(q+1):numpar]
  u = qnorm(prob)
  
  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  # MLE
  PredictorVariablesMean      <- paste("X[,", 2:q,"]", sep="")
  if(m==1)
    Formula <- formula(paste("y ~ ", paste(PredictorVariablesMean, collapse=" + ")))
  else{
    PredictorVariablesPrecision <- paste("Z[,", 2:m,"]", sep="")
    Formula <- formula(paste("y ~ ", paste(PredictorVariablesMean, collapse=" + "),"|", paste(PredictorVariablesPrecision, collapse=" + ")))
  }  
  
  fit<-betareg::betareg(Formula, link.phi="log")
  if(fit$optim$convergence != 0){
    cat(paste("\n# r = ", r, " # no convergence in betareg", as.character(Sys.time()),"\n"), file = filenamefail, append = TRUE, sep = " ", fill = 100, labels = NULL)
    FAIL=TRUE
    return(FAIL)
  }
  thetahat<-as.vector(fit$optim$par)
  stderrors <- sqrt(diag(fit$vcov))
  
  if(is.nan(stderrors) || is.infinite(stderrors)){
    cat(paste("\n# r = ", r, " # stderrors error", as.character(Sys.time()),"\n"), file = filenamefail, append = TRUE, sep = " ", fill = 100, labels = NULL)
    FAIL=TRUE
    return(FAIL)
  }
  
  MLEU<-thetahat[psi]-(-u)*stderrors[psi]
  MLEL<-thetahat[psi]-( u)*stderrors[psi]
  
  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  # MBR
  MBR<-rep(0,numpar)
  for(i in 1:numpar){
    
    up<- qnorm(0.5)
    
    if(i<=q)# psi<=q => parameter of interest is beta_i
    { 
      # put the covariate corresponding to the parameter of interest (in the beta parameter vector) in the first column of X
      s_mXpsi <- cbind(s_mX[,i],s_mX[,-i])
      
      p<-c(thetahat[i],thetahat[-i])  # initial values: MLE for psi and for lambda
      ans<-try(nleqslv(p, scoreQBRbeta, jac=NULL, method="Newton",global="cline",control=list(trace=0,maxit=500,allowSingular=TRUE),y,s_mXpsi,Z,beta_score,beta_cumulants,link.mu="logit", link.phi="log",nu=NULL,up),TRUE)
      if(m==1){
        if(inherits(ans, "try-error") || ans$termcd != 1 || ans$x[numpar]<0){#error handling code, maybe just skip this iteration using
          cat(paste("\n# r = ", r, " # MBRbeta error", as.character(Sys.time()),"\n"), file = filenamefail, append = TRUE, sep = " ", fill = 100, labels = NULL)
          FAIL=TRUE
          return(FAIL)
        }
      }
      else{ # m>1
        if(inherits(ans, "try-error") || ans$termcd != 1){#error handling code, maybe just skip this iteration using
          cat(paste("\n# r = ", r, " # MBRbeta error", as.character(Sys.time()),"\n"), file = filenamefail, append = TRUE, sep = " ", fill = 100, labels = NULL)
          FAIL=TRUE
          return(FAIL)
        }
      }
      MBR[i] <- ans$x[1]
    }
    else# i>q => parameter of interest is gamma_i
    { 
      # put the covariate corresponding to the parameter of interest (in the gamma parameter vector) in the first column of Z
      if(m==1) Zpsi <- iota else Zpsi <- cbind(Z[,(i-q)],Z[,-(i-q)])
      
      p<-c(thetahat[1:q],thetahat[i],thetahat[-c(1:q,i)]) # initial values: MLE for psi and for lambda
      ans<-try(nleqslv(p, scoreQBRgamma, jac=NULL, method="Newton",global="cline",control=list(trace=0,maxit=500,allowSingular=TRUE),y,X,Zpsi,beta_score,beta_cumulants,link.mu="logit", link.phi="log",nu=NULL,up),TRUE)
      if(m==1){
        if(inherits(ans, "try-error") || ans$termcd != 1 || ans$x[numpar]<0){#error handling code, maybe just skip this iteration using
          cat(paste("\n# r = ", r, " # MBRgamma error", as.character(Sys.time()),"\n"), file = filenamefail, append = TRUE, sep = " ", fill = 100, labels = NULL)
          FAIL=TRUE
          return(FAIL)
        }
      }
      else{ # m>1
        if(inherits(ans, "try-error") || ans$termcd != 1){#error handling code, maybe just skip this iteration using
          cat(paste("\n# r = ", r, " # MBRgamma error", as.character(Sys.time()),"\n"), file = filenamefail, append = TRUE, sep = " ", fill = 100, labels = NULL)
          FAIL=TRUE
          return(FAIL)
        }
      }
      MBR[i] <- ans$x[(q+1)]    
    }
  }
  
  etatil <- s_mX%*%MBR[1:q];
  mutil <- exp(etatil) / (1+exp(etatil));
  
  deltatil <- Z%*%MBR[(q+1):numpar]
  phitil<-exp(deltatil)
  
  T <- as.vector( exp(etatil) / (1.0 + exp(etatil))^2 ) # d mu/d eta = 1/g'(mu)=mu(1-mu)
  H <- as.vector( phitil )  # d phi/d delta = phi 
  
  K22 = as.vector( phitil^2*(         psigamma(mutil*phitil,1) +             psigamma((1-mutil)*phitil,1) ) ) 
  K21 = as.vector( phitil  *( mutil  *psigamma(mutil*phitil,1) - (1-mutil)*  psigamma((1-mutil)*phitil,1) ) ) 
  K20 = as.vector(            mutil^2*psigamma(mutil*phitil,1) + (1-mutil)^2*psigamma((1-mutil)*phitil,1) - psigamma(phitil,1) )
  
  fisher = matrix(0,numpar,numpar)                               # kappa_(theta_i,theta_j)
  fisher[1:q,1:q]                  = crossprod(X,T*T*K22*X)      # kappa_betarbetas
  fisher[1:q,(q+1):numpar]         = crossprod(X,T*H*K21*Z)      # kappa_betargammas
  fisher[(q+1):numpar,1:q]         = t(fisher[1:q,(q+1):numpar]) # kappa_gammarbetas
  fisher[(q+1):numpar,(q+1):numpar]= crossprod(Z,H*H*K20*Z)      # kappa_gammargammas
  
  stderrors <- sqrt(diag(solve(fisher)));
  
  MBRU<-MBR[psi]-(-u)*stderrors[psi]
  MBRL<-MBR[psi]-( u)*stderrors[psi]
  
  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  # QBR
  if(psi<=q){# psi<=q => parameter of interest is beta_i
    
    # put the covariate corresponding to the parameter of interest (in the beta parameter vector) in the first column of X
    s_mXpsi <- cbind(s_mX[,psi],s_mX[,-psi])
    
    # CI upper limit QBRE
    up<- -u
    p<-c(thetahat[psi],thetahat[-psi])  # initial values: MLE for psi and for lambda
    ans<-try(nleqslv(p, scoreQBRbeta, jac=NULL, method="Newton",global="cline",control=list(trace=0,maxit=500,allowSingular=TRUE),y,s_mXpsi,Z,beta_score,beta_cumulants,link.mu="logit", link.phi="log",nu=NULL,up),TRUE)
    ans
    if(m==1){
      if(inherits(ans, "try-error") || ans$termcd != 1 || ans$x[numpar]<0){#error handling code, maybe just skip this iteration using
        cat(paste("\n# r = ", r, " # QBRUbeta error", as.character(Sys.time()),"\n"), file = filenamefail, append = TRUE, sep = " ", fill = 100, labels = NULL)
        FAIL=TRUE
        return(FAIL)
      }
    }else{ # m>1
      if(inherits(ans, "try-error") || ans$termcd != 1){#error handling code, maybe just skip this iteration using
        cat(paste("\n# r = ", r, " # QBRUbeta error", as.character(Sys.time()),"\n"), file = filenamefail, append = TRUE, sep = " ", fill = 100, labels = NULL)
        FAIL=TRUE
        return(FAIL)
      }
    }
    QBRU <- ans$x[1]
    
    # CI lower limit QBRE
    up<- u
    p<-c(thetahat[psi],thetahat[-psi])  # initial values: MLE for psi and for lambda
    ans<-try(nleqslv(p, scoreQBRbeta, jac=NULL, method="Newton",global="cline",control=list(trace=0,maxit=500,allowSingular=TRUE),y,s_mXpsi,Z,beta_score,beta_cumulants,link.mu="logit", link.phi="log",nu=NULL,up),TRUE)
    if(m==1){
      if(inherits(ans, "try-error") || ans$termcd != 1 || ans$x[numpar]<0){#error handling code, maybe just skip this iteration using
        cat(paste("\n# r = ", r, " # QBRLbeta error", as.character(Sys.time()),"\n"), file = filenamefail, append = TRUE, sep = " ", fill = 100, labels = NULL)
        FAIL=TRUE
        return(FAIL)
      }
    }else{ # m>1
      if(inherits(ans, "try-error") || ans$termcd != 1){#error handling code, maybe just skip this iteration using
        cat(paste("\n# r = ", r, " # QBRLbeta error", as.character(Sys.time()),"\n"), file = filenamefail, append = TRUE, sep = " ", fill = 100, labels = NULL)
        FAIL=TRUE
        return(FAIL)
      }
    }
    QBRL <- ans$x[1]
  }else{ # psi>q => parameter of interest is gamma_i
    
    # put the covariate corresponding to the parameter of interest (in the gamma parameter vector) in the first column of Z
    if(m==1)Zpsi <- iota else Zpsi <- cbind(Z[,(psi-q)],Z[,-(psi-q)])
    
    # CI upper limit QBRE
    up<- -u
    p<-c(thetahat[1:q],thetahat[psi],thetahat[-c(1:q,psi)]) # initial values: MLE for psi and for lambda
    ans<-try(nleqslv(p, scoreQBRgamma, jac=NULL, method="Newton",global="cline",control=list(trace=0,maxit=500,allowSingular=TRUE),y,X,Zpsi,beta_score,beta_cumulants,link.mu="logit", link.phi="log",nu=NULL,up),TRUE)
    if(m==1){
      if(inherits(ans, "try-error") || ans$termcd != 1 || ans$x[numpar]<0){#error handling code, maybe just skip this iteration using
        cat(paste("\n# r = ", r, " # QBRUgamma error", as.character(Sys.time()),"\n"), file = filenamefail, append = TRUE, sep = " ", fill = 100, labels = NULL)
        FAIL=TRUE
        return(FAIL)
      }
    }
    else{ # m>1
      if(inherits(ans, "try-error") || ans$termcd != 1){#error handling code, maybe just skip this iteration using
        cat(paste("\n# r = ", r, " # QBRUgamma error", as.character(Sys.time()),"\n"), file = filenamefail, append = TRUE, sep = " ", fill = 100, labels = NULL)
        FAIL=TRUE
        return(FAIL)
      }
    }
    QBRU <- ans$x[(q+1)]    
    
    # CI lower limit QBRE    
    up<- u
    p<-c(thetahat[1:q],thetahat[psi],thetahat[-c(1:q,psi)]) # initial values: MLE for psi and for lambda
    ans<-try(nleqslv(p, scoreQBRgamma, jac=NULL, method="Newton",global="cline",control=list(trace=0,maxit=500,allowSingular=TRUE),y,X,Zpsi,beta_score,beta_cumulants,link.mu="logit", link.phi="log",nu=NULL,up),TRUE)
    if(m==1){
      if(inherits(ans, "try-error") || ans$termcd != 1 || ans$x[numpar]<0){#error handling code, maybe just skip this iteration using
        cat(paste("\n# r = ", r, " # QBRLgamma error", as.character(Sys.time()),"\n"), file = filenamefail, append = TRUE, sep = " ", fill = 100, labels = NULL)
        FAIL=TRUE
        return(FAIL)
      }
    }
    else{ # m>1
      if(inherits(ans, "try-error") || ans$termcd != 1){#error handling code, maybe just skip this iteration using
        cat(paste("\n# r = ", r, " # QBRLgamma error", as.character(Sys.time()),"\n"), file = filenamefail, append = TRUE, sep = " ", fill = 100, labels = NULL)
        FAIL=TRUE
        return(FAIL)
      }
    }
    QBRL <- ans$x[(q+1)]    
  }
  
  coverageEMV_L = as.integer(MLEL <= parameter[psi]);
  coverageMBR_L = as.integer(MBRL <= parameter[psi]);
  coverageQBR_L = as.integer(QBRL <= parameter[psi]);
  coverageEMV_U = as.integer(parameter[psi] <= MLEU);
  coverageMBR_U = as.integer(parameter[psi] <= MBRU);
  coverageQBR_U = as.integer(parameter[psi] <= QBRU);
  coverageEMV_B = as.integer(MLEL <= parameter[psi])&&(parameter[psi] <= MLEU);
  coverageMBR_B = as.integer(MBRL <= parameter[psi])&&(parameter[psi] <= MBRU);
  coverageQBR_B = as.integer(QBRL <= parameter[psi])&&(parameter[psi] <= QBRU);
  
  lengthEMV =	MLEU - MLEL;
  lengthMBR =	MBRU - MBRL;
  lengthQBR =	QBRU - QBRL;
  
  return(c(FAIL,1-prob, coverageEMV_L, coverageMBR_L, coverageQBR_L, coverageEMV_U, coverageMBR_U, coverageQBR_U, coverageEMV_B, coverageMBR_B, coverageQBR_B, lengthEMV, lengthMBR, lengthQBR, MLEL, MBRL, QBRL, MLEU, MBRU, QBRU))
  
}

# ----Defining the file names---------------------------------------------------
filenameestimates<-paste("estimates_psi",psi,"_prob",prob*1000,"_Rep",Rep,".txt",sep="")
filenamecoverage<-paste("coverage_psi",psi,"_prob",prob*1000,"_Rep",Rep,".txt",sep="")
filenamemensagem <-paste("message_psi",psi,"_prob",prob*1000,"_Rep",Rep,".txt",sep="")
filenamefail <-paste("fail_psi",psi,"_prob",prob*1000,"_Rep",Rep,".txt",sep="")

## ----Genereting samples--------------------------------------------------------
theta<-c(beta,gamma)
numpar=length(theta)

setRNG(kind="Mersenne-Twister",normal.kind="Inversion",seed=777)
s_mX<-cbind( rep(1, nobs),matrix(runif(nobs*(numpar-2))-0.5,nobs,(length(beta)-1)) ) # for any dimension of beta
s_mZ<-cbind( rep(1, nobs),matrix(runif(nobs*(numpar-2))+1,nobs,(length(gamma)-1)) ) # for any dimension of gamma
eta = s_mX %*% beta
mu = exp(eta) / (1.0 + exp(eta)); #/* column vector nobs x 1 */
delta = s_mZ %*% gamma
phi = exp(delta) ; #/* column vector nobs x 1 */
samples = matrix(0,Rep,nobs);		#/* Inicializing matrix of samples s_vY */
for(i in 1:Rep){
  for(j in 1:nobs){
    #/* Beta (a,b) */
    samples[i,j] = rbeta(1, mu[j]*phi[j], (1-mu[j])*phi[j]);
    if(samples[i,j]==1.0000)samples[i,j]=0.9999999999;
    if(samples[i,j]==0.0000)samples[i,j]=0.0000000001;
  }
}

start.time <- Sys.time()

## ----Alocating Processors------------------------------------------------------
closeAllConnections()
ncores <- detectCores()
cat("\n ncores = ",ncores,"\n\n")
cl <- makeCluster(ncores)
registerDoParallel(cl, cores=ncores)
on.exit(stopCluster(cl))

## ----Obtaining the coverage----------------------------------------------------
tab <- foreach(r = 1:Rep, .combine = rbind, .packages=c('doParallel'), .export = ls(globalenv())) %dopar% {
  
  if(r == 1)cat(paste("\n# r = ",r, " # Inicio: ", as.character(Sys.time()),"\n"), file = filenamemensagem, append = TRUE, sep = " ", fill = 100, labels = NULL)
  
  s_vY<-samples[r,]
  result<-simula_hoaci(s_vY = s_vY, s_mX = s_mX, s_mZ = s_mZ, parameter = theta, psi = psi, prob = prob, r = r)

  if(r%%500==0)cat(paste("\n# r = ",r,             as.character(Sys.time()),"\n"), file = filenamemensagem, append = TRUE, sep = " ", fill = 100, labels = NULL)
  if(r  == Rep)cat(paste("\n# r = ",r, " # Fim: ", as.character(Sys.time()),"\n"), file = filenamemensagem, append = TRUE, sep = " ", fill = 100, labels = NULL)
  
  if(result[1]){
    return(5)
  }else{
    return(result[-1])
  }
}
fails<-tab[tab[,1]==5,]
validos<-tab[tab[,1]!=5,]
estimates<-validos[,14:19]
coverage<-c((1-prob),mean(validos[,2])*100,mean(validos[,3])*100,mean(validos[,4])*100,mean(validos[,5])*100,mean(validos[,6])*100,mean(validos[,7])*100,mean(validos[,8])*100,mean(validos[,9])*100,mean(validos[,10])*100,mean(validos[,11]),mean(validos[,12]),mean(validos[,13]))

# ----Saving the results-----------------------------------------------------------
cat(paste("\n# fails : ",nrow(rbind(fails,c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)))-1), file = filenamefail, append = TRUE, sep = " ", fill = 100, labels = NULL)
cat(estimates, file = filenameestimates, append = TRUE, sep = " ", fill = 100, labels = NULL)
cat(c("       |----------------- coverage  --------------| |--------- length ---------|"), file = filenamecoverage, append = TRUE, sep = " ", fill = 100, labels = NULL)
cat(c("1-prob MLEL MBRL QBRL MLEU MBRU QBRU MLEB MBRB QBRB    MLE       MBR       QBR"), file = filenamecoverage, append = TRUE, sep = " ", fill = 100, labels = NULL)
cat(coverage, file = filenamecoverage, append = TRUE, sep = "  ", fill = 100, labels = NULL)

## ----Conclusion-----------------------------------------------------------------
end.time <- Sys.time()
time.taken <- end.time - start.time
cat("\n# Conclusion: ", as.character(Sys.time()),"\n# Duration: ", format(time.taken), "\n", file = filenamemensagem, append = TRUE)
stopCluster(cl)
