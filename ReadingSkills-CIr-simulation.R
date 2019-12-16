# Pinheiro, Ferrari and Medeiros (2019) - Higher-order approximate confidence intervals
#
# Monte Carlo simulation experiment to evaluate the finite-sample performance of confidence 
# intervals obtained from the signed log-likelihood ratio statistic (r) and its third-order 
# modification, r* (Barndorff-Nielsen, 1986), using the beta regression model with n = 44 and 
# the values of the parameters taken as the MLEs computed with the reading skills data set.
# Define: 
# psi - index of the interest parameter theta=(beta,gamma)
# Rep - number of replicates
psi = 1 # psi=1 for beta_0, psi=2 for beta_1, psi=3 for beta_2, psi=4 for beta_3, psi=5 for gamma_0, psi=6 for gamma_1, psi=7 for gamma_2
Rep = 10000

# See the working directory
getwd() 
# Set the working directory on the terminal where the HOACI repository is located
setwd("put here the path where the HOACI repository is located") 

# Install packages :
# install.packages("likelihoodAsy", repos="https://cran.rstudio.com")
# library(likelihoodAsy)
# install.packages("knitr", repos="https://cran.rstudio.com")
# library(knitr)

require(knitr)
require(doParallel)
require(setRNG)

## ----Defining functions-------------------------------------------------------
opts_chunk$set(comment="", message=FALSE, warning=FALSE,
               tidy.opts=list(keep.blank.line=TRUE, width.cutoff=180),
               fig.width=8, fig.height=6,out.width='1\\textwidth',
               options(width=180),fig.show='hold',fig.align='center',cache=TRUE)
## ----Log likelihood for a beta regression model--------------------------------
loglik.beta<-function(parameter,data){
  numpar<-length(parameter)
  y<-data$y
  X<-data$X
  Z<-data$Z
  q<-ncol(X)
  m<-ncol(Z)
  beta <-parameter[1:q]
  gamma<-parameter[(q+1):numpar]
  
  nobs<-length(y)
  iota=rep(1,nobs)
  
  eta <- X%*%beta 
  mu <- exp(eta) / (1.0 + exp(eta)) #link.mu=logit
  
  delta <- Z%*%gamma
  phi<-exp(delta) # link.phi=log
  
  ystar <- log( y / (iota - y) );
  ydag <- log(iota-y);
  
  #  c <- log( gamma(phi) / ( gamma(mu*phi) * gamma((iota - mu)*phi) )  )
  c <- -lbeta((iota - mu)*phi,mu*phi)
  l = crossprod(iota, ((mu*phi-iota)*ystar + (phi-2*iota)*ydag + c ) ) 
  return(l)
  
}
## ----Data generation for a beta regression model-------------------------------
gendat.beta <- function(parameter, data){  
  numpar<-length(parameter)
  y<-data$y
  X<-data$X
  Z<-data$Z
  q<-ncol(X)
  m<-ncol(Z)
  beta <-parameter[1:q]
  gamma<-parameter[(q+1):numpar]
  
  nobs<-length(y)
  iota=rep(1,nobs)
  
  eta <- X%*%beta 
  mu <- exp(eta) / (1.0 + exp(eta)) #link.mu=logit
  
  delta <- Z%*%gamma
  phi<-exp(delta) # link.phi=log
  
  data$y <- rbeta(nobs, mu*phi, (1-mu)*phi)
  data$y[data$y==1.0000]=0.9999999999
  data$y[data$y==0.0000]=0.0000000001
  return(data)  
}
## ----Score function of the beta regression with link.mu=logit and link.phi=log-
grad.beta<-function(parameter,data){
  numpar<-length(parameter)
  y<-data$y
  X<-data$X
  Z<-data$Z
  q<-ncol(X)
  m<-ncol(Z)
  beta <-parameter[1:q]
  gamma<-parameter[(q+1):numpar]
  
  nobs<-length(y)
  iota=rep(1,nobs)
  
  eta <- X%*%beta
  mu <- exp(eta) / (1.0 + exp(eta))  #link.mu=logit
  T <- exp(eta) / (1.0 + exp(eta))^2 # d mu/d eta = 1/g'(mu)=mu(1-mu) 
  
  delta <- Z%*%gamma
  phi<-exp(delta)  # link.phi=log
  H <- phi         # d phi/d delta = phi
  
  ystar = log( y / (iota - y) );
  mustar = psigamma(mu*phi, 0) - psigamma((iota - mu)*phi, 0);
  ydag = log(iota-y);
  mudag =  -psigamma(phi,0)+psigamma((iota-mu)*phi,0);
  
  U <- cbind( phi*(ystar-mustar) , (mu*(ystar-mustar)+ydag-mudag) )
  
  c(crossprod(X,T*U[,1]),crossprod(Z,H*U[,2] ))
}
## ----
psifcn.beta <- function(theta){
  return(theta[psi])
}
## ----

## ----Defining file names-------------------------------------------------------
filenamecoverage<-paste("ReadingSkills_CIr_coverage_psi",psi,"_rep",Rep,".txt",sep="")
filenamemessage <- paste("ReadingSkills_CIr_message_psi",psi,"_rep",Rep,".txt",sep="")
filenamefail <-paste("ReadingSkills_CIr_fail_psi",psi,"_rep",Rep,".txt",sep="")

beta<-c(1.1232,      -0.7416,       0.4864,      -0.5813)
gamma<-c(3.304,        1.747,        1.229)
theta<-c(beta,gamma)
numpar=length(theta)

# values of the parameters taken as the MLEs computed with the reading skills data set
dyslexia<-c(-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1)
iq<-c(0.827,  0.590,  0.471,  1.144, -0.676, -0.795, -0.281, -0.914, -0.043,  0.907,  0.511,  1.223,  0.590,  1.856, -0.399,  0.590, -0.043,  1.738,  0.471,  1.619,  1.144, -0.201,
      -0.281,  0.590,  1.777, -0.083, -0.162, -0.795, -0.281, -0.874,  0.313,  0.709,  1.223, -1.230, -0.162, -0.993, -1.191, -1.745, -1.745, -0.439, -1.666, -1.507, -0.518, -1.270)
nobs=length(dyslexia)
s_mX<-cbind( rep(1, nobs),dyslexia,iq, dyslexia*iq) # para qualquer dimens?o de beta
s_mZ<-cbind( rep(1, nobs),dyslexia,iq) # para qualquer dimens?o de gamma

eta = s_mX %*% beta
mu = exp(eta) / (1.0 + exp(eta)); #/* column vector nobs x 1 */
delta = s_mZ %*% gamma
phi = exp(delta) ; #/* column vector nobs x 1 */

# parameters beta and gamma (Example 7)
# nobs<-25
# beta<-c(1,1)
# gamma<-c(1,2)
# # gamma<-c(50)
# theta<-c(beta,gamma)
# numpar=length(theta)
# 
# s_mX<-cbind( rep(1, nobs),matrix(runif(nobs*(numpar-2))-0.5,nobs,(length(beta)-1)) ) 
# s_mZ<-cbind( rep(1, nobs),matrix(runif(nobs*(numpar-2))+1,nobs,(length(gamma)-1)) )
# eta = s_mX %*% beta
# mu = exp(eta) / (1.0 + exp(eta)); #/* column vector nobs x 1 */
# delta = s_mZ %*% gamma
# phi = exp(delta) ; #/* column vector nobs x 1 */

setRNG(kind="Mersenne-Twister",normal.kind="Inversion",seed=777)
samples = matrix(0,Rep,nobs);		#/* Inicializing matrix of samples s_vY */
for(i in 1:Rep){
  for(j in 1:nobs){
    #/* Beta (a,b) */
    samples[i,j] = rbeta(1, mu[j]*phi[j], (1-mu[j])*phi[j]);
    if(samples[i,j]==1.0000)samples[i,j]=0.9999999999;
    if(samples[i,j]==0.0000)samples[i,j]=0.0000000001;
  }
}

simula_CIr<- function(s_vY, s_mX, s_mZ, parameter, psi, r){
  require(likelihoodAsy)
  
  FAIL=FALSE
  
  data.gen<-list(X = s_mX, Z = s_mZ, y = s_vY)
  
  rs.int<-try(likelihoodAsy::rstar.ci(data=data.gen, thetainit = parameter, floglik = loglik.beta,
                                   fpsi = psifcn.beta, fscore = grad.beta, datagen = gendat.beta,
                                   trace=FALSE, seed=1223, psidesc=NULL) ,TRUE)
  if(inherits(rs.int, "try-error")){#error handling code, maybe just skip this iteration using
    cat(paste("\n# r = ", r, " # error", as.character(Sys.time()),"\n"), file = filenamefail, append = TRUE, sep = " ", fill = 100, labels = NULL)
    FAIL=TRUE
    return(FAIL)
  }
  
  coverage_CIr<-array(0,3)
  coverage_CIrs<-array(0,3)
  for(i in 1:3){
    coverage_CIr[i] <-as.integer((rs.int$CIr[i,1]  <= parameter[psi])&&(parameter[psi] <= rs.int$CIr[i,2]))
    coverage_CIrs[i]<-as.integer((rs.int$CIrs[i,1] <= parameter[psi])&&(parameter[psi] <= rs.int$CIrs[i,2]))
  }
  return(c(FAIL,coverage_CIr[1],coverage_CIr[2],coverage_CIr[3],coverage_CIrs[1],coverage_CIrs[2],coverage_CIrs[3]))
}

start.time <- Sys.time()
## ----Alocating Processors------------------------------------------------------
closeAllConnections()
ncores <- detectCores()
cat("\n ncores = ",ncores,"\n\n")
cl <- makeCluster(ncores)
registerDoParallel(cl, cores=ncores)
# on.exit(stopCluster(cl))
## ----Obtaining the coverage----------------------------------------------------
# tab <- foreach(r = 1:(Rep+fail), .combine = rbind, .packages=c('doParallel'), .export = ls(globalenv())) %dopar% {
# tab <- foreach(r = 1:Rep, .combine = rbind, .packages=c('doParallel')) %dopar% {
tab <- foreach(r = 1:Rep, .combine = rbind, .packages=c('doParallel'), .export = ls(globalenv())) %dopar% {
  if(r == 1)cat(paste("\n# r = ",r, " # Inicio: ", as.character(Sys.time()),"\n"), file = filenamemessage, append = TRUE, sep = " ", fill = 100, labels = NULL)
  
  s_vY<-samples[r,]
  result<-simula_CIr(s_vY = s_vY, s_mX = s_mX, s_mZ = s_mZ, parameter = theta, psi = psi, r = r)
  
  if(r%%500==0)cat(paste("\n# r = ",r,             as.character(Sys.time()),"\n"), file = filenamemessage, append = TRUE, sep = " ", fill = 100, labels = NULL)
  if(r  == Rep)cat(paste("\n# r = ",r, " # Fim: ", as.character(Sys.time()),"\n"), file = filenamemessage, append = TRUE, sep = " ", fill = 100, labels = NULL)
  
  if(result[1]){
    return(5)
    cat(paste("\n# fail in r = ",r,                as.character(Sys.time()),"\n"), file = filenamefail, append = TRUE, sep = " ", fill = 100, labels = NULL)
  }
  return(result[-1])
}
fails<-tab[tab[,1]==5,]
validos<-tab[tab[,1]!=5,]
coverage<-c(mean(validos[,1])*100,mean(validos[,2])*100,mean(validos[,3])*100,mean(validos[,4])*100,mean(validos[,5])*100,mean(validos[,6])*100)

# ----Saving the results--------------------------------------------------------
cat(paste("\n# fails : ",nrow(rbind(fails,c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)))-1), file = filenamefail, append = TRUE, sep = " ", fill = 100, labels = NULL)
cat("     CIr            CIRs   ",file=filenamecoverage, append = TRUE, sep = "  ", fill = 100, labels = NULL)
cat("90%  95%  99%  90%  95%  99%",file=filenamecoverage, append = TRUE, sep = "  ", fill = 100, labels = NULL)
cat(coverage,file=filenamecoverage, append = TRUE, sep = "  ", fill = 100, labels = NULL)

# ----Conclusion-----------------------------------------------------------------
duration = Sys.time() - start.time
cat("\n# Conclusion: ", as.character(Sys.time()),"\n# Duration: ", format(duration), "\n", file = filenamemessage, append = TRUE)
stopCluster(cl)
