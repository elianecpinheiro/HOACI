# Pinheiro, Ferrari and Medeiros (2019) - Higher-order approximate confidence intervals
#
# Functions needed to HOACI-examples.R, ReadingSkills-simulation.R and HOACI-simulation-beta.R. 

options(warn=-1)# turn off warnings
# options(warn=0)# turn on warnings

require(nleqslv)
require(betareg)
require(ssym)
require(doParallel)
require(tensorA)

###########################################  
# functions for all distributions
###########################################  

regressioncumulants<-function(parameter,v,Nx,Nz,cumulants,link.mu,link.phi,nu){
  
  y<-v
  X<-Nx
  Z<-Nz
  q<-ncol(X)
  m<-ncol(Z)
  numpar<-q+m

  LINK<-link.matrix(parameter,X,Z,link.mu,link.phi)
  T<-as.vector(LINK[[1]])
  H<-as.vector(LINK[[2]])
  S<-as.vector(LINK[[3]])
  Q<-as.vector(LINK[[4]])
  
  K <- cumulants(parameter,v,Nx,Nz,link.phi,nu)
  
  K22 = as.vector(K[[1]]);K21 = as.vector(K[[2]]);K20 = as.vector(K[[3]]);         #K22=K_mu,mu        K21=K_mu,phi        K20=K_phi,phi
  K33 = as.vector(K[[4]]);K32 = as.vector(K[[5]]);K31 = as.vector(K[[6]]); 
  K30 = as.vector(K[[7]]);                                                         #K33=K_mu,mu,mu     K32=K_mu,mu,phi     K31=K_mu,phi,phi     K30=K_phi,phi,phi
  K44 = as.vector(K[[8]]);K43 = as.vector(K[[9]]);K42 = as.vector(K[[10]]);
  K41 = as.vector(K[[11]]);K40 = as.vector(K[[12]])                                #K44=K_mu,mu,mu,mu  K43=K_mu,mu,mu,phi  K42=K_mu,mu,phi,phi  K41=K_mu,phi,phi,phi K40=K_phi,phi,phi,phi
  Kmu2  = as.vector(K[[13]]);Kmu1  = as.vector(K[[14]]);Kmu0  = as.vector(K[[15]]) #Kmu2=K_mu,mumu     Kmu1=K_mu,muphi     Kmu0=K_mu,phiphi
  Kphi2 = as.vector(K[[16]]);Kphi1 = as.vector(K[[17]]);Kphi0 = as.vector(K[[18]]) #Kphi2=K_phi,mumu   Kphi1=K_phi,muphi   Kphi0=K_phi,phiphi
  
  #///////////////////////////////////////////////////////////////////////////
  #///////////////////////////////////////////////////////////////////////////
  K2 = matrix(NA_real_,numpar,numpar)   # kappa_(theta_i,theta_j)
  K2[1:q,1:q]                  = crossprod(X,T*T*K22*X)  # kappa_betarbetas
  K2[1:q,(q+1):numpar]         = crossprod(X,T*H*K21*Z)  # kappa_betargammas
  K2[(q+1):numpar,1:q]         = t(K2[1:q,(q+1):numpar]) # kappa_gammarbetas
  K2[(q+1):numpar,(q+1):numpar]= crossprod(Z,H*H*K20*Z)  # kappa_gammargammas
  #///////////////////////////////////////////////////////////////////////////
  #///////////////////////////////////////////////////////////////////////////
  K3 = array(NA_real_,dim=c(numpar,numpar,numpar)) # kappa_(theta_i,theta_j,theta_k)
  Kab = array(NA_real_,dim=c(numpar,numpar,numpar)) # kappa_(psi,ab) psi,a,b={beta_1,...,beta_q,gamma_1,....gamma_m}
  
  for (t in 1:q ) {
    xt <- X[,t]
    if (q > 0L) {
      bb1 <- crossprod(X,(T*T*T*Kmu2+T*T*S*K22)*xt*X)
      bb2 <- crossprod(X,T*T*T*K33*xt*X)
    } else bb1 <- bb2 <- crossprod(X)
    
    if ((q > 0L) & (m > 0L)) {
      bg1 <- crossprod(X,(T*T*H*Kphi2+T*H*S*K21)*xt*Z)
      bg2 <- crossprod(X,T*T*H*K32*xt*Z)
    } else  bg1 <- bg2 <- crossprod(X ,Z)
    
    if (m > 0L) {
      gg1  <- crossprod(Z,(T*H*H*Kmu0+T*H*Q*K21)*xt*Z)
      bg12 <- crossprod(X,T*T*H*Kmu1*xt*Z)
      gg2  <- crossprod(Z,T*H*H*K31*xt*Z)
    }  else gg1 <- gg2 <- crossprod(Z)
    
    K3[t,,] <- rbind(cbind(bb2,bg2),cbind(t(bg2),gg2))
    Kab[t,,]<- rbind(cbind(bb1,bg1),cbind(t(bg12),gg1))
  }
  
  for (t in 1:m ) {
    zt <- Z[,t]
    if (q > 0L) {
      bb1 <- crossprod(X,T*T*H*Kmu1*zt*X)
      bb2 <- crossprod(X,T*T*H*K32*zt*X)
    } else bb1 <- bb2 <- crossprod(X)
    
    if ((q > 0L) & (m > 0L)) {
      bg1 <- crossprod(X,T*H*H*Kphi1*zt*Z)
      bg2 <- crossprod(X,T*H*H*K31*zt*Z)
    } else  bg1 <- bg2 <- crossprod(X ,Z)
    
    if (m > 0L) {
      gg1  <- crossprod(Z,(H*H*H*Kphi0+H*H*Q*K20)*zt*Z)
      bg12 <- crossprod(X,(T*H*H*Kmu0+T*H*Q*K21)*zt*Z)
      gg2  <- crossprod(Z,H*H*H*K30*zt*Z)
    }  else gg1 <- gg2 <- crossprod(Z)
    
    K3[(q+t),,] <- rbind(cbind(bb2,bg2),cbind(t(bg2),gg2))
    Kab[(q+t),,]<- rbind(cbind(bb1,bg1),cbind(t(bg12),gg1))
  }
  #///////////////////////////////////////////////////////////////////////////
  #///////////////////////////////////////////////////////////////////////////
  K4=array(NA_real_,dim=c(numpar,numpar,numpar,numpar)) # kappa_(theta_i,theta_j,theta_k,theta_l)
  
  for (t in 1:q ) {
    xt <- X[,t]
    for (s in 1:q ) {
      xs <- X[,s]
      
      if (q > 0L) {
        bb2 <- crossprod(X,T*T*T*T*K44*xt*xs*X)
      } else bb1 <- bb2 <- crossprod(X)
      
      if ((q > 0L) & (m > 0L)) {
        bg2 <- crossprod(X,T*T*T*H*K43*xt*xs*Z)
      } else  bg1 <- bg2 <- crossprod(X ,Z)
      
      if (m > 0L) {
        gg2 <- crossprod(Z,T*T*H*H*K42*xt*xs*Z)
      }  else gg1 <- gg2 <- crossprod(Z)
      
      K4[t,s,,] <- rbind(cbind(bb2,bg2),cbind(t(bg2),gg2))
    }
    
    for (s in 1:m ) {
      zs <- Z[,s]
      
      if (q > 0L) {
        bb2 <- crossprod(X,T*T*T*H*K43*xt*zs*X)
      } else bb1 <- bb2 <- crossprod(X)
      
      if ((q > 0L) & (m > 0L)) {
        bg2 <- crossprod(X,T*T*H*H*K42*xt*zs*Z)
      } else  bg1 <- bg2 <- crossprod(X ,Z)
      
      if (m > 0L) {
        gg2 <- crossprod(Z,T*H*H*H*K41*xt*zs*Z)
      }  else gg1 <- gg2 <- crossprod(Z)
      
      K4[t,(q+s),,] <- rbind(cbind(bb2,bg2),cbind(t(bg2),gg2))
    }
  }
  
  for (t in 1:m ) {
    zt <- Z[,t]
    for (s in 1:q ) {
      xs <- X[,s]
      if (q > 0L) {
        bb2 <- crossprod(X,T*T*T*H*K43*zt*xs*X)
      } else bb1 <- bb2 <- crossprod(X)
      
      if ((q > 0L) & (m > 0L)) {
        bg2 <- crossprod(X,T*T*H*H*K42*zt*xs*Z)
      } else  bg1 <- bg2 <- crossprod(X ,Z)
      
      if (m > 0L) {
        gg2 <- crossprod(Z,T*H*H*H*K41*zt*xs*Z)
      }  else gg1 <- gg2 <- crossprod(Z)
      
      K4[(q+t),s,,] <- rbind(cbind(bb2,bg2),cbind(t(bg2),gg2))
    }
    
    for (s in 1:m ) {
      zs <- Z[,s]
      
      if (q > 0L) {
        bb2 <- crossprod(X,T*T*H*H*K42*zt*zs*X)
      } else bb1 <- bb2 <- crossprod(X)
      
      if ((q > 0L) & (m > 0L)) {
        bg2 <- crossprod(X,T*H*H*H*K41*zt*zs*Z)
      } else  bg1 <- bg2 <- crossprod(X ,Z)
      
      if (m > 0L) {
        gg2 <- crossprod(Z,H*H*H*H*K40*zt*zs*Z)
      }  else gg1 <- gg2 <- crossprod(Z)
      
      K4[(q+t),(q+s),,] <- rbind(cbind(bb2,bg2),cbind(t(bg2),gg2))
    }
  }
  return(list(K2, K3, K4, Kab))
}

profilecumulantsbeta<-function(parameter,v,Nx,Nz,cumulants,link.mu,link.phi,nu){
  #parameter=(beta_0,beta_1,..,beta_(q-1),gamma_0,gamma_1,...,gamma_(m-1))
  #obtaining the cumulants of the profile likelihood of beta0 
  require(tensorA)
  
  q<-ncol(Nx)
  m<-ncol(Nz)
  numpar<-q+m
  
  c<-regressioncumulants(parameter,v,Nx,Nz,cumulants,link.mu,link.phi,nu)
  
  K2  = c[[1]] # kappa_(parameter_i,parameter_j)
  K3  = c[[2]] # kappa_(parameter_i,parameter_j,parameter_k)
  K4  = c[[3]] # kappa_(parameter_i,parameter_j,parameter_k,parameter_l)
  Kab = c[[4]] # kappa_(parameter_i,parameter_jparameter_k)
  #///////////////////////////////////////////////////////////////////////////
  #///////////////////////////////////////////////////////////////////////////
  Beta = rep(NA_real_,numpar-1) # beta_beta0^a  a={beta_1,...,beta_q,gamma_1,....gamma_m}
  K2inv = solve(K2[-1,-1],tol=1e-100)
  Beta=K2[1,-1]%*%K2inv
  #///////////////////////////////////////////////////////////////////////////
  #///////////////////////////////////////////////////////////////////////////
  Kabtensor <- tensorA::to.tensor(array(Kab[-1,-1,-1]),c(A=(numpar-1),B=(numpar-1),C=(numpar-1)))
  K3tensor  <- tensorA::to.tensor(array(K3[-1,-1,-1]),c(A=(numpar-1),B=(numpar-1),C=(numpar-1)))
  K4tensor3 <- to.tensor(array(K4[-1,-1,-1, 1]),c(A=(numpar-1),B=(numpar-1),C=(numpar-1)))
  K4tensor  <- to.tensor(array(K4[-1,-1,-1,-1]),c(A=(numpar-1),B=(numpar-1),C=(numpar-1),D=(numpar-1)))
  A <- to.tensor(array(Beta),c(A=(numpar-1),E=1))
  B <- to.tensor(array(Beta),c(B=(numpar-1),E=1))
  C <- to.tensor(array(Beta),c(C=(numpar-1),E=1))
  D <- to.tensor(array(Beta),c(D=(numpar-1),E=1))
  
  if(q==1 & m==1){
    k1<- 0;
    k1<- - 1/2 *sum( K2inv*( Kab[-1,-1,1] + K3[-1,-1,1] - Kab[-1,-1,-1]*Beta - K3 [-1,-1,-1]*Beta  )  )
  }else{
    k1<- 0;
    k1<- - 1/2 *sum( K2inv*( Kab[-1,-1,1] + K3[-1,-1,1] - apply((sweep(Kab[-1,-1,-1],3,Beta,FUN="*")),c(1,2),sum)
                             - apply((sweep(K3 [-1,-1,-1],3,Beta,FUN="*")),c(1,2),sum)  )  )
  }
  k2<-K2[1,1]- tcrossprod(Beta%*%K2[-1,-1],Beta)
  k3<-K3[1,1,1] - 3%*%Beta%*%K3[-1,1,1] + 3%*%tcrossprod(Beta%*%K3[-1,-1,1],Beta) - sum(K3tensor %e% A %e% B %e% C)
  k4<-K4[1,1,1,1] - 4%*%Beta%*%K4[-1,1,1,1] + 6%*%tcrossprod(Beta%*%K4[-1,-1,1,1],Beta) - 4%*%sum(K4tensor3 %e% A %e% B %e% C) + sum(K4tensor %e% A %e% B %e% C %e% D)
  return(k=c(k1,k2,k3,k4))
}

Mbeta<-function(parameter,v,Nx,Nz,cumulants,link.mu,link.phi,nu,u){
  k<-profilecumulantsbeta(parameter,v,Nx,Nz,cumulants,link.mu,link.phi,nu)
  up<-u
  return(- k[1] - sqrt(k[2])*up - k[3]/k[2]/6*(up^2-1) - k[4]/k[2]^(3/2)/24*(up^3-3*up)+ k[3]^2/k[2]^(5/2)/36*(2*up^3 - 5*up))
}

scoreQBRbeta<-function(parameter,v,Nx,Nz,score,cumulants,link.mu,link.phi,nu,u){
  numpar<-length(parameter)
  y<-array(v)
  nobs<-length(y)
  y<-matrix(y,nobs,1)
  X<-Nx
  Z<-Nz
  up<-u
  
  q<-ncol(X)
  m<-ncol(Z)
  
  LINK<-link.matrix(parameter,X,Z,link.mu,link.phi)
  T<-as.vector(LINK[[1]])
  H<-as.vector(LINK[[2]])
  
  U <- score(parameter, y, X, Z, link.phi, nu)
  
  if(q==1 & m==1){
    c(crossprod(X[,1]  ,T*U[,1]) + Mbeta(parameter, y, X, Z, cumulants, link.mu, link.phi, nu, up),
      crossprod(Z      ,H*U[,2])
    )
  }else{
    c(crossprod(X[,1]  ,T*U[,1]) + Mbeta(parameter, y, X, Z, cumulants, link.mu, link.phi, nu, up),
      crossprod(X[,2:q],T*U[,1]),
      crossprod(Z      ,H*U[,2])
    )
  }
}

profilecumulantsgamma<-function(parameter,v,Nx,Nz,cumulants,link.mu,link.phi,nu){
  #parameter=(beta_0,beta_1,..,beta_q,gamma_0,gamma_1,...,gamma_m)
  #obtaining the cumulants of the profile likelihood of gamma0
  require(tensorA)
  
  q<-ncol(Nx)
  m<-ncol(Nz)
  numpar<-q+m

  c<-regressioncumulants(parameter,v,Nx,Nz,cumulants,link.mu,link.phi,nu)
  
  K2  = c[[1]] # kappa_(parameter_i,parameter_j)
  K3  = c[[2]] # kappa_(parameter_i,parameter_j,parameter_k)
  K4  = c[[3]] # kappa_(parameter_i,parameter_j,parameter_k,parameter_l)
  Kab = c[[4]] # kappa_(parameter_i,parameter_jparameter_k)
  #///////////////////////////////////////////////////////////////////////////
  #///////////////////////////////////////////////////////////////////////////
  Beta = rep(0,numpar-1) # beta_beta0^a  a={beta_1,...,beta_q,gamma_1,....gamma_m}
  K2inv = solve(K2[-(q+1),-(q+1)],tol=1e-100)
  Beta=K2[(q+1),-(q+1)]%*%K2inv
  #///////////////////////////////////////////////////////////////////////////
  #///////////////////////////////////////////////////////////////////////////
  K3tensor <- to.tensor(array(K3[-(q+1),-(q+1),-(q+1)]),c(A=(numpar-1),B=(numpar-1),C=(numpar-1)))
  K4tensor3 <- to.tensor(array(K4[-(q+1),-(q+1),-(q+1), (q+1)]),c(A=(numpar-1),B=(numpar-1),C=(numpar-1)))
  K4tensor  <- to.tensor(array(K4[-(q+1),-(q+1),-(q+1),-(q+1)]),c(A=(numpar-1),B=(numpar-1),C=(numpar-1),D=(numpar-1)))
  A <- to.tensor(array(Beta),c(A=(numpar-1),E=1))
  B <- to.tensor(array(Beta),c(B=(numpar-1),E=1))
  C <- to.tensor(array(Beta),c(C=(numpar-1),E=1))
  D <- to.tensor(array(Beta),c(D=(numpar-1),E=1))

  k1<- 0;
  if(q==1 & m==1){
    k1<- - 1/2 *sum( K2inv*( Kab[-(q+1),-(q+1),(q+1)] + K3[-(q+1),-(q+1),(q+1)] - Kab[-(q+1),-(q+1),-(q+1)]*Beta - K3 [-(q+1),-(q+1),-(q+1)]*Beta  )  )
  }else{
    k1=- 1/2 *sum( K2inv*( Kab[-(q+1),-(q+1),(q+1)] + K3[-(q+1),-(q+1),(q+1)]
                           - apply((sweep(Kab[-(q+1),-(q+1),-(q+1)],3,Beta,FUN="*")),c(1,2),sum)
                           - apply((sweep(K3 [-(q+1),-(q+1),-(q+1)],3,Beta,FUN="*")),c(1,2),sum) )  )
  }
  k2<-K2[(q+1),(q+1)]- tcrossprod(Beta%*%K2[-(q+1),-(q+1)],Beta)
  k3<-K3[(q+1),(q+1),(q+1)] - 3%*%Beta%*%K3[-(q+1),(q+1),(q+1)] + 3%*%tcrossprod(Beta%*%K3[-(q+1),-(q+1),(q+1)],Beta) - sum(K3tensor %e% A %e% B %e% C)
  k4<-K4[(q+1),(q+1),(q+1),(q+1)] - 4%*%Beta%*%K4[-(q+1),(q+1),(q+1),(q+1)] + 6%*%tcrossprod(Beta%*%K4[-(q+1),-(q+1),(q+1),(q+1)],Beta) - 4%*%sum(K4tensor3 %e% A %e% B %e% C) + sum(K4tensor %e% A %e% B %e% C %e% D)
  
  return(c(k1,k2,k3,k4))
}

Mgamma<-function(parameter,v,Nx,Nz,cumulants,link.mu,link.phi,nu,u){ #recebe os parâmetros, calcula os cumulantes para cada parâmetro e calcula M
  k<-profilecumulantsgamma(parameter,v,Nx,Nz,cumulants,link.mu,link.phi,nu)
  up<-u
  return(- k[1] - sqrt(k[2])*up - k[3]/k[2]/6*(up^2-1) - k[4]/k[2]^(3/2)/24*(up^3-3*up)+ k[3]^2/k[2]^(5/2)/36*(2*up^3 - 5*up))
}

scoreQBRgamma<-function(parameter,v,Nx,Nz,score,cumulants,link.mu,link.phi,nu,u){

  y<-array(v)
  nobs<-length(y)
  y<-matrix(y,nobs,1)
  X<-Nx
  Z<-Nz
  up<-u
  
  q<-ncol(X)
  m<-ncol(Z)
  numpar<-q+m
  
  LINK<-link.matrix(parameter,X,Z,link.mu,link.phi)
  T<-as.vector(LINK[[1]])
  H<-as.vector(LINK[[2]])
  
  U <- score(parameter, y, X, Z, link.phi, nu)

  if(m==1){
    c(crossprod(X    ,T*U[,1]),
      crossprod(Z[,1],H*U[,2]) + Mgamma(parameter, y, X, Z,cumulants,link.mu,link.phi,nu,up)
    )
  }
  else{
    c(crossprod(X      ,T*U[,1]),
      (crossprod(Z[,1]  ,H*U[,2]) + Mgamma(parameter, y, X, Z,cumulants,link.mu,link.phi,nu,up)),
      crossprod(Z[,2:m],H*U[,2])
    )
  }
}

QBRE<-function(parameter,v,Nx,Nz,score,cumulants,link.mu,link.phi,nu,u){
  require(nleqslv)
  
  y<-v
  nobs<-length(y)
  y<-matrix(y,nobs,1)
  X<-Nx
  Z<-Nz
  
  q<-ncol(X)
  m<-ncol(Z)
  numpar<-q+m
  
  beta <-parameter[1:q]
  gamma<-parameter[(q+1):numpar]
  nobs<-length(y)
  iota =  matrix(rep(1,nobs),nobs,1);
  
  if(psi<=q){# beta
    Xpsi <- cbind(X[,psi],X[,-psi])
    
    if(q==1 & m==1){
      # MBRE
      p<-c(parameter[psi],parameter[-psi])  # initial guess : MLE
      ans<-try(nleqslv(p, scoreQBRbeta, jac=NULL, method="Newton",control=list(trace=0,maxit=500),y,Xpsi,Z,score,cumulants,link.mu,link.phi,nu,qnorm(0.5)),TRUE)
      if(inherits(ans, "try-error") || ans$termcd != 1){
        cat(paste0("\nQBRM beta_",(psi-1)," convergence termination code : ",ans$termcd),file="warning.txt",append=T)
        QBRM <- NaN
      }else{QBRM <- ans$x[1]}
      
      # CI upper limit QBRE
      up<- -u
      p<-c(parameter[psi],parameter[-psi])  # initial guess : MLE
      ans<-try(nleqslv(p, scoreQBRbeta, jac=NULL, method="Newton",control=list(trace=0,maxit=50,allowSingular=T),y,Xpsi,Z,score,cumulants,link.mu,link.phi,nu,up),TRUE)
      if(inherits(ans, "try-error") || ans$termcd != 1){
        p<-c(parameter[psi],parameter[-psi]) # initial guess : MLE
        ans<-try(nleqslv(p, scoreQBRbeta, jac=NULL, method="Newton",global="hook",control=list(trace=0,maxit=100,allowSingular=T),y,Xpsi,Z,score,cumulants,link.mu,link.phi,nu,up),TRUE)
      }
      if(inherits(ans, "try-error") || (ans$termcd != 1 && ans$termcd != 3)){
        cat(paste0("\nQBRU beta_",(psi-1)," convergence termination code : ",ans$termcd),file="warning.txt",append=T)
        QBRU <- NaN
      }else{QBRU <- ans$x[1]}
      if(ans$termcd == 3)cat(paste0("\nQBRU beta_",(psi-1),"  : No better point found (algorithm has stalled)"),file="warning.txt",append=T)
      
      # CI lower limit QBRE
      up<- u
      p<-c(parameter[psi],parameter[-psi])  # initial guess : MLE
      ans<-try(nleqslv(p, scoreQBRbeta, jac=NULL, method="Newton",control=list(trace=0,maxit=50,allowSingular=T),y,Xpsi,Z,score,cumulants,link.mu,link.phi,nu,up),TRUE)
      if(inherits(ans, "try-error") || ans$termcd != 1){
        p<-c(parameter[psi],parameter[-psi])  # initial guess : MLE
        ans<-try(nleqslv(p, scoreQBRbeta, jac=NULL, method="Newton",global="hook",control=list(trace=0,maxit=100,allowSingular=T),y,Xpsi,Z,score,cumulants,link.mu,link.phi,nu,up),TRUE)
      }
      if(inherits(ans, "try-error") || (ans$termcd != 1 && ans$termcd != 3)){
        cat(paste0("\nQBRL beta_",(psi-1)," convergence termination code : ",ans$termcd),file="warning.txt",append=T)
        QBRL <- NaN
      }else{QBRL <- ans$x[1]}
      if(ans$termcd == 3)cat(paste0("\nQBRL beta_",(psi-1),"  : No better point found (algorithm has stalled)"),file="warning.txt",append=T)
    }else{
      # MBRE
      p<-c(parameter[psi],parameter[-psi])  # initial guess : MLE
      ans<-try(nleqslv(p, scoreQBRbeta, jac=NULL, method="Newton",control=list(trace=0,maxit=500),y,Xpsi,Z,score,cumulants,link.mu,link.phi,nu,qnorm(0.5)),TRUE)
      if(inherits(ans, "try-error") || ans$termcd != 1){
        cat(paste0("\nQBRM beta_",(psi-1)," convergence termination code : ",ans$termcd),file="warning.txt",append=T)
        QBRM <- NaN
      }else{QBRM <- ans$x[1]}
      
      # CI upper limit QBRE
      up<- -u
      p<-c(parameter[psi],parameter[-psi])  # initial guess : MLE
      ans<-try(nleqslv(p, scoreQBRbeta, jac=NULL, method="Newton",control=list(trace=0,maxit=50,allowSingular=T),y,Xpsi,Z,score,cumulants,link.mu,link.phi,nu,up),TRUE)
      if(inherits(ans, "try-error") || ans$termcd != 1){
        p<-c(parameter[psi],parameter[-psi]) # initial guess : MLE
        ans<-try(nleqslv(p, scoreQBRbeta, jac=NULL, method="Newton",global="hook",control=list(trace=0,maxit=100,allowSingular=T),y,Xpsi,Z,score,cumulants,link.mu,link.phi,nu,up),TRUE)
      }
      if(inherits(ans, "try-error") || (ans$termcd != 1 && ans$termcd != 3)){
        cat(paste0("\nQBRU beta_",(psi-1)," convergence termination code : ",ans$termcd),file="warning.txt",append=T)
        QBRU <- NaN
      }else{QBRU <- ans$x[1]}
      if(ans$termcd == 3)cat(paste0("\nQBRU beta_",(psi-1),"  : No better point found (algorithm has stalled)"),file="warning.txt",append=T)
      
      # CI lower limit QBRE
      up<- u
      p<-c(parameter[psi],parameter[-psi])  # initial guess : MLE
      ans<-try(nleqslv(p, scoreQBRbeta, jac=NULL, method="Newton",control=list(trace=0,maxit=50,allowSingular=T),y,Xpsi,Z,score,cumulants,link.mu,link.phi,nu,up),TRUE)
      if(inherits(ans, "try-error") || ans$termcd != 1){
        p<-c(parameter[psi],parameter[-psi])  # initial guess : MLE
        ans<-try(nleqslv(p, scoreQBRbeta, jac=NULL, method="Newton",global="hook",control=list(trace=0,maxit=100,allowSingular=T),y,Xpsi,Z,score,cumulants,link.mu,link.phi,nu,up),TRUE)
      }
      if(inherits(ans, "try-error") || (ans$termcd != 1 && ans$termcd != 3)){
        cat(paste0("\nQBRL beta_",(psi-1)," convergence termination code : ",ans$termcd),file="warning.txt",append=T)
        QBRL <- NaN
      }else{QBRL <- ans$x[1]}
      if(ans$termcd == 3)cat(paste0("\nQBRL beta_",(psi-1),"  : No better point found (algorithm has stalled)"),file="warning.txt",append=T)
    }    
    
  }else{# (psi>q) gamma
    
    if(m==1)Zpsi <- iota else Zpsi <- cbind(Z[,(psi-q)],Z[,-(psi-q)])

    if(q==1 & m==1){
      # MBRE
      p<-parameter  # initial guess : MLE
      ans<-try(nleqslv(p, scoreQBRgamma, jac=NULL, method="Newton",control=list(trace=0,maxit=500),y,X,Zpsi,score,cumulants,link.mu,link.phi,nu,qnorm(0.5)),TRUE)
      if(inherits(ans, "try-error") || ans$termcd != 1){
        cat(paste0("\nQBRM beta_",(psi-1)," convergence termination code : ",ans$termcd),file="warning.txt",append=T)
        QBRM <- NaN
      }else{QBRM <- ans$x[2]}
      # CI upper limit QBRE
      up<- -u
      p<-parameter  # initial guess : MLE
      ans<-try(nleqslv(p, scoreQBRgamma, jac=NULL, method="Newton",control=list(trace=0,maxit=50,allowSingular=T),y,X,Zpsi,score,cumulants,link.mu,link.phi,nu,up),TRUE)
      if(inherits(ans, "try-error") || ans$termcd != 1){
        p<-parameter # initial guess : MLE
        ans<-try(nleqslv(p, scoreQBRgamma, jac=NULL, method="Newton",global="hook",control=list(trace=0,maxit=100,allowSingular=T),y,Xpsi,Z,score,cumulants,link.mu,link.phi,nu,up),TRUE)
      }
      if(inherits(ans, "try-error") || (ans$termcd != 1 && ans$termcd != 3)){
        cat(paste0("\nQBRU beta_",(psi-1)," convergence termination code : ",ans$termcd),file="warning.txt",append=T)
        QBRU <- NaN
      }else{QBRU <- ans$x[2]}
      if(ans$termcd == 3)cat(paste0("\nQBRU beta_",(psi-1),"  : No better point found (algorithm has stalled)"),file="warning.txt",append=T)
      
      # CI lower limit QBRE
      up<- u
      p<-parameter  # initial guess : MLE
      ans<-try(nleqslv(p, scoreQBRgamma, jac=NULL, method="Newton",control=list(trace=0,maxit=50,allowSingular=T),y,X,Zpsi,score,cumulants,link.mu,link.phi,nu,up),TRUE)
      if(inherits(ans, "try-error") || ans$termcd != 1){
        p<-parameter  # initial guess : MLE
        ans<-try(nleqslv(p, scoreQBRgamma, jac=NULL, method="Newton",global="hook",control=list(trace=0,maxit=100,allowSingular=T),y,X,Zpsi,score,cumulants,link.mu,link.phi,nu,up),TRUE)
      }
      if(inherits(ans, "try-error") || (ans$termcd != 1 && ans$termcd != 3)){
        cat(paste0("\nQBRL beta_",(psi-1)," convergence termination code : ",ans$termcd),file="warning.txt",append=T)
        QBRL <- NaN
      }else{QBRL <- ans$x[2]}
      if(ans$termcd == 3)cat(paste0("\nQBRL beta_",(psi-1),"  : No better point found (algorithm has stalled)"),file="warning.txt",append=T)
    }else{
      # MBRE
      p<-c(beta,gamma[(psi-q)],gamma[-(psi-q)]) # initial guess : MLE
      ans<-try(nleqslv(p, scoreQBRgamma, jac=NULL, method="Newton",control=list(trace=0,maxit=500),y,X,Zpsi,score,cumulants,link.mu,link.phi,nu,qnorm(0.5)),TRUE)
      if(inherits(ans, "try-error") || ans$termcd != 1){
        if(m==1)cat(paste("\nQBRM phi error convergence termination code : ",ans$termcd),file="warning.txt",append=T)
        else    cat(paste0("\nQBRM gamma_",((psi-q)-1)," convergence termination code : ",ans$termcd),file="warnings.txt",append=T)
        QBRM <- NaN
      }else{QBRM <- ans$x[(q+1)]}
      
      # CI upper limit QBRE
      up<- -u
      p<-c(beta,gamma[(psi-q)],gamma[-(psi-q)]) # initial guess : MLE
      ans<-try(nleqslv(p, scoreQBRgamma, jac=NULL, method="Newton",control=list(trace=0,maxit=50),y,X,Zpsi,score,cumulants,link.mu,link.phi,nu,up),TRUE)
      if(inherits(ans, "try-error") || ans$termcd != 1){
        p<-c(beta,gamma[(psi-q)],gamma[-(psi-q)]) # initial guess : MLE
        ans<-try(nleqslv(p, scoreQBRgamma, jac=NULL, method="Newton",global="hook",control=list(trace=0,maxit=100,allowSingular=T),y,X,Zpsi,score,cumulants,link.mu,link.phi,nu,up),TRUE)
      }
      if(inherits(ans, "try-error") || (ans$termcd != 1 && ans$termcd != 3)){
        if(m==1)cat(paste("\nQBRU phi error convergence termination code : ",ans$termcd),file="warning.txt",append=T)
        else    cat(paste0("\nQBRU gamma_",((psi-q)-1)," convergence termination code : ",ans$termcd),file="warning.txt",append=T)
        QBRU <- NaN
      }else{QBRU <- ans$x[(q+1)]}
      if(ans$termcd == 3)cat(paste0("\nQBRU gamma_",((psi-q)-1)," : No better point found (algorithm has stalled)"),file="warning.txt",append=T)
      
      # CI lower limit QBRE
      up<- u
      p<-c(beta,gamma[(psi-q)],gamma[-(psi-q)]) # initial guess : MLE
      ans<-try(nleqslv(p, scoreQBRgamma, jac=NULL, method="Newton",control=list(trace=0,maxit=50),y,X,Zpsi,score,cumulants,link.mu,link.phi,nu,up),TRUE)
      if(inherits(ans, "try-error") || ans$termcd != 1){
        p<-c(beta,gamma[(psi-q)],gamma[-(psi-q)]) # initial guess : MLE
        ans<-try(nleqslv(p, scoreQBRgamma, jac=NULL, method="Newton",global="hook",control=list(trace=0,maxit=100,allowSingular=T),y,X,Zpsi,score,cumulants,link.mu,link.phi,nu,up),TRUE)
      }
      if(inherits(ans, "try-error") || (ans$termcd != 1 && ans$termcd != 3)){
        if(m==1)cat(paste("\nQBRL phi error convergence termination code : ",ans$termcd),file="warning.txt",append=T)
        else    cat(paste0("\nQBRL gamma_",((psi-q)-1)," convergence termination code : ",ans$termcd),file="warning.txt",append=T)
        QBRL <- NaN
      }else{QBRL <- ans$x[(q+1)]}
      if(ans$termcd == 3)
        cat(paste0("\nQBRL gamma_",((psi-q)-1)," : No better point found (algorithm has stalled)"),file="warning.txt",append=T)
    }
  }
  return(c(QBRM,QBRL,QBRU))
}

link.matrix<-function(parameter,X,Z,link.mu,link.phi){
  
  numpar<-length(parameter)
  q<-ncol(X)
  m<-ncol(Z)
  nobs<-nrow(X)
  
  iota =  matrix(rep(1,nobs),nobs,1)
  
  eta <- X%*%parameter[1:q]
  if(link.mu=="logit"){
    mu <- exp(eta) / (1.0 + exp(eta))
    T <-  exp(eta) / (1.0 + exp(eta))^2   # d mu/d eta = 1/g'(mu)=mu(1-mu)
    S <- iota-2*mu	                       # d2 mu/d eta d mu = - g''(mu)/g'(mu)^2 = 1-2mu
  }
  if(link.mu=="identity"){
    mu <- eta
    T <- iota         # d mu/d eta = 1
    S <- rep(0,nobs)  # d2 mu/d eta d mu = 0
  }
  
  delta <- Z%*%parameter[(q+1):numpar]
  if(link.phi=="log"){
    phi<-exp(delta)
    H <- phi  # d phi/d delta = phi
    Q <- iota # d2 phi/d delta d phi = 1
  }
  if(link.phi=="identity"){
    phi<-delta;
    H <- iota         # d phi/d delta = 1
    Q <- rep(0,nobs)  # d2 phi/d delta d phi = 0
  }
  
  return(list(T,H,S,Q))
}

print_result<-function(estimates,intervals,X,Z,q,m,link.phi,CL){
  colnames(estimates)<-c("ML","MBR")
  colnames(intervals)<-c(paste0(((1-CL)/2)*100,"%"),paste0((CL+(1-CL)/2)*100,"%"),
                         paste0(((1-CL)/2)*100,"%"),paste0((CL+(1-CL)/2)*100,"%"),
                         paste0(((1-CL)/2)*100,"%"),paste0((CL+(1-CL)/2)*100,"%"))
  options(width=300,digits=4)
  
  if(m==1){
    if(link.phi=="identity"){
      rownames(estimates)<-c(colnames(X),'\u03c6')
      rownames(intervals)<-c(colnames(X),'\u03c6')
    }
    if(link.phi=="log"){
      rownames(estimates)<-c(colnames(X),paste0("log(",'\u03c6',")"))
      rownames(intervals)<-c(colnames(X),paste0("log(",'\u03c6',")"))
    }
    cat("\n")
    cat("\n=========================  Point estimates   =========================")
    cat("\n")
    estimates<-as.table(estimates)
    print(estimates, na.print = "-" , quote = FALSE)
    cat("\n=========================  Confidence limits  =========================")
    cat(paste("\n       ",(CL+(1-CL)/2)*100,"% CL - One-sided    and   ",CL*100,"% CL - Two-sided") )
    cat("\n")
    cat("                    ML           MBR         QBR")
    cat("\n")
    cat("=======================================================================")
    cat("\n")
    intervals<-as.table(intervals)
    print(intervals, na.print = "-" , quote = FALSE)
  }
  else{
    estimatesmean<-estimates[1:q,]
    intervalsmean<-intervals[1:q,]
    rownames(estimatesmean)<-c(colnames(X))
    rownames(intervalsmean)<-c(colnames(X))
    colnames(estimatesmean)<-c("ML","MBR")
    colnames(intervalsmean)<-c(paste0(((1-CL)/2)*100,"%"),paste0((CL+(1-CL)/2)*100,"%"),
                               paste0(((1-CL)/2)*100,"%"),paste0((CL+(1-CL)/2)*100,"%"),
                               paste0(((1-CL)/2)*100,"%"),paste0((CL+(1-CL)/2)*100,"%"))
    
    
    
    estimatesprecision<-estimates[(q+1):(q+m),]
    intervalsprecision<-intervals[(q+1):(q+m),]
    rownames(estimatesprecision)<-c(colnames(Z))
    rownames(intervalsprecision)<-c(colnames(Z))
    colnames(estimatesprecision)<-c("ML","MBR")
    colnames(intervalsprecision)<-c(paste0(((1-CL)/2)*100,"%"),paste0((CL+(1-CL)/2)*100,"%"),
                                    paste0(((1-CL)/2)*100,"%"),paste0((CL+(1-CL)/2)*100,"%"),
                                    paste0(((1-CL)/2)*100,"%"),paste0((CL+(1-CL)/2)*100,"%"))
    
    cat("\nCoefficients (mean model with logit link):\n")
    
    cat("\n=========================  Point estimates   =========================")
    cat("\n")
    estimatesmean<-as.table(estimatesmean)
    print(estimatesmean, na.print = "-" , quote = FALSE)
    cat("\n=========================  Confidence limits  =========================")
    cat(paste("\n       ",(CL+(1-CL)/2)*100,"% CL - One-sided    and   ",CL*100,"% CL - Two-sided") )
    cat("\n")
    cat("                    ML           MBR         QBR")
    cat("\n")
    cat("=======================================================================")
    cat("\n")
    intervalsmean<-as.table(intervalsmean)
    print(intervalsmean, na.print = "-" , quote = FALSE)
    
    cat("\nCoefficients (precision model with log link):\n")
    
    cat("\n=========================  Point estimates   =========================")
    cat("\n")
    estimatesprecision<-as.table(estimatesprecision)
    print(estimatesprecision, na.print = "-" , quote = FALSE)
    cat("\n=========================  Confidence limits  =========================")
    cat(paste("\n       ",(CL+(1-CL)/2)*100,"% CL - One-sided    and   ",CL*100,"% CL - Two-sided") )
    cat("\n")
    cat("                    ML           MBR         QBR")
    cat("\n")
    cat("=======================================================================")
    cat("\n")
    intervalsprecision<-as.table(intervalsprecision)
    print(intervalsprecision, na.print = "-" , quote = FALSE)
  }
  file.show("warning.txt", delete.file = T, pager = "console", encoding = "")
}  

###########################################  
# functions for beta distribution
###########################################  

beta_score<-function(parameter,v,Nx,Nz,link.phi,nu){
  numpar<-length(parameter)
  y<-v
  X<-Nx
  Z<-Nz
  q<-ncol(X)
  m<-ncol(Z)
  numpar<-q+m
  beta <-parameter[1:q]
  gamma<-parameter[(q+1):numpar]
  
  nobs<-length(y)
  iota=rep(1,nobs)
  
  eta <- X%*%beta 
  mu <- exp(eta) / (1.0 + exp(eta))
  
  delta <- Z%*%gamma
  if(link.phi=="identity")phi<-delta
  else phi<-exp(delta);
  
  ystar = log( y / (iota - y) );
  mustar = psigamma(mu*phi, 0) - psigamma((iota - mu)*phi, 0);
  ydag = log(iota-y);
  mudag =  -psigamma(phi,0)+psigamma((iota-mu)*phi,0);
  
  cbind(phi*(ystar-mustar), mu*(ystar-mustar)+ydag-mudag)
}

beta_cumulants<-function(parameter,v,Nx,Nz,link.phi,nu){
  numpar<-length(parameter)
  y<-v
  X<-Nx
  Z<-Nz
  q<-ncol(X)
  m<-ncol(Z)
  beta <-parameter[1:q]
  gamma<-parameter[(q+1):numpar]
  
  nobs<-length(y)
  iota=rep(1,nobs)
  
  eta <- X%*%beta 
  mu <- exp(eta) / (1.0 + exp(eta))
  
  delta <- Z%*%gamma
  if(link.phi=="identity")phi<-delta
  else phi<-exp(delta)
  
  K22   = phi^2*(     psigamma(mu*phi,1) +          psigamma((1-mu)*phi,1) )                 #K22=K_mu,mu
  K21   = phi  *(mu  *psigamma(mu*phi,1) - (1-mu)  *psigamma((1-mu)*phi,1) )                 #K21=K_mu,phi
  K20   =        mu^2*psigamma(mu*phi,1) + (1-mu)^2*psigamma((1-mu)*phi,1) - psigamma(phi,1) #K20=K_phi,phi
  
  K33   = phi^3*(     psigamma(mu*phi,2) -          psigamma((1-mu)*phi,2) )                 #K33=K_mu,mu,mu
  K32   = phi^2*(mu  *psigamma(mu*phi,2) + (1-mu)  *psigamma((1-mu)*phi,2) )                 #K32=K_mu,mu,phi
  K31   = phi*  (mu^2*psigamma(mu*phi,2) - (1-mu)^2*psigamma((1-mu)*phi,2) )                 #K31=K_mu,phi,phi
  K30   =        mu^3*psigamma(mu*phi,2) + (1-mu)^3*psigamma((1-mu)*phi,2) - psigamma(phi,2) #K30=K_phi,phi,phi
  
  K44   = phi^4*(     psigamma(mu*phi,3) +          psigamma((1-mu)*phi,3) )                 #K44=K_mu,mu,mu,mu   
  K43   = phi^3*(mu  *psigamma(mu*phi,3) - (1-mu)  *psigamma((1-mu)*phi,3) )                 #K43=K_mu,mu,mu,phi
  K42   = phi^2*(mu^2*psigamma(mu*phi,3) + (1-mu)^2*psigamma((1-mu)*phi,3) )                 #K42=K_mu,mu,phi,phi
  K41   = phi  *(mu^3*psigamma(mu*phi,3) - (1-mu)^3*psigamma((1-mu)*phi,3) )                 #K41=K_mu,phi,phi,phi  
  K40   =        mu^4*psigamma(mu*phi,3) + (1-mu)^4*psigamma((1-mu)*phi,3) - psigamma(phi,3) #K40=K_phi,phi,phi,phi
  
  Kmu2  = rep(0,nobs)                                                       #Kmu2=K_mu,mumu
  Kmu1  = phi  *(    psigamma(mu*phi,1) +          psigamma((1-mu)*phi,1) ) #Kmu1=K_mu,muphi
  Kmu0  = rep(0,nobs)                                                       #Kmu0=K_mu,phiph
  Kphi2 = rep(0,nobs)                                                       #Kphi2=K_phi,mumu
  Kphi1 =         mu  *psigamma(mu*phi,1) - (1-mu)  *psigamma((1-mu)*phi,1) #Kphi1=K_phi,muphi
  Kphi0 = rep(0,nobs)                                                       #Kphi0=K_phi,phiphi
  
  return(list(K22,K21,K20,K33,K32,K31,K30,K44,K43,K42,K41,K40,Kmu2,Kmu1,Kmu0,Kphi2,Kphi1,Kphi0))
}

beta_mle<-function(formula,data,link.mu,link.phi,nu){
  
  r<-betareg(formula, data, link.phi=link.phi, x = TRUE, y = TRUE)
  
  if(r$optim$convergence != 0){
    print("no convergence in betareg")
  }
  thetahat<-as.vector(r$optim$par)
  y<-r$y
  X<-r$x$mean
  Z<-r$x$precision
  
  return(list(thetahat,y,X,Z))
}

###########################################  
# functions for Student-t distribution
###########################################  

t_score<-function(parameter,v,Nx,Nz,link.phi,nu){
  
  y<-v
  X<-Nx
  Z<-Nz
  q<-ncol(X)
  m<-ncol(Z)
  numpar<-q+m
  
  beta <-parameter[1:q]
  gamma<-parameter[(q+1):numpar]
  
  nobs<-length(y)
  iota=matrix(rep(1,nobs),nobs,1)
  
  eta <- X%*%beta 
  mu <- eta
  
  delta <- Z%*%gamma
  if(link.phi=="log"){
    phi<-exp(delta)
  }
  if(link.phi=="identity"){
    phi<-delta;
  }
  
  z = ( y - mu ) / phi
  W = (nu + 1)/(nu + z^2)
  
  cbind((W*z)/phi,(W*(z^2) - iota)/phi)
}

t_cumulants<-function(parameter,v,Nx,Nz,link.phi,nu){
  numpar<-length(parameter)
  y<-v
  X<-Nx
  Z<-Nz
  
  q<-ncol(X)
  # m<-ncol(Z)
  # beta <-parameter[1:q]
  gamma<-parameter[(q+1):numpar]
  
  nobs<-length(y)
  iota=matrix(rep(1,nobs),nobs,1)
  
  # eta <- X%*%beta 
  # mu <- eta
  
  delta<-Z%*%gamma
  if(link.phi=="log"){
    phi<-exp(delta)
  }
  if(link.phi=="identity"){
    phi<-delta;
  }
  
  del00101 = (6*(nu+1))/((nu+3)*(nu+5))
  del00103 = (6*(3*nu-5))/((nu+3)*(nu+5))
  del01000 = -(nu+1)/(nu+3)
  del01002 = (3-nu)/(nu+3)
  del11001 = ((nu+1)*(nu-1))/((nu+3)*(nu+5))
  del11003 = 3*(nu+1)*(nu-5)/((nu+3)*(nu+5))
  del20000 = (nu+1)/(nu+3)
  del20002 = 3*(nu+1)/(nu+3)
  del21000 = -(nu+1)^3*(nu+2)/(nu*(nu+3)*(nu+5)*(nu+7))
  del30001 = -3*(nu+1)^2/((nu+3)*(nu+5))
  del30003 = -(15*(nu+1)^2)/((5+nu)*(3+nu))
  del40000 = -3*del21000
  del40002 = 15*(nu+1)^3/((nu+3)*(nu+5)*(nu+7)) 
  del40004 = (105*(nu+1)^3)/((7+nu)*(5+nu)*(3+nu))
  
  K22 = (del20000/phi^2)*iota     #K22=K_mu,mu
  K21 = 0*iota                    #K21=K_mu,phi
  K20 = ((del20002-1)/phi^2)*iota #K20=K_phi,phi
  
  K33 = 0*iota                 #K33=K_mu,mu,mu 
  K32 = (2*del11001/phi^3)     #K32=K_mu,mu,phi
  K31 = 0*iota                 #K31=K_mu,phi,phi
  K30 = (2*(del11003+1)/phi^3) #K30=K_phi,phi,phi 
  
  K44 = ((del40000-3*del20000^2)/phi^4) #K44=K_mu,mu,mu,mu
  
  # del20001= # this quantity was not obtained since K_mu,phi,phi,phi is not necessary because mu and phi are orthogonal
  # del30000= # this quantity was not obtained since K_mu,mu,mu,phi   is not necessary because mu and phi are orthogonal
  # del40001= # this quantity was not obtained since K_mu,mu,mu,phi   is not necessary because mu and phi are orthogonal
  # K43 = (del30000+del40001)*iota                        #K43=K_mu,mu,mu,phi
  # K42 = ((2*del30001+del40002-del01000*del01002)/phi^4) #K42=K_mu,mu,phi,phi
  # K41 = (del40001+3*del30001+3*del20001)*iota           #K41=K_mu,phi,phi,phi
  
  K43 = 1000*iota # this quantity is not necessary because mu and phi are orthogonal
  K42 = 1000*iota # this quantity is not necessary because mu and phi are orthogonal
  K41 = 1000*iota # this quantity is not necessary because mu and phi are orthogonal
  
  K40 = ((del40004+4*del30003+12*del20002-3*del20002^2-6)/phi^4) #K40=K_phi,phi,phi,phi
  
  Kmu2  = 0*iota                          #Kmu2=K_mu,mumu
  Kmu1  = (del01000-del11001)/phi^3       #Kmu1=K_mu,muphi
  Kmu0  = 0*iota                          #Kmu0=K_mu,phiphi     
  Kphi2 = (del00101/phi^3)*iota           #Kphi2=K_phi,mumu   
  Kphi1 = 0*iota                          #Kphi1=K_phi,muphi
  Kphi0 = ((4*del01002+del00103-2)/phi^3) #Kphi0=K_phi,phiphi
  
  return(list(K22,K21,K20,K33,K32,K31,K30,K44,K43,K42,K41,K40,Kmu2,Kmu1,Kmu0,Kphi2,Kphi1,Kphi0))
}

t_mle<-function(formula,data,link.mu,link.phi,nu){

  r <- ssym.l(formula , family='Student', xi=nu, link.phi="log", data=data)
  
  betahat <- coefficients(r)$mu
  if(link.phi=="identity"){
    gammahat <- sqrt( exp(coefficients(r)$phi) ) #adjust due to different parametrization
    thetahat <- c(betahat, gammahat)
  }
  if(link.phi=="log"){
    gammahat<-coefficients(r)$phi/2 #adjust due to different parametrization
    thetahat <- c(betahat, gammahat)
  }
  
  y<-r$y
  X<-r$model.matrix.mu
  Z<-r$model.matrix.phi 
  
  return(list(thetahat,y,X,Z))
}


###########################################  
# main functions
###########################################  

main<-function (formula,data,CL,mle,score,cumulants,link.mu,link.phi,nu){
  start.time <- Sys.time()
  
  prob<-CL+(1-CL)/2
  u<-qnorm(prob)
  
  MLE <- mle(formula,data,link.mu,link.phi,nu)
  
  thetahat <- MLE[[1]]
  y        <- MLE[[2]]
  X        <- MLE[[3]]
  Z        <- MLE[[4]]
  
  LINK<-link.matrix(thetahat,X,Z,link.mu,link.phi)
  T<-as.vector(LINK[[1]])
  H<-as.vector(LINK[[2]])
  
  q<-ncol(X)
  m<-ncol(Z)
  numpar<-q+m
  
  K <- cumulants(thetahat,y,X,Z,link.phi,nu)
  K22 = as.vector(K[[1]])
  K21 = as.vector(K[[2]])
  K20 = as.vector(K[[3]])
  
  fisher = matrix(0,numpar,numpar)   # kappa_(theta_i,theta_j)
  fisher[1:q,1:q]                  = crossprod(X,T*T*K22*X)     # kappa_betarbetas
  fisher[1:q,(q+1):numpar]         = crossprod(X,T*H*K21*Z)     # kappa_betargammas
  fisher[(q+1):numpar,1:q]         = t(fisher[1:q,(q+1):numpar])# kappa_gammarbetas
  fisher[(q+1):numpar,(q+1):numpar]= crossprod(Z,H*H*K20*Z)     # kappa_gammargammas
  
  stderrors <- sqrt(diag(solve(fisher))); 
  
  if(is.nan(stderrors) || is.infinite(stderrors)){
    print("stderrors error")
  }
  else{
    MLEM<-thetahat
    Aux <- matrix(rep(u,numpar),nrow=1,ncol=numpar,byrow=F)
    Aux <-sweep(Aux,MARGIN=2,stderrors,`*`)
    MLEL<-sweep(-Aux,MARGIN=2,thetahat,`+`) 
    MLEU<-sweep(+Aux,MARGIN=2,thetahat,`+`) 
    
    lengthEMV <-	MLEU - MLEL;
  }
  #---------------------------------------------------
  # Quantile bias reduction estimate 
  #---------------------------------------------------
  if(q==1 & m==1){ # iid sample
    tab<-matrix(NA_real_,numpar,3)
    for(i in 1:numpar){
      psi<<-i
      tab[i,]<-QBRE(parameter = thetahat, v = y, Nx = X, Nz = Z, score = score, cumulants = cumulants ,link.mu = link.mu, link.phi = link.phi, nu = nu,  u = u)  
    }
  }else{
    # Alocar Processadores ----
    closeAllConnections()
    ncores <- detectCores()
    cl <- makeCluster(ncores)
    registerDoParallel(cl, cores=ncores)
    cat("\ncores = ",getDoParWorkers(),"\n")
    on.exit(stopCluster(cl))
    
    tab <- foreach(psi = 1:numpar, .combine = rbind, .packages=c('doParallel'), .export = ls(globalenv())) %dopar% {
      result<-QBRE(parameter = thetahat, v = y, Nx = X, Nz = Z, score = score, cumulants = cumulants ,link.mu = link.mu, link.phi = link.phi, nu = nu,  u = u)
      return(result)
    }
  }
  QBRL <- tab[,2]
  QBRU <- tab[,3]
  lengthQBR <-	QBRU - QBRL
  
  QBRM <- tab[,1]
  thetatil<-QBRM
  
  LINK<-link.matrix(thetatil,X,Z,link.mu,link.phi)
  T<-as.vector(LINK[[1]])
  H<-as.vector(LINK[[2]])
  
  K <- cumulants(thetatil,y,X,Z,link.phi,nu)
  K22 = as.vector(K[[1]])
  K21 = as.vector(K[[2]])
  K20 = as.vector(K[[3]])
  
  fisher = matrix(0,numpar,numpar)   # kappa_(theta_i,theta_j)
  fisher[1:q,1:q]                  = crossprod(X,T*T*K22*X)     # kappa_betarbetas
  fisher[1:q,(q+1):numpar]         = crossprod(X,T*H*K21*Z)     # kappa_betargammas
  fisher[(q+1):numpar,1:q]         = t(fisher[1:q,(q+1):numpar])# kappa_gammarbetas
  fisher[(q+1):numpar,(q+1):numpar]= crossprod(Z,H*H*K20*Z)     # kappa_gammargammas
  
  stderrorstil <- sqrt(diag(solve(fisher))); 
  
  if(is.nan(stderrorstil) || is.infinite(stderrorstil)){
    print("stderrorstil error")
    MBRL<-matrix(NaN,1,numpar) 
    MBRU<-matrix(NaN,1,numpar)  
    lengthMBR<-matrix(NaN,1,numpar) 
  }else{
    Aux <- matrix(rep(u,numpar),nrow=1,ncol=numpar,byrow=F)
    Aux <- sweep(Aux,MARGIN=2,stderrorstil,`*`)
    MBRL<-sweep(-Aux,MARGIN=2,thetatil,`+`) 
    MBRU<-sweep(+Aux,MARGIN=2,thetatil,`+`)
    lengthMBR<-MBRU - MBRL
  }  
  
  estimates<-t(rbind(MLEM, QBRM))
  intervals<-t(rbind(MLEL, MLEU, MBRL, MBRU, QBRL, QBRU))
  
  print_result(estimates,intervals,X,Z,q,m,link.phi,CL)
  
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  cat("\n")
  print(time.taken)
  
} 

hoaci<-function (formula, data, CL, subset, na.action, weights, offset,
                 mle,
                 score,
                 cumulants,
                 # link.mu = c("logit", "probit", "cloglog", "cauchit", "log", "loglog","identity"),
                 # link.phi = c("identity", "log", "sqrt"),
                 link.mu = c("identity","logit"),
                 link.phi = c("identity", "log"),
                 # type = c("exp", "gamma", "beta", "Student", "normal"),
                 type = c("beta", "Student"),
                 nu = nu,
                 control = betareg.control(...), 
                 model = TRUE, 
                 y = TRUE, 
                 x = FALSE,
                 ...){
  cl <- match.call()
  if (missing(data)) data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action", "weights", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  
  oformula <- as.formula(formula)
  formula <- as.Formula(formula)
  
  if (length(formula)[2L] < 2L) {
    formula <- as.Formula(formula(formula), ~1)
    simple_formula <- TRUE
  }
  else {
    if (length(formula)[2L] > 2L) {
      formula <- Formula(formula(formula, rhs = 1:2))
      warning("formula must not have more than two RHS parts")
    }
    simple_formula <- FALSE
  }
  
  mf$formula <- formula
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- terms(formula, data = data)
  mtX <- terms(formula, data = data, rhs = 1L)
  mtZ <- delete.response(terms(formula, data = data, rhs = 2L))
  Y <- model.response(mf, "numeric")
  X <- model.matrix(mtX, mf)
  Z <- model.matrix(mtZ, mf)
  
  if (length(Y) < 1)  stop("empty model")
  
  if(type=="beta")
    if (!(min(Y) > 0 & max(Y) < 1)) stop("invalid dependent variable, all observations must be in (0, 1)")
  
  
  n <- length(Y)
  # if (is.character(link))link <- match.arg(link)
  if (is.character(link.mu))link.mu <- match.arg(link.mu)
  
  if (is.null(link.phi)) link.phi <- if (simple_formula) "identity"  else "log"
  if (is.character(link.phi)) link.phi <- match.arg(link.phi, c("identity", "log", "sqrt"))
  weights <- model.weights(mf)
  if (is.null(weights)) weights <- 1
  if (length(weights) == 1) weights <- rep.int(weights, n)
  weights <- as.vector(weights)
  names(weights) <- rownames(mf)
  
  expand_offset <- function(offset) {
    if (is.null(offset)) 
      offset <- 0
    if (length(offset) == 1) 
      offset <- rep.int(offset, n)
    as.vector(offset)
  }
  offsetX <- expand_offset(model.offset(model.part(formula, data = mf, rhs = 1L, terms = TRUE)))
  offsetZ <- expand_offset(model.offset(model.part(formula, data = mf, rhs = 2L, terms = TRUE)))
  if (!is.null(cl$offset))  offsetX <- offsetX + expand_offset(mf[, "(offset)"])
  offset <- list(mean = offsetX, precision = offsetZ)
  
  cat(paste("\nWait....") )
  
  main(formula, data, CL, mle, score, cumulants, link.mu, link.phi, nu)
}
