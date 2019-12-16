# Pinheiro, Ferrari and Medeiros (2019) - Higher-order approximate confidence intervals
#
# Applications using the hoaci function to obtain confidence intervals of the paper 
# for the symmetric Student-t and beta distributions (examples 6 and 7).
# ML, MBR, and QBR denote the usual Wald-type CIs, i.e. CIs that use the asymptotic normality of MLE, 
# the adjusted Wald-type CIs (MLEs are replaced by median bias reduced estimates), and the
# modified CIs proposed in this paper, respectively.

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

source("HOACI-functions.R") 

#***********************************************************************
#***********************************************************************
#***********************************************************************
# symmetric regression example

# Mirhosseini, H., Tan, C.P. (2010). 
# Discrimination of orange beverage emulsions with different formulations using multivariate analysis. 
# Journal of the Science of Food and Agriculture 90, 1308-1316.
x<-scan("orange.txt")
orange=matrix(x,length(x)/3,3)
orange
density<-orange[,1]
arabic_gum<-orange[,2]
orange_oil <- orange[,3]

density    <- density * 1000
arabic_gum <- arabic_gum / 100
orange_oil <- orange_oil / 100
orange<-data.frame(orange)
names(orange)<-c("density","arabic_gum","orange_oil")
orange
length(density)

# link.phi="identity"
hoaci(density ~ arabic_gum + orange_oil,            data=orange, CL=0.95, mle=t_mle, score=t_score, cumulants=t_cumulants, link.mu="identity", link.phi="identity", type="Student", nu=3)
# link.phi="log"
hoaci(density ~ arabic_gum + orange_oil,            data=orange, CL=0.95, mle=t_mle, score=t_score, cumulants=t_cumulants, link.mu="identity", link.phi="log", type="Student", nu=3)
hoaci(density ~ arabic_gum + orange_oil|arabic_gum, data=orange, CL=0.95, mle=t_mle, score=t_score, cumulants=t_cumulants, link.mu="identity", link.phi="log", type="Student", nu=3)
hoaci(density ~ arabic_gum + orange_oil|orange_oil, data=orange, CL=0.95, mle=t_mle, score=t_score, cumulants=t_cumulants, link.mu="identity", link.phi="log", type="Student", nu=3)

# iid
hoaci(density ~ 1,                                  data=orange, CL=0.95, mle=t_mle, score=t_score, cumulants=t_cumulants, link.mu="identity", link.phi="identity", type="Student", nu=3)

#***********************************************************************
#***********************************************************************
#***********************************************************************
# beta regression examples
require(betareg)
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
data("ReadingSkills", package = "betareg")
r<-betareg(accuracy ~ dyslexia * iq | dyslexia + iq, data = ReadingSkills, x=TRUE)
summary(r)
hoaci(accuracy ~ dyslexia * iq | dyslexia + iq, data = ReadingSkills, CL=0.95, mle=beta_mle, score=beta_score, cumulants=beta_cumulants, link.mu="logit", link.phi="log", type="beta", nu=NULL)

# iid
r<-betareg(accuracy ~ 1, data = ReadingSkills, x=TRUE)
summary(r)
hoaci(accuracy ~ 1, data = ReadingSkills, CL=0.95, mle=beta_mle, score=beta_score, cumulants=beta_cumulants, link.mu="logit", link.phi="identity", type="beta", nu=NULL)

#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
data("MockJurors", package = "betareg")
hoaci(confidence ~ verdict + conflict | verdict + conflict, data = MockJurors, CL=0.95, mle=beta_mle, score=beta_score, cumulants=beta_cumulants, link.mu="logit", link.phi="log", type="beta", nu=NULL)

#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
data("StressAnxiety", package = "betareg")
betareg(anxiety ~ stress | stress, data = StressAnxiety, hessian = TRUE)
hoaci(anxiety ~ stress | stress, data = StressAnxiety, CL=0.95, mle=beta_mle, score=beta_score, cumulants=beta_cumulants, link.mu="logit", link.phi="log", type="beta", nu=NULL)

#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
data("FoodExpenditure", package = "betareg")
# Bayer & Cribari-Neto (2013) - JSPI - Bartlett corrections in beta regression models - food expenditure data
hoaci(I(food/income) ~ income + persons, data = FoodExpenditure, CL=0.95, mle=beta_mle, score=beta_score, cumulants=beta_cumulants, link.mu="logit", link.phi="log", type="beta", nu=NULL)

#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
#------------------------------------------------------------------------
# Ferrari, S. L. P. (2017). 
# Beta Regression. 
# Wiley StatsRef: Statistics Reference Online, 1-5. 1st ed.: Wiley.
x<-scan("gas.txt")
gas=matrix(x,length(x)/2,2,byrow=T)
y<-gas[,1]
x1<-gas[,2] # natural logarithm of computed power

hoaci(y~x1, CL=0.95, mle=beta_mle, score=beta_score, cumulants=beta_cumulants, link.mu="logit", link.phi="log", type="beta", nu=NULL)
hoaci(y~x1|x1, CL=0.95, mle=beta_mle, score=beta_score, cumulants=beta_cumulants, link.mu="logit", link.phi="log", type="beta", nu=NULL)

#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
# Draper, N. R. & Smith, H. (1998). 
# Applied Regression Analysis. 
# 3rd edition. New York: Wiley
# pages: 101-102
# York: Wiley.
dados<-read.table("draper.txt")
dados<-matrix(dados$V1,34,2)
y<-dados[,1]
x1<-dados[,2]

hoaci(y~x1, CL=0.95, mle=beta_mle, score=beta_score, cumulants=beta_cumulants, link.mu="logit", link.phi="log", type="beta", nu=NULL)

#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
# Brownlee, K. A. (1965). 
# Statistical Theory and Methodology in Science and Engineering. 
# 2nd edition. London: John Wiley & Sons
# page 454
dados<-read.table("brownlee.dat")
y <-dados[,1]
x1<-dados[,2]#air
x2<-dados[,3]#water temperature
x3<-dados[,4]#acid concentration

hoaci(y~x1+x2, CL=0.95, mle=beta_mle, score=beta_score, cumulants=beta_cumulants, link.mu="logit", link.phi="log", type="beta", nu=NULL)

#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
# Kieschnick, R. & McCullough, B. D.(2003)
# Regression analysis of variates observed on (0,1): percentages, proportions and fractions
# Statistical Modelling, 3, 193-213
dados<-read.table("voting.dat")
dados<-matrix(dados[,1],51,6)

y <-dados[,1] # % of Votes for Bush - per_bush
x1<-dados[,2] # % of Civilian Labor Force Unemployed - per_clfu
x2<-dados[,3] # Unemployment Rate - unemprate
x3<-dados[,4] # Males per 100 females - mp100f
x4<-dados[,5] # Males > Age of 18 per 100 females - mgt80p100f
x5<-dados[,6] # % of Population Over Age of 65 - pergt65

hoaci(y~x1+x3+x5, CL=0.95, mle=beta_mle, score=beta_score, cumulants=beta_cumulants, link.mu="logit", link.phi="log", type="beta", nu=NULL)
