/**************************************************************************************************
 * Monte Carlo simulation to evaluate the finite-sample performance of confidence intervals based on
 * the alpha-quantile modified score (modified CIs).
 * For comparison, we included results for Wald-type CIs, i.e. CIs that use the asymptotic normality of
 * MLE (usual CIs), and those that are constructed similarly to Wald-type CIs with median bias reduced
 * estimates of all the parameters in place of MLEs (adjusted CIs).
 *
 *  Gamma(mu,phi) iid samples with density function
 *	f(x;mu,phi)= 1 / Gamma[phi] (phi/mu)^phi x^(phi - 1) Exp[- phi x / mu]
 *  E(X)=mu
 *  GLM parameterization - mu and phi are orthogonal
 *
 *
 * Download Ox at https://www.doornik.com/download.html
 *
 * AUTHOR: Eliane Cantinho Pinheiro.
 *
 * DATE: 1/12/2019.
 *
 **************************************************************************************************/

 #include <oxstd.h>
 #include <oxprob.h>
 #import <maximize>
 #include <oxdraw.h>
 #include <oxfloat.h>
 #import <maxsqp> 
 #import <solvenle>

/* Global variables */
static decl Y;            /* Random sample */
static decl upi;
static decl up;
 
/* Constants */
const decl nobs = 15;     /* Sample size */
const decl iRep = 100000; /* Number of replicates for MC */

const decl mu = 10;       /* true value of mu */
const decl phi = 3;       /* true value of phi*/
 
const decl prob = <0.005, 0.010, 0.020, 0.025, 0.030, 0.040, 0.050, 0.060, 0.070, 0.080, 0.090, 0.100, 0.150, 0.200, 0.250, 0.300, 0.400,
                    0.500,
					0.600, 0.700, 0.750, 0.800, 0.850, 0.900, 0.910, 0.920, 0.930, 0.940, 0.950, 0.960, 0.970, 0.975, 0.980, 0.990, 0.995>; 

/* Log-likelihood funtion */
floglik(const vtheta, const adFunc, const avScore, const amHess)
{
  decl mu = vtheta[0];
  decl phi = vtheta[1];

  decl iota =  ones(nobs,1);

  adFunc[0] = double( iota'*( (-log(gammafact(phi)) + phi*log(phi/mu) )*iota + (phi-1)*log(Y) - phi*Y/mu  )  );
                                                                                                              /* log-likelihood function */

   if(avScore)
   {
      (avScore[0])[0] = double( iota'*( -(phi/mu)*iota + (phi/mu^2)*Y ) );
      (avScore[0])[1] = double( iota'*(  (1 + log(phi/mu) - polygamma(phi,0) )*iota - Y/mu + log(Y)  ) );
   }

   if( isnan(adFunc[0]) || isdotinf(adFunc[0]) )
       return 0;
   else
       return 1; /* 1 is success */
 }

QBREmu(const avF, const vX)
{
  decl mu=vX[0];
  decl phi= vX[1];
  decl up=upi;

  decl k2mu = nobs*(phi/mu^2);
  decl k3mu = nobs*2*(phi/mu^3);
  decl k4mu = nobs*6*(phi/mu^4);

  decl M =
  //- k1mu
  - sqrt(k2mu)*up - k3mu/k2mu/6*(up^2-1)- k4mu/k2mu^(3/2)/24*(up^3-3*up)
  + k3mu^2/k2mu^(5/2)/36*(2*up^3 - 5*up);

  decl iota =  ones(nobs,1);
  decl Umu = double( iota'*( -(phi/mu)*iota + (phi/mu^2)*Y ) );
  decl Uphi = double( iota'*(  (1 + log(phi/mu) - polygamma(phi,0) )*iota - Y/mu + log(Y)  ) );

  /* Solving the equation : U_p(mu) + M_mu = 0   */
  avF[0] = (Umu  +  M)|(Uphi);

  return 1;
}

QBREphi(const avF, const vX)
{
  decl mu = vX[0];
  decl phi= vX[1];
  decl up=upi;

  decl k1phi = 1/(2*phi);
  decl k2phi = nobs*(-1/phi   + polygamma(phi,1));
  decl k3phi = nobs* (1/phi^2 + polygamma(phi,2));
  decl k4phi = nobs*(-2/phi^3 + polygamma(phi,3));

  decl M = - k1phi - sqrt(k2phi)*up - k3phi/k2phi/6*(up^2-1)- k4phi/k2phi^(3/2)/24*(up^3-3*up)
  + k3phi^2/k2phi^(5/2)/36*(2*up^3 - 5*up);

  decl iota =  ones(nobs,1);
  decl Umu  = double( iota'*( -(phi/mu)*iota + (phi/mu^2)*Y ) );
  decl Uphi = double( iota'*(  (1 + log(phi/mu) - polygamma(phi,0) )*iota - Y/mu + log(Y)  ) );

  /* Solving the equation : U_p(phi) + M_phi = 0   */
  avF[0] = (Umu)|(Uphi + M);

  return 1;
}

/*  Starting the program */
main()
{

  /* Declaring variables */

  up=quann(prob);

  /* System variables */
  decl dExecTime, fail, ci, i;
  
  decl failMLE, failstderrors, failQBR_L, failQBR_U, failQBR_phi;
  
  /* Model variables */
  decl iota, phiinicial, muinicial, Hess, converge;

  decl vthetahat, muhat, phihat, dfunchat;
  decl Imumu, Imuphi, Iphiphi, fisher, fisherinv, stderrors;

  decl muhatMBR, phihatMBR, stderrorsMBR;

  decl vthetaQBR, muhatQBR;
  decl convergeQBR;

  /* Result variables */
  decl vEMV, coverageEMV_L, coverageEMV_U, coverageEMV_B, lengthEMV;
  decl vMBR, coverageMBR_L, coverageMBR_U, coverageMBR_B, lengthMBR;
  decl vQBR, coverageQBR_L, coverageQBR_U, coverageQBR_B, lengthQBR;
  decl vsort;

  decl lprob = columns(prob);
  decl coverage = zeros(columns(prob), 3*lprob);

  /* Start the clock */
  dExecTime = timer();

  /* Initializing the variables */
  fail = 0;            
  failMLE = 0;      
  failstderrors = 0;  
  failQBR_L = 0;
  failQBR_U = 0;
  failQBR_phi = 0;
  
  vEMV = zeros(iRep, lprob);
  coverageEMV_L = zeros(iRep,idiv(lprob,2)+1);
  coverageEMV_U = zeros(iRep,idiv(lprob,2)+1);
  coverageEMV_B = zeros(iRep,idiv(lprob,2)+1);
  lengthEMV = zeros(iRep,idiv(lprob,2)+1);

  vMBR = zeros(iRep, lprob);
  coverageMBR_L = zeros(iRep,idiv(lprob,2)+1);
  coverageMBR_U = zeros(iRep,idiv(lprob,2)+1);
  coverageMBR_B = zeros(iRep,idiv(lprob,2)+1);
  lengthMBR = zeros(iRep,idiv(lprob,2)+1);

  vQBR = zeros(iRep, lprob);
  coverageQBR_L = zeros(iRep,idiv(lprob,2)+1);
  coverageQBR_U = zeros(iRep,idiv(lprob,2)+1);
  coverageQBR_B = zeros(iRep,idiv(lprob,2)+1);
  lengthQBR = zeros(iRep,idiv(lprob,2)+1);

  /* Defining the random numbers generator */
  ranseed("GM");
  /* Defining the seed */
  ranseed({1965, 2001});

  /* Monte Carlo loop */
  for(ci = 0; ci < iRep + fail; ci++)
  {
    Y = zeros(nobs, 1);
	decl r=phi;
	decl a=phi/mu;
	Y = rangamma(nobs,1,r,a); /* Gamma Distribution Gamma(r,a) */

  /* ############################################################################################################## */
  /* ############################################################################################################## */
  /* ############################################################################################################## */
  /* ############################################################################################################## */
  /* Wald-type CIs, i.e. CIs that use the asymptotic normality of MLE (usual CIs) */
  
  /* Obtaining the moment estimates sugested in Ferrari e Cribari-Neto (2004)*/
  iota = ones(nobs,1);
  muinicial = meanc(Y);
  phiinicial =  muinicial^2/varc(Y);
  vthetahat = muinicial|phiinicial;
  Hess=unit(rows(vthetahat));
  converge = MaxSQP(floglik, &vthetahat, &dfunchat, &Hess, 0, 0, 0, <0;0>, <>);   /* analytics derivatives */
  /* Checking the convergence */
  if(converge != MAX_CONV && converge != MAX_WEAK_CONV)
  {  
    failMLE++;
    print("\n mle do not converge");
    fail++;
    continue;       /* quit the MC loop */
  }

  muhat=vthetahat[0];
  phihat=vthetahat[1];

  Imumu = nobs*phihat/muhat^2;
  Imuphi = 0;
  Iphiphi = nobs*(-1/phihat+polygamma(phihat,1));

  fisher = ( (Imumu~Imuphi)|(Imuphi~Iphiphi) );
  fisherinv = invert(fisher);
  stderrors = sqrt(diagonal(fisherinv));

  if(isnan(stderrors))
  {
    failstderrors++;
    print("\n mle isnan(stderrors)");
    fail++;
    continue;       /* quit the MC loop */
  }

  vEMV[ci-fail][] = phihat - up*stderrors[1];

  /* ############################################################################################################## */
  /* ############################################################################################################## */
  /* ############################################################################################################## */
  /* ############################################################################################################## */
  /* Wald-type CIs with median bias reduced	estimates of all the parameters in place of MLEs (adjusted CIs) */
  /* Kenne Pagui, E. C., Salvan, A. & Sartori, N.(2017). Median bias reduction of maximum likelihood estimates.*/
  /* Biometrika, 104, 923-938.*/

  upi=0;
  MaxControl( 50, 0 );									  
  vthetaQBR = muhat|phihat;
  convergeQBR=SolveNLE(QBREmu, &vthetaQBR);
  /* Checking the convergence */
  if(convergeQBR != MAX_CONV && convergeQBR != MAX_WEAK_CONV)
  {
    print("\n MBR mu do not converge");
    fail++;
    continue;       /* quit the MC loop */
  }
  muhatMBR=vthetaQBR[0];
  if(vthetaQBR[1]<0)
  {
    failQBR_phi++;
    print("\n phihatQBR < 0");
    fail++;
    continue;       /* quit the MC loop */
  }

  vthetaQBR = muhat|phihat;
  convergeQBR=SolveNLE(QBREphi, &vthetaQBR);
  /* Checking the convergence */
  if(convergeQBR != MAX_CONV && convergeQBR != MAX_WEAK_CONV)
  {
    print("\n MBR phi do not converge");
	fail++;
    continue;       /* quit the MC loop */
  }
  phihatMBR=vthetaQBR[1];
  if(vthetaQBR[1]<0)
  {
    failQBR_phi++;
    print("\n phihatQBR < 0");
    fail++;
    continue;       /* quit the MC loop */
  }

  Imumu = nobs*phihatMBR/muhatMBR^2;
  Imuphi = 0;
  Iphiphi = nobs*(-1/phihatMBR+polygamma(phihatMBR,1));

  fisher = ( (Imumu~Imuphi)|(Imuphi~Iphiphi) );
  fisherinv = invert(fisher);
  stderrorsMBR = sqrt(diagonal(fisherinv));
  
  vMBR[ci-fail][] = phihatMBR - up*stderrorsMBR[1];

  /* ############################################################################################################## */
  /* ############################################################################################################## */
  /* ############################################################################################################## */
  /* ############################################################################################################## */
  /* Confidence intervals based on the alpha-quantile modified score (modified CIs) */
  /* Quantile bias reduction - QBR */
  for(i = 0; i < lprob ; i++)
  {
    upi=up[i];
	MaxControl( 50, 0 );
	vthetaQBR = muhat|phihat;
	convergeQBR=SolveNLE(QBREphi, &vthetaQBR);
	/* Checking convergence */
	if(convergeQBR != MAX_CONV && convergeQBR != MAX_WEAK_CONV)
	{
	  if(prob[i]<0.5){failQBR_U++;print("\n QBR_U do not converge");}
	  else {failQBR_L++;print("\n QBR_L do not converge");}
	  fail++;
	  continue;       /* quit the MC loop */
	}

	if(vthetaQBR[1]<0)
	{
	  failQBR_phi++;
	  print("\n phihatQBR_U<0");
	  fail++;
	  continue;       /* quit the MC loop */
	}

	vQBR[ci-fail][i] = vthetaQBR[1];
  }

  /* Coverage  */
  for(i = 0; i < idiv(lprob,2)+1 ; i++)	/* prob < 0.5 */
  {
    coverageEMV_L[ci-fail][i] = (vEMV[ci-fail][lprob-i-1] .<= phi);
	coverageMBR_L[ci-fail][i] = (vMBR[ci-fail][lprob-i-1] .<= phi);
	coverageQBR_L[ci-fail][i] = (vQBR[ci-fail][lprob-i-1] .<= phi);

	coverageEMV_U[ci-fail][i] = (phi .<= vEMV[ci-fail][i]);
	coverageMBR_U[ci-fail][i] = (phi .<= vMBR[ci-fail][i]);
	coverageQBR_U[ci-fail][i] = (phi .<= vQBR[ci-fail][i]);

	coverageEMV_B[ci-fail][i] = (vEMV[ci-fail][lprob-i-1] .<= phi)&&(phi .<= vEMV[ci-fail][i]);
	coverageMBR_B[ci-fail][i] = (vMBR[ci-fail][lprob-i-1] .<= phi)&&(phi .<= vMBR[ci-fail][i]);
	coverageQBR_B[ci-fail][i] = (vQBR[ci-fail][lprob-i-1] .<= phi)&&(phi .<= vQBR[ci-fail][i]);

	lengthEMV[ci-fail][i] =	vEMV[ci-fail][i] - vEMV[ci-fail][lprob-i-1];
	lengthMBR[ci-fail][i] =	vMBR[ci-fail][i] - vMBR[ci-fail][lprob-i-1];
	lengthQBR[ci-fail][i] =	vQBR[ci-fail][i] - vQBR[ci-fail][lprob-i-1];
  }

  if(ci-1000*idiv(ci,1000)==0)print("\n ci : ",ci);

} /* Monte Carlo loop */

coverage =  prob[0:idiv(lprob,2)]' ~ ((sumc(coverageEMV_L)/iRep)'*100) ~ ((sumc(coverageMBR_L)/iRep)'*100) ~ ((sumc(coverageQBR_L)/iRep)'*100)
                                   ~ ((sumc(coverageEMV_U)/iRep)'*100) ~ ((sumc(coverageMBR_U)/iRep)'*100) ~ ((sumc(coverageQBR_U)/iRep)'*100)
								   ~ ((sumc(coverageEMV_B)/iRep)'*100) ~ ((sumc(coverageMBR_B)/iRep)'*100) ~ ((sumc(coverageQBR_B)/iRep)'*100)
								   ~ (sumc(lengthEMV)./sumc(coverageEMV_B))' ~ (sumc(lengthMBR)./sumc(coverageMBR_B))' ~ (sumc(lengthQBR)./sumc(coverageQBR_B))';
decl pos=vecindex(isdotnan(coverage[idiv(lprob,2)][]),1)';
coverage[idiv(lprob,2)][pos]= 0;

/* Printing results */
print( "\n OX PROGRAM: ", oxfilename(0) );
print( "\n OX VERSION: ", oxversion() );
print( "\n DISTRIBUTION: Gamma(",mu,",",phi,") " );	
print( "\n SAMPLE SIZE: ", nobs );
print( "\n MC REPLICATES         : ", iRep );
print( "\n MC fails              : ", ci-iRep );
print( "\n\t MLE fails         : ", failMLE );
print( "\n\t stderrors fails   : ", failstderrors );
print( "\n\t QBR_L fails       : ", failQBR_L );
print( "\n\t QBR_U fails       : ", failQBR_U );
print( "\n\t QBR_sigma<0       : ", failQBR_phi );

print("\n\n------------------------------ COVERAGE FOR THE PARAMETER MU ------------------------------");
print("\n\n          CL        EMV_L        MBR_L        QBR_L        EMV_U        MBR_U        QBR_U");
print("\n", (1-coverage[][0])*100 ~ coverage[][1:6]);

print("\n\n          CL        EMV_B        MBR_B        QBR_B   EMV_length   MBR_length   QBR_length");
print("\n",  (2*(1-coverage[][0])-1)*100 ~ coverage[][7:9] ~ coverage[][10:]);

/* Saving coverage */
savemat("coverage.txt", coverage);
	   
/* Date, Time and Execution time */
print( "\nDATE: ", date() );
print( "\nTIME: ", time(), "\n" );
print( "\n\n\t\t EXECUTION TIME: ", timespan(dExecTime) );
print( "\n" );

}/* end of main */

