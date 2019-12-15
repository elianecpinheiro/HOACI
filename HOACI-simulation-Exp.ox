/**************************************************************************************************
 * Monte Carlo simulation to evaluate the finite-sample performance of confidence intervals based on
 * the alpha-quantile modified score (modified CIs).
 * For comparison, we included results for Wald-type CIs, i.e. CIs that use the asymptotic normality of
 * MLE (usual CIs), and those that are constructed similarly to Wald-type CIs with median bias reduced
 * estimates of all the parameters in place of MLEs (adjusted CIs).
 *
 *  Exp(lambda) iid samples with density function
 *	f(x)= lambda exp(-lambda x)
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
const decl nobs = 5;      /* Sample size */
const decl iRep = 100000; /* Number of replicates for MC */

const decl lambda = 1;	  /* true value of lambda */ 

const decl prob = <0.005, 0.010, 0.020, 0.025, 0.030, 0.040, 0.050, 0.060, 0.070, 0.080, 0.090, 0.100, 0.150, 0.200, 0.250, 0.300, 0.400,
                    0.500,
					0.600, 0.700, 0.750, 0.800, 0.850, 0.900, 0.910, 0.920, 0.930, 0.940, 0.950, 0.960, 0.970, 0.975, 0.980, 0.990, 0.995>; 

/*  Starting the program */
main()
{

  /* Declaring variables */

  up=quann(prob);
 
  /* System variables */
  decl dExecTime, fail, ci, i;
  
  decl failMLE, failstderrors, failQBR_L, failQBR_U, failQBR_phi;
  
  /* Model variables */
  decl lambdahat;
  decl Ilambdalambda, stderrors;

  decl lambdatil;
  decl Ilambdalambdatil, stderrorstil;

  decl lambdahatQBR;
  decl convergeQBR;

  /* Result variables */
  decl vEMV, coverageEMV_L, coverageEMV_U, coverageEMV_B, lengthEMV;
  decl vMBR, coverageMBR_L, coverageMBR_U, coverageMBR_B, lengthMBR;
  decl vQBR, coverageQBR_L, coverageQBR_U, coverageQBR_B, lengthQBR;
//  decl vsort;

  decl lprob = columns(prob);
  decl coverage = zeros(columns(prob), 3*lprob);

  decl C; 

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
  
//  vsort =  zeros(iRep,2);
  
  /* Defining the random numbers generator */
  ranseed("PM");
  /* Defining the seed */
  ranseed(5);

  /* Monte Carlo loop */
  for(ci = 0; ci < iRep + fail; ci++)
  {
    Y = zeros(nobs, 1);
	Y = -log( 1 - ranu(nobs,1) )/lambda;	/* exponential(lambda) */
 
  /* ############################################################################################################## */
  /* ############################################################################################################## */
  /* ############################################################################################################## */
  /* ############################################################################################################## */
  /* Wald-type CIs, i.e. CIs that use the asymptotic normality of MLE (usual CIs) */
  
  lambdahat=1/meanc(Y);
  Ilambdalambda = nobs/lambdahat^2;
  stderrors = sqrt(1/Ilambdalambda);
  if(isnan(stderrors))
  {
    failstderrors++;
    fail++;
    continue;       /* quit the MC loop */
  }

  vEMV[ci-fail][] = lambdahat - up*stderrors;

  /* ############################################################################################################## */
  /* ############################################################################################################## */
  /* ############################################################################################################## */
  /* ############################################################################################################## */
  /* Wald-type CIs with median bias reduced	estimates of all the parameters in place of MLEs (adjusted CIs) */
  /* Kenne Pagui, E. C., Salvan, A. & Sartori, N.(2017). Median bias reduction of maximum likelihood estimates.*/
  /* Biometrika, 104, 923-938.*/

  upi=0;
  C = upi - (upi^2-1)/3/sqrt(nobs) + (upi^3-7*upi)/36/nobs;
  lambdatil =  lambdahat*(1-C/sqrt(nobs));
  Ilambdalambdatil = nobs/lambdatil^2;
  stderrorstil = sqrt(1/Ilambdalambdatil);

  vMBR[ci-fail][] = lambdatil - up*stderrorstil;
  
  /* ############################################################################################################## */
  /* ############################################################################################################## */
  /* ############################################################################################################## */
  /* ############################################################################################################## */
  /* Confidence intervals based on the alpha-quantile modified score (modified CIs) */
  /* Quantile bias reduction - QBR */
  for(i = 0; i < lprob ; i++)
  {
    upi=up[i];
	C = upi - (upi^2-1)/3/sqrt(nobs) + (upi^3-7*upi)/36/nobs;
    vQBR[ci-fail][i] = lambdahat*(1-C/sqrt(nobs));
  }

  /* Coverage  */
  for(i = 0; i < idiv(lprob,2)+1 ; i++)	/* prob < 0.5 */
  {
    coverageEMV_L[ci-fail][i] = (vEMV[ci-fail][lprob-i-1] .<= lambda);
	coverageMBR_L[ci-fail][i] = (vMBR[ci-fail][lprob-i-1] .<= lambda);
	coverageQBR_L[ci-fail][i] = (vQBR[ci-fail][lprob-i-1] .<= lambda);

	coverageEMV_U[ci-fail][i] = (lambda .<= vEMV[ci-fail][i]);
	coverageMBR_U[ci-fail][i] = (lambda .<= vMBR[ci-fail][i]);
	coverageQBR_U[ci-fail][i] = (lambda .<= vQBR[ci-fail][i]);

	coverageEMV_B[ci-fail][i] = (vEMV[ci-fail][lprob-i-1] .<= lambda)&&(lambda .<= vEMV[ci-fail][i]);
	coverageMBR_B[ci-fail][i] = (vMBR[ci-fail][lprob-i-1] .<= lambda)&&(lambda .<= vMBR[ci-fail][i]);
	coverageQBR_B[ci-fail][i] = (vQBR[ci-fail][lprob-i-1] .<= lambda)&&(lambda .<= vQBR[ci-fail][i]);

	lengthEMV[ci-fail][i] =	vEMV[ci-fail][i] - vEMV[ci-fail][lprob-i-1];
	lengthMBR[ci-fail][i] =	vMBR[ci-fail][i] - vMBR[ci-fail][lprob-i-1];
	lengthQBR[ci-fail][i] =	vQBR[ci-fail][i] - vQBR[ci-fail][lprob-i-1];
  }

  if(ci-10000*idiv(ci,10000)==0)print("\n ci : ",ci);

} /* Monte Carlo loop */

coverage =  (prob[0:idiv(lprob,2)]') ~ (((sumc(coverageEMV_L)/iRep)') ~ ((sumc(coverageMBR_L)/iRep)') ~ ((sumc(coverageQBR_L)/iRep)') ~ ((sumc(coverageEMV_U)/iRep)') ~ ((sumc(coverageMBR_U)/iRep)') ~ ((sumc(coverageQBR_U)/iRep)') ~ ((sumc(coverageEMV_B)/iRep)') ~ ((sumc(coverageMBR_B)/iRep)') ~ ((sumc(coverageQBR_B)/iRep)'))*100 ~ (sumc(lengthEMV)./sumc(coverageEMV_B))' ~ (sumc(lengthMBR)./sumc(coverageMBR_B))' ~ (sumc(lengthQBR)./sumc(coverageQBR_B))';
decl pos=vecindex(isdotnan(coverage[idiv(lprob,2)][]),1)';
coverage[idiv(lprob,2)][pos]= 0;

/* Printing results */
print( "\n OX PROGRAM: ", oxfilename(0) );
print( "\n OX VERSION: ", oxversion() );
print( "\n DISTRIBUTION: Exp(",lambda,") " );	
print( "\n SAMPLE SIZE: ", nobs );
print( "\n MC REPLICATES         : ", iRep );
print( "\n MC fails              : ", ci-iRep );
print( "\n\t MLE fails         : ", failMLE );
print( "\n\t stderrors fails   : ", failstderrors );
print( "\n\t QBR_L fails       : ", failQBR_L );
print( "\n\t QBR_U fails       : ", failQBR_U );
print( "\n\t QBR_sigma<0       : ", failQBR_phi );

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
