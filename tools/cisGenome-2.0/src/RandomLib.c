/* ----------------------------------------------------------------------- */
/*  RandomLib.c : implementation of the random generator library           */
/*  Author : Ji HongKai ; Time: 2001.12                                    */
/* ----------------------------------------------------------------------- */

#include "stdlib.h"
#include "stdio.h"
#include "math.h"
#include "limits.h"
#include "time.h"

#include "RandomLib.h"
#include "MatrixLib.h"

/* ----------------------------------------------------------------------- */ 
/*                            Global Variables                             */
/* ----------------------------------------------------------------------- */ 
struct tagRand_U_Parameter RAND_U_PARA_1;
struct tagRand_U_Parameter RAND_U_PARA_2;
struct tagRand_U_Parameter RAND_U_PARA_3;
struct tagRand_U_Mixture_Parameter RAND_U_PARA_MIX;
double pi;
double rnorm_status = 0;
double rnorm_registor;

/* ----------------------------------------------------------------------- */
/*                             Implementation                              */
/* ----------------------------------------------------------------------- */

/* ----------------------------------------------------------------------- */
/*                     Init Uniform Rand Generator                         */
/* ----------------------------------------------------------------------- */
void rand_u_init(long ninitstep)
{
	/*rand_u_1_init(ninitstep);*/
	rand_u_mixture_init(rand_u_1, rand_u_1_init, ninitstep, rand_u_2, rand_u_2_init, ninitstep, 12);
	pi = 4.0*atan(1.0);
}

/* ----------------------------------------------------------------------- */
/*            Create Random Number By Uniform Rand Generator               */
/* ----------------------------------------------------------------------- */
double rand_u()
{
	/*return rand_u_interval(rand_u_1, 1);*/
	/*return rand_u_1();*/
	return rand_u_mixture();
}

/* ----------------------------------------------------------------------- */
/*                     Init Uniform Rand Generator 1                       */
/* ----------------------------------------------------------------------- */
void rand_u_1_init(long ninitstep)
{
	long nj;
	long nstep;

	nstep = ninitstep+2;

	RAND_U_PARA_1.dA = 179.0;
	RAND_U_PARA_1.dM = pow(2.0,35.0);
	RAND_U_PARA_1.dC = 0.0;
	RAND_U_PARA_1.dX = 11.0;

	RAND_U_PARA_1.dX = (RAND_U_PARA_1.dA*RAND_U_PARA_1.dX+RAND_U_PARA_1.dC)/RAND_U_PARA_1.dM;
	RAND_U_PARA_1.dX = RAND_U_PARA_1.dX-(int)RAND_U_PARA_1.dX;
	RAND_U_PARA_1.dC /= RAND_U_PARA_1.dM;

	for(nj=0; nj<nstep; nj++)
	{
		RAND_U_PARA_1.dX = RAND_U_PARA_1.dA*RAND_U_PARA_1.dX+RAND_U_PARA_1.dC;
		RAND_U_PARA_1.dX = RAND_U_PARA_1.dX-(int)RAND_U_PARA_1.dX;
	}
}

/* ----------------------------------------------------------------------- */
/*            Create Random Number By Uniform Rand Generator 1             */
/* ----------------------------------------------------------------------- */
double rand_u_1()
{
	RAND_U_PARA_1.dX = RAND_U_PARA_1.dA*RAND_U_PARA_1.dX+RAND_U_PARA_1.dC;
	RAND_U_PARA_1.dX = RAND_U_PARA_1.dX-(int)RAND_U_PARA_1.dX;

	return RAND_U_PARA_1.dX;
}

/* ----------------------------------------------------------------------- */
/*                     Init Uniform Rand Generator 2                       */
/* ----------------------------------------------------------------------- */
void rand_u_2_init(long ninitstep)
{
	long nj;
	long nstep;

	nstep = ninitstep+2;

	RAND_U_PARA_2.dA = pow(2.0,7.0)+1.0;
	RAND_U_PARA_2.dM = pow(2.0,35.0);
	RAND_U_PARA_2.dC = pow(2.0,35.0)*(0.5+sqrt(3.0)/6.0);
	RAND_U_PARA_2.dX = 17.0;

	RAND_U_PARA_2.dX = (RAND_U_PARA_2.dA*RAND_U_PARA_2.dX+RAND_U_PARA_2.dC)/RAND_U_PARA_2.dM;
	RAND_U_PARA_2.dX = RAND_U_PARA_2.dX-(int)RAND_U_PARA_2.dX;
	RAND_U_PARA_2.dC /= RAND_U_PARA_2.dM;

	for(nj=0; nj<nstep; nj++)
	{
		RAND_U_PARA_2.dX = RAND_U_PARA_2.dA*RAND_U_PARA_2.dX+RAND_U_PARA_2.dC;
		RAND_U_PARA_2.dX = RAND_U_PARA_2.dX-(int)RAND_U_PARA_2.dX;
	}
}

/* ----------------------------------------------------------------------- */
/*            Create Random Number By Uniform Rand Generator 2             */
/* ----------------------------------------------------------------------- */
double rand_u_2()
{
	RAND_U_PARA_2.dX = RAND_U_PARA_2.dA*RAND_U_PARA_2.dX+RAND_U_PARA_2.dC;
	RAND_U_PARA_2.dX = RAND_U_PARA_2.dX-(int)RAND_U_PARA_2.dX;

	return RAND_U_PARA_2.dX;
}

/* ----------------------------------------------------------------------- */
/*                     Init Uniform Rand Generator 3                       */
/* ----------------------------------------------------------------------- */
void rand_u_3_init(long ninitstep)
{
	long nj;
	long nstep;

	nstep = ninitstep+2;

	RAND_U_PARA_3.dA = pow(5.0,5.0);
	RAND_U_PARA_3.dM = pow(2.0,35.0)-31.0;
	RAND_U_PARA_3.dC = 0.0;
	RAND_U_PARA_3.dX = pow(2.0,31.0)-1.0;

	RAND_U_PARA_3.dX = (RAND_U_PARA_3.dA*RAND_U_PARA_3.dX+RAND_U_PARA_3.dC)/RAND_U_PARA_3.dM;
	RAND_U_PARA_3.dX = RAND_U_PARA_3.dX-(int)RAND_U_PARA_3.dX;
	RAND_U_PARA_3.dC /= RAND_U_PARA_3.dM;

	for(nj=0; nj<nstep; nj++)
	{
		RAND_U_PARA_3.dX = RAND_U_PARA_3.dA*RAND_U_PARA_3.dX+RAND_U_PARA_3.dC;
		RAND_U_PARA_3.dX = RAND_U_PARA_3.dX-(int)RAND_U_PARA_3.dX;
	}
}

/* ----------------------------------------------------------------------- */
/*            Create Random Number By Uniform Rand Generator 3             */
/* ----------------------------------------------------------------------- */
double rand_u_3()
{
	RAND_U_PARA_3.dX = RAND_U_PARA_3.dA*RAND_U_PARA_3.dX+RAND_U_PARA_3.dC;
	RAND_U_PARA_3.dX = RAND_U_PARA_3.dX-(int)RAND_U_PARA_3.dX;

	return RAND_U_PARA_3.dX;
}

/* ----------------------------------------------------------------------- */
/*                     Init Uniform Rand Generator (Mixture)               */
/* ----------------------------------------------------------------------- */
void rand_u_mixture_init(double (* randfunc1)(), void (* initfunc1)(long), long ninitstep1, double (* randfunc2)(), void (* initfunc2)(long), long ninitstep2, int npoolsize)
{
	/* local variables */
	int nj;

	/* check parameter */
	if((npoolsize > RAND_MIXTURE_POOL_MAX) || (npoolsize < 2))
	{
		printf("Can't create mixture rand generator!\n");
		return;
	}

	/* init */
	RAND_U_PARA_MIX.nPoolSize = npoolsize;
	RAND_U_PARA_MIX.randfunc1 = randfunc1;
	RAND_U_PARA_MIX.randfunc2 = randfunc2;
	(*initfunc1)(ninitstep1);
	(*initfunc2)(ninitstep2);

	for(nj=0; nj<npoolsize; nj++)
		(RAND_U_PARA_MIX.dPool)[nj] = (*randfunc1)();
}

/* ----------------------------------------------------------------------- */
/*          Create Random Number By Uniform Rand Generator (Mixture)       */
/* ----------------------------------------------------------------------- */
double rand_u_mixture()
{
	/* local variables */
	int nIndex;
	double d_rand;

	/* create rand */
	d_rand = (*(RAND_U_PARA_MIX.randfunc2))();
	nIndex = (int)(d_rand*RAND_U_PARA_MIX.nPoolSize);
	if(nIndex == RAND_U_PARA_MIX.nPoolSize)
		nIndex--;

	d_rand = (RAND_U_PARA_MIX.dPool)[nIndex];
	(RAND_U_PARA_MIX.dPool)[nIndex] = (*(RAND_U_PARA_MIX.randfunc1))();

	/* return */
	return d_rand;
}

/* ----------------------------------------------------------------------- */
/*          Create Random Number By Uniform Rand Generator (Interval)      */
/* ----------------------------------------------------------------------- */
double rand_u_interval(double (* randfunc)(), long nInterval)
{
	long ni;
	double d_rand;

	/* check parameter */
	if(nInterval <= 0)
	{
		printf("Can't create rand number by interval generator!\n");
		return 0.0;
	}

	/* create rand */
	for(ni=0; ni<nInterval; ni++)
	{
		d_rand = (*randfunc)();
	}

	/* return */
	return d_rand;
}

/* ----------------------------------------------------------------------- */
/*                  Perform Uniform Test For Rand_UNIFORM                  */
/* nK is region number, nN is sample number, randfunc is rand generator.   */
/* return Chi statistics = sum{[n(j)-N/K]^2/(N/K)} (j=1-K)                 */ 
/* ----------------------------------------------------------------------- */
double rand_u_uniform_test(void (* initfunc)(long), long ninitstep, double (* randfunc)(), long nK, long nN)
{
	/* pCountMatrix is used for counting random number. */
	struct DOUBLEMATRIX *pCountMatrix;
	/* dChi is the statistics to be returned */
	double dChi;
	/* nj is used for counting */
	long nj;
	/* d_rand is used for recording rand number */
	double d_rand;
	/* nrandindex is the index of the region corresponded to rand number */
	long nrandindex;
	/* dCount is used for processing count */
	double dCount;
	double *pElement;

	/* check parameter */
	if((nK <= 0) || (nN <= 0))
	{
		printf("Can't perform uniform test!\n");
		return -1;
	}

	/* create matrix for counting */
	pCountMatrix = NULL;
	pCountMatrix = CreateDoubleMatrix(1, nK);
	if(pCountMatrix == NULL)
	{
		printf("Can't perform uniform test!\n");
		return -1.0;
	}

	/* init rand generator */
	(*initfunc)(ninitstep);

	/* counting */
	for(nj=0; nj<nN; nj++)
	{
		d_rand = (*randfunc)();
		nrandindex = (long)(d_rand*nK);
		if(nrandindex == nK)
			nrandindex--;
		dCount = DMGETAT(pCountMatrix, 0, nrandindex);
		dCount++;
		DMSETAT(pCountMatrix, 0, nrandindex, dCount);
	}

	/* calculate chi statistics */
	pElement = pCountMatrix->pMatElement;
	dCount = (double)nN/(double)nK;
	dChi = 0.0;
	for(nj=0; nj<nK; nj++)
	{
		dChi += (*pElement-dCount)*(*pElement-dCount)/dCount;
		pElement++;
	}

	/* return */
	return dChi;
}

/* ----------------------------------------------------------------------- */
/*                Perform Independence Test for Rand_Uniform               */
/* nJ is delay level, nN is sample number, randfunc is rand generator.     */
/* return Mu statistics = sqrt(N-J)*rho;                                   */
/* rho = sum[U(i)*U(i+J)-(1/2)*(1/2)]/[(N-J)*S^2];                         */
/* S^2 = sum[U(i)-1/2]^2/(N-1).                                            */ 
/* ----------------------------------------------------------------------- */
double rand_u_independence_test(void (* initfunc)(long), long ninitstep, double (* randfunc)(), long nJ, long nN)
{
	/* dMu is the statistics to be returned */
	double dMu;
	/* ni, nj is used for counting */
	long ni,nj;
	/* d_rand is used for recording rand number */
	double *d_rand;
	double *d_rand_j;
	double d_current;
	/* dS2 is used for calculating S^2. */
	double dS2;
	/* dMean is used for calculating mean. */
	double dMean;
	/* dRho is used for calculating rho. */
	double dRho;

	/* check parameter */
	if((nJ <= 0) || (nN <= nJ))
	{
		printf("Can't perform uniform test!\n");
		return -1;
	}

	/* init */
	dMu = 0.0;
	dS2 = 0.0;
	dMean = 0.0;
	dRho = 0.0;
	d_rand = (double *)calloc(nJ, sizeof(double));
	if(d_rand == NULL)
	{
		printf("Can't perform uniform test!\n");
		return -1;
	}

	/* init rand generator */
	(*initfunc)(ninitstep);
	
	/* init d_rand */
	d_rand_j = d_rand;
	for(ni=0; ni<nJ; ni++)
	{
		*d_rand_j = (*randfunc)();
		dMean += *d_rand_j;
		/*d_current = *d_rand_j-0.5;*/
		dS2 += (*d_rand_j)*(*d_rand_j);
		d_rand_j++;
	}

	/* calculate statistics */
	d_rand_j = d_rand;
	nj = 0;
	for(; ni<nN; ni++)
	{
		d_current = (*randfunc)();
		dRho += (*d_rand_j)*d_current;
		*d_rand_j = d_current;
		d_rand_j++;
		nj++;
		if(nj == nJ)
		{
			nj = 0;
			d_rand_j = d_rand;
		}
		dMean += d_current;
		/*d_current -= 0.5;*/
		dS2 += d_current*d_current;
	}

	dMean = dMean/(double)nN;
	dS2 = (dS2-nN*dMean*dMean)/(double)(nN-1);
	if(dS2 < ZERO_BOUND)
		return 1.0/ZERO_BOUND;

	dRho = dRho/(double)(nN-nJ);
	dRho = (dRho-dMean*dMean)/dS2;
	dMu = fabs(dRho*sqrt((double)(nN-nJ)));

	/* free memory */
	free(d_rand);

	/* return */
	return dMu;
}

/* ----------------------------------------------------------------------- */
/*          Generate exponential random number                             */
/* pdf(x) = exp(-x/mu)/mu                                                  */
/* ----------------------------------------------------------------------- */
double exprnd(double mu)
{
	double du;

	if(mu <= 0.0)
	{
		printf("Error: exprnd, mu must be greater than 0!\n");
		exit(EXIT_FAILURE);
	}

	du = rand_u();
	return -log(du)*mu;
}

/* ----------------------------------------------------------------------- */
/*          Generate gamma random number                                   */ 
/* pdf(x) = x^(alpha-1)*exp(-x/beta)/[gamma(alpha)*beta^alpha]             */
/* ----------------------------------------------------------------------- */
double gamrnd(double alpha, double beta)
{
	double dgamma,du,dv,W,X,Y,Z,b,c;
	/* double e,V,dr,a1,a2; */
	int flag;

	if(alpha<=0.0 || beta<=0.0)
	{
		printf("Error: gamrnd, alpha and beta must be greater than 0!\n");
		exit(EXIT_FAILURE);
	}

	if(alpha > 1.0)
	{
		/* We implement Best's Rejection algorithm(1978), Devroye(page 410) for alpha >= 1 */
	    b = alpha - 1.0; 
		c = 3.0*alpha - 0.75;
		flag = 1;
		while(flag == 1) 
		{
			du = rand_u();
			dv = rand_u();
			W = du*(1.0-du);
			Y = sqrt(c/W)*(du-0.5); 
			X = b+Y;
		    
			if(X >= 0)
			{
				Z = 64.0*pow(W,3.0)*pow(dv,2.0);
				if(Z <= (1.0-(2.0*Y*Y)/X))
				{
					flag = 0;
				}
				else 
				{
					if( log(Z) <= (2.0*( b*log(X/b)-Y )) )
						flag = 0;
				}
			}
		}
		dgamma = X;
	}
	else if(alpha < 1.0)
	{
		/* We implement Johnk's gamma generator, Devroye(page 418) for alpha < 1 */
	    flag = 1;
		while(flag == 1) 
		{
			du = rand_u(); 
			dv = rand_u();
			X = pow(du, 1.0/alpha); 
			Y = pow(dv, 1.0/(1.0-alpha));
			if((X + Y) <= 1.0)
				flag = 0;
		}
		Z = rand_u(); 
		Z = -log(Z);
		dgamma = (Z*X)/(X+Y);

		/* From Longe Book, Excercise. Lowest acceptance rate = 0.7196. */
		/* e = exp(1.0);
		a1 = e/(e+alpha);
		a2 = 1-a1;
		W = 1/alpha;
		V = alpha-1.0;
		do {
			if( rand_u() < a1 )
			{
				du = rand_u();
				dgamma = pow(du, W);
				dr = exp(-dgamma);
			}
			else
			{
				du = rand_u();
				dgamma = 1-log(du);
				dr = pow(dgamma, V);
			}
		} while( rand_u() > dr ); */
	}
	else if(alpha == 1.0)
	{
		/* Generate exponential random number */
		du = rand_u();
		dgamma = -log(du);
	}
	
	/* return */
	dgamma = dgamma*beta;
	return dgamma;
}

/* ----------------------------------------------------------------------- */
/*          Generate dirichlet random number                               */ 
/* pdf(x) = Product[xi^(alphai-1)]/normalizing constant                    */
/* ----------------------------------------------------------------------- */
double *dirichletrnd(double alpha[], int n, double output[])
{
	int ni;
	double dSum;

	dSum = 0.0;
	for(ni=0; ni<n; ni++)
	{
		output[ni] = gamrnd(alpha[ni], 1.0);
		dSum += output[ni];
	}
	for(ni=0; ni<n; ni++)
	{
		output[ni] /= dSum;
	}

	return output;
}


/* ----------------------------------------------------------------------- */
/*          Generate normal random number                                  */
/* pdf(x) = exp{-(x-mu)^2/(2*sigma^2)}/[sqrt(2*pi)*sigma]                  */
/* ----------------------------------------------------------------------- */
double normrnd(double mu, double sigma)
{
	double du1,du2;
	double dnorm;

	if(sigma < 0.0)
	{
		printf("Error: normrnd, sigma must be no less than 0!\n");
		exit(EXIT_FAILURE);
	}

	du1 = rand_u();
	du2 = rand_u();
	dnorm = sqrt(-2.0*log(du1))*cos(2.0*pi*du2);
	dnorm = mu+sigma*dnorm;

	return dnorm;
}

/* ----------------------------------------------------------------------- */
/*          Create normal random number                                    */
/* pdf(x) = exp{-(x-mu)^2/(2*sigma^2)}/[sqrt(2*pi)*sigma]                  */
/* Note: normrnd is preferred to rnorm, but normrnd has lower efficiency.  */
/* ----------------------------------------------------------------------- */
double rnorm(double mu, double sigma)
{
	double v1,v2,rsq,fac;

	if(sigma < 0.0)
	{
		printf("Error: rnorm, sigma must be no less than 0!\n");
		exit(EXIT_FAILURE);
	}

	if(rnorm_status == 0)
	{
		/* generate two normal random numbers */
		do {
			v1 = 2.0*rand_u()-1.0;
			v2 = 2.0*rand_u()-1.0;
			rsq = v1*v1+v2*v2; 
		} while( (rsq>=1.0) || (rsq==0.0) );

		fac = sqrt(-2.0*log(rsq)/rsq);
		rnorm_registor = v1*fac;
		rnorm_status = 1;
		return v2*fac*sigma+mu;
	}
	else
	{
		rnorm_status = 0;
		return rnorm_registor*sigma+mu;
	}
}

/* ----------------------------------------------------------------------- */
/*          Get normal density                                             */
/* ----------------------------------------------------------------------- */
double dnorm(double dx, double mu, double sigma)
{
	double df;
	
	if(sigma <= 0.0)
	{
		printf("rnorm error: sigma=0!\n");
		exit(EXIT_FAILURE);
	}

	df = exp(-pow((dx-mu), 2.0)/(2.0*sigma*sigma))/(sqrt(2.0*pi)*sigma);

	return df;
}






/* ----------------------------------------------------------------------- */
/*                Temporary Functions                                      */
/* ----------------------------------------------------------------------- */

/* ----------------------------------------------------------------------- */
/*                Temporary functions for creating Gamma                   */
/* ----------------------------------------------------------------------- */
double rand_gamma(double dalpha, /* the shape parameter */ double dbeta  /* the scale parameter */)
{
	int flag;
	double b, c, W, Y, X, U, V, Z, E;
	double dgamma;

	/* Longe: a low efficiency algorithm for alpha >1.0 */
	/*if(alpha > 1.0)
	{
		e = exp(1.0);
		a1 = alpha-1.0;
		a2 = alpha+1.0;
		ku = pow((a1/e), (a1/2.0));
		kv = pow((a2/e), (a2*2.0));
		do {
		du = rand_u();
		dv = rand_u();
		dgamma = kv*dv/(ku*du);
		W = dgamma/a1;
		} while( (2.0*log(du)/a1-1-log(W)+W) > 0.0);
	}
	/* From Longe Book, Excercise. Lowest acceptance rate = 0.7196. */
	/*else if(alpha < 1.0)
	{
		e = exp(1.0);
		a1 = e/(e+alpha);
		a2 = 1-a1;
		W = 1/alpha;
		V = alpha-1.0;
		do {
			if( rand_u() < a1 )
			{
				du = rand_u();
				dgamma = pow(du, W);
				dr = exp(-dgamma);
			}
			else
			{
				du = rand_u();
				dgamma = 1-log(du);
				dr = pow(dgamma, V);
			}
		} while( rand_u() > dr );
	}*/

	if(dalpha > 1.0)
	{
		/* We implement Best's Rejection algorithm(1978), Devroye(page 410) for alpha >= 1 */
	    b = dalpha - 1.0; c = 3.0 * dalpha - 0.75;
		flag = 1;
		while(flag == 1) 
		{
			U = rand_u();
			V = rand_u();
			W = U * (1 - U);
			Y = sqrt(c / W) * (U - 0.5); 
			X = b + Y;
		    
			if(X >= 0)
			{
				Z = 64.0 * pow(W, 3) * pow(V, 2);
				if(Z <= (1.0 - (2.0 * Y * Y) / X))
				{
					flag = 0;
				}
				else 
				{
					if(log(Z) <= (2.0 * (b * log(X / b) - Y)))
						flag = 0;
				}
			}
		}
		dgamma = X / dbeta;
	}
	else if (dalpha < 1.0)
	{
		/* We implement Johnk's gamma generator, Devroye(page 418) for alpha < 1 */
	    flag = 1;
		while(flag == 1) 
		{
			U = rand_u(); 
			V = rand_u();
			X = pow(U, 1.0 / dalpha); 
			Y = pow(V, 1.0 / (1.0 - dalpha));
			if((X + Y) <= 1.0)
				flag = 0;
		}
		E = rand_u(); 
		E = -log(E);
		dgamma = (E * X) / ((X + Y) * dbeta);
	}
	else if (dalpha == 1.0)
	{
		U = rand_u();
		dgamma = -log(U) / dbeta;
	}
	else
	{
		printf("Error: rand_gamma() parameter error!\n");
	}

	/* return */
	return dgamma;
}

#define TWO_PI 6.283185307179586
#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

double runif(long *idum, double a, double b)
{
        static int inext,inextp;
        static long ma[56];
        static int iff=0;
        long mj,mk;
        int i,ii,k;

        if (*idum < 0 || iff == 0) {
                iff=1;
                mj=labs(MSEED-labs(*idum));
                mj %= MBIG;
                ma[55]=mj;
                mk=1;
                for (i=1;i<=54;i++) {
                        ii=(21*i) % 55;
                        ma[ii]=mk;
                        mk=mj-mk;
                        if (mk < MZ) mk += MBIG;
                        mj=ma[ii];
                }
                for (k=1;k<=4;k++)
                        for (i=1;i<=55;i++) {
                                ma[i] -= ma[1+(i+30) % 55];
                                if (ma[i] < MZ) ma[i] += MBIG;
                        }
                inext=0;
                inextp=31;
                *idum=1;
        }
        if (++inext == 56) inext=1;
        if (++inextp == 56) inextp=1;
        mj=ma[inext]-ma[inextp];
        if (mj < MZ) mj += MBIG;
        ma[inext]=mj;
        return a + (b - a) * mj * FAC;
}

#undef MBIG
#undef MSEED
#undef MZ
#undef FAC


void rgamma(long *seed,   /* the random number generator seed */
            double alpha, /* the shape parameter */
            double beta,  /* the scale parameter */
            int n,        /* the sample size */
            double *x)
{

  int i, flag;
  double b, c, W, Y, X, U, V, Z, E;

  if(alpha >= 1) {
     /* We implement Best's Rejection algorithm(1978), Devroye(page
       410) for alpha >= 1 */

    b = alpha - 1.0; c = 3.0 * alpha - 0.75;
    for(i = 0; i < n; i++) {
      flag = 1;
      while(flag == 1) {
        U = runif(seed, 0.0, 1.0); V = runif(seed, 0.0, 1.0);
        W = U * (1 - U); Y = sqrt(c / W) * (U - 0.5); X = b + Y;
        if(X >= 0) {
          Z = 64.0 * pow(W, 3) * pow(V, 2);
          if(Z <= (1.0 - (2.0 * Y * Y) / X)) {
            flag = 0;
          }
          else {
            if(log(Z) <= (2.0 * (b * log(X / b) - Y)))
              flag = 0;
          }
        }
      }
      x[i] = X / beta;
    }
  }
  else {
    /* We implement Johnk's gamma generator, Devroye(page 418) for
       alpha < 1 */

    for(i = 0; i < n; i++) {
      flag = 1;
      while(flag == 1) {
        U = runif(seed, 0.0, 1.0); V = runif(seed, 0.0, 1.0);
        X = pow(U, 1.0 / alpha); Y = pow(V, 1.0 / (1.0 - alpha));
        if((X + Y) <= 1.0)
          flag = 0;
      }
      E = runif(seed, 0.0, 1.0); E = -log(E);
      x[i] = (E * X) / ((X + Y) * beta);
    }
  }

}
