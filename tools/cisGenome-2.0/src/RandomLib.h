/* ----------------------------------------------------------------------- */
/*  Random.h : interface of the random generator library                   */
/*  Author : Ji HongKai ; Time: 2001.12                                    */
/* ----------------------------------------------------------------------- */

/* ----------------------------------------------------------------------- */
/*                              Declaration                                */
/* ----------------------------------------------------------------------- */
#define RAND_MIXTURE_POOL_MAX 24

/* ----------------------------------------------------------------------- */ 
/*                                 Struct                                  */
/* ----------------------------------------------------------------------- */ 
/* tagRand_U_Parameter is used for creating uniformly distributed */
/* pseudo random number: X(k+1) = A*X(k)+C (mod M) */
struct tagRand_U_Parameter
{
	double dA;
	double dM;
	double dC;
	double dX;
};

/* tagRand_U_Mixture_Parameter is used for creating uniformly distributed */
/* pseudo random number by mixture rand generator. */
struct tagRand_U_Mixture_Parameter
{
	double dPool[RAND_MIXTURE_POOL_MAX];
	int nPoolSize;
	double (* randfunc1)();
	double (* randfunc2)();
};

/* ----------------------------------------------------------------------- */
/*                              Declaration                                */
/* ----------------------------------------------------------------------- */

/* ----------------------------------------------------------------------- */
/*                     Init Uniform Rand Generator                         */
/* ----------------------------------------------------------------------- */
void rand_u_init(long ninitstep);

/* ----------------------------------------------------------------------- */
/*            Create Random Number By Uniform Rand Generator               */
/* ----------------------------------------------------------------------- */
double rand_u();

/* ----------------------------------------------------------------------- */
/*                     Init Uniform Rand Generator 1                       */
/* ----------------------------------------------------------------------- */
void rand_u_1_init(long ninitstep);

/* ----------------------------------------------------------------------- */
/*            Create Random Number By Uniform Rand Generator 1             */
/* ----------------------------------------------------------------------- */
double rand_u_1();

/* ----------------------------------------------------------------------- */
/*                     Init Uniform Rand Generator 2                       */
/* ----------------------------------------------------------------------- */
void rand_u_2_init(long ninitstep);

/* ----------------------------------------------------------------------- */
/*            Create Random Number By Uniform Rand Generator 2             */
/* ----------------------------------------------------------------------- */
double rand_u_2();

/* ----------------------------------------------------------------------- */
/*                     Init Uniform Rand Generator 3                       */
/* ----------------------------------------------------------------------- */
void rand_u_3_init(long ninitstep);

/* ----------------------------------------------------------------------- */
/*            Create Random Number By Uniform Rand Generator 3             */
/* ----------------------------------------------------------------------- */
double rand_u_3();

/* ----------------------------------------------------------------------- */
/*                     Init Uniform Rand Generator (Mixture)               */
/* ----------------------------------------------------------------------- */
void rand_u_mixture_init(double (* randfunc1)(), void (* initfunc1)(long), long ninitstep1, double (* randfunc2)(), void (* initfunc2)(long), long ninitstep2, int npoolsize);

/* ----------------------------------------------------------------------- */
/*          Create Random Number By Uniform Rand Generator (Mixture)       */
/* ----------------------------------------------------------------------- */
double rand_u_mixture();

/* ----------------------------------------------------------------------- */
/*          Create Random Number By Uniform Rand Generator (Interval)      */
/* ----------------------------------------------------------------------- */
double rand_u_interval(double (* randfunc)(), long nInterval);

/* ----------------------------------------------------------------------- */
/*                  Perform Uniform Test For Rand_UNIFORM                  */
/* nK is region number, nN is sample number, randfunc is rand generator.   */
/* return Chi statistics = sum{[n(j)-N/K]^2/(N/K)} (j=1-K)                 */ 
/* ----------------------------------------------------------------------- */
double rand_u_uniform_test(void (* initfunc)(long), long ninitstep, double (* randfunc)(), long nK, long nN);

/* ----------------------------------------------------------------------- */
/*                Perform Independence Test for Rand_Uniform               */
/* nJ is delay level, nN is sample number, randfunc is rand generator.     */
/* return Mu statistics = sqrt(N-J)*rho;                                   */
/* rho = sum[U(i)*U(i+J)-(1/2)*(1/2)]/[(N-J)*S^2];                         */
/* S^2 = sum[U(i)-1/2]^2/(N-1).                                            */ 
/* ----------------------------------------------------------------------- */
double rand_u_independence_test(void (* initfunc)(long), long ninitstep, double (* randfunc)(), long nJ, long nN);

/* ----------------------------------------------------------------------- */
/*          Generate exponential random number                             */
/* pdf(x) = exp(-x/mu)/mu                                                  */
/* ----------------------------------------------------------------------- */
double exprnd(double mu);

/* ----------------------------------------------------------------------- */
/*          Generate gamma random number                                   */ 
/* pdf(x) = x^(alpha-1)*exp(-x/beta)/[gamma(alpha)*beta^alpha]             */
/* ----------------------------------------------------------------------- */
double gamrnd(double alpha, double beta);

/* ----------------------------------------------------------------------- */
/*          Generate dirichlet random number                               */ 
/* pdf(x) = Product[xi^(alphai-1)]/normalizing constant                    */
/* ----------------------------------------------------------------------- */
double *dirichletrnd(double alpha[], int n, double output[]);

/* ----------------------------------------------------------------------- */
/*          Generate normal random number                                  */
/* pdf(x) = exp{-(x-mu)^2/(2*sigma^2)}/[sqrt(2*pi)*sigma]                  */
/* ----------------------------------------------------------------------- */
double normrnd(double mu, double sigma);

/* ----------------------------------------------------------------------- */
/*          Create normal random number                                    */
/* pdf(x) = exp{-(x-mu)^2/(2*sigma^2)}/[sqrt(2*pi)*sigma]                  */
/* Note: normrnd is preferred to rnorm, but normrnd has lower efficiency.  */
/* ----------------------------------------------------------------------- */
double rnorm(double mu, double sigma);

/* ----------------------------------------------------------------------- */
/*          Get normal density                                             */
/* ----------------------------------------------------------------------- */
double dnorm(double dx, double mu, double sigma);


/* ----------------------------------------------------------------------- */
/*                Temporary Gamma                                          */ 
/* ----------------------------------------------------------------------- */
double rand_gamma(double dalpha, /* the shape parameter */ double dbeta  /* the scale parameter */);

double runif(long *idum, double a, double b);

void rgamma(long *seed,   /* the random number generator seed */
            double alpha, /* the shape parameter */
            double beta,  /* the scale parameter */
            int n,        /* the sample size */
            double *x);