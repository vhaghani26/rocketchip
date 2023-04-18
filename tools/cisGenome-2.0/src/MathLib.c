/* ----------------------------------------------------------------------- */
/*  MathLib.c : implementation of the mathematics library                  */
/*  Author : Ji HongKai ; Time: 2004.07                                    */
/* ----------------------------------------------------------------------- */

#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include "math.h"
#include "limits.h"
#include "time.h"

#include "MathLib.h"

#define ITMAX 100000
#define EPS 3.0e-7
#define FPMIN 1.0e-30

/* ----------------------------------------------------------------------- */
/*                            lnGamma Function                             */
/* Return ln(gamma(dx)).                                                   */
/* ----------------------------------------------------------------------- */
double gammaln(double dx)
{
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,
						-86.50532032941677,
						24.01409824083091,
						-1.231739572450155,
						0.1208650973866179e-2,
						-0.5395239384953e-5};
	int j;

	if(dx<0.0)
	{
		printf("Error: gammaln(dx), dx must be greater than 0!\n");
		exit(EXIT_FAILURE);
	}
	y = dx;
	x = dx;
	tmp = x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser = 1.000000000190015;
	for(j=0; j<=5; j++)
		ser += cof[j]/++y;

	return -tmp+log(2.5066282746310005*ser/x);
}

/* ----------------------------------------------------------------------- */
/*                              Gamma Function                             */
/* Return gamma(dx).                                                       */
/* ----------------------------------------------------------------------- */
double gamma(double dx)
{
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,
						-86.50532032941677,
						24.01409824083091,
						-1.231739572450155,
						0.1208650973866179e-2,
						-0.5395239384953e-5};
	int j;

	if(dx<0.0)
	{
		printf("Error: gammaln(dx), dx must be greater than 0!\n");
		exit(EXIT_FAILURE);
	}
	y = dx;
	x = dx;
	tmp = x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser = 1.000000000190015;
	for(j=0; j<=5; j++)
		ser += cof[j]/++y;

	return exp(-tmp+log(2.5066282746310005*ser/x));
}

/* ----------------------------------------------------------------------- */
/*                           Factorial Function                            */
/* Return nx! = 1*2*...*nx.                                                */
/* ----------------------------------------------------------------------- */
double factorial(int nx)
{
	static int ntop=4;
	static double a[33] = {1.0,1.0,2.0,6.0,24.0};
	int j;

	if(nx<0)
	{
		printf("Error: factorial(nx), nx must be nonnegative!\n");
		exit(EXIT_FAILURE);
	}
	if(nx>32)
		return gamma(nx+1.0);

	while(ntop<nx)
	{
		j=ntop++;
		a[ntop]=a[j]*ntop;
	}
	return a[nx];
}

/* ----------------------------------------------------------------------- */
/*                               Beta Function                             */
/* Return beta(dx).                                                        */
/* ----------------------------------------------------------------------- */
double beta(double da, double db)
{
	return exp(gammaln(da)+gammaln(db)-gammaln(da+db));
}

/* ----------------------------------------------------------------------- */
/*                            Normal CDF Function                          */
/* Return normcdf(dx).                                                     */
/* ----------------------------------------------------------------------- */
double normcdf(double mu, double sigma, double dx)
{
	double dy;

	dy = (dx-mu)/sigma;
	if(dy>=0.0)
	{
		if(dy>7.5)
			return (1.0+erf(7.5/sqrt(2.0)))/2.0;
		else
			return (1.0+erf(dy/sqrt(2.0)))/2.0;
	}
	else
	{
		if(dy<-7.5)
			return (1.0-erf(7.5/sqrt(2.0)))/2.0;
		else
			return (1.0-erf(-dy/sqrt(2.0)))/2.0;
	}
}

/* ----------------------------------------------------------------------- */
/*                               Erf Function                              */
/* Return erf(dx).                                                         */
/* ----------------------------------------------------------------------- */
/* double erf(double dx)
{
	double y1=0.0;
	double y2=0.0;
	double y3=0.0;
	double y4=0.0;
	double y5=0.0;
	double x1 = -sqrt(245.0+14.0*sqrt(70.0))/21.0;
	double x2 = -sqrt(245.0-14.0*sqrt(70.0))/21.0;
	double x3=0.0;
	double x4=-x2;
	double x5=-x1;
	double w1=(322.0-13.0*sqrt(70.0))/900.0;
	double w2=(322.0+13.0*sqrt(70.0))/900.0;
	double w3=128.0/225.0;
	double w4=w2;
	double w5=w1;
	int n=4800;
	int i=0;
	double s=0.0;
	double h=dx/n;
	/* double coef = 2.0/sqrt(3.14159265358979323846264338327950); */
	/* double coef = 2.0/sqrt(4.0*atan(1.0));
	
	for(i=0; i<=n; i++)
	{
		y1=h*x1/2.0+(h+2.0*i*h)/2.0;
		y2=h*x2/2.0+(h+2.0*i*h)/2.0;
		y3=h*x3/2.0+(h+2.0*i*h)/2.0;
		y4=h*x4/2.0+(h+2.0*i*h)/2.0;
		y5=h*x5/2.0+(h+2.0*i*h)/2.0;
		s=s+h*coef*(w1*exp(-pow(y1,2.0))+w2*exp(-pow(y2,2.0))+w3*exp(-pow(y3,2.0))+w4*exp(-pow(y4,2.0))+w5*exp(-pow(y5,2.0)))/2.0;
	}

	
	return s;
} */

/* ----------------------------------------------------------------------- */
/*               Erf Function (Numerical Recipe Version)                   */
/* Return erf(dx).                                                         */
/* ----------------------------------------------------------------------- */
double erf(double dx)
{
	if(dx < 0.0)
		return -gammp(0.5,dx*dx);
	else
		return gammp(0.5,dx*dx);
}

/* ----------------------------------------------------------------------- */
/*                      Incomplete gamma function                          */
/* Return the incomplete gamma function P(a,x)                             */
/* ----------------------------------------------------------------------- */
double gammp(double da, double dx)
{
	double gamser = 0.0;
	double gammcf = 0.0;
	double gln = 0.0;

	if( (dx<0.0) || (da<=0.0) )
	{
		printf("Error: invalid arguments in routine gammp\n");
		exit(EXIT_FAILURE);
	}
	if( dx < (da+1.0) )
	{
		gser(&gamser,da,dx,&gln);
		return gamser;
	}
	else
	{
		gcf(&gammcf,da,dx,&gln);
		return 1.0-gammcf;
	}
}

/* ----------------------------------------------------------------------- */
/*                       Incomplete Gamma Function                         */
/* Returns the incomplete gamma function P(a,x) evaluated by its series    */
/* representation as gamser. Also returns lnGamma(a) as gln.               */     
/* ----------------------------------------------------------------------- */
void gser(double *gamser, double da, double dx, double *gln)
{
	int n;
	double sum,del,ap;

	*gln = gammaln(da);
	if(dx<=0.0)
	{
		if(dx<0.0)
			printf("x less than 0 in routine gser");
		*gamser=0.0;
		return;
	}
	else
	{
		ap = da;
		del = 1.0/da;
		sum = del;
		for(n=1; n<=ITMAX; n++)
		{
			++ap;
			del *= dx/ap;
			sum += del;
			if( fabs(del) < fabs(sum)*EPS)
			{
				*gamser = sum*exp(-dx+da*log(dx)-(*gln));
				return;
			}
		}
		
		printf("gamser: a too large, ITMAX too small in routine gser");
		*gamser = sum*exp(-dx+da*log(dx)-(*gln));
		return;
	}
}

/* ----------------------------------------------------------------------- */
/*                       Incomplete Gamma Function                         */
/* Returns the incomplete gamma function Q(a,x) evaluated by its continued */
/* fraction representation as gammcf. Also returns lnGamma(a) as gln.      */     
/* ----------------------------------------------------------------------- */
void gcf(double *gammcf, double da, double dx, double *gln)
{
	int i;
	double an,b,c,d,del,h;

	*gln = gammaln(da);
	b = dx+1.0-da;
	c = 1.0/FPMIN;
	d = 1.0/b;
	h = d;
	for(i=1; i<=ITMAX; i++)
	{
		an = -i*(i-da);
		b += 2.0;
		d = an*d+b;
		if(fabs(d) < FPMIN)
			d = FPMIN;
		c = b+an/c;
		if(fabs(c) < FPMIN)
			c = FPMIN;
		d = 1.0/d;
		del = d*c;
		h *= del;
		if(fabs(del-1.0) < EPS) 
			break;
	}
	if(i > ITMAX)
	{
		printf("gcf: a too large, ITMAX too small in gcf!\n");
	}
	*gammcf = exp(-dx+da*log(dx)-(*gln))*h;
	return;
}

/* ----------------------------------------------------------------------- */
/*                       Inverse Normal CDF Function                       */
/* Inverse of normal cumulative function.                                  */
/* ----------------------------------------------------------------------- */
double norminv(double mu, double sigma, double dx)
{
	double dy = -20;
	double dz = 20;
	double du = 0.0;
	double dv;
	int i;
	int maxIter = 10000;

	dv = normcdf(0, 1, dy);
	if(dv>dx)
		return mu+dy*sigma;
	dv = normcdf(0, 1, dz);
	if(dv<dx)
		return mu+dz*sigma;


	dv = normcdf(0, 1, du);
	i = 0;
	while( (fabs(dx-dv) > 1e-10) || (fabs(dz-dy) > 1e-10) )
	{
		if(dv>dx)
		{
			dz = du;
		}
		else
		{
			dy = du;
		}

		du = (dz+dy)/2.0;
		dv = normcdf(0, 1, du);
		i++;

		if(i>= maxIter)
		{
			printf("maxIter not sufficient to achieve the precision in the routine norminv!\n");
			break;
		}
	}

	du = mu+du*sigma;
	return du;
}