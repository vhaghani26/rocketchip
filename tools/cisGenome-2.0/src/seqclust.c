#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include "math.h"
#include "limits.h"
#include "time.h"

#include "StringLib.h"
#include "MatrixLib.h"
#include "RandomLib.h"
#include "MathLib.h"
#include "MotifLib.h"
#include "SequenceLib.h"
#include "GenomeLib.h"
#include "MicroarrayLib.h"
#include "AffyLib.h"
#include "TilingArrayLib.h"
#include "HTSequencingLib.h"

int menu_seqclust(int argv, char **argc);

int main(int argv, char **argc)
{
	int nLen;
	int nseed;

	/* init rand */
	srand( (unsigned)time( NULL ) );
	if(strcmp(OS_SYSTEM, "WINDOWS") == 0)
	{
		nseed = (int)(rand()*1000/RAND_MAX);
	}
	else
	{
		nseed = rand()%1000;
	}
	rand_u_init(nseed);


	/* ---- */
	/* menu */
	/* ---- */
	menu_seqclust(argv, argc);

	/* exit */
	exit(EXIT_SUCCESS);
}

int menu_seqclust(int argv, char **argc)
{
	char strInputPath[LINE_LENGTH];
	char strOutputPath[LINE_LENGTH];
	char strOutputFile[LINE_LENGTH];
	int nKmin = 2;
	int nKmax = 20;
	int nBmin = 1;
	int nBmax = 1;
	int nKr = 5;
	int nMethod = 0;
	int nTransform = 0;
	double dTL = 1.0;
	int nRowStandardize = 0;
	double dCut = 0.0;
	int nSeed = 13;
	int nSkipCol = 0;
	int nMaxIter = 1000;
	double dTol = 1e-4;

	int iOK;
	int oOK;
	int dOK;
	int ni;
	int nResult;

	/* ------------------------------- */
	/*     menu_seqclust               */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    seqclust                    \n");
		printf(" -i input \n");
		printf(" -d output folder  \n");
		printf(" -o output file  \n");
		printf(" -kmin minimal cluster number, default = 2 \n");
		printf(" -kmax maximal cluster number, default = 20 \n");
		printf(" -kr number of trials (using different seeds) for each cluster number, default = 5 \n");
		printf(" -m method for clustering (0: (default) k normals; 1: k dirichlet; \n");
		printf("     2: k multivariate normals (MVN) with user specified group ID. The IDs should be given in \n");
		printf("     the data file as a line starting with GROUPID, samples within a group are modeled as MVN \n");
		printf("     samples from different groups are assumed to be independent. \n");
		printf("     3: k multivariate normals (MVN) with automatically determined group ID.\n");
		printf("     4: k multivariate normals with local factor dimension reduction \n"); 
		printf(" -bmin minimal block/factor level. Required if -m 3 or -m -4 is used. \n");
		printf("     if -m 3, then Block level = 1 means all samples are assumed to be independent. \n");
		printf("     if -m 4, then bmin should be > 0; it is the minimal number of factors used. \n");
		printf("           (note samples are independent only when number of factors = 0. \n");
		printf("     Block level = sample size means all samples are modeled jointly by a MVN. \n");
		printf(" -bmax maximal block level. \n");
		printf(" -t how to transform the data before clustering \n");
		printf("    0: no transform (default)       \n");
		printf("    1: log2 transform               \n");
		printf(" -tl truncation lower bound (default = 1.0), used to avoid log(0) when -t 1 is set \n");
		printf(" -srow 1: standardize rows (x-mean)/sd; 0: no standardization (default) \n");
		printf(" -c cutoff for reporting the clusters, default = 0.0 (i.e. no cut, classify genes to the most likely cluster) \n");
		printf(" -seed seed for generating random numbers, default = 13. \n");
		printf(" -skipcol how many columns to skip before the data start, default = 0 \n");
		printf(" -maxiter maximal iteration for each run, default = 1000 \n");
		printf(" -tol convergence criteria, default tol = 1e-4\n");
		printf(" Example: \n");
		printf("    seqclust -i inputfile.txt -d /workingfolder -o mycluster\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	iOK = 0;
	oOK = 0;
	dOK = 0;
	
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strInputPath, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strOutputPath, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputFile, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-kmin") == 0)
		{
			ni++;
			nKmin = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-kmax") == 0)
		{
			ni++;
			nKmax = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-kr") == 0)
		{
			ni++;
			nKr = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-m") == 0)
		{
			ni++;
			nMethod = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-bmin") == 0)
		{
			ni++;
			nBmin = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-bmax") == 0)
		{
			ni++;
			nBmax = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-t") == 0)
		{
			ni++;
			nTransform = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-tl") == 0)
		{
			ni++;
			dTL = atof(argc[ni]);
		}
		else if(strcmp(argc[ni], "-srow") == 0)
		{
			ni++;
			nRowStandardize = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-c") == 0)
		{
			ni++;
			dCut = atof(argc[ni]);
		}
		else if(strcmp(argc[ni], "-seed") == 0)
		{
			ni++;
			nSeed = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-skipcol") == 0)
		{
			ni++;
			nSkipCol = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-maxiter") == 0)
		{
			ni++;
			nMaxIter = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-tol") == 0)
		{
			ni++;
			dTol = atof(argc[ni]);
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	if((iOK == 0) || (oOK == 0) || (dOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = SeqClust_Main(strInputPath, strOutputPath, strOutputFile, nSkipCol,
			nKmin, nKmax, nKr, nMethod, nBmin, nBmax, nSeed,
			nTransform, dTL, nRowStandardize, dCut, 
			nMaxIter, dTol);
	}

	return nResult;
}