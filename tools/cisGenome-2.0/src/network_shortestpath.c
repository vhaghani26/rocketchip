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
#include "PhysicalNetworkLib.h"
#include "WorkLib.h"

int menu_network_shortestpath(int argv, char **argc);

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
	menu_network_shortestpath(argv, argc);

	/* exit */
	exit(EXIT_SUCCESS);
}

int menu_network_shortestpath(int argv, char **argc)
{
	/* define */
	char strNetworkPath[LINE_LENGTH];
	char strAnnotPath[LINE_LENGTH];
	char strSourcePath[LINE_LENGTH];	
	char strDestPath[LINE_LENGTH];	
	char strOutPath[LINE_LENGTH];
	int nMaxIteration = 100;
	
	int nResult;
	int nOK;
	int aOK;
	int sOK;
	int dOK;
	int oOK;
	int mOK;
	int ni;

	/* ------------------------------- */
	/*  network_shortestpath           */
	/* -n network path                 */
	/* -a annotation path              */
	/* -s source node path             */
	/* -d dest node path               */
	/* -o output path                  */
	/* -m maximum iteration            */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf(" network_shortestpath       \n");
		printf(" -n network path            \n");
		printf(" -a annotation path         \n");
		printf(" -s source node path        \n");
		printf(" -d dest node path          \n");
		printf(" -o path for output         \n");
		printf(" -m maximum iteration       \n");
		printf(" example: \n");
		printf("    network_shortestpath -n net.sif -a netannotation.txt -s src.txt -d dst.txt -o out.txt -m 10\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	nOK = 0;
	aOK = 0;
	sOK = 0;
	dOK = 0;
	oOK = 0;
	mOK = 0;
	
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-n") == 0)
		{
			ni++;
			strcpy(strNetworkPath, argc[ni]);
			nOK = 1;
		}
		else if(strcmp(argc[ni], "-a") == 0)
		{
			ni++;
			strcpy(strAnnotPath, argc[ni]);
			aOK = 1;
		}
		else if(strcmp(argc[ni], "-s") == 0)
		{
			ni++;
			strcpy(strSourcePath, argc[ni]);
			sOK = 1;
		}
		else if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strDestPath, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutPath, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-m") == 0)
		{
			ni++;
			nMaxIteration = atoi(argc[ni]);
			mOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	if(nMaxIteration <= 0)
	{
		printf("Error: Max iter <= 0!\n");
		exit(EXIT_FAILURE);
	}

	if((nOK == 0) || (aOK == 0) || (sOK == 0) || (dOK == 0) || (oOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = Network_FindShortestPath_Main(strNetworkPath, strAnnotPath,
								  strSourcePath, strDestPath, 
								  strOutPath, nMaxIteration);
	}

	/* return */
	return nResult;
}