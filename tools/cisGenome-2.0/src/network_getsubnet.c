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

int menu_network_getsubnet(int argv, char **argc);

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
	menu_network_getsubnet(argv, argc);

	/* exit */
	exit(EXIT_SUCCESS);
}

int menu_network_getsubnet(int argv, char **argc)
{
	/* define */
	char strNetworkPath[LINE_LENGTH];
	char strAnnotationPath[LINE_LENGTH];
	char strInPath[LINE_LENGTH];
	char strOutPath[LINE_LENGTH];
	
	int nResult;
	int nOK;
	int aOK;
	int iOK;
	int oOK;
	int ni;

	/* ------------------------------- */
	/*  network_getsubnet              */
	/* -n network path                 */
	/* -a annotation path              */
	/* -i target node path             */
	/* -o output path                  */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf(" network_getsubnet          \n");
		printf(" -n network path            \n");
		printf(" -a annotation path         \n");
		printf(" -i target node path        \n");
		printf(" -o path for output         \n");
		printf(" example: \n");
		printf("    network_getsubnet -n net.sif -h annotation.txt -i targetnode.txt -o subnet.sif\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	nOK = 0;
	aOK = 0;
	iOK = 0;
	oOK = 0;
	strcpy(strAnnotationPath, "");

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
			strcpy(strAnnotationPath, argc[ni]);
			aOK = 1;
		}
		else if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strInPath, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutPath, argc[ni]);
			oOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	if((nOK == 0) || (iOK == 0) || (oOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = Network_GetSubNet_Main(strNetworkPath, strInPath, 
						   aOK, strAnnotationPath, strOutPath);
	}
	
	/* return */
	return nResult;
}