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

int menu_network_createortholognet(int argv, char **argc);

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
	menu_network_createortholognet(argv, argc);

	/* exit */
	exit(EXIT_SUCCESS);
}

int menu_network_createortholognet(int argv, char **argc)
{
	/* define */
	char strNetworkPath[LINE_LENGTH];
	char strHomoloPath[LINE_LENGTH];
	int nSrcSpecies,nDestSpecies;
	char strOutPath[LINE_LENGTH];
	
	int nResult;
	int nOK;
	int hOK;
	int sOK;
	int dOK;
	int oOK;
	int ni;

	/* ------------------------------- */
	/*  network_createortholognet      */
	/* -n network path                 */
	/* -h homologene path              */
	/* -s source species ID            */
	/* -d dest species ID              */
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
		printf(" network_createortholognet  \n");
		printf(" -n network path            \n");
		printf(" -h homologene path         \n");
		printf(" -s source species ID       \n");
		printf(" -d dest species ID         \n");
		printf(" -o path for output         \n");
		printf(" example: \n");
		printf("    network_createortholognet -n humannet.sif -h homologene.data -s 9606 -d 10090 -o mousenet.sif\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	nOK = 0;
	hOK = 0;
	sOK = 0;
	dOK = 0;
	oOK = 0;
	
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-n") == 0)
		{
			ni++;
			strcpy(strNetworkPath, argc[ni]);
			nOK = 1;
		}
		else if(strcmp(argc[ni], "-h") == 0)
		{
			ni++;
			strcpy(strHomoloPath, argc[ni]);
			hOK = 1;
		}
		else if(strcmp(argc[ni], "-s") == 0)
		{
			ni++;
			nSrcSpecies = atoi(argc[ni]);
			sOK = 1;
		}
		else if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			nDestSpecies = atoi(argc[ni]);
			dOK = 1;
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

	if((nOK == 0) || (hOK == 0) || (sOK == 0) || (dOK == 0) || (oOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = Network_CreateOrthogloNet_Main(strNetworkPath, strHomoloPath,
								  nSrcSpecies, nDestSpecies, strOutPath);
	}
	
	/* return */
	return nResult;
}