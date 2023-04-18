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

int menu_network_bind2entrez(int argv, char **argc);

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
	menu_network_bind2entrez(argv, argc);

	/* exit */
	exit(EXIT_SUCCESS);
}

int menu_network_bind2entrez(int argv, char **argc)
{
	/* define */
	char strBindPath[LINE_LENGTH];
	char strEntrezPath[LINE_LENGTH];
	char strOutPath[LINE_LENGTH];
	
	int nResult;
	int bOK;
	int eOK;
	int oOK;
	int ni;

	/* ------------------------------- */
	/*  network_bind2entrez            */
	/* -b BIND path                    */
	/* -e Entrez path                  */
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
		printf(" network_bind2entrez        \n");
		printf(" -b BIND path               \n");
		printf(" -e Entrez path             \n");
		printf(" -o output path             \n");
		printf(" example: \n");
		printf("    network_bind2entrez -b bind.dat -e gene2accession.txt -o bind.sif\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	bOK = 0;
	eOK = 0;
	oOK = 0;
	
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-b") == 0)
		{
			ni++;
			strcpy(strBindPath, argc[ni]);
			bOK = 1;
		}
		else if(strcmp(argc[ni], "-e") == 0)
		{
			ni++;
			strcpy(strEntrezPath, argc[ni]);
			eOK = 1;
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

	if((bOK == 0) || (eOK == 0) || (oOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = Network_LinkBINDEntrez_Main(strBindPath, strEntrezPath, 
						   strOutPath);
	}
	
	
	/* return */
	return nResult;
}