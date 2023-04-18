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
#include "WorkLib.h"

int menu_refgene_getmultiortholog1way(int argv, char **argc);

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
	menu_refgene_getmultiortholog1way(argv, argc);

	/* exit */
	exit(EXIT_SUCCESS);
}

int menu_refgene_getmultiortholog1way(int argv, char **argc)
{
	/* define */
	char strTargetPath[LINE_LENGTH];
	char strOutPath[LINE_LENGTH];
	char strParamPath[LINE_LENGTH];	

	int nResult;
	int iOK;
	int dOK;
	int oOK;
	int ni;

	/* ------------------------------- */
	/*    refgene_getmultiortholog1way */
	/* -i target list                  */
	/* -d database info                */
	/* -o output                       */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    refgene_getmultiortholog1way   \n");
		printf(" -i path of target list           \n");
		printf(" -d database infomation           \n");
		printf(" -o path for saving retrieved orthologs \n");
		printf(" example: \n");
		printf("    refgene_getmultiortholog1way -i testrefid.msomap -d orthologsetting.txt -o testortho\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	iOK = 0;
	dOK = 0;
	oOK = 0;
	
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strTargetPath, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strParamPath, argc[ni]);
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

	

	if((iOK == 0) || (dOK == 0) || (oOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = RefGene_GetMultiOrtholog1way_Main(strTargetPath, strParamPath, strOutPath);
	}

	/* return */
	return nResult;
}