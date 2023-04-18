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

int menu_refgene_getaffy(int argv, char **argc);

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
	menu_refgene_getaffy(argv, argc);

	/* exit */
	exit(EXIT_SUCCESS);
}

int menu_refgene_getaffy(int argv, char **argc)
{
	/* define */
	char strDatabasePath[LINE_LENGTH];
	char strInputPath[LINE_LENGTH];
	char strOutputPath[LINE_LENGTH];
	int nColumn = 0;

	int nResult;
	int dOK;
	int iOK;
	int oOK;
	int cOK;
	int ni;

	/* ------------------------------- */
	/*      refgene_getaffy            */
	/* -d affy-refgene map             */
	/* -i input file                   */
	/* -o output file                  */
	/* -c starting column (0-based)    */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("        refgene_getaffy            \n");
		printf(" -d affy-refgene map               \n");
		printf(" -i input file                     \n");
		printf(" -o output file                    \n");
		printf(" -c starting column (0-based)      \n");
		printf(" example: \n");
		printf("    refgene_getaffy -d Mouse430_2_affy2refid.txt -i target.txt -c 1 -o target_affy.txt\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	dOK = 0;
	iOK = 0;
	oOK = 0;
	cOK = 0;
	
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strDatabasePath, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strInputPath, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputPath, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-c") == 0)
		{
			ni++;
			nColumn = atoi(argc[ni]);
			cOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	

	if((dOK == 0) || (iOK == 0) || (oOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = RefGene_GetAffy_Main(strDatabasePath, 
			strInputPath, nColumn, strOutputPath);
	}

	/* return */
	return nResult;
}