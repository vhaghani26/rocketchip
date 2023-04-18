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

int menu_expression_geneselection(int argv, char **argc);

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
	menu_expression_geneselection(argv, argc);

	/* exit */
	exit(EXIT_SUCCESS);
}

int menu_expression_geneselection(int argv, char **argc)
{
	int nResult;
	int dOK;
	int aOK;
	int cOK;
	int oOK;
	int ni;
	char strDataPath[LINE_LENGTH];
	char strAnnotationPath[LINE_LENGTH];
	char strCompInfoPath[LINE_LENGTH];
	char strOutPath[LINE_LENGTH];

	/* ------------------------------- */
	/*        powexpress               */
	/* -d data file                    */
	/* -a annotation file              */
	/* -c parameter file               */
	/* -o output file                  */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("          powexpress            \n");
		printf(" -d data file      \n");
		printf(" -a annotation file      \n");
		printf(" -c parameter file      \n");
		printf(" -o output file      \n");
		printf(" example: \n");
		printf("    powexpress -d shhdata1.txt -a 3mousechipsinfo.txt -c shhcompinfo1.txt -o shh_8som_pos\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	dOK = 0;
	aOK = 0;
	cOK = 0;
	oOK = 0;
	
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strDataPath, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-a") == 0)
		{
			ni++;
			strcpy(strAnnotationPath, argc[ni]);
			aOK = 1;
		}
		else if(strcmp(argc[ni], "-c") == 0)
		{
			ni++;
			strcpy(strCompInfoPath, argc[ni]);
			cOK = 1;
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

	

	if((dOK == 0) || (aOK == 0) || (cOK == 0) || (oOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = Expression_GeneSelection_Main(strDataPath, strAnnotationPath, strCompInfoPath, strOutPath);
	}

	/* return */
	return nResult;
}
