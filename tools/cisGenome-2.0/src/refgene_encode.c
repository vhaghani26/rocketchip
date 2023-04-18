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

int menu_refgene_sort(int argv, char **argc);

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
	menu_refgene_sort(argv, argc);

	/* exit */
	exit(EXIT_SUCCESS);
}

int menu_refgene_sort(int argv, char **argc)
{
	/* define */
	char strRefGenePath[LINE_LENGTH];
	char strOutputPath[LINE_LENGTH];
	char strSpecies[LINE_LENGTH];
	int nChrNum;
	int nResult;
	int dOK;
	int oOK;
	int sOK;
	int numOK;
	int ni;

	/* ------------------------------- */
	/*        refgene_sort             */
	/* -d database                     */
	/* -o output                       */
	/* -s species                      */
	/* -n number of chromosome         */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("          refgene_encode          \n");
		printf(" -d path of refgene database      \n");
		printf(" -o path for saving coded and sorted refgene \n");
		printf(" -s species \n");
		printf(" -n number of chromosome \n");
		printf(" example: \n");
		printf("    refgene_encode -d /data/genomes/human/b35_hg17/annotation/refGene.txt -o /data/genomes/human/b35_hg17/annotation/refGene_sorted.txt -s human -n 24\n\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	dOK = 0;
	oOK = 0;
	sOK = 0;
	numOK = 0;

	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strRefGenePath, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputPath, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-s") == 0)
		{
			ni++;
			strcpy(strSpecies, argc[ni]);
			sOK = 1;
		}
		else if(strcmp(argc[ni], "-n") == 0)
		{
			ni++;
			nChrNum = atoi(argc[ni]);
			numOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	

	if((dOK == 0) || (oOK == 0) || (sOK == 0) || (numOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = RefGene_Convert_UCSC_To_Lab(strRefGenePath, strOutputPath, strSpecies, nChrNum);
	}

	return nResult;
}