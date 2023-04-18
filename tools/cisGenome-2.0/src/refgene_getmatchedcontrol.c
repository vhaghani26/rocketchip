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

int menu_refgene_getmatchedcontrol(int argv, char **argc);

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
	menu_refgene_getmatchedcontrol(argv, argc);

	/* exit */
	exit(EXIT_SUCCESS);
}

int menu_refgene_getmatchedcontrol(int argv, char **argc)
{
	/* define */
	char strDatabasePath[LINE_LENGTH];
	int nDatabaseType = 0;
	char strInputPath[LINE_LENGTH];
	char strOutputPath[LINE_LENGTH];
	char strSpecies[LINE_LENGTH];
	char strChrLenPath[LINE_LENGTH];
	int nRepNum = 1;
	int nRegionLen = 2000;
	int nNR = 1;
	
	int nResult;
	int dOK;
	int dtOK;
	int sOK;
	int cOK;
	int iOK;
	int oOK;
	int nOK;
	int lOK;
	int nrOK;
	int ni;

	/* ------------------------------- */
	/*      refgene_getmatchedcontrol  */
	/* -d database                     */
	/* -s species                      */
	/* -dt database type 0: refGene    */
	/*     1: refFlat                  */
	/*     2: refLocus                 */
	/* -c chromosome length            */
	/* -i input coordinates            */
	/* -o output                       */
	/* -n number of replications       */
	/* -l region length                */
	/* -nr remove redundancy           */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("     refgene_getmatchedcontrol     \n");
		printf(" -d path of refgene database      \n");
		printf(" -dt database type                \n");
		printf("     0: UCSC refGene format(default) \n");
		printf("     1: UCSC refFlat format          \n");
		printf("     2: refLocus format              \n");
		printf(" -s species (current support: human, mouse, dog, cow, chicken, zebrafish, D_melanogaster, C_elegans, yeast, arabidopsis) \n");
		printf(" -c path of chromosome length \n");
		printf(" -i input coordinates         \n");
		printf(" -o file for saving results   \n");
		printf(" -n number of replications when selecting controls \n");
		printf(" -l control region length     \n");
		printf(" -nr remove redundancy or not \n");
		printf("     0: no                    \n");
		printf("     1: yes (default)         \n");
		printf(" example: \n");
		printf("    refgene_getmatchedcontrol -d /data/genomes/human/b35_hg17/annotation/refFlat_sorted.txt -dt 1 -s human -c /data/genomes/human/b35_hg17/chrlen.txt -i target.txt -o target_ct.txt -n 3 -l 2000 -nr 1\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	dOK = 0;
	dtOK = 0;
	sOK = 0;
	cOK = 0;
	iOK = 0;
	oOK = 0;
	nOK = 0;
	lOK = 0;
	nrOK = 0;
	
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strDatabasePath, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-dt") == 0)
		{
			ni++;
			nDatabaseType = atoi(argc[ni]);
			dtOK = 1;
		}
		else if(strcmp(argc[ni], "-s") == 0)
		{
			ni++;
			strcpy(strSpecies, argc[ni]);
			sOK = 1;
		}
		else if(strcmp(argc[ni], "-c") == 0)
		{
			ni++;
			strcpy(strChrLenPath, argc[ni]);
			cOK = 1;
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
		else if(strcmp(argc[ni], "-n") == 0)
		{
			ni++;
			nRepNum = atoi(argc[ni]);
			nOK = 1;

			if(nRepNum <= 0)
			{
				printf("Error: -n must >0! \n");
				exit(EXIT_FAILURE);
			}
		}
		else if(strcmp(argc[ni], "-l") == 0)
		{
			ni++;
			nRegionLen = atoi(argc[ni]);
			lOK = 1;

			if(nRegionLen <= 0)
			{
				printf("Error: -l must >0! \n");
				exit(EXIT_FAILURE);
			}
		}
		else if(strcmp(argc[ni], "-nr") == 0)
		{
			ni++;
			nNR = atoi(argc[ni]);
			nrOK = 1;

			if( (nNR < 0) || (nNR > 1) )
			{
				printf("Error: -nr must be 0 or 1! \n");
				exit(EXIT_FAILURE);
			}
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	if((dOK == 0) || (dtOK == 0) || (sOK == 0) || (cOK == 0) || (iOK == 0) || (oOK == 0) )
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = RefGene_GetMatchedControl_Main(strDatabasePath, nDatabaseType,
			strSpecies, strChrLenPath, strInputPath, strOutputPath,
			nRepNum, nRegionLen, nNR);
	}

	/* return */
	return nResult;
}