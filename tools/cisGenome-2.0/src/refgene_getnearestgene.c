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

int menu_refgene_getnearestgene(int argv, char **argc);

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
	menu_refgene_getnearestgene(argv, argc);

	/* exit */
	exit(EXIT_SUCCESS);
}

int menu_refgene_getnearestgene(int argv, char **argc)
{
	/* define */
	char strDatabasePath[LINE_LENGTH];
	int nDatabaseType = 0;
	char strInputPath[LINE_LENGTH];
	char strOutputPath[LINE_LENGTH];
	char strSpecies[LINE_LENGTH];
	int nRefType = 0;
	int nUP;
	int nDOWN;

	int nResult;
	int dOK;
	int dtOK;
	int sOK;
	int iOK;
	int oOK;
	int rOK;
	int upOK;
	int downOK;
	int ni;

	/* ------------------------------- */
	/*      refgene_getnearestgene     */
	/* -d database                     */
	/* -s species                      */
	/* -dt database type 0: refGene    */
	/*     1: refFlat                  */
	/*     2: refLocus                 */
	/* -i input coordinates            */
	/* -o output                       */
	/* -r reference type               */
	/*    0: TSS-up, TES-down          */
	/*    1: TSS-up, TSS-down          */
	/*    2: TES-up, TES-down          */
	/*    3: CDSS-up, CDSE-down        */
	/*    4: CDSS-up, CDSS-down        */
	/*    5: CDSE-up, CDSE-down        */
	/* -up up distance                 */
	/* -down down distance down        */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("        refgene_getnearestgene     \n");
		printf(" -d path of refgene database      \n");
		printf(" -dt database type                \n");
		printf("     0: UCSC refGene format(default) \n");
		printf("     1: UCSC refFlat format          \n");
		printf("     2: refLocus format              \n");
		printf(" -s species (current support: human, mouse, dog, cow, chicken, zebrafish, D_melanogaster, C_elegans, yeast, arabidopsis) \n");
		printf(" -i input coordinates              \n");
		printf(" -o path for saving results \n");
		printf(" -r reference type \n");
		printf("    0: TSS-up, TES-down (default)    \n");
		printf("    1: TSS-up, TSS-down     \n");
		printf("    2: TES-up, TES-down     \n");
		printf("    3: CDSS-up, CDSE-down   \n");
		printf("    4: CDSS-up, CDSS-down   \n");
		printf("    5: CDSE-up, CDSE-down   \n");
		printf(" -up up distance limit\n");
		printf(" -down down distance limit \n");
		printf(" example: \n");
		printf("    refgene_getnearestgene -d /data/genomes/human/b35_hg17/annotation/refFlat_sorted.txt -dt 1 -s human -i target.txt -o target_gene.txt -r 0 -up 5000 -down 1000\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	dOK = 0;
	dtOK = 0;
	sOK = 0;
	iOK = 0;
	oOK = 0;
	rOK = 0;
	upOK = 0;
	downOK = 0;
	nUP = 5000;
	nDOWN = 5000;
	
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
		else if(strcmp(argc[ni], "-r") == 0)
		{
			ni++;
			nRefType = atoi(argc[ni]);
			rOK = 1;
		}
		else if(strcmp(argc[ni], "-up") == 0)
		{
			ni++;
			nUP = atoi(argc[ni]);
			if(nUP < 0)
			{
				printf("Error: -up must >=0! \n");
				exit(EXIT_FAILURE);
			}

			upOK = 1;
		}
		else if(strcmp(argc[ni], "-down") == 0)
		{
			ni++;
			nDOWN = atoi(argc[ni]);
			if(nDOWN < 0)
			{
				printf("Error: -down must >=0! \n");
				exit(EXIT_FAILURE);
			}
			downOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	

	if((dOK == 0) || (sOK == 0) || (iOK == 0) || (oOK == 0) || (upOK == 0) || (downOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = RefGene_GetNearestGene_Main(strDatabasePath, nDatabaseType,
			strSpecies, strInputPath, strOutputPath,
			nRefType, nUP, nDOWN);
	}

	/* return */
	return nResult;
}
