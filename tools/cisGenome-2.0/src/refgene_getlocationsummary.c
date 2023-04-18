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

int menu_refgene_getlocationsummary(int argv, char **argc);

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
	menu_refgene_getlocationsummary(argv, argc);

	/* exit */
	exit(EXIT_SUCCESS);
}

int menu_refgene_getlocationsummary(int argv, char **argc)
{
	/* define */
	char strDatabasePath[LINE_LENGTH];
	int nDatabaseType = 0;
	char strInputPath[LINE_LENGTH];
	char strOutputPath[LINE_LENGTH];
	char strSpecies[LINE_LENGTH];
	int nInputType = 0;
	
	int nResult;
	int dOK;
	int dtOK;
	int sOK;
	int iOK;
	int oOK;
	int rOK;
	int ni;

	/* ------------------------------- */
	/*   refgene_getlocationsummary    */
	/* -d database                     */
	/* -s species                      */
	/* -dt database type 0: refGene    */
	/*     1: refFlat                  */
	/*     2: refLocus                 */
	/* -i input coordinates            */
	/* -o output                       */
	/* -r input type (0: cod, 1: bed   */
	/*      2: codp, 3: bedp)          */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    refgene_getlocationsummary    \n");
		printf(" -d path of refgene database      \n");
		printf(" -dt database type                \n");
		printf("     0: UCSC refGene format(default) \n");
		printf("     1: UCSC refFlat format          \n");
		printf("     2: refLocus format              \n");
		printf(" -s species (current support: human, mouse, dog, cow, chicken, zebrafish, D_melanogaster, C_elegans, yeast, arabidopsis) \n");
		printf(" -i input coordinates              \n");
		printf(" -o path for saving results \n");
		printf(" -r input type \n");
		printf("    0: cod; 1: bed; 2: codp; 3: bedp \n");
		printf(" example: \n");
		printf("    refgene_getlocationsummary -d /data/genomes/human/b35_hg17/annotation/refFlat_sorted.txt -dt 1 -s human -i target.txt -o target_gene.txt -r 0\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	dOK = 0;
	dtOK = 0;
	sOK = 0;
	iOK = 0;
	oOK = 0;
	rOK = 0;
	
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
			nInputType = atoi(argc[ni]);
			rOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	

	if((dOK == 0) || (sOK == 0) || (iOK == 0) || (oOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = RefGene_GetLocationSummary_Main(strDatabasePath, nDatabaseType,
			strSpecies, strInputPath, nInputType, strOutputPath);
	}

	/* return */
	return nResult;
}