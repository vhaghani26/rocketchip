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
#include "TilingArrayLib.h"
#include "HTSequencingLib.h"

int menu_hts_countreads4refgene(int argv, char **argc);

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
	menu_hts_countreads4refgene(argv, argc);

	/* exit */
	exit(EXIT_SUCCESS);
}

int menu_hts_countreads4refgene(int argv, char **argc)
{
	/* define */
	char strDatabasePath[LINE_LENGTH];
	int nDatabaseType = 1;
	char strInputPath[LINE_LENGTH];
	char strOutputPath[LINE_LENGTH];
	char strSpecies[LINE_LENGTH];
	int nRefType = 0;
	int nUP = 0;
	int nDOWN = 0;
	int nInputType = 0; 
	int nStandardizebyT = 1;
	int nStandardizebyL = 1;

	int nResult;
	int dOK;
	int dtOK;
	int sOK;
	int iOK;
	int itOK;
	int oOK;
	int rOK;
	int upOK;
	int downOK;
	int ni;

	/* ------------------------------- */
	/*      hts_countreads4refgene     */
	/* -d database                     */
	/* -s species                      */
	/* -dt database type 0: refGene    */
	/*     1: refFlat (default)        */
	/*     2: refLocus                 */
	/* -i input bar file               */
	/* -it input type. 0: bar file     */
	/*     1: a txt file listing input bar files */
	/* -o output text file             */
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
		printf("        hts_countreads4refgene     \n");
		printf(" -d path of refgene database      \n");
		printf(" -dt database type                \n");
		printf("     0: UCSC refGene format       \n");
		printf("     1: UCSC refFlat format (default) \n");
		printf("     2: refLocus format              \n");
		printf(" -s species (current support: human, mouse, dog, cow, chicken, zebrafish, D_melanogaster, C_elegans, yeast, arabidopsis) \n");
		printf(" -i input bar file              \n");
		printf(" -it input type                 \n");
		printf("	0: bar file (default) \n");
		printf("    1: a txt file listing all input bar files. \n");
		printf(" -o output text file \n");
		printf(" -r reference type \n");
		printf("    0: TSS-up, TES-down (default)    \n");
		printf("    1: TSS-up, TSS-down     \n");
		printf("    2: TES-up, TES-down     \n");
		printf("    3: CDSS-up, CDSE-down   \n");
		printf("    4: CDSS-up, CDSS-down   \n");
		printf("    5: CDSE-up, CDSE-down   \n");
		printf(" -up up distance limit (default = 0)\n");
		printf(" -down down distance limit (default = 0)\n");
		printf(" -normt normalize by total count (1: yes (default), 0: no) \n");
		printf(" -norml normalize by total gene length (1: yes (default), 0: no) \n");
		printf(" example: \n");
		printf("    hts_countreads4refgene -d refFlat_sorted.txt -dt 1 -s human -i K36.bar -o geneK36.txt -r 0 -up 5000 -down 1000 -normt 0\n");
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
	nUP = 0;
	nDOWN = 0;
	itOK = 0;
	
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
		else if(strcmp(argc[ni], "-it") == 0)
		{
			ni++;
			nInputType = atoi(argc[ni]);
			itOK = 1;
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
		else if(strcmp(argc[ni], "-normt") == 0)
		{
			ni++;
			nStandardizebyT = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-norml") == 0)
		{
			ni++;
			nStandardizebyL = atoi(argc[ni]);
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	

	if((dOK == 0) || (sOK == 0) || (iOK == 0) || (oOK == 0) )
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = HTS_CountReads4RefGene_Main(strDatabasePath, nDatabaseType,
			strSpecies, strInputPath, nInputType, strOutputPath,
			nRefType, nUP, nDOWN, nStandardizebyT, nStandardizebyL);
	}

	/* return */
	return nResult;
}