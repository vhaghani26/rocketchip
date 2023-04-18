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

int menu_codingCDS(int argv, char **argc);

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
	menu_codingCDS(argv, argc);

	/* exit */
	exit(EXIT_SUCCESS);
}

int menu_codingCDS(int argv, char **argc)
{
	/* define */
	char strRefGenePath[LINE_LENGTH];
	char strOutputPath[LINE_LENGTH];
	char strGenomePath[LINE_LENGTH];
	char strSpecies[LINE_LENGTH];
	int nResult;
	int dOK;
	int oOK;
	int gOK;
	int gtOK;
	int sOK;
	int nGType;
	int ni;

	/* ------------------------------- */
	/*        genome_codingCDS         */
	/* -d genome database              */
	/* -g path of refgene database     */
	/* -gt refgene database type       */
	/* -s species                      */
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
		printf("    genome_codingCDS             \n");
		printf(" -d path of genome database      \n");
		printf(" -g path of refgene database     \n");
		printf(" -gt refgene database type       \n");
		printf("     0: UCSC refGene format(default) \n");
		printf("     1: UCSC refFlat format          \n");
		printf(" -s species                      \n");
		printf(" -o path for saving coded vector (*.cds files) \n");
		printf(" example: \n");
		printf("    genome_codingCDS -d /data/genomes/human/b35_hg17/ -g /data/genomes/human/b35_hg17/annotation/refFlat_sorted.txt -gt 1 -s human -o /data/genomes/human/b35_hg17/cds/\n\n");
		printf(" [note]                           \n");
		printf("    Before coding, there must be a chrlist.txt file to list all chromosomes that need to be coded and a chrlen.txt file for chromosome lengths \n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	dOK = 0;
	oOK = 0;
	gOK = 0;
	gtOK = 0;
	sOK = 0;
	nGType = 0;

	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strGenomePath, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-g") == 0)
		{
			ni++;
			strcpy(strRefGenePath, argc[ni]);
			gOK = 1;
		}
		else if(strcmp(argc[ni], "-gt") == 0)
		{
			ni++;
			nGType = atoi(argc[ni]);
			gtOK = 1;
		}
		else if(strcmp(argc[ni], "-s") == 0)
		{
			ni++;
			strcpy(strSpecies, argc[ni]);
			sOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputPath, argc[ni]);
			oOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	
	if((dOK == 0) || (gOK == 0) || (gtOK == 0) || (sOK == 0) || (oOK == 0) )
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = Genome_CodeCDS_FromRefGene_Main(strGenomePath, strRefGenePath, nGType, strSpecies, strOutputPath);
	}

	return nResult;
}

