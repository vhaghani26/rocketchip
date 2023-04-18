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

int menu_refgene_gettargettssaround(int argv, char **argc);

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
	menu_refgene_gettargettssaround(argv, argc);

	/* exit */
	exit(EXIT_SUCCESS);
}

int menu_refgene_gettargettssaround(int argv, char **argc)
{
	/* define */
	char strRefGenePath[LINE_LENGTH];
	char strTargetPath[LINE_LENGTH];
	char strOutputPath[LINE_LENGTH];
	char strSpecies[LINE_LENGTH];
	char strChrLen[LINE_LENGTH];
	int nTSSUP;
	int nTSSDOWN;

	int nResult;
	int dOK;
	int tOK;
	int oOK;
	int sOK;
	int upOK;
	int downOK;
	int cOK;
	int ni;

	/* ------------------------------- */
	/*        refgene_gettssaround     */
	/* -d database                     */
	/* -t target list                  */
	/* -o output                       */
	/* -s species                      */
	/* -up TSS up                      */
	/* -down TSS down                  */
	/* -c chromosome length            */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("          refgene_gettssaround    \n");
		printf(" -d path of refgene database      \n");
		printf(" -t path of target list           \n");
		printf(" -o path for saving retrieved target coordinates \n");
		printf(" -s species (current support: human, mouse, dog, cow, chicken, zebrafish, D_melanogaster, C_elegans, yeast, arabidopsis) \n");
		printf(" -up TSS up \n");
		printf(" -down TSS down \n");
		printf(" -c chromosome length\n");
		printf(" example: \n");
		printf("    refgene_gettssaround -d /data/genomes/human/b35_hg17/annotation/xenoRefGene_sorted.txt -t /data/genomes/human/b35_hg17/annotation/testrefid.txt -o humcod.txt -s human -up 5000 -down 1000 -c /data/genomes/human/b35_hg17/chrlen.txt\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	dOK = 0;
	tOK = 0;
	oOK = 0;
	sOK = 0;
	upOK = 0;
	downOK = 0;
	cOK = 0;

	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strRefGenePath, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-t") == 0)
		{
			ni++;
			strcpy(strTargetPath, argc[ni]);
			tOK = 1;
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
		else if(strcmp(argc[ni], "-up") == 0)
		{
			ni++;
			nTSSUP = atoi(argc[ni]);
			upOK = 1;
		}
		else if(strcmp(argc[ni], "-down") == 0)
		{
			ni++;
			nTSSDOWN = atoi(argc[ni]);
			downOK = 1;
		}
		else if(strcmp(argc[ni], "-c") == 0)
		{
			ni++;
			strcpy(strChrLen, argc[ni]);
			cOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	

	if((dOK == 0) || (tOK == 0) || (oOK == 0) || (sOK == 0) || (upOK == 0) || (downOK == 0) || (cOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = RefGene_GetTargetTSSAround(strRefGenePath, strTargetPath,
			nTSSUP, nTSSDOWN, strSpecies, strChrLen, strOutputPath);
	}

	/* return */
	return nResult;
}