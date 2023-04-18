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

int menu_reflocus_getneighborgenes(int argv, char **argc);

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
	menu_reflocus_getneighborgenes(argv, argc);

	/* exit */
	exit(EXIT_SUCCESS);
}

int menu_reflocus_getneighborgenes(int argv, char **argc)
{
	/* define */
	char strDatabasePath[LINE_LENGTH];
	char strInputPath[LINE_LENGTH];
	char strOutputPath[LINE_LENGTH];
	char strSpecies[LINE_LENGTH];
	char strAnnotationPath[LINE_LENGTH];
	int nGap = 0;
	int nUP = 1;
	int nDOWN = 1;

	int nResult;
	int dOK;
	int sOK;
	int iOK;
	int oOK;
	int aOK;
	int gOK;
	int upOK;
	int downOK;
	int ni;

	/* ------------------------------- */
	/*    reflocus_getneighborgenes    */
	/* -d database                     */
	/* -s species                      */
	/* -i input coordinates            */
	/* -o output                       */
	/* -a annotation                   */
	/* -g distance upper limit         */
	/* -up no. of upstream genes       */
	/* -down no. of downstream genes   */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("      reflocus_getneighborgenes     \n");
		printf(" -d path of refgene database      \n");
		printf(" -s species (current support: human, mouse, dog, cow, chicken, zebrafish, D_melanogaster, C_elegans, yeast, arabidopsis) \n");
		printf(" -i input coordinates              \n");
		printf(" -o path for saving results \n");
		printf(" -a annotation  \n");
		printf(" -g distance upper limit  \n");
		printf(" -up no. of upstream genes \n");
		printf(" -down no. of downstream genes \n");
		printf(" example: \n");
		printf("    reflocus_getneighborgenes -d /data/genomes/human/b35_hg17/annotation/refLocus_sorted.txt -s human -i target.txt -o target_gene.txt -a expression.txt -g 1000000 -up 3 -down 3\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	dOK = 0;
	sOK = 0;
	iOK = 0;
	oOK = 0;
	aOK = 0;
	gOK = 0;
	upOK = 0;
	downOK = 0;
	nUP = 1;
	nDOWN = 1;
	nGap = 1000000000;
	strcpy(strAnnotationPath, "NULL");

	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strDatabasePath, argc[ni]);
			dOK = 1;
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
		else if(strcmp(argc[ni], "-a") == 0)
		{
			ni++;
			strcpy(strAnnotationPath, argc[ni]);
			aOK = 1;
		}
		else if(strcmp(argc[ni], "-g") == 0)
		{
			ni++;
			nGap = atoi(argc[ni]);
			gOK = 1;
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

	

	if((dOK == 0) || (sOK == 0) || (iOK == 0) || (oOK == 0) )
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = RefGene_GetNeighborGenes_Main(strDatabasePath, 2,
			strSpecies, strInputPath, strOutputPath,
			strAnnotationPath, nUP, nDOWN, nGap);
	}

	/* return */
	return nResult;
}