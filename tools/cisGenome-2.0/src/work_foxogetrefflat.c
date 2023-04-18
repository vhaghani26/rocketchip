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

int ESClustGetConsRegion(int argv, char **argc);

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
	ESClustGetConsRegion(argv, argc);

	/* exit */
	exit(EXIT_SUCCESS);
}

int ESClustGetConsRegion(int argv, char **argc)
{
	/* define */
	char strInFile[LINE_LENGTH];
	char strOutFile[LINE_LENGTH];
	FILE *fpIn;
	FILE *fpOut;
	char strLine[LONG_LINE_LENGTH];

	char strRefDatabase[LINE_LENGTH];
	char strSpecies[LINE_LENGTH];
	int nSourceRefNum;
	struct tagRefGene **vSourceRefGene;
	int ni;


	/* init */
	if(argv != 3)
	{
		printf("Error: param wrong!\n");
		exit(EXIT_FAILURE);
	}
	strcpy(strInFile, argc[1]);
	strcpy(strOutFile, argc[2]);
	
	strcpy(strSpecies, "mouse");
	strcpy(strRefDatabase, "/Volumes/FTP/genomes/mouse/mm6/annotation/refFlat_sorted.txt");

	vSourceRefGene = NULL;
	vSourceRefGene = RefGene_LoadDatabase(strRefDatabase, 1, strSpecies, &nSourceRefNum);
	if(vSourceRefGene == NULL)
	{
		printf("Error: cannot open file!\n");
		exit(EXIT_FAILURE);
	}


	fpOut = NULL;
	fpOut = fopen(strOutFile, "w");
	if(fpOut == NULL)
	{
		printf("Error: cannot open file!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<nSourceRefNum; ni++)
	{
		fpIn = NULL;
		fpIn = fopen(strInFile, "r");
		if(fpIn == NULL)
		{
			printf("Error: cannot open file!\n");
			exit(EXIT_FAILURE);
		}

		while(fgets(strLine, LINE_LENGTH, fpIn) != NULL)
		{
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strLine[0] == '\0')
				continue;

			if(strcmp(strLine, vSourceRefGene[ni]->strName) == 0)
			{
				RefFlatWrite(vSourceRefGene[ni], fpOut);
				break;
			}
		}

		fclose(fpIn);
	}	

	RefGene_ClearDatabase(&vSourceRefGene, nSourceRefNum);
	

	fclose(fpOut);

	/* return */
	return PROC_SUCCESS;
}
