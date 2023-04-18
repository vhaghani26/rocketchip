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

int menu_cs_getdistn(int argv, char **argc);

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
	menu_cs_getdistn(argv, argc);

	/* exit */
	exit(EXIT_SUCCESS);
}

int menu_cs_getdistn(int argv, char **argc)
{
	/* define */
	char strDataPath[LINE_LENGTH];
	char strSpecies[LINE_LENGTH];
	char strChrListFile[LINE_LENGTH];
	char strChrLenFile[LINE_LENGTH];
	char strOutputFile[LINE_LENGTH];
	
	int ni;
	int dOK,sOK,iOK,oOK,lOK;
	int nResult;

	/* ------------------------------- */
	/*       genome_csgetdistn         */
	/* -d database                     */
	/* -s species                      */
	/* -l chrmosome length file        */
	/* -i target region                */
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
		printf("        genome_csgetdistn         \n");
		printf(" -d path of conservation score database \n");
		printf(" -s species name                 \n");
		printf("		currently supporting: human, mouse, dog, D_melanogaster\n");
		printf(" -i target chromosome \n");
		printf(" -o file for saving the output statistics\n");
		printf(" -l file for chromosome length        \n");
		printf(" example: \n");
		printf("    genome_csgetdistn -d /data/genomes/human/b33_hg15/conservation/cs/ -s human -l chrlen.txt -i chrlist.txt -o csstat.txt \n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}


	dOK = 0;
	sOK = 0;
	iOK = 0;
	oOK = 0;
	lOK = 0;

	/* set default */
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strDataPath, argc[ni]);
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
			strcpy(strChrListFile, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputFile, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-l") == 0)
		{
			ni++;
			strcpy(strChrLenFile, argc[ni]);
			lOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	if( (dOK == 0) || (sOK == 0) || (iOK == 0) || (oOK == 0) || (lOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = Genome_CS_Get_Distribution(strDataPath, strChrListFile, strChrLenFile, strOutputFile, strSpecies);
	}
		
	/* int nCount;

	nCount = Genome_CS_Get_Distribution("C:\\Projects\\research_harvard\\genomelab_project\\genomes\\crossspecies\\hg16Mm3Rn3\\hg_cs\\",
			"chrlist.txt", "C:\\Projects\\research_harvard\\genomelab_project\\genomes\\crossspecies\\hg16Mm3Rn3\\chrlen.txt",
			"C:\\Projects\\research_harvard\\genomelab_project\\genomes\\crossspecies\\hg16Mm3Rn3\\hg_cs\\csstat.txt", "human");

	return nCount; */

	return nResult;
}