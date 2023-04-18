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

int menu_codingphastcons(int argv, char **argc);

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
	menu_codingphastcons(argv, argc);

	/* exit */
	exit(EXIT_SUCCESS);
}

int menu_codingphastcons(int argv, char **argc)
{
	/* define */
	char strPhastPath[LINE_LENGTH];
	char strOutputPath[LINE_LENGTH];
	char strGenomePath[LINE_LENGTH];
	int nResult;
	int dOK;
	int oOK;
	int cOK;
	int ni;
	char strExt[LINE_LENGTH];

	/* ------------------------------- */
	/*        codingphastcons          */
	/* -d genome database              */
	/* -c conservation database        */
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
		printf("       genome_codephastcons        \n");
		printf(" -d path of genome database      \n");
		printf(" -c path of phastCons score database \n");
		printf(" -o path for saving coded score (*.cs files) \n");
		printf(" -e original phastCons file extension (e.g. \"-e .pp\" means chr1.pp is the original score file for chr1). Default = no extension. \n");
		printf(" example: \n");
		printf("    codephastcons -d /data/genomes/human/b35_hg17/ -c /data/genomes/human/b35_hg17/conservation/phastcons/ -o /data/genomes/human/b35_hg17/conservation/phastcons/\n\n");
		printf(" [note]                           \n");
		printf("    Before coding, there must be a chrlist.txt file to list all chromosomes that need to be coded and a chrlen.txt file for chromosome lengths \n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	dOK = 0;
	oOK = 0;
	cOK = 0;
	strcpy(strExt, "");

	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strGenomePath, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-c") == 0)
		{
			ni++;
			strcpy(strPhastPath, argc[ni]);
			cOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputPath, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-e") == 0)
		{
			ni++;
			strcpy(strExt, argc[ni]);
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	
	if((dOK == 0) || (cOK == 0) || (oOK == 0) )
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		AdjustDirectoryPath(strGenomePath);
		AdjustDirectoryPath(strPhastPath);
		AdjustDirectoryPath(strOutputPath);
		nResult = Genome_PhastCons_To_Code_8bit_Main(strGenomePath, strPhastPath, strOutputPath, strExt);
	}

	return nResult;
}